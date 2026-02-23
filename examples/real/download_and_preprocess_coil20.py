import os, re, sys, zipfile, math
import urllib.request
import numpy as np

# ----------------------------
# Config
# ----------------------------
COIL20_PROC_URLS = [
    # Primary: Columbia CAVE (referenced from COIL-20 page). :contentReference[oaicite:2]{index=2}
    "https://www.cs.columbia.edu/CAVE/databases/SLAM_coil-20_coil-100/coil-20/coil-20-proc.zip",

    # Optional fallback mirrors (may or may not be accessible in your environment).
    # Academic Torrents has a COIL-20 entry. :contentReference[oaicite:3]{index=3}
    # NOTE: this is not the same URL as above; keep as a last resort.
]

DEFAULT_IMG_SIZE = 32
DEFAULT_PCA_DIM = 30
DEFAULT_N_PER_OBJECT = 50  # 20 * 50 = 1000

# ----------------------------
# Small utilities
# ----------------------------
def _mkdirp(p: str) -> None:
    os.makedirs(p, exist_ok=True)

def _download(url: str, dst_path: str) -> None:
    print(f"[download] {url}")
    _mkdirp(os.path.dirname(dst_path))
    with urllib.request.urlopen(url) as r, open(dst_path, "wb") as f:
        total = r.length if hasattr(r, "length") else None
        done = 0
        chunk = 1 << 20
        while True:
            b = r.read(chunk)
            if not b:
                break
            f.write(b)
            done += len(b)
            if total:
                pct = 100.0 * done / total
                print(f"\r  {done/1e6:8.1f} MB / {total/1e6:8.1f} MB  ({pct:5.1f}%)", end="")
        print()
    print(f"[ok] saved to {dst_path}")

def _extract_zip(zip_path: str, out_dir: str) -> None:
    print(f"[extract] {zip_path} -> {out_dir}")
    _mkdirp(out_dir)
    with zipfile.ZipFile(zip_path, "r") as z:
        z.extractall(out_dir)
    print("[ok] extracted")

def _load_image_gray_resize(path: str, size: int) -> np.ndarray:
    # Try Pillow first (most reliable).
    try:
        from PIL import Image
        im = Image.open(path).convert("L").resize((size, size))
        arr = np.asarray(im, dtype=np.float32)
        return arr
    except Exception:
        pass

    # Fallback: imageio
    try:
        import imageio.v2 as imageio
        im = imageio.imread(path)
        if im.ndim == 3:
            im = im[..., 0]
        # naive resize if Pillow missing: use numpy slicing / nearest neighbor
        im = im.astype(np.float32)
        h, w = im.shape[:2]
        ys = (np.linspace(0, h - 1, size)).astype(int)
        xs = (np.linspace(0, w - 1, size)).astype(int)
        return im[np.ix_(ys, xs)]
    except Exception as e:
        raise RuntimeError(
            "Could not read images. Install pillow or imageio.\n"
            "Try: pip install pillow\n"
        ) from e

def _pca_svd(X: np.ndarray, p: int) -> np.ndarray:
    """
    PCA via SVD on centered data.
    Returns X projected to p dims (n x p).
    """
    X = X.astype(np.float64)
    mu = X.mean(axis=0, keepdims=True)
    Xc = X - mu
    # economy SVD
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    V = Vt[:p].T  # d x p
    return (Xc @ V).astype(np.float64)

def _pick_evenly_spaced(indices, k):
    if k >= len(indices):
        return list(indices)
    # evenly spaced selection across sorted indices
    pos = np.linspace(0, len(indices) - 1, k)
    pos = np.round(pos).astype(int)
    return [indices[i] for i in pos]

# ----------------------------
# Main
# ----------------------------
def main():
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    out_dir = os.path.join(repo_root, "data", "coil20")
    raw_dir = os.path.join(out_dir, "raw")
    zip_path = os.path.join(raw_dir, "coil-20-proc.zip")

    img_size = int(os.environ.get("COIL20_IMG_SIZE", DEFAULT_IMG_SIZE))
    pca_dim = int(os.environ.get("COIL20_PCA_DIM", DEFAULT_PCA_DIM))
    n_per_obj = int(os.environ.get("COIL20_PER_OBJECT", DEFAULT_N_PER_OBJECT))

    _mkdirp(out_dir)
    _mkdirp(raw_dir)

    if not os.path.exists(zip_path):
        ok = False
        for url in COIL20_PROC_URLS:
            try:
                _download(url, zip_path)
                ok = True
                break
            except Exception as e:
                print(f"[warn] failed: {e}")
        if not ok:
            raise RuntimeError("All download URLs failed. Try downloading manually and place the zip at:\n"
                               f"  {zip_path}\n")
    else:
        print(f"[skip] zip already exists: {zip_path}")

    extract_dir = os.path.join(raw_dir, "coil-20-proc")
    if not os.path.isdir(extract_dir):
        _extract_zip(zip_path, raw_dir)

    # Find images: accept common COIL naming patterns
    # Typical: obj<id>__<angle>.png or obj<id>__<angle>.pgm/ppm
    img_paths = []
    for root, _, files in os.walk(raw_dir):
        for fn in files:
            if fn.lower().endswith((".png", ".pgm", ".ppm", ".jpg", ".jpeg", ".bmp", ".tif", ".tiff")):
                img_paths.append(os.path.join(root, fn))

    if not img_paths:
        raise RuntimeError(f"No images found under {raw_dir}. Check zip contents.")

    pat = re.compile(r"obj(\d+)[_]+(\d+)", re.IGNORECASE)

    by_obj = {i: [] for i in range(1, 21)}
    unparsed = 0
    for pth in img_paths:
        m = pat.search(os.path.basename(pth))
        if not m:
            unparsed += 1
            continue
        oid = int(m.group(1))
        ang = int(m.group(2))
        if 1 <= oid <= 20:
            by_obj[oid].append((ang, pth))

    if sum(len(v) for v in by_obj.values()) == 0:
        raise RuntimeError("Could not parse COIL filenames. Please inspect names and adjust regex in script.")

    if unparsed:
        print(f"[note] {unparsed} image files did not match obj<id>__<angle> pattern; ignored.")

    # Sort and select n_per_obj evenly across angles
    selected = []
    for oid in range(1, 21):
        lst = sorted(by_obj[oid], key=lambda t: t[0])
        if not lst:
            raise RuntimeError(f"Missing images for object {oid}.")
        idxs = list(range(len(lst)))
        pick = _pick_evenly_spaced(idxs, n_per_obj)
        for j in pick:
            ang, pth = lst[j]
            selected.append((oid, ang, pth))

    # Stable ordering: object, angle
    selected.sort(key=lambda t: (t[0], t[1]))
    n = len(selected)
    print(f"[select] {n} images total ({n_per_obj} per object)")

    # Build raw feature matrix
    X = np.zeros((n, img_size * img_size), dtype=np.float32)
    meta = []
    for idx, (oid, ang, pth) in enumerate(selected):
        arr = _load_image_gray_resize(pth, img_size)
        X[idx, :] = arr.reshape(-1)
        meta.append((idx, oid, ang, os.path.relpath(pth, repo_root)))

    print(f"[matrix] X = {X.shape}")

    # PCA
    if pca_dim > 0 and pca_dim < X.shape[1]:
        Xp = _pca_svd(X, pca_dim)
        print(f"[pca] projected to {Xp.shape}")
    else:
        Xp = X.astype(np.float64)

    # Save
    X_csv = os.path.join(out_dir, "coil20_n1000.csv")
    M_csv = os.path.join(out_dir, "coil20_n1000_meta.csv")

    np.savetxt(X_csv, Xp, delimiter=",")
    with open(M_csv, "w", encoding="utf-8") as f:
        f.write("idx,obj_id,angle,relpath\n")
        for row in meta:
            f.write("{},{},{},{}\n".format(*row))

    print(f"[write] {X_csv}")
    print(f"[write] {M_csv}")
    print("[done]")

if __name__ == "__main__":
    main()

