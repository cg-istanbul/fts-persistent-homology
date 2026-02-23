# scripts/report/barcodes_io.sage

import os, json

def _jsonable(x):
    if x is None or isinstance(x, (bool, int, float, str)):
        return x

    try:
        if hasattr(x, "is_integer") and x.is_integer():
            return int(x)
    except Exception:
        pass

    try:
        if isinstance(x, int):
            return int(x)
    except Exception:
        pass

    try:
        return int(x)
    except Exception:
        pass

    try:
        return float(x)
    except Exception:
        pass

    if isinstance(x, dict):
        return {str(_jsonable(k)): _jsonable(v) for k, v in x.items()}

    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]

    if isinstance(x, (set, frozenset)):
        arr = [_jsonable(v) for v in x]
        try:
            return sorted(arr)
        except Exception:
            return arr

    try:
        if hasattr(x, "nrows") and hasattr(x, "ncols") and hasattr(x, "rows"):
            return [_jsonable(r) for r in x.rows()]
    except Exception:
        pass

    try:
        if hasattr(x, "list"):
            return _jsonable(x.list())
    except Exception:
        pass

    return str(x)

def write_bars_jsonl(bars, outfile):
    outdir = os.path.dirname(outfile)
    if outdir:
        os.makedirs(outdir, exist_ok=True)

    with open(outfile, "w", encoding="utf-8") as f:
        for b in bars:
            f.write(json.dumps(_jsonable(b), ensure_ascii=False) + "\n")

    return outfile


def barcode_outfile(name, dim, outdir="examples/barcodes"):
    os.makedirs(outdir, exist_ok=True)
    return os.path.join(outdir, f"{name}_beta{int(dim)}.svg")

def barcode_jsonl_path(name, dim, outdir="examples/barcodes"):
    os.makedirs(outdir, exist_ok=True)
    return os.path.join(outdir, f"{name}_beta{int(dim)}.jsonl")

def barcode_title(name, dim):
    return f"{name}: β{int(dim)} barcode"
