# scripts/report/barcodes_svg.sage

import os


def _safe_int(x, default=0):
    try:
        return int(x)
    except Exception:
        return int(default)

def _as_float(x, default=None):
    try:
        if x is None:
            return default
        return float(x)
    except Exception:
        return default

def _format_ratio(a, b):
    a = _safe_int(a, 0)
    b = _safe_int(b, 0)
    return "{}{}".format(a, "" if b == 0 else "/{}".format(b))

def _format_pct(k, n):
    k = _safe_int(k, 0)
    n = _safe_int(n, 0)
    if n <= 0:
        return str(k)
    pct = 100.0 * float(k) / float(n)
    if pct >= 10:
        return "{} ({:.0f}%)".format(k, pct)
    return "{} ({:.1f}%)".format(k, pct)

def write_enhanced_barcode_svg(
    bars,
    nX,
    outfile,
    title="Enhanced barcode",
    max_bars=200,
    x_padding=60,
    y_padding=40,
    bar_h=10,
    bar_gap=10,
    width=1200,
):
    if outfile is None:
        raise ValueError("outfile must be provided (use barcode_svg_path(name, dim))")

    outdir = os.path.dirname(outfile)
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        
    nX = _safe_int(nX, default=0)

    finite_ts = []
    for b in bars:
        bt = _as_float(b.get("birth_time", None), None)
        dt = _as_float(b.get("death_time", None), None)
        if bt is not None: finite_ts.append(bt)
        if dt is not None: finite_ts.append(dt)

    has_inf = False
    end_taus = []
    for b in bars:
        dt = _as_float(b.get("death_time", None), None)
        if dt is None:
            has_inf = True
            meta = b.get("meta", {})
            mend = meta.get("end", {}) if isinstance(meta, dict) else {}
            et = _as_float(mend.get("tau", None), None)
            if et is not None:
                end_taus.append(et)

    if has_inf and len(end_taus) > 0:
        finite_ts.append(max(end_taus))

    if len(finite_ts) == 0:
        tmin, tmax = 0.0, 1.0
    else:
        tmin_raw, tmax_raw = min(finite_ts), max(finite_ts)
        if tmax_raw <= tmin_raw:
            tmax_raw = tmin_raw + 1.0

        pad = 0.03 * (tmax_raw - tmin_raw)
        if pad <= 0:
            pad = 1.0
        tmin, tmax = tmin_raw, tmax_raw + pad

    t_inf = tmax

    max_bars_i = int(max_bars) if max_bars is not None else None
    bars2 = list(bars)[:int(max_bars_i)]
    W = int(width)
    RIGHT_GUTTER = 420
    W = int(W) + int(RIGHT_GUTTER)

    dim0 = int(bars2[0].get("dim", 0)) if bars2 else 0
    if dim0 == 1:
        bar_h = max(bar_h, 8)
        bar_gap = max(bar_gap, 14)


    def x_of(t):
        left = float(x_padding)
        right = float(W - 350)
        return left + (float(t) - tmin) * (right - left) / (tmax - tmin)

    lines = []
    lines.append('<?xml version="1.0" encoding="UTF-8"?>')
    H = y_padding*2 + len(bars2)*(bar_h + bar_gap) + 30

    lines.append('<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}">'.format(W, H))
    lines.append('<rect x="0" y="0" width="{}" height="{}" fill="white"/>'.format(W, H))
    lines.append('<text x="{}" y="{}" font-family="monospace" font-size="16">{}</text>'.format(
        x_padding, y_padding-15, title.replace("&","&amp;")
    ))

    y_axis = H - y_padding + 10
    lines.append('<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="black" stroke-width="1"/>'.format(
        x_of(tmin), y_axis, x_of(tmax), y_axis
    ))
    lines.append('<text x="{}" y="{}" font-family="monospace" font-size="12">{:.3f}</text>'.format(
        x_of(tmin), y_axis+15, tmin
    ))
    lines.append('<text x="{}" y="{}" font-family="monospace" font-size="12">{:.3f}</text>'.format(
        x_of(tmax)-30, y_axis+15, tmax
    ))

    y0 = y_padding + 10
    for idx, b in enumerate(bars2):
        bt = _as_float(b.get("birth_time", None), tmin)
        dt = _as_float(b.get("death_time", None), None)
        dt_draw = t_inf if dt is None else dt

        if dt is not None and abs(float(dt_draw) - float(t_inf)) < 1e-12:
            try:
                print("[svg][WARN] finite bar reaches t_inf boundary:", b)
            except Exception:
                pass

        x1 = x_of(bt)
        x2 = x_of(dt_draw)
        y = y0 + idx*(bar_h + bar_gap)

        if dt is None:
            lines.append('<line x1="{:.2f}" y1="{:.2f}" x2="{:.2f}" y2="{:.2f}" '
                         'stroke="black" stroke-width="3"/>'.format(
                             x1, y, x2, y
                         ))

            ah = 5.0  
            aw = 8.0  
            x_tip = x2
            x_base = x2 - aw
            lines.append('<polygon points="{:.2f},{:.2f} {:.2f},{:.2f} {:.2f},{:.2f}" fill="black"/>'.format(
                x_tip, y, x_base, y-ah, x_base, y+ah
            ))
        else:
            lines.append('<line x1="{:.2f}" y1="{:.2f}" x2="{:.2f}" y2="{:.2f}" stroke="black" stroke-width="3"/>'.format(
                x1, y, x2, y
            ))

            cap = 4.0
            lines.append('<line x1="{:.2f}" y1="{:.2f}" x2="{:.2f}" y2="{:.2f}" stroke="black" stroke-width="2"/>'.format(
                x2, y-cap, x2, y+cap
            ))

        bs = b.get("birth_step", "")
        ds = b.get("death_step", None)
        ds_s = "∞" if ds is None else str(ds)
        lines.append('<text x="{}" y="{}" font-family="monospace" font-size="11">{}</text>'.format(
            5, y+4, "{}→{}".format(bs, ds_s)
        ))

        meta = b.get("meta", {})
        if isinstance(meta, dict) and ("core" in meta):
            core = meta.get("core", {})
            diam = core.get("suppX_diameter", None)
            ratio = core.get("suppX_diameter_ratio", None)
            supp = core.get("suppX_size", None)
            loc = core.get("locality", None)
            diam_s = "-" if diam is None else "{:.2f}".format(float(diam))
            ratio_s = "-" if ratio is None else "{:.2f}".format(float(ratio))
            supp_s = "-" if supp is None else _format_pct(supp, nX)
            loc_s = "-" if loc is None else str(loc)
            lbl = "core_diam={}  ({}g)  X={}  {}".format(diam_s, ratio_s, supp_s, loc_s)
        else:
            mb = meta.get("birth", {}) if isinstance(meta, dict) else {}

            dim_bar = _safe_int(b.get("dim", 0), 0)
            is_beta1 = (dim_bar == 1)
            
            if is_beta1:
                q = _safe_int(mb.get("q", 0))
                btype = mb.get("birth_type", "")
                h1i = mb.get("h1_basis_index", None)
                
                suppC = mb.get("cycle_support_classes", [])
                suppE = mb.get("cycle_edges", [])
                suppP = mb.get("cycle_support_points_prefix", [])

                nC = len(suppC) if isinstance(suppC, list) else 0
                nE = len(suppE) if isinstance(suppE, list) else 0
                nP = len(suppP) if isinstance(suppP, list) else 0
                
                parts = []
                parts.append(f"q={q}")
                if btype:
                    parts.append(str(btype))
                if h1i is not None:
                    parts.append(f"i={_safe_int(h1i)}")
                parts.append(f"|C|={nC}")
                parts.append(f"|E|={nE}")
                parts.append(f"|P|≈{nP}")

                lbl = "  ".join(parts)

            else:
                q = _safe_int(mb.get("q", 0))
                C = _safe_int(mb.get("comp_num_classes", 0))
                supp = _safe_int(mb.get("suppX_size", 0))
                rep_cls = _safe_int(mb.get("rep_class_size", 0))
                lbl = "C/q={}  X={}  repcls={}".format(_format_ratio(C, q), _format_pct(supp, nX), rep_cls)
        lines.append('<text x="{}" y="{}" font-family="monospace" font-size="12">{}</text>'.format(
            int(W-340), y+4, lbl.replace("&","&amp;")
        ))

    lines.append('</svg>')

    with open(outfile, "w") as f:
        f.write("\n".join(lines))

    return outfile
