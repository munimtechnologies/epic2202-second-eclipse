# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""Pixel-level vetting of flagged TESS dips (pre-registered in AMENDMENT-2026-11-25-pixvet.md).

Pure-stdlib helpers on small images (lists of lists, row-major [y][x]); the job (job_pixvet.py)
does the I/O with numpy/astropy and calls these for every decision, so the rules are unit-tested.
"""
import math

DEPTH_MIN = 0.03          # below this the flagged dip is not reproduced
CORR_SYS = 0.3            # PSF correlation below this -> SYSTEMATIC
CORR_ON = 0.7             # PSF correlation at least this -> can be ON_TARGET
RATIO_ON = 0.75           # peak-pixel / aperture depth ratio needed for ON_TARGET
CTRL_SYS = 0.5            # median control depth >= CTRL_SYS * target depth -> SYSTEMATIC
CTRL_ON = 0.25            # |median control depth| <= CTRL_ON * target depth needed for ON_TARGET
MIN_CONTROLS = 2


def median(xs):
    v = sorted(x for x in xs if x is not None and not (isinstance(x, float) and math.isnan(x)))
    n = len(v)
    if n == 0:
        return float("nan")
    m = n // 2
    return float(v[m]) if n % 2 else 0.5 * (v[m - 1] + v[m])


def percentile(xs, p):
    """Linear interpolation between closest ranks (numpy's default)."""
    v = sorted(x for x in xs if x is not None and not (isinstance(x, float) and math.isnan(x)))
    if not v:
        return float("nan")
    r = (len(v) - 1) * p / 100.0
    lo = int(math.floor(r))
    hi = min(lo + 1, len(v) - 1)
    return v[lo] + (v[hi] - v[lo]) * (r - lo)


def pearson(a, b):
    pairs = [(x, y) for x, y in zip(a, b) if not (math.isnan(x) or math.isnan(y))]
    n = len(pairs)
    if n < 3:
        return float("nan")
    ma = sum(p[0] for p in pairs) / n
    mb = sum(p[1] for p in pairs) / n
    sab = sum((p[0] - ma) * (p[1] - mb) for p in pairs)
    saa = sum((p[0] - ma) ** 2 for p in pairs)
    sbb = sum((p[1] - mb) ** 2 for p in pairs)
    if saa <= 0 or sbb <= 0:
        return float("nan")
    return sab / math.sqrt(saa * sbb)


def window_idx(times, t_mid, width, out_lo=0.15, out_hi=0.50):
    """In-dip indices |t - t_mid| <= width/2; comparison indices out_lo <= |t - t_mid| <= out_hi."""
    ins, outs = [], []
    for i, t in enumerate(times):
        d = abs(t - t_mid)
        if d <= width / 2.0:
            ins.append(i)
        elif out_lo <= d <= out_hi:
            outs.append(i)
    return ins, outs


def find_target(img):
    """Brightest pixel of the central 3x3 (the cutout is centred on the target)."""
    ny, nx = len(img), len(img[0])
    cy, cx = ny // 2, nx // 2
    best = None
    for y in range(cy - 1, cy + 2):
        for x in range(cx - 1, cx + 2):
            v = img[y][x]
            if v is None or math.isnan(v):
                continue
            if best is None or v > best[0]:
                best = (v, y, x)
    return (best[1], best[2]) if best else (cy, cx)


def bkg_mask(img, sy, sx, ry, rx, half=5, min_dist=4.0, pct=60.0):
    """Background pixels as in io_tess.photometry: inside the (2*half+1)^2 region centred on (ry, rx)
    (clipped to the image), at least min_dist pixels from the star (sy, sx), and below the pct-th
    percentile of the region's median-image values."""
    ny, nx = len(img), len(img[0])
    y0, y1 = max(0, ry - half), min(ny - 1, ry + half)
    x0, x1 = max(0, rx - half), min(nx - 1, rx + half)
    region = [img[y][x] for y in range(y0, y1 + 1) for x in range(x0, x1 + 1)]
    thr = percentile(region, pct)
    out = []
    for y in range(y0, y1 + 1):
        for x in range(x0, x1 + 1):
            v = img[y][x]
            if v is None or math.isnan(v):
                continue
            if math.hypot(y - sy, x - sx) >= min_dist and v < thr:
                out.append((y, x))
    return out


def ap_sum(img, y, x, r=1):
    s = 0.0
    for yy in range(y - r, y + r + 1):
        for xx in range(x - r, x + r + 1):
            if 0 <= yy < len(img) and 0 <= xx < len(img[0]):
                v = img[yy][xx]
                if v is not None and not math.isnan(v):
                    s += v
    return s


def local_maxima(img, ty, tx, min_dist=4.0, edge=2, min_ratio=0.05, max_n=8):
    """Control stars: pixels brighter than all 8 neighbours, >= min_dist from the target, >= edge
    pixels from the border, with 3x3 sum >= min_ratio x the target's 3x3 sum. Brightest first."""
    ny, nx = len(img), len(img[0])
    ref = ap_sum(img, ty, tx)
    found = []
    for y in range(edge, ny - edge):
        for x in range(edge, nx - edge):
            v = img[y][x]
            if v is None or math.isnan(v) or math.hypot(y - ty, x - tx) < min_dist:
                continue
            nb = [img[y + dy][x + dx] for dy in (-1, 0, 1) for dx in (-1, 0, 1) if dy or dx]
            if all(n is not None and not math.isnan(n) and v > n for n in nb):
                s = ap_sum(img, y, x)
                if ref > 0 and s >= min_ratio * ref:
                    found.append((s, y, x))
    found.sort(reverse=True)
    return [(y, x) for _, y, x in found[:max_n]]


def depth(in_vals, out_vals):
    mo = median(out_vals)
    if not mo or math.isnan(mo):
        return float("nan")
    return 1.0 - median(in_vals) / mo


def box(img, y, x, r=2):
    """Flattened (2r+1)^2 box around (y, x); NaN outside the image."""
    out = []
    for yy in range(y - r, y + r + 1):
        for xx in range(x - r, x + r + 1):
            if 0 <= yy < len(img) and 0 <= xx < len(img[0]) and img[yy][xx] is not None:
                out.append(float(img[yy][xx]))
            else:
                out.append(float("nan"))
    return out


def label_event(d_t, corr_psf, ratio_peak, control_depths):
    """Pre-registered label of one flagged dip."""
    if d_t is None or math.isnan(d_t) or d_t < DEPTH_MIN:
        return "REPRO_FAIL"
    ctrl = [c for c in control_depths if c is not None and not math.isnan(c)]
    mc = median(ctrl) if len(ctrl) >= MIN_CONTROLS else float("nan")
    if not math.isnan(mc) and mc >= CTRL_SYS * d_t:
        return "SYSTEMATIC"
    if not math.isnan(corr_psf) and corr_psf < CORR_SYS:
        return "SYSTEMATIC"
    if (not math.isnan(corr_psf) and corr_psf >= CORR_ON and not math.isnan(ratio_peak)
            and ratio_peak >= RATIO_ON and (math.isnan(mc) or abs(mc) <= CTRL_ON * d_t)):
        return "ON_TARGET"
    return "UNCLEAR"


def control_gate(pos_label, neg_label):
    """The injected dip must come out ON_TARGET and the same window without injection REPRO_FAIL."""
    return "CONTROLS_OK" if (pos_label == "ON_TARGET" and neg_label == "REPRO_FAIL") else "CONTROL_FAIL"


def overall(gate, labels):
    """labels: {name: label} for the flagged events. CANDIDATE_ECLIPSE only if the gate passed and
    some event is ON_TARGET; ALL_SYSTEMATIC if every event is SYSTEMATIC; else MIXED."""
    if gate != "CONTROLS_OK":
        return "UNINTERPRETED"
    vals = list(labels.values())
    if any(v == "ON_TARGET" for v in vals):
        return "CANDIDATE_ECLIPSE"
    if vals and all(v == "SYSTEMATIC" for v in vals):
        return "ALL_SYSTEMATIC"
    return "MIXED"
