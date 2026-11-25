# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""Pixel-level vetting job (AMENDMENT-2026-11-25-pixvet.md).

python job_pixvet.py --data-dir /data --out-dir /out
Needs numpy and astropy. All decisions are made by pixvet.py (stdlib, unit-tested).
"""
import argparse
import glob
import io
import json
import math
import os
import sys
import time
import urllib.request
import zipfile

import pixvet

RA, DEC = 17.731554, 0.314032
CUT_URL = ("https://mast.stsci.edu/tesscut/api/v0.1/astrocut?ra={ra:.6f}&dec={dec:.6f}"
           "&y={n}&x={n}&units=px&sector={s}")
BAD_BITS = 1 | 2 | 4 | 8 | 32 | 128   # same as io_tess.py
NCUT = 21
HALF_WIN = 0.5

# Flagged events from results/v0.1.0/results.json (t_mid, width of the box detection).
EVENTS = [
    {"name": "S30_2459127.00", "sector": 30, "t_mid": 2459126.9979289854, "width": 0.12},
    {"name": "S30_2459127.53", "sector": 30, "t_mid": 2459127.5279289735, "width": 0.12},
    {"name": "S42_2459456.12", "sector": 42, "t_mid": 2459456.121552408, "width": 0.08},
]
# Controls: Sector 97 at the P235 n=15 prediction (flat in the light curve), with and without an
# injected 10% dip of the target's own image.
CONTROLS = [
    {"name": "S97_negative", "sector": 97, "t_mid": 2460942.805, "width": 0.12, "inject": 0.0},
    {"name": "S97_positive", "sector": 97, "t_mid": 2460942.805, "width": 0.12, "inject": 0.10},
]


def log(msg):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), msg), flush=True)


def get(url, tries=4, timeout=900):
    last = None
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                return r.read()
        except Exception as exc:  # noqa: BLE001
            last = repr(exc)
            log("GET try %d failed: %s" % (k + 1, last))
            time.sleep(30 * (k + 1))
    raise RuntimeError(last)


def fetch(sector, data_dir):
    d = os.path.join(data_dir, "tesscut21")
    os.makedirs(d, exist_ok=True)
    have = glob.glob(os.path.join(d, "*s%04d*.fits" % sector))
    if have:
        return have[0]
    blob = get(CUT_URL.format(ra=RA, dec=DEC, n=NCUT, s=sector))
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        names = [x for x in z.namelist() if x.endswith(".fits")]
        z.extractall(d, members=names)
    log("sector %d cutout: %s" % (sector, names))
    return os.path.join(d, names[0])


def tolist2(a):
    return [[float(v) for v in row] for row in a]


def analyse(cube_t, cube_f, med_all, ev):
    import numpy as np
    sel = np.abs(cube_t - ev["t_mid"]) <= HALF_WIN
    t = cube_t[sel]
    f = cube_f[sel].astype(float).copy()
    ins, outs = pixvet.window_idx(t.tolist(), ev["t_mid"], ev["width"])
    res = {"name": ev["name"], "sector": ev["sector"], "t_mid": ev["t_mid"], "width": ev["width"],
           "n_in": len(ins), "n_out": len(outs)}
    if len(ins) < 3 or len(outs) < 10:
        res["label"] = "REPRO_FAIL"
        res["note"] = "too few cadences"
        return res, None
    med_l = tolist2(med_all)
    ty, tx = pixvet.find_target(med_l)
    ny, nx = med_all.shape
    tmask = pixvet.bkg_mask(med_l, ty, tx, ny // 2, nx // 2)
    ys = np.array([p[0] for p in tmask])
    xs = np.array([p[1] for p in tmask])

    def bkg_series(arr, mys, mxs):
        return np.nanmedian(arr[:, mys, mxs], axis=1)

    inject = float(ev.get("inject", 0.0))
    if inject > 0:
        b0 = bkg_series(f[outs], ys, xs)
        m_out = np.nanmedian(f[outs], axis=0) - float(np.nanmedian(b0))
        psf = np.zeros_like(m_out)
        psf[ty - 2:ty + 3, tx - 2:tx + 3] = np.clip(m_out[ty - 2:ty + 3, tx - 2:tx + 3], 0, None)
        for k in ins:
            f[k] = f[k] - inject * psf
        res["injected_depth_of_5x5"] = inject

    b = bkg_series(f, ys, xs)
    f3 = np.nansum(f[:, ty - 1:ty + 2, tx - 1:tx + 2], axis=(1, 2)) - 9.0 * b
    f1 = f[:, ty, tx] - b
    d_t = pixvet.depth([float(f3[k]) for k in ins], [float(f3[k]) for k in outs])
    d_1 = pixvet.depth([float(f1[k]) for k in ins], [float(f1[k]) for k in outs])
    ratio = d_1 / d_t if d_t and not math.isnan(d_t) and d_t != 0 else float("nan")

    raw_in = np.nanmedian(f[ins], axis=0)
    raw_out = np.nanmedian(f[outs], axis=0)
    b_in = float(np.nanmedian(b[ins]))
    b_out = float(np.nanmedian(b[outs]))
    dimg = raw_in - raw_out - (b_in - b_out)
    mimg = raw_out - b_out
    corr = -pixvet.pearson(pixvet.box(tolist2(dimg), ty, tx), pixvet.box(tolist2(mimg), ty, tx))
    f3_out = float(np.nanmedian(f3[outs]))
    bkg_equiv = 9.0 * (b_in - b_out) / f3_out if f3_out else float("nan")

    ctrls = []
    for (cy, cx) in pixvet.local_maxima(tolist2(mimg), ty, tx):
        cm = pixvet.bkg_mask(med_l, cy, cx, cy, cx)
        if len(cm) < 5:
            continue
        cys = np.array([p[0] for p in cm])
        cxs = np.array([p[1] for p in cm])
        cb = bkg_series(f, cys, cxs)
        c3 = np.nansum(f[:, cy - 1:cy + 2, cx - 1:cx + 2], axis=(1, 2)) - 9.0 * cb
        cd = pixvet.depth([float(c3[k]) for k in ins], [float(c3[k]) for k in outs])
        ctrls.append({"y": int(cy), "x": int(cx), "flux3_ratio": float(np.nanmedian(c3[outs]) / f3_out),
                      "depth": cd})
    label = pixvet.label_event(d_t, corr, ratio, [c["depth"] for c in ctrls])
    res.update({"target_pixel": [int(ty), int(tx)], "n_bkg_pixels": len(tmask),
                "depth_aperture": d_t, "depth_peak": d_1, "ratio_peak": ratio,
                "corr_psf": corr, "bkg_change_equiv": bkg_equiv,
                "controls": ctrls, "control_median_depth": pixvet.median([c["depth"] for c in ctrls]),
                "label": label})
    lc = [(float(t[k]), float(f3[k] / f3_out), float(b[k])) for k in range(len(t))]
    return res, {"diff": tolist2(dimg), "psf": tolist2(mimg), "lc": lc}


def main():
    import numpy as np
    from astropy.io import fits
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    out = {"version": "pixvet-0.1.0", "amendment": "AMENDMENT-2026-11-25-pixvet.md", "cases": {}}
    cubes = {}
    for ev in EVENTS + CONTROLS:
        s = ev["sector"]
        if s not in cubes:
            path = fetch(s, a.data_dir)
            with fits.open(path, memmap=False) as h:
                d = h[1].data
                t = np.array(d["TIME"], dtype=float) + 2457000.0
                fl = np.array(d["FLUX"], dtype=np.float32)
                q = np.array(d["QUALITY"], dtype=int)
            ny, nx = fl.shape[1], fl.shape[2]
            good = np.isfinite(t) & ((q & BAD_BITS) == 0) & np.isfinite(fl[:, ny // 2, nx // 2])
            t, fl = t[good], fl[good]
            med_all = np.nanmedian(fl, axis=0)
            near = np.zeros(len(t), dtype=bool)
            for e in EVENTS + CONTROLS:
                if e["sector"] == s:
                    near |= np.abs(t - e["t_mid"]) <= HALF_WIN
            cubes[s] = (t[near], fl[near], med_all)
            del fl
            log("sector %d: %d good cadences kept near events" % (s, int(near.sum())))
        t, fl, med_all = cubes[s]
        res, extra = analyse(t, fl, med_all, ev)
        out["cases"][ev["name"]] = res
        log("%s: %s" % (ev["name"], json.dumps({k: res.get(k) for k in (
            "label", "depth_aperture", "depth_peak", "ratio_peak", "corr_psf", "bkg_change_equiv",
            "control_median_depth")})))
        if extra:
            with open(os.path.join(a.out_dir, "pix_%s.json" % ev["name"]), "w") as fh:
                json.dump(extra, fh)
    gate = pixvet.control_gate(out["cases"]["S97_positive"]["label"],
                               out["cases"]["S97_negative"]["label"])
    labels = {e["name"]: out["cases"][e["name"]]["label"] for e in EVENTS}
    out["gate"] = gate
    out["overall"] = pixvet.overall(gate, labels)
    with open(os.path.join(a.out_dir, "results_pixvet.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=lambda o: None)
    log("OUTCOME gate=%s overall=%s labels=%s" % (gate, out["overall"], labels))
    return 0


if __name__ == "__main__":
    sys.exit(main())
