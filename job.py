# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""EPIC 220208795 second-eclipse test (pre-registered 2026-11-24, PREREGISTRATION.md).

python job.py --data-dir /data --out-dir /out
Writes /out/results.json (rewritten after every stage, key "stage"), CSVs and summary.md.
"""
import argparse
import csv
import json
import math
import os
import sys
import time
import traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import ecl  # noqa: E402
import io_ztf  # noqa: E402

RA, DEC = 17.731554, 0.314032  # 2MASS J01105556+0018507 = EPIC 220208795
START = time.time()


def log(msg):
    print("[%7.1fs] %s" % (time.time() - START, msg), flush=True)


def write_json(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=1, default=str)
    os.replace(tmp, path)


def write_csv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for r in rows:
            w.writerow(r)


def get_ztf(data_dir, res):
    path = os.path.join(data_dir, "ztf_lc.csv")
    info = io_ztf.fetch(RA, DEC, path, log=log)
    res["data"]["ztf"] = {"fetch": info}
    if info["status"] != "ok":
        return []
    with open(path) as fh:
        pts, meta = io_ztf.parse(fh.read())
    res["data"]["ztf"]["meta"] = meta
    res["data"]["ztf"]["n_used"] = len(pts)
    log("ZTF points used: %d" % len(pts))
    return pts


def get_asassn(data_dir, res):
    import asassn_io
    out = {}
    res["data"]["asassn"] = out
    try:
        out["client"] = asassn_io.ensure_client()
        path = os.path.join(data_dir, "asassn_lc.csv")
        info = None
        for k in range(3):
            try:
                info = asassn_io.fetch(RA, DEC, path, log=log)
                if info.get("status") == "ok":
                    break
            except Exception as exc:
                info = {"status": "error", "error": repr(exc)}
                log("ASAS-SN try %d: %r" % (k + 1, exc))
            time.sleep(60)
        out["fetch"] = info
        if not info or info.get("status") != "ok":
            return []
        by_band, lmeta = asassn_io.load_points(path)
        pts, nmeta = asassn_io.normalise(by_band)
        out["load"] = lmeta
        out["groups"] = nmeta
        out["n_used"] = len(pts)
        log("ASAS-SN points used: %d" % len(pts))
        return pts
    except Exception as exc:
        out["error"] = repr(exc)
        out["traceback"] = traceback.format_exc()[-2000:]
        log("ASAS-SN failed: %r" % exc)
        return []


def get_tess(data_dir, out_dir, res):
    import io_tess
    info = {"sectors": {}}
    res["data"]["tess"] = info
    veto_pts = []
    events = []
    try:
        secs = io_tess.sectors(RA, DEC, log=log)
    except Exception as exc:
        info["error"] = "sector list: %r" % exc
        return veto_pts, events
    info["sector_list"] = secs
    for s in secs:
        si = {}
        info["sectors"][str(s)] = si
        try:
            path = io_tess.fetch_cutout(RA, DEC, s, os.path.join(data_dir, "tesscut"), log=log)
            ph = io_tess.photometry(path)
        except Exception as exc:
            si["error"] = repr(exc)
            log("TESS sector %d failed: %r" % (s, exc))
            continue
        t = ph["bjd"]
        r3 = ecl.normalise_binned(t, ph["f3"])
        r1 = ecl.normalise_binned(t, ph["f1"])
        keep = [i for i in range(len(t)) if r3[i] == r3[i] and r1[i] == r1[i]]
        t = [t[i] for i in keep]
        r3 = [r3[i] for i in keep]
        r1 = [r1[i] for i in keep]
        cad = ph["cadence_days"]
        s3 = ecl.robust_scatter(r3)
        si.update({k: ph[k] for k in ("cadence_days", "n_good", "n_all", "aperture_centre",
                                      "n_bkg_pixels", "frac_aperture_in_peak")})
        si.update({"t_first": t[0] if t else None, "t_last": t[-1] if t else None,
                   "scatter_r3": s3, "n_used": len(t)})
        if len(t) < 50:
            si["note"] = "too few points"
            continue
        best, dets = ecl.box_search(t, r3, cad)
        si["best_box"] = best
        confirmed = []
        for d in dets:
            # same box on the single-pixel light curve
            w, c = d["width"], d["t_mid"]
            inbox = [r1[i] for i in range(len(t)) if abs(t[i] - c) <= w / 2]
            d1 = 1.0 - sum(inbox) / len(inbox) if inbox else float("nan")
            d["depth_peak_pixel"] = d1
            d["confirmed"] = bool(d1 == d1 and d1 >= 0.5 * d["depth"])
            # time of minimum: lowest 30-min bin inside +-0.2 d
            win = [(t[i], r3[i]) for i in range(len(t)) if abs(t[i] - c) <= 0.2]
            b = ecl.bin_series([x[0] for x in win], [x[1] for x in win], 1.0 / 48)
            d["t_min"] = min(b, key=lambda x: x[1])[0] if b else c
            d["sector"] = s
            confirmed.append(d)
            events.append({"t_mid": d["t_min"], "n": d["n_in"], "depth_max": d["depth"],
                           "sources": ["tess"], "strong": d["confirmed"], "sector": s,
                           "snr": d["snr"]})
        si["detections"] = confirmed
        # injection-recovery: template at half depth, 10 times evenly spaced inside the sector
        inj = []
        span = t[-1] - t[0]
        for k in range(10):
            T = t[0] + span * (k + 0.5) / 10.0
            ri = [r3[i] * (1.0 - 0.5 * ecl.template_depth(t[i] - T)) for i in range(len(t))]
            have = sum(1 for x in t if abs(x - (T + ecl.T_MIN_OFFSET)) <= 0.08)
            if have < 0.7 * 0.16 / cad:
                inj.append({"T": T, "covered": False})
                continue
            _, dd = ecl.box_search(t, ri, cad)
            ok = any(abs(x["t_mid"] - (T + ecl.T_MIN_OFFSET)) <= 0.15 for x in dd)
            inj.append({"T": T, "covered": True, "recovered": ok})
        cov = [x for x in inj if x["covered"]]
        si["injection"] = {"n_covered": len(cov), "n_recovered": sum(1 for x in cov if x["recovered"])}
        # 30-min bins as veto points
        b = ecl.bin_series(t, r3, 1.0 / 48)
        for tb, rb, nb in b:
            veto_pts.append((tb, rb, s3 / math.sqrt(nb), "tess"))
        write_csv(os.path.join(out_dir, "tess_s%03d_30min.csv" % s), ["bjd", "r_3x3", "n"],
                  [("%.5f" % tb, "%.5f" % rb, nb) for tb, rb, nb in b])
        log("TESS sector %d: %d pts, scatter %.4f, best box %s, detections %d, injection %s"
            % (s, len(t), s3, None if best is None else "%.3f@%.2f snr %.1f" % (best["depth"], best["t_mid"], best["snr"]),
               len(confirmed), si["injection"]))
    return veto_pts, events


def window_test(pts, P0, half, step, t_lo, t_hi):
    n = int(round(2 * half / step))
    grid = [P0 - half + k * step for k in range(n + 1)]
    flags = []
    for P in grid:
        ex, inf = ecl.trial(P, pts, t_lo, t_hi)
        flags.append((P, ex, inf))
    nex = sum(1 for _, e, _ in flags if e)
    return flags, nex / len(flags)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    a = ap.parse_args()
    os.makedirs(a.data_dir, exist_ok=True)
    os.makedirs(a.out_dir, exist_ok=True)
    rpath = os.path.join(a.out_dir, "results.json")
    res = {"version": ecl.VERSION, "target": "EPIC 220208795", "ra": RA, "dec": DEC,
           "stage": "start", "data": {}, "preregistration": "PREREGISTRATION.md (2026-11-24)"}
    write_json(rpath, res)

    ztf = get_ztf(a.data_dir, res)
    res["stage"] = "ztf"
    write_json(rpath, res)
    asn = get_asassn(a.data_dir, res)
    res["stage"] = "asassn"
    write_json(rpath, res)
    tess_pts, tess_events = get_tess(a.data_dir, a.out_dir, res)
    res["stage"] = "tess"
    write_json(rpath, res)

    rows = [(p[0], p[1], p[2], "ztf") for p in ztf] + [(p[0], p[1], p[2], "asassn") for p in asn] + tess_pts
    write_csv(os.path.join(a.out_dir, "ground_points.csv"), ["bjd", "r", "sigma", "survey", "group", "filter"],
              [("%.5f" % p[0], "%.5f" % p[1], "%.5f" % p[2], "ztf", p[4], p[5]) for p in ztf]
              + [("%.5f" % p[0], "%.5f" % p[1], "%.5f" % p[2], "asassn", p[5], p[4]) for p in asn])
    pts = ecl.Points(rows)
    if len(pts) == 0:
        res["stage"] = "done"
        res["outcome"] = "NO_DATA"
        write_json(rpath, res)
        log("OUTCOME NO_DATA")
        return
    t_lo, t_hi = pts.t[0], pts.t[-1]
    res["span"] = {"t_lo": t_lo, "t_hi": t_hi, "date_lo": ecl.bjd_to_date(t_lo), "date_hi": ecl.bjd_to_date(t_hi),
                   "n_points": len(pts), "n_veto_capable": len(pts.vt),
                   "by_source": {c: pts.src.count(c) for c in ("ztf", "asassn", "tess")}}
    pre = pts.subset(lambda t, c: t < ecl.AUTHORS_CUTOFF_BJD)

    # ground dip points and events
    dips = ecl.dip_points(pts)
    gev = ecl.group_events(dips)
    write_csv(os.path.join(a.out_dir, "dip_points_ground.csv"), ["bjd", "date", "r", "sigma", "survey"],
              [("%.5f" % d[0], ecl.bjd_to_date(d[0]), "%.4f" % d[1], "%.4f" % d[2], d[3]) for d in dips])
    events = gev + tess_events
    strong = [e for e in events if e["strong"]]
    res["events"] = {"ground_dip_points": len(dips), "ground_events": gev, "tess_events": tess_events,
                     "n_strong": len(strong)}
    res["stage"] = "events"
    write_json(rpath, res)
    log("dip points %d, ground events %d (strong %d), tess events %d (strong %d)"
        % (len(dips), len(gev), sum(1 for e in gev if e["strong"]), len(tess_events),
           sum(1 for e in tess_events if e["strong"])))

    # pre-registered candidate windows
    cand = {}
    rows_c = []
    for name, P0 in ecl.CANDIDATES.items():
        flags, frac = window_test(pts, P0, ecl.LOCAL_HALF, ecl.LOCAL_STEP, t_lo, t_hi)
        flags_pre, frac_pre = window_test(pre, P0, ecl.LOCAL_HALF, ecl.LOCAL_STEP, t_lo, t_hi) if len(pre) else ([], float("nan"))
        _, frac_wide = window_test(pts, P0, ecl.WIDE_HALF, 0.001, t_lo, t_hi)
        allowed = [P for P, e, _ in flags if not e]
        support = []
        for e in strong:
            for P in allowed:
                n = round((e["t_mid"] - ecl.T_MIN_OFFSET - ecl.T_REF) / P)
                if n != 0 and abs(e["t_mid"] - (ecl.T_REF + ecl.T_MIN_OFFSET + n * P)) <= ecl.SUPPORT_DAYS:
                    support.append({"event_t": e["t_mid"], "date": ecl.bjd_to_date(e["t_mid"]), "n": n,
                                    "P": P, "sources": e["sources"]})
                    break
        verdict = ecl.candidate_verdict(frac, bool(support))
        # which predicted eclipse of P0 itself is vetoed, and by what
        ex0, inf0 = ecl.trial(P0, pts, t_lo, t_hi)
        nxt = []
        for P in ([P0] + ([min(allowed), max(allowed)] if allowed else [])):
            for n, T in ecl.eclipse_refs(P, 2461369.5, 2461771.5):  # 2026-12-01 .. 2027-12-31
                nxt.append({"P": P, "n": n, "min_bjd": T + ecl.T_MIN_OFFSET,
                            "date": ecl.bjd_to_date(T + ecl.T_MIN_OFFSET)})
        cand[name] = {"P0": P0, "window_half": ecl.LOCAL_HALF, "excluded_frac": frac,
                      "excluded_frac_data_to_2021_02_12": frac_pre, "excluded_frac_wide_half_0.5": frac_wide,
                      "n_allowed": len(allowed), "allowed_min": min(allowed) if allowed else None,
                      "allowed_max": max(allowed) if allowed else None,
                      "P0_excluded": ex0, "P0_veto": inf0, "support": support, "verdict": verdict,
                      "next_eclipses": nxt}
        for (P, e, inf), (_, epre, _) in zip(flags, flags_pre or [(None, None, None)] * len(flags)):
            rows_c.append((name, "%.4f" % P, int(e), "" if epre is None else int(epre),
                           "" if not inf else inf["n"], "" if not inf else "%.4f" % inf["T"],
                           "" if not inf else ecl.bjd_to_date(inf["T"]), "" if not inf else "+".join(inf["by"])))
        log("%s P0=%.3f: excluded %.3f (to 2021-02: %s) verdict %s" % (name, P0, frac, frac_pre, verdict))
    res["candidates"] = cand
    write_csv(os.path.join(a.out_dir, "candidate_windows.csv"),
              ["candidate", "P", "excluded", "excluded_data_to_2021_02_12", "veto_n", "veto_T_bjd", "veto_date", "veto_by"],
              rows_c)
    res["stage"] = "candidates"
    write_json(rpath, res)

    # implied periods from strong events
    implied = []
    for e in strong:
        for n, P in ecl.periods_from_event(e["t_mid"]):
            ex, inf = ecl.trial(P, pts, t_lo, t_hi)
            implied.append({"event_t": e["t_mid"], "date": ecl.bjd_to_date(e["t_mid"]), "n": n, "P": P,
                            "excluded": ex, "veto": inf, "sources": e["sources"]})
    res["implied_periods"] = implied
    write_csv(os.path.join(a.out_dir, "implied_periods.csv"), ["event_date", "n", "P", "excluded"],
              [(x["date"], x["n"], "%.4f" % x["P"], int(x["excluded"])) for x in implied])

    # blind scan 60-1500 d
    grid = ecl.period_grid(ecl.P_MIN, ecl.P_MAX, t_lo, t_hi)
    res["scan"] = {"n_trials": len(grid), "tol_days": ecl.GRID_TOL}
    write_json(rpath, res)
    log("scan: %d trials" % len(grid))
    edges = [60.0 + 20.0 * k for k in range(73)]  # 60..1500
    for label, p, mv in (("primary", pts, ecl.MIN_VETO), ("one_point_veto", pts, 1),
                         ("data_to_2021_02_12", pre, ecl.MIN_VETO)):
        if len(p) == 0:
            continue
        n_ex, runs = ecl.scan(p, grid, t_lo, t_hi, mv)
        fr = ecl.binned_fraction(grid, runs, edges)
        res["scan"][label] = {"n_excluded": n_ex, "frac_excluded": n_ex / len(grid), "n_allowed_runs": len(runs)}
        write_csv(os.path.join(a.out_dir, "allowed_periods_%s.csv" % label), ["P_first", "P_last", "n_trials"],
                  [("%.5f" % x[0], "%.5f" % x[1], x[2]) for x in runs])
        write_csv(os.path.join(a.out_dir, "chart_allowed_fraction_%s.csv" % label),
                  ["P_lo", "P_hi", "n_trials", "n_allowed", "frac_allowed"],
                  [("%.0f" % x[0], "%.0f" % x[1], x[2], x[3], "%.4f" % x[4]) for x in fr])
        log("scan %s: excluded %d / %d (%.3f), %d allowed runs" % (label, n_ex, len(grid), n_ex / len(grid), len(runs)))
        write_json(rpath, res)

    # outcome
    out = []
    if any(e["sources"] == ["tess"] for e in strong):
        out.append("NEW_ECLIPSE_TESS")
    elif strong:
        out.append("NEW_ECLIPSE_GROUND")
    out += ["%s_%s" % (k, v["verdict"]) for k, v in cand.items()]
    res["outcome"] = out
    res["stage"] = "done"
    res["runtime_s"] = time.time() - START
    write_json(rpath, res)
    with open(os.path.join(a.out_dir, "summary.md"), "w") as fh:
        fh.write("# EPIC 220208795 second-eclipse test v%s\n\n" % ecl.VERSION)
        fh.write("Outcome: %s\n\n" % ", ".join(out))
        for k, v in cand.items():
            fh.write("- %s (P0 = %.3f d): excluded fraction of P0 +- %.2f d = %.4f (data to 2021-02-12: %s); verdict %s\n"
                     % (k, v["P0"], ecl.LOCAL_HALF, v["excluded_frac"], v["excluded_frac_data_to_2021_02_12"], v["verdict"]))
        fh.write("\nStrong events: %d. Scan 60-1500 d: %s\n" % (len(strong), json.dumps(res["scan"], default=str)))
    log("OUTCOME " + " ".join(out))


if __name__ == "__main__":
    main()
