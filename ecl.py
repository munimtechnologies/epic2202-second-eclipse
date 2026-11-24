# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""Core of the EPIC 220208795 second-eclipse test (stdlib only; imported by the tests).

The K2 eclipse template below is restated (rounded to 4 decimals) from the file
"data 2202/eclipse_lightcurve2202.txt" of github.com/lizvdkamp/EPIC2202_disk (BSD-2-Clause,
van der Kamp et al. 2022, A&A 658, A38; arXiv:2110.15086). Time is BJD - 2457409.0, flux is
normalised K2 flux. Outside the listed range the flux is 1.

Everything here is fixed by PREREGISTRATION.md (2026-11-24) before any post-2021 data were
analysed by this code.
"""
import bisect
import math
import statistics

VERSION = "0.1.0"

T_REF = 2457409.0  # BJD; zero of the template time axis (the K2 eclipse)

TEMPLATE = [
    (-0.2089, 1.0000),
    (-0.1885, 1.0006),
    (-0.1680, 1.0005),
    (-0.1476, 0.9928),
    (-0.1272, 0.9826),
    (-0.1067, 0.9517),
    (-0.0863, 0.8920),
    (-0.0659, 0.8139),
    (-0.0454, 0.7447),
    (-0.0250, 0.7166),
    (-0.0046, 0.7421),
    (0.0159, 0.7719),
    (0.0363, 0.8087),
    (0.0567, 0.8600),
    (0.0771, 0.8970),
    (0.0976, 0.9419),
    (0.1180, 0.9679),
    (0.1384, 0.9868),
    (0.1589, 0.9976),
    (0.1793, 1.0000),
    (0.1997, 0.9995),
    (0.2202, 1.0000),
]
TEMPLATE_T = [p[0] for p in TEMPLATE]
T_MIN_OFFSET = -0.0250  # template time of minimum flux

# Veto zone: template depth >= CORE_DEPTH. Its edges (by linear interpolation of the
# template) are CORE_LO..CORE_HI; points outside this offset range can never veto.
CORE_DEPTH = 0.15
CORE_LO = -0.0760
CORE_HI = 0.0530

# Pre-registered thresholds (PREREGISTRATION.md section 3)
VETO_MAX_DIP = {"ztf": 0.05, "asassn": 0.05, "tess": 0.03}
VETO_NSIGMA = 4.0
MIN_VETO = 2                 # vetoing points needed inside one predicted eclipse
DIP_MIN_DEPTH = 0.10         # ground dip point
DIP_NSIGMA = 5.0
DIP_GROUP_DAYS = 0.30
SUPPORT_DAYS = 0.30          # strong dip event within this of a predicted eclipse minimum
TESS_BOX_WIDTHS = (0.08, 0.12, 0.16)
TESS_MIN_DEPTH = 0.04
TESS_MIN_SNR = 10.0
P_MIN, P_MAX = 60.0, 1500.0
GRID_TOL = 0.02              # max drift (days) of the farthest predicted eclipse between grid steps
CANDIDATES = {"P290": 290.230, "P235": 235.587}
LOCAL_HALF = 0.10            # pre-registered window: P0 +- 0.10 d
WIDE_HALF = 0.50             # reported for context only
LOCAL_STEP = 0.0002
AUTHORS_CUTOFF_BJD = 2459257.5  # 2021-02-12, last ground data used by van der Kamp et al.


def template_flux(u):
    """Template flux at offset u (days from T_REF-aligned eclipse reference)."""
    if u <= TEMPLATE[0][0] or u >= TEMPLATE[-1][0]:
        return 1.0
    i = bisect.bisect_right(TEMPLATE_T, u)
    t0, f0 = TEMPLATE[i - 1]
    t1, f1 = TEMPLATE[i]
    return f0 + (f1 - f0) * (u - t0) / (t1 - t0)


def template_depth(u):
    return 1.0 - template_flux(u)


def robust_scatter(values):
    """1.4826 x median absolute deviation."""
    v = [x for x in values if x == x]
    if len(v) < 2:
        return float("nan")
    m = statistics.median(v)
    return 1.4826 * statistics.median([abs(x - m) for x in v])


class Points:
    """Normalised photometry: parallel lists sorted by time. r = flux / out-of-eclipse level."""

    def __init__(self, rows):
        rows = sorted(rows, key=lambda x: x[0])
        self.t = [x[0] for x in rows]
        self.r = [x[1] for x in rows]
        self.s = [x[2] for x in rows]
        self.src = [x[3] for x in rows]
        # veto-capable subset (flux close to normal)
        vc = [x for x in rows if (1.0 - x[1]) < VETO_MAX_DIP.get(x[3], 0.05)]
        self.vt = [x[0] for x in vc]
        self.vr = [x[1] for x in vc]
        self.vs = [x[2] for x in vc]
        self.vsrc = [x[3] for x in vc]

    def __len__(self):
        return len(self.t)

    def subset(self, keep):
        rows = [(t, r, s, c) for t, r, s, c in zip(self.t, self.r, self.s, self.src) if keep(t, c)]
        return Points(rows)


def vetoes(pts, T, limit=None):
    """Points that veto an eclipse whose template zero is at T. Returns list of (t, source)."""
    out = []
    i = bisect.bisect_left(pts.vt, T + CORE_LO)
    n = len(pts.vt)
    while i < n and pts.vt[i] <= T + CORE_HI:
        u = pts.vt[i] - T
        f = template_flux(u)
        if 1.0 - f >= CORE_DEPTH:
            r, s, c = pts.vr[i], pts.vs[i], pts.vsrc[i]
            if (1.0 - r) < VETO_MAX_DIP.get(c, 0.05) and (r - f) > VETO_NSIGMA * s:
                out.append((pts.vt[i], c))
                if limit is not None and len(out) >= limit:
                    return out
        i += 1
    return out


def eclipse_refs(P, t_lo, t_hi):
    """Template zeros T_REF + n P inside [t_lo, t_hi] (n != 0)."""
    n0 = int(math.ceil((t_lo - T_REF) / P))
    n1 = int(math.floor((t_hi - T_REF) / P))
    return [(n, T_REF + n * P) for n in range(n0, n1 + 1) if n != 0]


def trial(P, pts, t_lo, t_hi, min_veto=MIN_VETO):
    """(excluded?, info). Excluded if one predicted eclipse has >= min_veto vetoing points."""
    for n, T in eclipse_refs(P, t_lo, t_hi):
        v = vetoes(pts, T, limit=min_veto)
        if len(v) >= min_veto:
            return True, {"n": n, "T": T, "by": sorted(set(c for _, c in v))}
    return False, None


def period_grid(p_min, p_max, span_lo, span_hi, tol=GRID_TOL):
    D = max(abs(span_hi - T_REF), abs(T_REF - span_lo), 1.0)
    fac = 1.0 + tol / D
    P = p_min
    out = []
    while P <= p_max:
        out.append(P)
        P *= fac
    return out


def scan(pts, grid, t_lo, t_hi, min_veto=MIN_VETO):
    """Returns (n_excluded, allowed_runs) with runs [P_first, P_last, n_trials]."""
    runs = []
    cur = None
    n_ex = 0
    for P in grid:
        ex, _ = trial(P, pts, t_lo, t_hi, min_veto)
        if ex:
            n_ex += 1
            if cur is not None:
                runs.append(cur)
                cur = None
        else:
            if cur is None:
                cur = [P, P, 1]
            else:
                cur[1] = P
                cur[2] += 1
    if cur is not None:
        runs.append(cur)
    return n_ex, runs


def binned_fraction(grid, runs, edges):
    """Fraction of grid trials allowed in each [edges[k], edges[k+1]) bin."""
    allowed = bytearray(len(grid))
    # mark allowed trials by index ranges (grid is sorted)
    for a, b, _ in runs:
        i = bisect.bisect_left(grid, a)
        j = bisect.bisect_right(grid, b)
        for x in range(i, j):
            allowed[x] = 1
    out = []
    for k in range(len(edges) - 1):
        i = bisect.bisect_left(grid, edges[k])
        j = bisect.bisect_left(grid, edges[k + 1])
        tot = j - i
        ok = sum(allowed[i:j])
        out.append((edges[k], edges[k + 1], tot, ok, (ok / tot) if tot else float("nan")))
    return out


def dip_points(pts, min_depth=DIP_MIN_DEPTH, nsig=DIP_NSIGMA, sources=("ztf", "asassn")):
    return [(t, r, s, c) for t, r, s, c in zip(pts.t, pts.r, pts.s, pts.src)
            if c in sources and (1.0 - r) >= min_depth and (1.0 - r) >= nsig * s]


def group_events(dips, gap=DIP_GROUP_DAYS):
    """Group dip points closer than gap; each event: dict(t_mid, n, depth_max, sources)."""
    ev = []
    for d in sorted(dips):
        if ev and d[0] - ev[-1]["t_last"] <= gap:
            e = ev[-1]
            e["pts"].append(d)
            e["t_last"] = d[0]
        else:
            ev.append({"pts": [d], "t_last": d[0]})
    out = []
    for e in ev:
        p = e["pts"]
        out.append({"t_mid": sum(x[0] for x in p) / len(p), "n": len(p),
                    "depth_max": max(1.0 - x[1] for x in p),
                    "sources": sorted(set(x[3] for x in p)),
                    "strong": len(p) >= 2})
    return out


def normalise_binned(t, f, bin_days=1.5):
    """Divide f by the median of its bin (bins of bin_days from the first time), with linear
    interpolation of bin medians between bin centres. Returns list of ratios."""
    if not t:
        return []
    t0 = t[0]
    bins = {}
    for ti, fi in zip(t, f):
        bins.setdefault(int((ti - t0) // bin_days), []).append((ti, fi))
    centres = []
    for k in sorted(bins):
        vals = bins[k]
        if len(vals) < 5:
            continue
        centres.append((statistics.median(v[0] for v in vals), statistics.median(v[1] for v in vals)))
    if not centres:
        return [float("nan")] * len(t)
    ct = [c[0] for c in centres]
    out = []
    for ti, fi in zip(t, f):
        j = bisect.bisect_left(ct, ti)
        if j <= 0:
            m = centres[0][1]
        elif j >= len(ct):
            m = centres[-1][1]
        else:
            (ta, ma), (tb, mb) = centres[j - 1], centres[j]
            m = ma + (mb - ma) * (ti - ta) / (tb - ta) if tb > ta else ma
        out.append(fi / m if m else float("nan"))
    return out


def bin_series(t, r, width):
    """Average into bins of `width` days. Returns [(t_mean, r_mean, n)]."""
    out = []
    if not t:
        return out
    t0 = t[0]
    cur_k, acc = None, []
    for ti, ri in zip(t, r):
        k = int((ti - t0) // width)
        if cur_k is not None and k != cur_k and acc:
            out.append((sum(a for a, _ in acc) / len(acc), sum(b for _, b in acc) / len(acc), len(acc)))
            acc = []
        cur_k = k
        acc.append((ti, ri))
    if acc:
        out.append((sum(a for a, _ in acc) / len(acc), sum(b for _, b in acc) / len(acc), len(acc)))
    return out


def box_search(t, r, cadence, widths=TESS_BOX_WIDTHS, step=0.01):
    """Box-dip search on a normalised series. Returns (best, detections). Each entry:
    dict(t_mid, width, depth, snr, n_in). Detection: depth >= TESS_MIN_DEPTH, snr >= TESS_MIN_SNR,
    >= 70% of expected in-box cadences present and >= 10 points on each side within 0.4 d with
    median ratio within 0.01 of 1."""
    n = len(t)
    if n < 20:
        return None, []
    s = robust_scatter(r)
    pre = [0.0]
    for x in r:
        pre.append(pre[-1] + x)
    best = None
    dets = []
    c = t[0]
    while c <= t[-1]:
        for w in widths:
            i = bisect.bisect_left(t, c - w / 2)
            j = bisect.bisect_right(t, c + w / 2)
            m = j - i
            if m < max(3, 0.7 * w / cadence):
                continue
            mean = (pre[j] - pre[i]) / m
            depth = 1.0 - mean
            snr = depth / (s / math.sqrt(m)) if s > 0 else 0.0
            cand = {"t_mid": c, "width": w, "depth": depth, "snr": snr, "n_in": m}
            if best is None or snr > best["snr"]:
                best = cand
            if depth >= TESS_MIN_DEPTH and snr >= TESS_MIN_SNR:
                li = bisect.bisect_left(t, c - 0.4)
                ri = bisect.bisect_right(t, c + 0.4)
                left = [r[k] for k in range(li, i) if t[k] < c - w / 2 - 0.05]
                right = [r[k] for k in range(j, ri) if t[k] > c + w / 2 + 0.05]
                if (len(left) >= 10 and len(right) >= 10 and abs(statistics.median(left) - 1) <= 0.01
                        and abs(statistics.median(right) - 1) <= 0.01):
                    dets.append(cand)
        c += step
    # merge detections within 0.3 d, keep highest snr
    merged = []
    for d in sorted(dets, key=lambda x: x["t_mid"]):
        if merged and d["t_mid"] - merged[-1]["t_mid"] <= 0.3:
            if d["snr"] > merged[-1]["snr"]:
                merged[-1] = d
        else:
            merged.append(d)
    return best, merged


def candidate_verdict(local_excluded_frac, supported):
    if supported:
        return "SUPPORTED"
    if local_excluded_frac >= 0.99:
        return "EXCLUDED"
    return "OPEN"


def periods_from_event(t_min, P_lo=P_MIN, P_hi=P_MAX):
    """Periods implied by an eclipse minimum at t_min together with the K2 eclipse."""
    T = t_min - T_MIN_OFFSET
    dt = T - T_REF
    out = []
    n = 1
    while True:
        P = abs(dt) / n
        if P < P_lo:
            break
        if P <= P_hi:
            out.append((n if dt > 0 else -n, P))
        n += 1
    return out


def bjd_to_date(bjd):
    """Calendar date (UTC, ignoring the BJD-UTC offset of ~1 min) as YYYY-MM-DD HH:MM."""
    import datetime
    d = datetime.datetime(2000, 1, 1, 12) + datetime.timedelta(days=bjd - 2451545.0)
    return d.strftime("%Y-%m-%d %H:%M")
