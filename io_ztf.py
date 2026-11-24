# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""ZTF public light curves from the IRSA API (fetch = network; parse/normalise = stdlib, tested)."""
import csv
import io
import time
import urllib.request

from ecl import robust_scatter

URL = ("https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?POS=CIRCLE%20{ra:.6f}%20{dec:.6f}"
       "%20{rad:.6f}&BANDNAME=g,r,i&FORMAT=csv&BAD_CATFLAGS_MASK=32768")


def fetch(ra, dec, out_path, radius_deg=0.0005, tries=4, log=print):
    url = URL.format(ra=ra, dec=dec, rad=radius_deg)
    last = None
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=300) as resp:
                data = resp.read()
            with open(out_path, "wb") as fh:
                fh.write(data)
            log("ZTF: %d bytes (try %d)" % (len(data), k + 1))
            return {"status": "ok", "url": url, "bytes": len(data)}
        except Exception as exc:  # network hiccup
            last = repr(exc)
            log("ZTF try %d failed: %s" % (k + 1, last))
            time.sleep(30 * (k + 1))
    return {"status": "failed", "url": url, "error": last}


def parse(text, max_magerr=0.05):
    """Rows -> per-oid normalised points [(bjd, r, sigma_eff, 'ztf', oid, filtercode)].
    Keeps catflags == 0 and magerr <= max_magerr. r = 10**(-0.4 (mag - median_oid)).
    sigma_eff = max(0.921 magerr, robust scatter of r for that oid)."""
    rd = csv.DictReader(io.StringIO(text))
    by = {}
    n_raw = 0
    for row in rd:
        n_raw += 1
        try:
            mag = float(row["mag"])
            err = float(row["magerr"])
            flags = int(float(row["catflags"]))
            hjd = float(row["hjd"])
        except (KeyError, ValueError, TypeError):
            continue
        if flags != 0 or err > max_magerr or mag != mag:
            continue
        t = hjd if hjd > 2400000 else hjd + 2400000.5
        by.setdefault((row.get("oid", "?"), row.get("filtercode", "?")), []).append((t, mag, err))
    out = []
    meta = {"n_raw": n_raw, "oids": {}}
    for (oid, fc), rows in by.items():
        if len(rows) < 10:
            meta["oids"][oid] = {"filter": fc, "n": len(rows), "used": False}
            continue
        mags = sorted(m for _, m, _ in rows)
        med = mags[len(mags) // 2] if len(mags) % 2 else 0.5 * (mags[len(mags) // 2 - 1] + mags[len(mags) // 2])
        rs = [10 ** (-0.4 * (m - med)) for _, m, _ in rows]
        sc = robust_scatter(rs)
        for (t, m, e), r in zip(rows, rs):
            out.append((t, r, max(0.921 * e, sc), "ztf", oid, fc))
        meta["oids"][oid] = {"filter": fc, "n": len(rows), "median_mag": med, "scatter": sc,
                             "used": True, "t_first": min(x[0] for x in rows),
                             "t_last": max(x[0] for x in rows)}
    return out, meta
