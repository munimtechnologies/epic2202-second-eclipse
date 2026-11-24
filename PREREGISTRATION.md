# Pre-registration: does the disk around EPIC 220208795's companion come back on either published period?

Frozen 2026-11-24 (UTC), before the job below has run. MV2, Munim, Inc. The work, code and text
were made with AI agents working in MV2's name; no human expert has checked them yet.
Licence: Apache-2.0 (code), CC BY 4.0 (text and results).

## 1. The near miss

van der Kamp, van Dam, Kenworthy, Mamajek & Pojmanski (2022, A&A 658, A38;
[arXiv:2110.15086](https://arxiv.org/abs/2110.15086)) found a single 7.2-hour, ~28%-deep,
asymmetric eclipse of the K dwarf EPIC 220208795 (2MASS J01105556+0018507, G = 14.26) in K2
Campaign 8 (BJD 2457409, 2016-01-16). They modelled it as a tilted, inclined opaque disk
(radius 1.14 R_sun) around an unseen companion. Folding ASAS, ASAS-SN and ATLAS photometry up to
2021-02 gave two "most likely" periods, **290.230 d** and **235.587 d** (eccentricity ~0.7), and
their Table 8 listed the next predicted eclipses (2021-2024). They wrote that the period was not
settled, because the eclipse is short and ground coverage is sparse. Mustill, Davies & Kenworthy
(2024, MNRAS; [arXiv:2404.12239](https://arxiv.org/abs/2404.12239)) took P = 290 d, e = 0.72 as
the system's orbit in their simulations of circumplanetary debris disks.

No later paper reports a second eclipse or a test of either period (searches of arXiv, Google,
OpenAlex citations on 2026-11-24; see notes/scoop-check-2026-11-24.md).

Since then, public data have covered many predicted dates:
- TESS full-frame images: Sectors 42 and 43 (Aug-Oct 2021), 70 (Sep-Oct 2023) and **97
  (2025-09-15 to 2025-11-09, four orbits)**, besides Sectors 3 and 30 already used by the authors.
  The 235.587-d period predicts an eclipse minimum near **2025-09-24 07 UT (n = 15)**, inside Sector 97.
- ZTF public light curves (g, r; 2018-2025), and ASAS-SN Sky Patrol v2 (V 2012-2018, g 2017-2026).

The wow if it works: the second-ever eclipse of a disk orbiting an unseen companion, fixing its
orbit, or a clean measurement that both published periods are wrong.

## 2. Data (all public, fetched by job.py)

| Source | Access | Use |
|---|---|---|
| ZTF | IRSA light-curve API, 1.8" cone, catflags = 0, magerr <= 0.05 | veto and dip points |
| ASAS-SN | Sky Patrol v2 via pyasassn, nearest source within 10" with rows, quality G, rel. error <= 0.05 | veto and dip points |
| TESS | MAST TESScut 11x11 FFI cutouts, every sector the sector service lists | box search, veto points (30-min bins) |

ATLAS forced photometry needs an account and is not used. Gaia DR4 epoch photometry (2014-2020)
will be added on/after 2026-12-02 as a separately dated amendment.

Disclosure: before freezing this file we read an automatic text summary of the ZTF API output
(about 1,000 points; it reported no point more than 0.12 mag below its light curve's median). The
thresholds below were set from the K2 template and survey noise, not from that summary.

## 3. Method (code: ecl.py, job.py)

**Template.** The K2 eclipse as published by the authors (their GitHub file
eclipse_lightcurve2202.txt, BSD-2-Clause), time t - 2457409.0, linearly interpolated. A trial
period P predicts eclipses with template zero at T_n = 2457409.0 + n P (n != 0); the minimum is at
T_n - 0.025 d. Assumption A1: later eclipses have the same timing and shape, and at least half the
depth (disk precession could change the depth; section 5).

**Normalisation.** ZTF: per object ID and filter, r = 10^(-0.4 (m - median)). ASAS-SN: per band
and camera, r = flux / median. TESS: 3x3-pixel aperture minus a per-cadence median background,
divided by 1.5-day running bin medians; also the single peak pixel as a check. sigma_eff = max(own
error, 1.4826 x MAD of that group); TESS 30-min bins use MAD/sqrt(n).

**Veto.** A point vetoes a predicted eclipse if it lies where the template is at least 15% deep
(offsets -0.075 to +0.053 d, about 3.1 h), the observed dip is below 5% (ZTF, ASAS-SN) or 3%
(TESS, to allow for blending in 21" pixels), and r - template > 4 sigma_eff. A trial period is
**excluded** if any single predicted eclipse in the data span has **at least 2 vetoing points**.
(Sensitivity analysis: 1 point.)

**Dips.** Ground dip point: depth >= 10% and >= 5 sigma_eff. Points within 0.3 d form an event;
an event is **strong** with >= 2 points. TESS box search (widths 0.08/0.12/0.16 d, 0.01-d steps):
detection if depth >= 4%, S/N >= 10, >= 70% of in-box cadences present and >= 10 points within
0.4 d on each side with median within 1% of 1; a TESS event is strong if the peak-pixel light
curve shows at least half the aperture depth. Injection test: the template at half depth at 10
evenly spaced times per sector.

**Candidate windows (primary).** For each published period P0, trial periods P0 +- 0.10 d in
0.0002-d steps (1,001 trials), all data. Also reported: the same with only data before
2021-02-12 (the authors' cutoff; our data, not theirs), and +- 0.50 d.

**Blind scan.** P = 60-1500 d (60 d = the K2 lower limit), geometric grid with at most 0.02 d drift
of the farthest predicted eclipse between steps. Output: excluded fraction, allowed intervals,
allowed fraction per 20-d bin.

**Implied periods.** For every strong event, P = (T - 2457409.0)/n with T = t_min + 0.025, each
checked against the vetoes.

## 4. Pre-registered outcomes

Per candidate (P290, P235):
- **SUPPORTED**: a strong event lies within 0.30 d of a predicted minimum for an allowed trial
  period in the window.
- **EXCLUDED**: not supported and at least 99% of the window's trial periods are excluded.
- **OPEN**: otherwise.

Overall: **NEW_ECLIPSE_TESS** if any strong TESS event, else **NEW_ECLIPSE_GROUND** if any strong
ground event; then the two candidate verdicts. A strong event that matches neither candidate is
reported with its implied periods; we would call it a candidate second eclipse, not a confirmed
one, unless it is a TESS event that passes the peak-pixel check.

## 5. Limits fixed in advance

- A1 can fail: a precessing or evolving disk could make later eclipses shallower or change
  their shape. Vetoes need a dip below 5% (3% TESS) where the template is >= 15% deep, so they
  tolerate eclipses down to about a third of the K2 depth in the core.
- Strong TESS blending (more than ~80% of the aperture light from other stars) could make TESS
  vetoes wrong; the job reports the peak-pixel to aperture ratio per sector.
- An excluded period means "no eclipse at that period in these data", not that the companion
  does not exist; very long periods (> 1500 d) are not tested.
- Times: ZTF HJD, ASAS-SN HJD, TESS BJD (TDB) are mixed; differences are minutes, far below the
  3-hour veto zone.

## 6. Deviations

Any change after the job runs is logged in a dated DEVIATION file and reported with the
results. A failed download is retried in a new job version; it does not change the rules.

Before any run (2026-11-24, logged for completeness): the first publication attempt (10:01 UTC)
was refused because one unit test failed. Its input times 0.5 and 0.51 sat on a floating-point
bin edge (0.5 // 0.1 == 4.0 in Python, 0.51 // 0.1 == 5.0), so the test, not the code, was wrong.
The test input was moved mid-bin (0.53, 0.55). No rule, threshold or analysis code changed, and
no post-2021 data had been analysed. The job is filed as epic2202-eclipse-v2.

## 7. Credit

Setup: L. van der Kamp, D. M. van Dam, M. A. Kenworthy, E. E. Mamajek, G. Pojmanski (2022).
Data: K2 (NASA), TESS (NASA; MAST TESScut, Brasseur et al. 2019), ZTF (Caltech/IPAC IRSA),
ASAS-SN Sky Patrol v2 (Hart et al. 2023; Kochanek et al. 2017).
