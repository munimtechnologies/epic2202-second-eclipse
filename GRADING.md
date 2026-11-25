# Grading: epic2202-eclipse-v2 (pre-registered 2026-11-24)

Graded 2026-11-25 by AI agents working in MV2's name (MV2, Munim, Inc.). No human expert has checked
this yet. Job: `epic2202-eclipse-v2`, run on 2026-11-24 (~21 min, exit codes 0). Code at commit
9dfda1bf (this repository). Raw outputs: `results/v0.1.0/`.

## 1. Pre-registered outcome (mechanical, as the rules in PREREGISTRATION.md section 4 give it)

**NEW_ECLIPSE_TESS, P290_OPEN, P235_EXCLUDED.**

| Candidate | P0 (d) | Excluded share of P0 +- 0.10 d | Same, data before 2021-02-12 only | Verdict |
|---|---|---|---|---|
| P290 | 290.230 | 0.191 | 0.015 | **OPEN** |
| P235 | 235.587 | **1.000** | 0.058 | **EXCLUDED** |

- The data: 9,192 points from 2012-10-16 to 2025-11-08 (ZTF 1,099; ASAS-SN 1,157; TESS 30-min bins
  6,936 from Sectors 3, 30, 42, 70 and 97; Sector 43 had too few good cadences).
- Injection tests (template at half depth, 10 times per sector): recovered 8/8, 10/10, 9/9, 6/7, 10/10.
- Blind scan 60-1500 d (576,065 trial periods): **52.9%** excluded (primary rule), 67.2% with the
  one-point veto, and 28.4% with only data before 2021-02-12 (the cutoff of the discovery paper).

## 2. The 235.587-day period is excluded

That period puts eclipse n = 15 at BJD 2460942.805 (2025-09-24 07:19 UTC), inside TESS Sector 97.
The background-subtracted 3x3-aperture light curve in 30-minute bins
(`results/v0.1.0/tess_s097_30min.csv`, 9 cadences per bin) is flat across the whole night: 43 bins
from BJD 2460942.405 to 2460943.280 all lie between 0.9944 and 1.0063. The K2 eclipse was about 25%
deep (discovery paper); under assumption A1 (at least half the depth) a 12% dip would stand out at
about 40 times the bin scatter. Every one of the 1,001 trial periods within +- 0.10 d of 235.587 d is excluded (at n = 15
that window spans +- 1.5 days around the prediction), and the core of the window is also vetoed by
ZTF points taken the same night (`candidate_windows.csv`, veto_by `tess+ztf`). With only the data
the discovery paper could use (before 2021-02-12), just 5.8% of the window was excluded, so this is
new information.

Robustness notes (checked after grading; not pre-registered):
- **Time zero.** Table 8 of the discovery paper lists predicted eclipses that back-extrapolate to
  T0 = BJD 2457408.80, 0.175 d earlier than the template minimum we used (2457408.975). At n = 15 this
  equals a trial-period shift of 0.012 d, well inside the fully excluded window. The verdict does not
  depend on which T0 is used.
- **Blending.** A Gaia DR3 cone of 108" (queried 2026-11-25 through the Gaia archive TAP service:
  `SELECT ... FROM gaiadr3.gaia_source WHERE CONTAINS(POINT(ra,dec),CIRCLE(17.731554,0.314032,0.03))=1`)
  returns 16 sources; the target (G = 14.26) is the brightest, and the only other sources brighter
  than G = 18 are a G = 15.65 star about 66" away and a G = 17.53 star about 55" away. A 3x3 TESS
  aperture (63" square) centred on the target is dominated by it, so dilution cannot hide a deep
  eclipse.
- **Assumption A1** can fail: a precessing disk could make a later transit shallower or miss the star.
  "Excluded" means no eclipse at that period in these data, under A1.

## 3. The 290.230-day period stays open

81% of its window is still allowed, including P0 itself (allowed range 290.1456-290.3272 d). Its
eclipses since 2021 fell outside TESS sectors and between sparse ground points. Next predicted
minima for P0 and for the ends of the allowed range (`results.json`, `next_eclipses`):

| n | P = 290.1456 d | P = 290.230 d | P = 290.3272 d |
|---|---|---|---|
| 14 | 2027-03-06 12:19 UTC | 2027-03-07 16:40 UTC | 2027-03-09 01:20 UTC |
| 15 | 2027-12-21 15:48 UTC | 2027-12-22 22:12 UTC | 2027-12-24 09:11 UTC |

(n = 13, around 2026-05-20, falls after the end of these data, 2025-11-08.) The December 2027
window, three days wide, is the cleanest test of the surviving period.

## 4. The "new eclipse" flags: not credible as eclipses (post-run inspection)

The rules flag any strong TESS event. Three events were flagged. On inspection, none looks like a
second eclipse. This inspection was not pre-registered; the flags stay in the record, and a
pre-registered pixel-level check (AMENDMENT-2026-11-25-pixvet.md) is filed to settle them.

1. **Sector 30, BJD 2459126.95 (2020-10-04), depth 10.9%, S/N 80, ~3.4 h.** The dip is followed
   within 0.65 d by 30-minute bins at 1.63, 0.70, 0.59 and 1.75 times the median, which no star does;
   it sits about 1.2 d before the mid-sector data gap. The discovery paper notes that Sector 30 had
   "significantly more scattered light" and reports no eclipse in it. The single-pixel light curve
   shows 0.058 depth, 53% of the aperture depth. A real dip on the target gives the same fractional
   depth in both (ratio 1); an error in the per-cadence background, spread over the 9 aperture
   pixels, gives a ratio near 1/(9 x 0.29) = 0.38 for this sector's peak-pixel share of 0.29. The
   observed 0.53 is closer to the background case. It is also shallower and half as long as the K2
   eclipse (about 25%, 7.2 h).
2. **Sector 30, BJD 2459127.56 (S/N 61, not strong).** Same noisy stretch (bins at 0.59 and 1.75).
3. **Sector 42, BJD 2459456.14 (2021-08-29), depth 6.1%, S/N 34, ~3 h.** Within +- 1.6 d the light
   curve has bins at 1.078, 1.072, 0.920, 1.088, 1.111, 1.225 and 0.887: the dip is no larger than
   the neighbouring excursions.
4. **ASAS-SN ground event, 2022-12-25 01:50 UTC.** Three g-band exposures from one camera within
   3.7 minutes read 0.976, 0.309 and 0.554 of the median. A disk eclipse lasting hours cannot change
   by 67 percentage points in two minutes; this is a bad-image artefact.

Implied periods from these events (`results/v0.1.0/implied_periods.csv`) are kept for completeness
but are not candidates.

## 5. Summary

- One of the two published periods of the circumsecondary disk around EPIC 220208795 is ruled out:
  no eclipse on 2025-09-24, when the 235.587-d period required one, in clean TESS photometry.
- The 290.230-d period survives (81% of its window allowed); the next good test is 21-24 Dec 2027.
- About half of all periods from 60 to 1500 d are now excluded, against 28% with data to early 2021.
- No credible second eclipse was found. Three flagged events are most likely scattered light and
  bad images; the pixel-level check is pre-registered and queued.

## Credit

Setup, template and period candidates: L. van der Kamp, D. M. van Dam, M. A. Kenworthy,
E. E. Mamajek, G. Pojmanski (2022, A&A 658, A38; arXiv:2110.15086). Orbit used in simulations by
Mustill, Davies & Kenworthy (2024, arXiv:2404.12239). Data: TESS (NASA; MAST TESScut), ZTF (IRSA),
ASAS-SN Sky Patrol v2, Gaia DR3 (ESA).
