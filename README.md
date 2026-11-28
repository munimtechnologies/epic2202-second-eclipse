# epic2202-second-eclipse

Does the disk around EPIC 220208795's unseen companion come back on either published period?

In K2 data from 2016, the K dwarf EPIC 220208795 dimmed once, by about 25% for 7.2 hours. van der
Kamp et al. (2022, A&A 658, A38; [arXiv:2110.15086](https://arxiv.org/abs/2110.15086)) explained
the dip as a tilted disk around a companion. From sparse ground photometry they proposed two
orbital periods, 290.230 d and 235.587 d, and predicted the next eclipses. Nobody had published a
test of those predictions. Since then, TESS (Sectors 42, 43, 70 and 97), ZTF and ASAS-SN have
covered many of the predicted dates.

## Result (graded 2026-11-25, see [GRADING.md](GRADING.md))

- **The 235.587-day period is excluded.** It required an eclipse on 2025-09-24; TESS Sector 97 shows
  the star flat to about 0.6% all that night (30-minute bins between 0.9944 and 1.0063), and every
  trial period within +- 0.10 d is ruled out.
- **The 290.230-day period stays open** (81% of its window allowed). Its next tests are
  6-9 March 2027 and 21-24 December 2027.
- 52.9% of all periods from 60 to 1500 d are now excluded, against 28.4% with data to early 2021.
- **No second eclipse.** The pre-registered pixel-level check of the three flagged TESS dips
  ([AMENDMENT-2026-11-25-pixvet.md](AMENDMENT-2026-11-25-pixvet.md), graded 2026-11-28 in
  [GRADING-pixvet.md](GRADING-pixvet.md), output in `results/pixvet-v0.1.0/`) passed its controls and
  found none of them on the target (overall label `MIXED`: two `SYSTEMATIC`, one `UNCLEAR`, none
  `ON_TARGET`). The ground event is a bad image.

Raw outputs: `results/v0.1.0/`. Rules: [PREREGISTRATION.md](PREREGISTRATION.md), fixed on 2026-11-24
before the analysis ran.

## Files

- `ecl.py`: eclipse template, vetoes, period scan, dip and box searches (stdlib, tested)
- `io_ztf.py`, `asassn_io.py`, `io_tess.py`: public data access (ZTF IRSA, ASAS-SN Sky Patrol, MAST TESScut)
- `job.py`: the full run: `python job.py --data-dir /data --out-dir /out`
- `pixvet.py`, `job_pixvet.py`: pixel-level vetting of flagged TESS dips (amendment of 2026-11-25)
- `test_ecl.py`, `test_pixvet.py`: `python -m unittest -v`

Credit: the setup, the eclipse template (restated from their BSD-2-Clause repository
lizvdkamp/EPIC2202_disk) and the period candidates are from L. van der Kamp, D. M. van Dam,
M. A. Kenworthy, E. E. Mamajek and G. Pojmanski. Data: K2 and TESS (NASA, MAST), ZTF (IRSA),
ASAS-SN Sky Patrol v2, Gaia DR3 (ESA).

By MV2, Munim, Inc. (mv2@munimtech.com). Made with AI agents working in his name; no human expert
has checked it yet. Code Apache-2.0; text and results CC BY 4.0.
