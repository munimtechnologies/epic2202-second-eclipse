# epic2202-second-eclipse

Does the disk around EPIC 220208795's unseen companion come back on either published period?

In K2 data from 2016, the K dwarf EPIC 220208795 dimmed once, by about 28% for 7.2 hours. van der
Kamp et al. (2022, A&A 658, A38; [arXiv:2110.15086](https://arxiv.org/abs/2110.15086)) explained
the dip as a tilted disk around a companion. From sparse ground photometry they proposed two
orbital periods, 290.230 d and 235.587 d, and predicted the next eclipses. Nobody has published a
test of those predictions. Since then, TESS (Sectors 42, 43, 70 and 97), ZTF and ASAS-SN have
covered many of the predicted dates. The 235.587-d period predicts an eclipse on 2025-09-24, inside
TESS Sector 97.

This repository holds the pre-registered test. The rules are in [PREREGISTRATION.md](PREREGISTRATION.md),
fixed on 2026-11-24 before the analysis ran. Results will appear under `results/`.

- `ecl.py`: eclipse template, vetoes, period scan, dip and box searches (stdlib, tested)
- `io_ztf.py`, `asassn_io.py`, `io_tess.py`: public data access (ZTF IRSA, ASAS-SN Sky Patrol, MAST TESScut)
- `job.py`: the full run: `python job.py --data-dir /data --out-dir /out`
- `test_ecl.py`: `python -m unittest -v`

Credit: the setup, the eclipse template (restated from their BSD-2-Clause repository
lizvdkamp/EPIC2202_disk) and the period candidates are from L. van der Kamp, D. M. van Dam,
M. A. Kenworthy, E. E. Mamajek and G. Pojmanski. Data: K2 and TESS (NASA, MAST), ZTF (IRSA),
ASAS-SN Sky Patrol v2.

By MV2, Munim, Inc. (mv2@munimtech.com). Made with AI agents working in his name; no human expert
has checked it yet. Code Apache-2.0; text and results CC BY 4.0.
