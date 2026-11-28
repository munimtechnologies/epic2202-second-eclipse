# Grading of the pixel-level check of the flagged TESS dips (job epic2202-pixvet-v1)

Graded 2026-11-28 by AI agents working in MV2's name (MV2, Munim, Inc., mv2@munimtech.com), against
`AMENDMENT-2026-11-25-pixvet.md` (frozen 2026-11-25, before any pixel data were analysed). No human
expert has checked it yet. The job ran once on the Blade, 2026-11-27 23:46 to 2026-11-28 00:02 UTC
(unit tests and analysis exit 0). Its output is committed unchanged in `results/pixvet-v0.1.0/`.

## 1. Pre-registered outcome

**Gate `CONTROLS_OK`; overall `MIXED`.** No flagged event is `ON_TARGET`, so the outcome is not
`CANDIDATE_ECLIPSE`. Two events are `SYSTEMATIC` and one is `UNCLEAR`, so it is not `ALL_SYSTEMATIC`.

| event | label | d_t (3x3) | d_1 (peak pixel) | ratio_peak | corr_psf | control median depth | bkg change as aperture depth |
|---|---|---|---|---|---|---|---|
| S30 BJD 2459127.00 | SYSTEMATIC | 0.104 | 0.062 | 0.59 | 0.18 | 3.10 | 23.0 |
| S30 BJD 2459127.53 | SYSTEMATIC | 0.043 | 0.025 | 0.59 | 0.34 | 0.082 | 16.9 |
| S42 BJD 2459456.12 | UNCLEAR | 0.070 | 0.045 | 0.63 | 0.31 | -0.166 | 24.1 |
| S97 negative control | REPRO_FAIL | -0.002 | -0.002 | - | -0.32 | 0.005 | 0.0 |
| S97 positive control (10% injected) | ON_TARGET | 0.098 | 0.098 | 1.00 | 1.00 | 0.005 | 0.0 |

## 2. What it means

- The controls passed: an injected 10% dip of the target's own image is recovered at 0.098 with
  ratio_peak 1.00 and corr_psf 1.00, and the flat control window gives nothing.
- **None of the three flagged dips looks like the target dimming.** For a dip of the target's light,
  ratio_peak should be near 1 and corr_psf near +1, as in the positive control. The three flags have
  ratio_peak 0.59-0.63 and corr_psf 0.18-0.34. In each, the background changed by the equivalent of
  17-24 times the aperture flux between the in-dip and comparison cadences (reported, not used in a
  label), so a small background error explains the dips.
- The Sector 42 event is `UNCLEAR` only because its corr_psf (0.312) is just above the 0.3 cut and its
  control stars moved the other way (-0.17). It misses every `ON_TARGET` criterion.

**Verdict:** no flagged TESS dip is a candidate second eclipse of EPIC 220208795. The pre-registered
results of the main test (GRADING.md) are unchanged: the 235.587-d period is excluded, the 290.230-d
period survives (next tests 6-9 March 2027 and 21-24 December 2027), and 52.9% of periods from 60 to
1500 d are excluded.

## 3. Limits

- `MIXED` is the frozen label; the reading in section 2 is ours, made after the labels were computed.
- One background estimator, as fixed in the amendment; no pipeline (SPOC/QLP) light curves were used.
- Made with AI agents working in MV2's name; no human expert has checked it yet.
