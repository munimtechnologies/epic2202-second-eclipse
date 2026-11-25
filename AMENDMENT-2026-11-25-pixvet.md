# Amendment 2026-11-25: pixel-level check of the flagged TESS dips

Frozen 2026-11-25 (UTC), after the main job (epic2202-eclipse-v2) was graded (GRADING.md) and
before any pixel-level data were analysed. MV2, Munim, Inc.; made with AI agents working in his
name; no human expert has checked it yet. Code: `pixvet.py` (rules, stdlib, tested in
`test_pixvet.py`) and `job_pixvet.py` (I/O). Job: `epic2202-pixvet-v1`.

## Why

The main job's rules flagged TESS dips in Sector 30 (BJD 2459127.00, 10.9% deep; 2459127.53) and
Sector 42 (BJD 2459456.12, 6.1%). Our inspection of the light curves (GRADING.md section 4) says they
look like scattered-light systematics, but that inspection was not pre-registered. If either dip is
on the target it would be the second eclipse ever seen of this system, so it gets a fixed test.
This amendment does not change any verdict of the main test.

## Data

MAST TESScut 21x21-pixel cutouts of Sectors 30, 42 and 97 around the target (same quality-bit mask
as io_tess.py). Only cadences within 0.5 d of each event are analysed.

## Measurements (per event)

- In-dip cadences: |t - t_mid| <= width/2 (t_mid and width from results/v0.1.0/results.json).
  Comparison cadences: 0.15 <= |t - t_mid| <= 0.50 d. All levels are medians.
- Target pixel: brightest of the central 3x3 of the whole-sector median image. Background: per-cadence
  median of the io_tess.py background pixels (central 11x11, >= 4 px from the target, below the
  60th percentile). This reproduces the main job's photometry.
- d_t: depth of the 3x3 aperture (background subtracted); d_1: depth of the target pixel;
  ratio_peak = d_1 / d_t (about 1 for a dip on the target; near 1/(9 x peak share) for a background
  error).
- corr_psf: minus the Pearson correlation, over the 5x5 pixels around the target, between the
  background-corrected difference image (in minus comparison) and the comparison image. A dip of the
  target's own light gives a value near +1.
- Control stars: up to 8 local maxima of the comparison image, >= 4 px from the target and >= 2 px
  from the edge, with 3x3 flux >= 5% of the target's; each gets its own local background (same rule)
  and depth. control_median_depth is their median.

## Labels (pixvet.label_event)

1. **REPRO_FAIL**: d_t < 0.03 (dip not reproduced).
2. **SYSTEMATIC**: at least 2 controls and control_median_depth >= 0.5 d_t, or corr_psf < 0.3.
3. **ON_TARGET**: corr_psf >= 0.7, ratio_peak >= 0.75, and (fewer than 2 controls, or
   |control_median_depth| <= 0.25 d_t).
4. **UNCLEAR**: anything else.

## Controls (gate)

Sector 97 at the 235.587-d prediction (BJD 2460942.805, width 0.12 d), where the light curve is flat:
- negative control, no change: must be **REPRO_FAIL**;
- positive control, a 10% dip of the target's own 5x5 image injected into the in-dip cadences: must be
  **ON_TARGET**.
Gate CONTROLS_OK only if both hold; otherwise the event labels are reported but not interpreted
(UNINTERPRETED).

## Overall

- **CANDIDATE_ECLIPSE**: gate passed and at least one event ON_TARGET. We would then report a
  candidate second eclipse (not a confirmed one) with its implied periods, and ask for independent
  photometry (SPOC/QLP/eleanor) before any stronger claim.
- **ALL_SYSTEMATIC**: gate passed and all three events SYSTEMATIC: the flags are closed as artefacts.
- **MIXED**: otherwise.

A failed download is retried as a new job version without changing these rules.
