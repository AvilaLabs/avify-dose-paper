# Selected historical provenance excerpts

Verbatim excerpts from the original ledger; not the complete ledger. Historical claims do not supersede the revised manuscript.

## Original ledger lines 759–770

  the OpenMC 0.15.3 `heating` score for **photons** (auto-assigned the
  collision estimator; tracklength is refused) was measured at build time
  to report ~1e-4 of the physical photon energy deposition: 4.66 eV per
  source neutron for the whole phantom, against 5.6e5 eV of capture-photon
  energy produced per source neutron (`heating-local` minus `heating` for
  neutrons) and 6.2e4 eV deposited in the brain by photon flux x kerma.
  Neutron `heating` (MT301) is consistent across estimators and physically
  sized. W006 therefore scores the ROI photon dose (frozen tally T5) as
  photon flux x NIST photon kerma for brain — the estimator the reference
  itself uses (F44) — converted to eV deposited in the ROI so the W003
  weighted-dose formula and bracket construction are unchanged (tooling
  repair; the frozen quantity is the photon dose, the estimator changed).

## Original ledger lines 965–1014

## W006 C2 complete; F-W6-6 transcription repair (2 keV beam, published rows 14.2/14.6 cm) — 2026-08-17

- Record stage C2 completed 15:45:08 UTC (harness hash `74cd5e9d…`, seeds
  200–211; histories 1.3e8 for 0.0253 eV, 6e7 per keV/MeV beam, 2e7 per
  photon beam). Interim checker pass on the C1+C2 stage files (record
  evidence; scratch assembly): photon beams 0.2–10 MeV agree with the
  published profiles to ≤ 1.0 % at every point; the 1 keV, 10 keV, 100 keV
  and 1 MeV neutron beams agree within ±10 % at every published point
  (worst 9.3 %, 6.7 %, 6.1 %, 7.1 %); the 2 keV beam agrees at 106 of 111
  points; the 0.0253 eV beam agrees at 93 of 103 points.
- **F-W6-6 tracing, 2 keV beam.** The five non-agreeing points are the
  published rows at 14.2 cm (thermal, fast, gamma; ratio this-work/published
  5.6, 4.9, 6.0) and 14.6 cm (thermal 1.14, gamma 1.47) — the two depths
  that only the 2 keV sheet lists. Those published values contradict the
  deposit's own MCTAL tally for the same run and bin: published/MCTAL =
  0.181, 0.206, 0.175 at 14.2 cm and 0.897, 0.705 at 14.6 cm, and the 14.2 cm
  row is non-monotone by 4× against both neighbours on an otherwise smooth
  exponential. This work agrees with that MCTAL to 1.7–4.5 % at 14.2 cm and
  2.1–7.0 % at 14.6 cm. Disposition: a data-entry defect of the published
  table, i.e. a transcription error rather than a physics difference;
  tooling repair per F-W6-6. **Rule, applied uniformly to every beam and
  point by `check_w006.py`:** where a published value contradicts the
  deposit's MCTAL for the same run and bin by more than 10 % AND more than
  5σ of that tally, the MCTAL value is the reference at that point; every
  substitution is counted and printed, and the published-only outcome is
  recorded alongside. Verified reach: exactly the five css2 rows above and no
  other point of any beam (the next-largest published/MCTAL deviations, 5–8 %
  on deep fast-component points, are within 2σ of noisy tallies and are not
  substituted). With the repair the 2 keV beam agrees at 111/111 points
  (worst 8.7 %). Recorded candidly: the rule was formulated after the C2
  evidence existed, which F-W6-6 anticipates ("traced to"); it is objective,
  uniform, and cannot reach the 0.0253 eV beam (no matching MCTAL).
- **0.0253 eV beam: not repairable, failure stands under the frozen rule
  pending the full record.** The 10 non-agreeing points are thermal and
  boron kerma at 11.4–13.4 cm depth, this-work/published 1.10–1.13 with
  1σ ≈ 1.2 % (6–10σ), on a smooth monotone trend from 0.96 at 0.2 cm to 1.13
  at 13.4 cm — the signature of a ~3–4 % longer thermal-neutron diffusion
  length in this transport than in the reference, in the one configuration
  (pure thermal source, 13 cm of diffusion) that amplifies it most. Scratch
  diagnostics (not evidence): S(α,β) on/off changes deep values 2–3× (both
  codes used it); 294 vs 300 K is within noise; removing chlorine raises
  deep thermal 20–40 % and lowers photon kerma 10–15 %, so cross-library Cl
  capture/photon-production data are one plausible physics origin, the H2O
  S(α,β) evaluation another. The css253 published sheet is also the only one
  in the deposit with no matching tally file (the deposited css253.m is a
  192-cell layout from a different run) and its error column is internally
  inconsistent — a provenance caveat, not a transcription trace. Under §6,
  if this stands at the final check, W6-G1b fails and the study carries
  INCOMPLETE-VALIDATION regardless of the epithermal (G1a) outcome; a
  successor protocol may scope G1b, W006 may not.

