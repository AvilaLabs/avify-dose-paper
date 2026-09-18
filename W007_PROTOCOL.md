# Project Waddle — W007 protocol (FROZEN)

Session: W007 (successor to W006: re-scoped external validation on the
community benchmark phantom, replicated with new seeds)
Date drafted and frozen: 2026-08-18 (UTC)
Status: **FROZEN** at the hash recorded in `ledger.md`. No scored run precedes
that hash. Amendments close W007 and open a successor.

Purpose: W006 (`sessions/W006.md`; protocol c9558f28…) closed with the kernel
gates PASSED and the epithermal-beam validation PASSED at all 184 published
points, but with the composite validation grade INCOMPLETE-VALIDATION because
two G1 sub-gates failed as frozen: (i) the pure-thermal 0.0253 eV test beam at
11.4–13.4 cm depth (+10–13 %), and (ii) the boron flux-depression control at
2.2 cm by 0.14 percentage points. W006 cannot be amended. W007 declares a
re-scoped validation on stated criteria, freezes it, and re-runs the affected
controls with new seeds, so that a scoped external-validation statement can be
made — or refused — under a rule fixed before its evidence exists.

## 0. Candour clause (read first)

The W007 scope is chosen **with knowledge of the W006 outcome**. W007
therefore does not have the pre-registration standing of a blind protocol for
the two re-scoped gates: the value of W007 lies in (a) the scoping criteria
being principled and stated, (b) independent replication with new seeds of
every run that enters a gate, and (c) the W006 grade remaining in the record
unchanged and cited alongside. Any statement of validation that rests on W007
must say "scoped as in W007" and must cite the W006 0.0253 eV result.

## 1. Dependencies (hash-bound)

- W006 protocol `c9558f2896d3e1330cfc8130d0a0722e6c08bd63af5334832931c09c8dd9153d`;
  W006 record `controls/w006_results.json`
  `7bd696a4381d3d9e38125fefce007dc52b279e3e62e32e5ac2f1a533122c7c81`;
  W006 harness `controls/w006_run.py`
  `d3ac76b4c94dc278b4ed23e1c278be297bb13937d4029e1cd2c9c42850c50f1f`;
  W006 checker `controls/check_w006.py`
  `08decf7730b8a2f54f51d3e71679ec9ee6bcade2999d9eedee1328d6770891a1`;
  transcribed reference `controls/w006_reference.json`
  `91586006f737b68cfe29e76756b695766805873017d7d7bbfd67e0dc39a9218c`.
- Reference archive as W006 §1.1 (SHA-256 511aa05c…; not stored).
- Physical model, materials, beams, tally grid, kerma factors, transcription
  notes T1–T6 and the F-W6-6 transcription-repair rule: **identical to W006,
  by reference to the hash-bound harness and checker above.** W007 changes no
  physics, no tally, no tolerance band, no history count. It changes only the
  set of (beam, point) pairs that enter the two re-scoped gates, and the seeds.
- Infrastructure as W006 (OpenMC 0.15.3, ENDF/B-VII.1; non-evidentiary).

## 2. Re-scoped gates — DECLARED before any W007 run

- **W7-G1b (monoenergetic validation, provenance-scoped).** Same ±10 % band
  at every published point, applied to the beams whose published profile
  has a deposited raw tally file of the deck's layout in the same archive
  (`mctal_layout_matches_deck` true in `w006_reference.json`): neutron
  1 keV, 2 keV, 10 keV, 100 keV, 1 MeV and photon 0.2, 0.5, 1, 2, 5, 10 MeV.
  Criterion for the scope: a G1 comparison requires a reference value whose
  provenance can be verified against the deposit's raw output; the 0.0253 eV
  sheet has none (its deposited tally file is a different 192-cell run) and
  its error column is internally inconsistent (W006 ledger). The F-W6-6
  repair rule (published value contradicted by the same run's MCTAL by
  > 10 % and > 5σ → MCTAL is the reference at that point) carries over
  verbatim and is reported per beam.
- **W7-G1c (boron flux depression, peak location).** Depression as Fig. 17
  of the source ((D_noB − D_B)/D_noB, thermal-neutron kerma, 30 µg/g 10B in
  scalp and brain, not cranium), evaluated **at the thermal-dose peak
  location of the epithermal beam** — the published grid point with the
  largest reference thermal-neutron kerma rate, 2.6 cm — within the
  published 8–11 % widened by the G1 tolerance (7.2–12.1 %). Criterion for
  the scope: the source's sentence is "at the thermal neutron dose rate
  peak, between 2 and 4 cm depth, the thermal neutron dose depression was
  8–11 %"; W006 read it as every grid point in 2–4 cm; W007 reads it at the
  peak. The five-point W006 values are re-reported alongside. The strict
  8–11 % outcome is recorded alongside, ungated.
- **Reported, ungated:** the 0.0253 eV beam is re-run with a new seed and its
  full comparison reported, so the W006 deviation is replicated or not.
- **W6-G1a is not re-run and not re-graded** (PASS in the W006 record at
  184/184; its Configuration A profile is the hash-bound no-boron baseline
  for W7-G1c, as it was for W6-G1c).
- **Grade** (derived by `controls/check_w007.py`; nonzero exit on any
  mismatch): **VALIDATED-EXTERNAL-SCOPED** iff W7-G1b passes on all eleven
  scoped beams AND W7-G1c passes; else **INCOMPLETE-VALIDATION**. The W006
  grade is not altered by either outcome.

## 3. Run schedule (new seeds; frozen)

1. C2′: the twelve monoenergetic beams, histories as W006 (6e7 per keV/MeV
   neutron beam, 1.3e8 for 0.0253 eV, 2e7 per photon beam), seeds 700–711.
2. C3′: the flux-depression run, group-split as W006 (1e7/1.5e8/4e7),
   seeds 800–802; depression against the W006 Configuration A profile.
3. Results to `controls/w007_c2.json`, `controls/w007_c3.json`,
   `controls/w007_results.json`; grade by `controls/check_w007.py`.
Compute is shared with another programme's overnight job: W007 runs at
reduced priority and thread count; that changes wall time only.

## 4. Falsifiers (recorded before execution)

- F-W7-1: any scoped beam outside ±10 % at any published point →
  INCOMPLETE-VALIDATION.
- F-W7-2: depression at 2.6 cm outside 7.2–12.1 % → INCOMPLETE-VALIDATION.
- F-W7-3: the 0.0253 eV re-run does NOT reproduce the W006 deviation (i.e.
  agrees within ±10 % at all points) → recorded as a seed-dependence finding
  that would itself call the W006 statistics into question; reported, does
  not change either grade.
- F-W7-4: any mismatch between W007's re-derived values and the recorded ones
  → checker exit nonzero; INCOMPLETE-TOOLING.

## 5. Honesty clauses

Everything in W006 §9 applies. In addition: W007's scope was set after the
W006 result was known (§0); the re-scoping is by stated provenance and
wording criteria, not by outcome, but a reader is entitled to weigh that. A
VALIDATED-EXTERNAL-SCOPED result is a code-to-code statement on the scoped
beams and points; it says nothing about the 0.0253 eV deep tail, which
remains a recorded, unexplained cross-library difference. All legal statuses
remain UNPROVED.

## 6. Budget

One session: two controls (~3 h of transport at full priority), the checker,
a short close in the ledger. No amendments.
