# Project Waddle — W006 protocol (FROZEN)

Session: W006 (reference study on the community benchmark phantom)
Date drafted: 2026-08-16 (UTC)
Date frozen: 2026-08-17 (UTC)
Status: **FROZEN.** Reviewed and adopted by the principal on 2026-08-17.
Section 2 transcription is complete and verified against the source paper's
own prose (§2.1). All declarations are populated with cited values and
recorded rationale, including DS-2 (§2.7), resolved from the full text.
No open placeholders remain.

Adoption asserts that the declared tumour geometry (§2.6), weighting set
(§2.8), declared uncertainty sets (§2.7) and ROI clinical criteria (§6.1)
are the principal's own choices, and that each judgment call flagged in
those sections — the tumour depth, the scalp CBE sitting near the low end
of its published range, the choice of escalation cohort for the brain
criterion, and the adverse finding recorded in §2.7 — has been read and
accepted.

**No scored run may precede the freeze hash recorded in `ledger.md`.**
Amendments to physics, gates, thresholds, declared sets, or criteria are
not permitted after this point; they close the session and open a successor
protocol. Repairs are append-only and tooling-scoped.

Purpose: run the verifier end-to-end on the BNCT community's own reference
phantom in order to (a) validate the instrument against published reference
dosimetry computed by an independent group with an independent transport
code, and (b) measure whether planning on a single assumed tumour-to-blood
ratio can return an action that certified bounds over the measured patient
range do not support.

Role: this is NOT a patentability or filing gate. W003 already closed the
value question. W006 is the product/publication deliverable — the artifact
that makes a pilot request credible and that can be submitted for
publication without any collaborator, patient data, or institutional
affiliation. All legal statuses remain `UNPROVED`.

## 1. Dependencies (hash-bound)

- Kernel: `theory/WADDLE_KERNEL.md`
  1dc8a70ee8195a12a6c184e2ae833a358fcfe30055bc8f91be1806c30daeccf7
- W003 close: **PASS** (soundness 9/9, gaps 11.4% / 3.6%). The bracket
  machinery under test here is the same machinery W003 validated; W006
  changes the geometry, the source, and the declared set, not the method.
- Infrastructure (declared, non-evidentiary): OpenMC 0.15.3 (conda-forge),
  ENDF/B-VII.1 HDF5, environment `~/.venvs/w003env`, data
  `~/nuclear-data/endfb-vii.1-hdf5`. Installed long before this draft;
  infrastructure setup is not evidence.

### 1.1 External reference (obtained; NOT redistributed)

Goorley JT, Kiger WS III, Zamenhof RG, "Reference dosimetry calculations
for neutron capture therapy with comparison of analytical and voxel
models," *Medical Physics* 29(2):145–156 (2002), DOI 10.1118/1.1428758.
AIP deposit E-MPHYA6-29-009201.

The publisher's supplementary archive was obtained 2026-08-16.

- File: `supplementary_material_1_1428758-sup-0001.zip`
- SHA-256:
  `511aa05c4253c77b2778a2bd7aa52334a1fd996cca3d08b7c1899f32db38fd62`
- Size: 14,658,000 bytes; 80 files; 108 MB extracted.

**This archive is deliberately NOT stored in this repository and NOT listed
in `MANIFEST.sha256`.** It is copyrighted supplementary material to a
paywalled article; the project holds it locally and binds it by hash and
DOI only. Any party with lawful access to the article can obtain the
identical file and verify the hash above, so reproducibility is preserved
without redistribution. Every transcribed value in §2 cites the archive
member it came from.

## 2. Frozen physical model

### 2.1 Phantom — transcribed from `MCNPinput/csse`

Modified Snyder head phantom, **analytic representation**: three nested
ellipsoids (MCNP `sq` quadrics), brain innermost, then cranium, then scalp.

| Region | Surface | `sq` coefficients (A, B, C) | Semi-axes x, y, z (cm) | Centre (cm) |
| --- | --- | --- | --- | --- |
| Brain (outer bound) | 10 | 0.0277778, 0.0123456, 0.0236686 | 6.0, 9.0, 6.5 | (0, 0, 1.0) |
| Cranium (outer bound) | 20 | 0.0216263, 0.0104123, 0.0145159 | 6.8, 9.8, 8.3 | (0, 0, 0) |
| Scalp (outer bound) | 30 | 0.0187652, 0.0094260, 0.0129133 | 7.3, 10.3, 8.8 | (0, 0, 0) |

Semi-axes are 1/√A, 1/√B, 1/√C and are recorded for readability; the `sq`
coefficients as transcribed are authoritative. Brain occupies the interior
of surface 10; cranium lies between 10 and 20; scalp between 20 and 30;
air outside 30.

**Transcription verified 2026-08-17 against the source paper's prose.**
Equations 1–3 of the paper give the boundaries independently of the MCNP
deck as (x/6)² + (y/9)² + ((z−1)/6.5)² = 1, (x/6.8)² + (y/9.8)² +
(z/8.3)² = 1, and (x/7.3)² + (y/10.3)² + (z/8.8)² = 1. These agree with the
deck-derived semi-axes above in every value. The paper further records that
the original Snyder phantom comprised two ellipsoids (cranium and brain)
and that the 5 mm scalp shell was a later addition — consistent with the
uniform 0.5 cm separation between surfaces 20 and 30 on all three axes.

Cross-model tolerances from the same source, recorded because they bound
what this project may claim: the 4 mm and 8 mm voxel models agree with the
analytic model within 2% and 4% respectively, while **the 16 mm voxel model
is reported to produce "unacceptably large discrepancies for all dose
components" and must not be used** for any cross-check. Updating kerma data
from ICRU 46 to ICRU 63 changes profiles by less than 2%, which is why the
G1 band is dominated by cross-code and cross-section-library differences
rather than by kerma provenance.

### 2.2 Materials — transcribed from `MCNPinput/csse` (ICRU 46)

Mass fractions, negative-density convention as in the source deck.

- **Brain** (ICRU 46 Adult Whole Brain), ρ = 1.040 g/cm³: H 0.107,
  C 0.145, N 0.022, O 0.712, Na 0.002, P 0.004, S 0.002, Cl 0.003, K 0.003.
- **Cranium** (ICRU 46 Skeleton, Whole Cranium), ρ = 1.610 g/cm³: H 0.050,
  C 0.212, N 0.040, O 0.435, Na 0.001, Mg 0.002, P 0.081, S 0.003,
  Ca 0.176.
- **Scalp** (ICRU 46 Adult Skin), ρ = 1.090 g/cm³: H 0.100, C 0.204,
  N 0.042, O 0.645, Na 0.002, P 0.001, S 0.002, Cl 0.003, K 0.001.
- **Air** (medical physics air), ρ = 1.293×10⁻³ g/cm³: C 0.00012,
  N 0.75527, O 0.23178, Ar 0.01283.
- **Thermal scattering:** `lwtr.01t` (light water S(α,β)) applied to
  hydrogen in all three tissue materials. The source paper examines
  alternative hydrogen scattering treatments; this protocol adopts the
  reference deck's choice and does not vary it.

**The reference deck contains boron tallies but NO boron in any material.**
Its own comment: *"For Boron Flux depression studies, Boron MUST be added
in the materials."* Boron loading is therefore this project's addition and
appears only in Configuration B and the C3 control.

### 2.3 Beam — transcribed from `MCNPinput/csse`

Generic epithermal neutron beam. Monodirectional disk source, plane
z = +15 cm, direction vector (0, 0, −1), radius 0–5 cm (10 cm diameter),
source area 78.53982 cm², coaxial with the phantom z-axis.

Energy composition: **10% thermal (< 0.5 eV), 89% epithermal, 1% fast
(> 100 keV)**, each component sampled with an equal-lethargy distribution
(100 uniform lethargy bins per component in the source deck).

Positive-control beams, same geometry, from the same suite:
monoenergetic neutrons at 0.0253 eV, 1, 2, 10, 100 and 1000 keV
(`css253`, `css1`, `css2`, `css10`, `css100`, `css1M`) and photons at 0.2,
0.5, 1, 2, 5 and 10 MeV (`cssp02` … `cssp10`).

### 2.4 Tally geometry — transcribed from the paper, §II

Kerma rates are scored in rectangular prisms of **16 × 16 mm cross-section
centred on the z-axis, 4 mm thick in z**. This project reproduces that tally
geometry exactly for G1; any deviation invalidates the comparison, because
kerma rates computed over different tally volumes are not comparable.

### 2.5 Kerma factors — from `Kermas/`

Neutron kerma factors for ICRU 46 adult brain based on ICRU 63 and
ENDF/B-VI cross sections and Q values, with S and Cl from JENDL-3.2.
Photon kerma factors from NIST mass energy-absorption coefficients. Boron
kerma supplied per 1 ppm B-10 and scaled linearly with concentration.

Note a deliberate asymmetry: the reference kermas derive from ENDF/B-VI
while this project's transport uses ENDF/B-VII.1. This is a genuine
cross-library comparison and is the principal reason the G1 tolerance is
set where it is. It is disclosed, not minimised.

### 2.6 Tumour region — DECLARED (adoption at freeze)

The Snyder phantom contains no tumour. One is declared here so that ROI
criteria can be evaluated. It is this project's extension to the benchmark
and is labelled as such wherever it appears; it is present in
Configuration B only.

- **Shape and size:** sphere, radius 1.5 cm.
- **Centre:** (0, 0, 3.8) cm, i.e. on the beam central axis.
- **Depth:** the scalp/air boundary on the +z axis is z = 8.8 cm (Eq. 3 of
  the source at x = y = 0), so the tumour centre lies at **5.0 cm depth**
  and the region spans 3.5–6.5 cm depth.
- **Material:** ICRU 46 adult whole brain, **identical to the surrounding
  brain material**. The tumour region differs from brain in its declared
  boron loading only.

Rationale, recorded before execution. Spherical tumours of this scale in
ellipsoidal head phantoms are the common construction in BNCT Monte Carlo
studies. The depth is the load-bearing choice: the source paper places the
thermal neutron dose peak at 2–4 cm and states that flux depression from
¹⁰B *"will reduce the ¹⁰B, thermal neutron, and induced photon dose
components at depth, often where the need to assure adequate dose to tumor
is most critical."* Siting the tumour just beyond the thermal peak puts it
where the coupling the kernel brackets is largest, which is the honest
place to test the instrument rather than the flattering one.

Giving the tumour the same tissue composition as brain is deliberate: it
avoids inventing a tumour composition absent from the benchmark, keeps the
validated material set intact, and leaves boron concentration as the only
variable distinguishing the ROIs — which is precisely the quantity under
declaration. The exact radius and depth are engineering judgment within
literature-typical values, not transcribed constants, and are labelled as
this project's choice.

### 2.7 Declared uncertainty sets

Two declared sets are run. Both are reported. Neither is selected after
seeing results.

- **DS-1 (continuity):** the W001/W003 bands — B ∈ [15, 25] ppm
  (nominal 20), ρ_T ∈ [2.5, 4.5] (nominal 3.5), ρ_M ∈ [0.8, 1.2]
  (nominal 1.0).
- **DS-2 (clinical, PET-derived):** blood boron and tumour-to-blood ratio
  over the ranges reported from ¹⁸F-BPA PET measurement, against the
  conventional planning assumptions of 25 ppm and a constant ratio of 3.5.

**RESOLVED 2026-08-17 from the full text.** Source: Kobayashi Y. et al.,
"Comparison of dose distribution with and without reflecting heterogeneous
boron distribution using ¹⁸F-BPA positron emission tomography in boron
neutron capture therapy," *Applied Radiation and Isotopes* 219 (2025)
111720. Retrospective study, **27 patients**, National Cancer Center
Hospital, Tokyo.

The parenthetical spread quoted in the abstract and body text is the
**1st–3rd quartile (interquartile range), not the full range** — the paper
states this explicitly: *"The median blood boron concentration (1st–3rd
quartile) and T/B ratio calculated from ¹⁸F-BPA were 25.57 (23.9–27.84) ppm
and 3.75 (2.54–4.59)."* Declaring DS-2 over those figures would therefore
have covered only the central half of the cohort while appearing to cover
all of it.

Table 3 of the same paper reports the **full observed ranges**, and those
are what DS-2 adopts:

| Quantity | Median | Interquartile range | **Full observed range (adopted)** |
| --- | --- | --- | --- |
| Blood boron concentration (aorta EBC) | 25.57 ppm | 23.90–27.84 | **19.36–31.55 ppm** |
| Tumour-to-blood ratio | 3.75 | 2.54–4.59 | **1.58–5.86** |

**DS-2 (declared):** blood boron ∈ [19.36, 31.55] ppm, tumour-to-blood
ratio ∈ [1.58, 5.86], against the conventional planning assumptions of a
uniform 25 ppm and a constant ratio of 3.5. Over these bounds the coverage
statement is exact: every patient in the cited cohort lies inside the
declared set. Normal-tissue ratios for brain and scalp remain declared from
cited literature.

Recorded consequence for DS-1: the W001/W003 band ρ_T ∈ [2.5, 4.5] does
**not** contain the observed patient range 1.58–5.86. DS-1 is retained
unchanged for continuity with the graded record, but it is not a
population-covering set and must never be described as one.

**Disclosed adverse finding from the same paper — recorded, not buried.**
The paper's own conclusion is that PET-derived blood boron and T/B ratio
were *"comparable to those in the conventional dose evaluations,"* that no
statistically significant differences in tumour dose indices were found
between conventional evaluation and evaluation using PET-derived blood
concentration and T/B ratio, and hence that *"inter-patient variations in
EBC in the blood would be negligible in the dose evaluation."* What the
paper did find significant (p < 0.01) was that reflecting **heterogeneous
boron distribution within the tumour** lowered all dose indices relative to
conventional evaluation.

This cuts against any naive framing that point-estimate planning is simply
wrong, and W006 must not adopt that framing. The defensible reading, which
this protocol adopts, is narrower and is exactly the distinction the
instrument exists to draw: the paper establishes **group-level statistical
adequacy of the median assumptions across a 27-patient cohort**, which is
not the same claim as **certified coverage for the individual plan in front
of the physicist**. A patient at a ratio of 1.58 is not the median patient,
the ratio cannot be measured at treatment time, and no patient-specific QA
exists to catch it. Whether that distinction produces an action divergence
on this phantom is precisely what W6-G4 tests, and a NULL-DIVERGENCE result
would be consistent with this paper's conclusion and must be reported as
such.

Normal-tissue boron ratios for brain and scalp are declared from cited
literature, not assumed.

### 2.8 Weighting — DECLARED (adoption at freeze)

The standard BPA weighting set for epithermal-beam BNCT, as established in
Coderre JA & Morris GM, "The radiation biology of boron neutron capture
therapy," *Radiation Research* 151:1–18 (1999), and used throughout the
clinical literature.

| Component | Factor | Applies to |
| --- | --- | --- |
| Photon dose | RBE 1.0 | all regions |
| Fast-neutron (recoil proton) dose | RBE 3.2 | all regions |
| ¹⁴N(n,p) nitrogen capture dose | RBE 3.2 | all regions |
| ¹⁰B(n,α) boron dose | CBE 3.8 | tumour ROI |
| ¹⁰B(n,α) boron dose | CBE 1.35 | brain |
| ¹⁰B(n,α) boron dose | CBE 2.5 | scalp |

Continuity check: W003 used (w_B, w_N, w_f, w_g) = (3.8, 3.2, 3.2, 1.0) for
its tumour ROI. The tumour weighting above is identical, so W006 changes
geometry, beam, and declared set without changing tumour weighting.

**Disclosed sensitivity, recorded before execution.** These are point
values drawn from measured distributions, not constants. Published brain
CBE for BPA spans roughly 1.31–1.48 from PET-derived normal-to-blood
ratios, bracketing the 1.35 adopted here. Published skin CBE spans roughly
2.4–3.7 depending on endpoint, and the adopted 2.5 sits near the **low**
end of that range — which understates biological dose to an
organ-at-risk and is therefore anti-conservative for the scalp criterion.
This is disclosed rather than minimised, and the scalp result is reported
with that caveat attached.

Noted for a successor study, not for W006: weighting factors are exactly
the kind of declared-interval quantity this instrument is built to bound.
Treating CBE as a declared set rather than a point value is a natural
extension and is deliberately **out of scope here**, so that W006 varies
boron alone and its result remains attributable.

## 3. Frozen tallies

Per ROI (declared tumour, brain, scalp): thermal flux (E < 0.5 eV);
¹⁰B(n,α) reaction rate; ¹⁴N(n,p) reaction rate; neutron heating above
10 keV; photon heating. Global: total ¹⁰B(n,α) and ¹H(n,γ) rates for the
photon decomposition. All tallies carry Monte Carlo standard errors.
Configuration A additionally scores central-axis depth-kerma profiles on
the §2.4 tally grid at the depths published in the reference.

## 4. Frozen bracket construction

Unchanged from W003 §"Frozen bracket construction" — mixed-corner boron
bracketing, nitrogen sandwich, fast-component applicability check with
declared remainder, two-source photon decomposition with non-negativity
check, statistical widening by 3σ. Any change to the method invalidates the
continuity claim with W003 and must be recorded as an amendment, not a
repair.

## 5. Frozen run schedule

1. **C1 setup control** (reduced histories): determinism under repeated
   seed; direction and zero-boron dominance checks as in W003 PC-B/PC-C.
2. **C2 monoenergetic control:** each monoenergetic beam of §2.3 on the
   unmodified phantom, compared against its published profile. Component-
   by-component validation before any composite spectrum is run.
3. **C3 flux-depression control:** the generic epithermal beam on the
   unmodified phantom with **30 µg/g ¹⁰B added to scalp and brain but not
   cranium**, compared against the run without boron, reproducing the
   published sensitivity study. See W6-G1c.
4. **Configuration A validation runs:** generic epithermal beam,
   unmodified phantom, histories sufficient for ≤1% relative error at each
   published comparison point.
5. **Configuration B scored block:** two designated extremal evaluations
   per declared set, plus interior evaluations at the box corners not used
   as extremal maps, the nominal point, and two declared interior points —
   matching the W003 interior schedule so soundness is tested identically.
6. **Preset-schedule comparator:** the customary practice — nominal point
   plus blood-concentration endpoints at nominal uptake ratio — evaluated
   directly, for G4.

Fixed seeds recorded per run. Results to `controls/w006_results.json`.
Pre-authorised amendment **AD-W6-1**: if any scored ROI tally exceeds 2.5%
relative error, double histories once for the affected run and re-run; both
runs logged, neither discarded.

## 6. Frozen gates

- **W6-G1a (external validation, epithermal; Configuration A).** Each
  computed central-axis quantity agrees with the published reference at
  every published depth within a pre-committed band of **10% for thermal
  neutron kerma rate and each individual dose component**. Rationale, fixed
  before execution: this is a cross-code and cross-library comparison
  (OpenMC with ENDF/B-VII.1 against MCNP4B with ENDF/B-VI and ICRU 63
  kermas), necessarily looser than the ≈4% the reference reports between
  its own analytic and voxel models, and tighter than the 5–10% the
  literature accepts for code-versus-measurement.
- **W6-G1b (external validation, monoenergetic).** Same tolerance, applied
  to the monoenergetic neutron and photon beams of C2.
- **W6-G1c (boron handling).** With 30 µg/g ¹⁰B in scalp and brain, the
  computed thermal neutron dose depression at the thermal peak (2–4 cm
  depth) falls in **8–11%**, the published result, within the G1 tolerance.
  This gate validates the project's boron handling specifically rather than
  its base transport, and is the closest available external check on the
  coupling mechanism the kernel depends upon.
- Failure of any G1 sub-gate → **INCOMPLETE-VALIDATION**. The study may not
  be described as validated against the benchmark, and no publication claim
  resting on that validation may be made.
- **W6-G2 (soundness).** For every interior evaluation, every ROI, and both
  declared sets, the directly computed weighted dose lies within
  [L_r − 3σ, U_r + 3σ]. Any violation → **FAIL-SOUNDNESS**. This would
  contradict W003 and is the most serious available outcome; report in
  full, do not repair.
- **W6-G3 (tightness).** Conservatism gaps computed as in W003.
  **PASS iff ≤ 0.35**, the same pre-committed threshold, for continuity.
  Exceeding → **FAIL-TIGHTNESS** on realistic geometry, which would
  materially weaken the product story even though W003 passed.
- **W6-G4 (clinical divergence; the study's reason for existing).**
  Under DS-2, does the preset-schedule comparator return PASS on a
  criterion for which the certified bounds do not establish PASS? Recorded
  per ROI and per criterion.
  - Divergence found → **DIVERGENCE-DEMONSTRATED**.
  - No divergence on any criterion → **NULL-DIVERGENCE**, reported as the
    headline result with equal prominence. A null is a real and publishable
    finding — it would say that on this geometry, point-estimate planning
    is adequate — and it must not be buried, re-run with different
    settings, or reframed. Criteria are declared before execution and not
    tuned afterwards.

### 6.1 ROI criteria for W6-G4 — DECLARED (adoption at freeze)

Fixed before execution and not tuned afterwards.

| ROI | Criterion | Value | Source |
| --- | --- | --- | --- |
| Tumour | minimum weighted dose | **≥ 20.0 Gy-w** | tumour coverage of ≥20 Gy-Eq mandated as an inclusion criterion in the multicentre dose-escalation phase 1 trial in recurrent high-grade glioma |
| Brain | maximum weighted dose | **≤ 11.0 Gy-w** | the middle cohort of the same trial's normal-brain D_max escalation (9, 11, 13 Gy-Eq) |
| Scalp | maximum weighted dose | **≤ 10.0 Gy-w** | reported clinical skin dose of 9.6 ± 1.4 Gy-Eq, rounded to the nearest whole unit above the mean |

Taking the tumour minimum and the brain maximum from the **same** trial is
deliberate: a criterion set drawn from one protocol is internally coherent
in a way that values mixed across institutions and eras are not. Normal-
brain constraints across the literature span roughly 7 Gy (Tsukuba phase I,
D2cc) to 13 Gy-Eq (top escalation cohort), and 11.0 sits mid-range; the
choice of cohort is this project's, and shifting it would shift the brain
criterion accordingly. The scalp value carries the additional CBE caveat of
§2.8.

These criteria exist to make W6-G4 decidable. They are not a treatment
prescription, not a recommendation, and not a standard; they are declared
thresholds for a phantom calculation.

- **Grade:** derived by `controls/check_w006.py` from
  `controls/w006_results.json` — recomputing brackets, gates, and grade
  from raw tallies, with nonzero exit on any mismatch. The grade is
  derived, never asserted.

## 7. Falsifiers (recorded before execution)

- **F-W6-1:** any interior evaluation escapes its bracket beyond combined
  3σ → FAIL-SOUNDNESS.
- **F-W6-2:** either conservatism gap exceeds 0.35 → FAIL-TIGHTNESS.
- **F-W6-3:** photon decomposition yields negative kernels beyond
  statistics → INCOMPLETE-TOOLING; recorded, no silent fallback.
- **F-W6-4:** no action divergence under DS-2 → NULL-DIVERGENCE, reported
  as the primary result.
- **F-W6-5 (RESOLVED, retained for the record):** the prior draft declared
  a fallback for the case where the reference data could not be obtained.
  The archive was obtained on 2026-08-16 and hash-bound in §1.1, so
  VALIDATED-EXTERNAL is available and the fallback does not apply. It is
  recorded rather than deleted because the falsifier existed before the
  data did.
- **F-W6-6:** if the transcribed model cannot reproduce the published
  profiles within the G1 band and the discrepancy is traced to a
  transcription error rather than a physics difference, the correction is
  a tooling repair, logged append-only; if it cannot be traced,
  INCOMPLETE-VALIDATION stands and is reported as such.

## 8. Prior-art note — read in full 2026-08-17 (recorded, not a gate)

**What §III.F of the source paper discloses.** Because of its very high
capture cross section, ¹⁰B in tissue alters neutron transport and causes
thermal neutron flux depression at clinically relevant concentrations. The
authors evaluated it by comparing thermal neutron depth-dose profiles in
the analytic model with and without **30 µg/g ¹⁰B added to skin and brain
but not bone** — stated to be about twice the normal-brain concentration
under BPA-F and about half the tumour concentration. Result: depression
increases with depth (neutrons thermalise deeper) and decreases with
increasing beam energy; at the thermal dose peak, **2–4 cm depth, the
depression was 8–11% for every beam evaluated**, consistent with earlier
published calculations the paper cites. The paper also explicitly notes the
effect reduces boron, thermal-neutron and induced-photon dose components
*at depth, where assuring adequate tumour dose is most critical.*

[Public release: internal patent-strategy discussion omitted. This copy is not the original frozen byte sequence; see PROVENANCE.md.]

## 9. Honesty clauses

The modified Snyder head phantom is a simplified three-shell ellipsoidal
model, not a patient geometry; the tumour region is this project's declared
extension and is not part of the published benchmark. G1 is a code-to-code
comparison, not validation against measurement — no measured dosimetry is
involved anywhere in W006, and the reference itself is a calculation. The
kerma libraries differ between reference and this project and that
difference is disclosed in §2.5 rather than minimised. Declared sets,
weighting factors, and clinical criteria are literature inputs, declared and
cited, not derived here. No patient data of any kind is used and none is
required. A divergence result demonstrates a mechanism on a phantom; it is
not a clinical finding, not a claim about any patient, plan, centre, or
planning system, and not evidence that any delivered treatment was
incorrect. Nothing in W006 is a medical device, informs a treatment
decision, or constitutes clinical advice. Protocol drafting and execution
occur in the same programme; hash ordering enforces protocol-before-outcome,
as in W001–W005. All legal statuses remain UNPROVED regardless of outcome.

## 10. Budget

Two bounded sessions: one to complete §2.6, §2.7, §2.8 and freeze; one to
execute. Repairs append-only and tooling-scoped. Amendments to physics,
gates, thresholds, or criteria are not permitted after the freeze — they
close the session and open a successor protocol.
