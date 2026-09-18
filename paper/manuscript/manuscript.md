# Two-evaluation dose envelopes for boron-uptake uncertainty in BNCT: a modified Snyder phantom feasibility study and reference-dosimetry comparison

AUTHOR: Connor Avila
AFFIL: Avila Labs, Oviedo, FL, USA
EMAIL: research@avilalabs.org
ORCID: 0009-0000-9957-9857

## Abstract

Boron neutron capture therapy (BNCT) dose calculations depend on uncertain tissue boron uptake. We investigate a two-evaluation dose-envelope construction and benchmark an independent OpenMC 0.15.3/ENDF-B-VII.1 model against the modified Snyder head-phantom reference dosimetry. Pure added absorption supports conditional bounds for boron and nitrogen components; the implemented total-dose envelope additionally uses fitted photon responses and a fast-component approximation, for which continuum-wide coverage is not established. All 184 epithermal-beam depth–kerma-rate points agreed within ±10 % (worst 6.0 %) at ≤ 1 % relative standard error; six photon beams agreed to ≤ 1.1 %. The full predeclared validation gate was not met. Thermal-neutron dose depression by 30 µg/g boron rose from 3 % at the surface to 29 % at 13 cm. For a declared tumour and two uncertainty sets, including a box spanning cohort-reported blood-boron and tumour-to-blood ranges, all 22 non-corner evaluations were compatible with the envelopes under the recorded 3σ comparison rule. Gaps to sampled dose extrema were 8.7 %/14.1 % for the tumour lower edge and 2.5–4.9 % for normal-tissue upper edges. For the cohort-derived box, the tumour lower edge was 0.44× the smallest dose in a fixed-ratio preset schedule. The two methods disagree for nominal tumour mean doses 1.2–2.8× an illustrative 20 Gy-w criterion. These phantom results motivate further investigation; they do not establish certified total-dose bounds, joint statistical coverage or clinical plan acceptability.

**Keywords:** boron neutron capture therapy; treatment plan verification; Monte Carlo; OpenMC; boron uptake uncertainty; boron-uptake sensitivity

## 1. Introduction

The weighted dose in BNCT combines boron capture, nitrogen capture, fast-neutron recoil and photon components with declared RBE/CBE factors (Coderre and Morris, 1999). Blood samples inform treatment-time boron estimates, while tissue concentrations additionally depend on uptake ratios and spatial distributions (Kim et al., 2026). Boron absorption also changes the neutron field, so dose is not generally affine in concentration. PET-informed uptake estimates and cohort variation therefore motivate evaluating more than a nominal boron map (Kobayashi et al., 2025).

Independent BNCT calculation is established: Hu et al. (2021) compared a clinical planning system against an independently developed PHITS model. Machine dosimetric quality assurance is also described by Hirose et al. (2022). The question here is how efficiently two evaluations can construct informative dose envelopes over a declared uptake set. Our preset comparator uses nominal T/B 3.5 and blood endpoints; it is an explicit study comparator, not a demonstrated universal clinical practice. A finite scenario list alone does not establish coverage between evaluated scenarios.

We (i) compare an independent transport model with community reference dosimetry (Goorley et al., 2002), (ii) evaluate a two-corner envelope on that phantom, distinguishing conditional component bounds from approximate total-dose construction, and (iii) compare its lower edge with the preset schedule. The protocol and subsequent pre-block declarations are described in Section 2.5. The present interpretation distinguishes empirical checks from mathematical or statistical certification; it does not change the historical protocol, tallies or grades.

## 2. Methods

### 2.1 Two-evaluation envelope and its scope

Let *c* be the boron concentration map and *C* the declared set {*c* : *c*_{r} = ρ_{r}·*B*, ρ_{r} ∈ [ρ_{r}^{−}, ρ_{r}^{+}], *B* ∈ [*B*^{−}, *B*^{+}]} for tissue compartments *r*. The extreme maps *c*^{−} and *c*^{+} are evaluated by Monte Carlo. Under a fixed-source, fixed-scattering transport model with boron entering only as added nonnegative absorption, the flux decreases pointwise as *c* increases. Consequently, the nitrogen component has endpoint bounds and the boron component has mixed-product bounds: ρ_{r}^{−}*B*^{−}·*K*_{r}(*c*^{+}) ≤ *D*_{B,r}(*c*) ≤ ρ_{r}^{+}*B*^{+}·*K*_{r}(*c*^{−}), where *K*_{r} is a positive per-unit-concentration dose response. Appendix A states the assumptions and argument. This is a statement about exact transport responses under that model, not noisy estimates or all components of the implemented dose.

The recorded photon construction fits two regional coefficients to endpoint photon doses using global hydrogen- and boron-capture totals, with a boron photon branching factor of 0.94. Positive fitted coefficients and agreement at the fitted endpoints do not establish a response valid for other boron maps: regional photon dose depends on source position and spectrum, and additional capture sources such as chlorine require treatment. Rescaling global boron captures by a total-boron-atom ratio also requires justification when compartment concentrations vary by different factors. We therefore treat this photon construction as an approximation. The fast component uses the endpoint mean and an ε_{f} = 2 % allowance; endpoint agreement checks this approximation only at the sampled corners.

The historical implementation combines these components into [*L*_{r}, *U*_{r}] and applies a 3σ expansion. We call the result a two-evaluation total-dose envelope, not a certified interval. Its uncertainty calculation does not establish simultaneous coverage: same-history tally covariances, uncertainty in fitted photon coefficients and deterministic approximation errors require separate treatment (Appendix A). A criterion met by an envelope edge is an illustrative envelope comparison, not proof of PASS or FAIL for all *c* ∈ *C*. Fig. 1 illustrates the construction for DS-2.

### 2.2 Reference phantom, beams, tallies and kerma factors

The analytic modified Snyder head phantom (three nested ellipsoids: brain, 6 × 9 × 6.5 cm centred 1 cm above the origin; cranium 6.8 × 9.8 × 8.3 cm; scalp 7.3 × 10.3 × 8.8 cm), the ICRU 46 compositions and densities (brain 1.040, cranium 1.610, skin 1.090 g cm^{−3}) (ICRU, 1992), the light-water S(α,β) treatment on hydrogen in all tissues, the beam (monodirectional 10 cm disk at *z* = +15 cm directed along −*z*), the axial tally grid (16 × 16 mm cross-section, 4 mm thick, 46 bins from 0.6 cm upstream of the scalp to 17.4 cm depth) and the kerma factors (^{10}B per ppm; ICRU 46 adult-brain neutron kerma from ICRU 63 / ENDF/B-VI with S and Cl from JENDL-3.2, extrapolated as 1/*v* to 10^{−4} eV (ICRU, 2000); NIST photon kerma for brain) were transcribed from the deposited MCNP input deck of Goorley et al. (2002) (Fig. 2). The generic epithermal beam is 10 % thermal (1 meV–0.5 eV), 89 % epithermal (0.5 eV–10 keV) and 1 % fast (10 keV–2 MeV), each group sampled with 100 equal-lethargy bins as in the deck; the reference run's energy-group source biasing was reproduced by running the three groups separately and combining them with their physical fractions. Monoenergetic control beams were the deck's 0.0253 eV, 1, 2, 10, 100 and 1000 keV neutron beams and 0.2, 0.5, 1, 2, 5 and 10 MeV photon beams. Transport used OpenMC 0.15.3 (Romano et al., 2015) with ENDF/B-VII.1 (Chadwick et al., 2011) in coupled neutron–photon mode on 8 threads; kerma-rate tallies were flux × dose function (log–log interpolation, endpoint rule as in MCNP), normalised to 10^{10} n cm^{−2} s^{−1} over the beam area as in the reference. Histories were 3 × 10^{7} / 6 × 10^{8} / 1.5 × 10^{8} for the three groups of the epithermal beam (≤ 1 % relative standard error at every published point), 6–13 × 10^{7} per monoenergetic neutron beam and 2 × 10^{7} per photon beam.

Two properties of the deposited reference were established during transcription and are reported: (a) for the epithermal beam and the 10 MeV photon beam the published table equals the deposited MCNP tally output exactly, whereas for the other monoenergetic beams the published table derives from a different run and differs from the deposited tally by up to ~4 % (0.2 cm depth) and typically < 0.5 % elsewhere; the published table is the comparison target throughout, with one exception: two rows of the 2 keV sheet (14.2 and 14.6 cm) contradict the deposit's own tally for the same run by a factor of up to 5.6 and are non-monotone by 4× against their neighbours, and were replaced by the deposited tally under a uniform rule (published value contradicted by > 10 % and > 5σ) that touched no other point of any beam; (b) the deposited tally file for the 0.0253 eV beam is from a different, 192-cell run and cannot be mapped to the published sheet, whose error column is internally inconsistent.

Photon dose estimator: the OpenMC `heating` score for photons was found in this build to report ~10^{−4} of the physical photon energy deposition (4.7 eV per source neutron for the whole phantom against 5.6 × 10^{5} eV of capture-photon energy produced); photon dose was therefore scored as photon flux × NIST kerma, the estimator the reference itself uses. Neutron heating (MT 301) was consistent across estimators.

### 2.3 Declared extensions: tumour, uncertainty sets, weights, criteria

The Snyder phantom contains no tumour; a 1.5 cm-radius sphere of brain composition centred on the beam axis at 5.0 cm depth (spanning 3.5–6.5 cm, just beyond the thermal peak) was declared, differing from brain only in boron loading. Two uncertainty sets were declared and both run: DS-1, the continuity bands of the programme's earlier surrogate benchmark (*B* ∈ [15, 25] ppm, ρ_{T} ∈ [2.5, 4.5], scalp ρ_{S} ∈ [0.8, 1.2]; nominal 20 ppm, 3.5, 1.0); DS-2, the full observed range of a 27-patient ^{18}F-BPA PET cohort (Kobayashi et al., 2025) (*B* ∈ [19.36, 31.55] ppm, ρ_{T} ∈ [1.58, 5.86], the same scalp band; nominal = the conventional planning assumptions 25 ppm, 3.5, 1.0). We note that the abstract of that paper quotes interquartile ranges (23.9–27.84 ppm, 2.54–4.59); its Table 3 gives the full ranges adopted here, so the box spans the reported marginal ranges. This Cartesian box may include combinations not observed together; it is neither a patient-specific confidence region nor a coverage guarantee for future patients or treatment-time tissue uptake. Brain follows blood at ratio 1.0. Weighting followed Coderre and Morris (1999): photon 1.0, fast and nitrogen 3.2, boron CBE 3.8 (tumour), 1.35 (brain), 2.5 (scalp; near the low end of published skin values, anti-conservative for that organ). Criteria for the divergence analysis: tumour ≥ 20.0 Gy-w, brain ≤ 11.0 Gy-w, scalp ≤ 10.0 Gy-w. The tumour and brain values are taken from the same multicentre phase 1 dose-escalation trial in recurrent high-grade glioma (Kim et al., 2026), whose eligibility required ≥ 90 % of the planning target volume to receive at least 20 Gy-Eq at an assumed T/B of ~3.5 and whose normal-brain maximum-dose cohorts were 9, 11 and 13 Gy-Eq (11 Gy-Eq selected as the recommended phase 2 dose); we use these numerical thresholds only as illustrative ROI-mean study criteria. A mean-dose requirement and a D90 requirement do not generally imply one another, and a whole-brain mean cannot substitute for the trial's maximum-dose constraint. The scalp value rounds up the mean maximum skin dose of 9.6 ± 1.4 Gy (weighted) reported for eight glioblastoma patients treated with BSH/BPA BNCT (Matsuda et al., 2009). ROI doses are whole-ROI mean weighted doses. The quantities *L*_{tumour} and *U*_{OAR} are constructed envelope edges, not spatial minimum/maximum doses or established extrema over the uptake set. Thresholds are literature-inspired illustrative inputs, not clinical acceptance criteria validated here.

### 2.4 Run schedule

C1 setup control (determinism; thermal flux and per-atom boron kernel lower at *c*^{+} than at *c*^{−} by > 3σ; zero-boron dominance; both sets, three ROIs); C2 twelve monoenergetic beams; Configuration A, the epithermal beam on the unmodified phantom (the external comparison); C3 flux-depression control (30 µg/g ^{10}B in scalp and brain, not cranium, as in the reference's own sensitivity study); Configuration B, the scored block: per set, two corner evaluations at 4 × 10^{7} histories, six mixed corners, the nominal point, two interior points and the two preset-comparator blood endpoints at 2 × 10^{7}, with a pre-authorised one-time doubling of any run whose worst ROI tally exceeded 2.5 % (nine such re-runs fired; both runs logged). A successor protocol (W007) re-ran C2 and C3 with independent seeds.

### 2.5 Pre-registration and independent checking

The protocol — physical model, tallies, bracket construction, run schedule, gates, tolerances, declared sets, weights and criteria — was frozen and its SHA-256 recorded before any run; results are checked by an independent script that re-derives every profile, bracket, dose and gate from raw tallies and derives the grade, exiting nonzero on any mismatch. Three quantities were found at build time to be absent from the frozen protocol (normal-tissue ratios; the plan normalisation that makes the Gy-w criteria decidable; the ROI-mean reading) and paused execution until they were adopted, before the scored block ran; the divergence result is therefore reported as a normalisation-free window rather than a binary. The composite pre-registered comparison gate — every published point of every beam within ±10 %, and the depression at each grid point in 2–4 cm within 7.2–12.1 % — was not met in W006 (the 0.0253 eV beam's deep tail; the 2.2 cm depression at 7.06 %) nor in W007 (five noise-limited deep fast-neutron points), so under the programme's own rule the study carries the internal label INCOMPLETE-VALIDATION; the per-gate facts are reported here in full and the label is not used as a claim. The run harness, checker and analysis scripts were written with the assistance of a large-language-model coding assistant (Claude, Anthropic) under the author's direction; every quantity reported here is re-derived from the raw tallies by the independent checker in the hashed record.

## 3. Results

Table 1 summarises the outcomes; Sections 3.1–3.4 give the detail.

TABLE1

### 3.1 External comparison

Fig. 3 shows the epithermal-beam profiles. All 184 published points agree within ±10 % (worst 6.0 %); this work's relative statistical error is ≤ 0.95 % at every point. Thermal-neutron and boron kerma agree within 2.5 % from the surface to 14 cm (ratio 0.984–1.023) and rise to +4–6 % in the deepest bins; fast-neutron kerma agrees within 3 % everywhere (0.996–1.026); photon kerma is +1 to +6 % throughout. Fig. 4 shows the monoenergetic beams for two independent seeds: the six photon beams agree to ≤ 1.1 % at every point; the 100 keV and 1 MeV neutron beams within 8 %; the 1, 2 and 10 keV beams within 10 % at every point with the first seed and at all but five deep fast-neutron points with the second (ratios 0.73–0.89, −1.5σ to −3.8σ of the combined statistics), where both codes carry 2–9 % statistics and the two seeds themselves differ by 3–26 %; the pure-thermal 0.0253 eV beam shows a smooth excess in thermal and boron kerma growing from 0.96 at the surface to 1.13 at 13.4 cm (replicated: 1.12 with the second seed), the signature of a ~3–4 % longer thermal diffusion length in this transport than in the reference. Scratch diagnostics attribute both this and the systematic photon excess plausibly to cross-library differences in chlorine capture/photon-production data and the H_{2}O S(α,β) evaluation (removing Cl from the tissues raises deep thermal kerma 20–40 % and lowers photon kerma 10–15 %); neither is resolved here.

### 3.2 Boron flux depression

Fig. 5: with 30 µg/g ^{10}B in scalp and brain the thermal-neutron dose is depressed by 3.2 % at 0.2 cm, 7.1–9.9 % across 2.2–3.8 cm (7.8 % at the 2.6 cm thermal-dose peak; 7.9 % with the second seed), 11 % at 4.6 cm and ~29 % at 13–14 cm, reproduced to within 0.15 percentage points by the second seed. The reference reports 8–11 % at the peak across all its beams and an increase with depth; the shape is reproduced and the magnitude agrees to about one percentage point.

### 3.3 Empirical envelope checks and sampled gaps

Fig. 6: for both sets and all three ROIs, every non-corner evaluation (nine interior points, two comparator points and the re-runs per set; 22 in total excluding re-runs) is compatible with its envelope under the recorded comparison rule, which allows both the envelope expansion and an additional 3σ uncertainty for the evaluated dose. The fitted photon coefficients are positive and the fast endpoint check holds; neither check establishes continuum-wide validity. Conservatism gaps — the distance from the unexpanded envelope edge to the extreme of the direct scan, relative to that extreme — are 8.7 % (DS-1) and 14.1 % (DS-2) on the tumour lower bound and 4.2 % / 4.9 % (brain), 2.5 % / 2.8 % (scalp) on the upper bounds. The smallest evaluated tumour dose lies at the all-low corner, with the mixed corner (low ρ_{T}, high ρ_{S}, low *B*) indistinguishable from it within statistics; the largest evaluated scalp dose lies at a mixed corner (low ρ_{T}, high ρ_{S}, high *B*), marginally above the all-high corner. At 10^{10} n cm^{−2} s^{−1} the DS-2 envelopes (with 3σ) are tumour [6.8, 41.7], brain [0.87, 1.23], scalp [0.48, 0.87] Gy-w min^{−1} against nominal 19.1, 1.04 and 0.65.

### 3.4 Divergence from the preset comparator

Fig. 7: under DS-2 the tumour lower envelope edge (with 3σ) is 0.436× the worst case of the preset schedule (nominal + blood endpoints at T/B 3.5); under DS-1, 0.685×. Because both comparisons scale with the plan, this fixes a window in the plan's nominal-assumption tumour ROI dose: for DS-2, any plan delivering 24.6–56.4 Gy-w at nominal assumptions — 1.2 to 2.8 times a 20 Gy-w minimum criterion — meets the preset check while its lower envelope edge is below the illustrative threshold; for DS-1, 25.3–36.9 Gy-w. The brain and scalp windows lie at 161–256 Gy-w nominal tumour dose: with these whole-ROI mean doses, the OAR thresholds do not bind within the reported tumour divergence windows. This does not establish hotspot control or clinical OAR acceptability.

## 4. Discussion

The transport model reproduces the community reference on the beam that matters to within a few per cent over four decades of kerma, with statistics better than the reference's own; the residual differences are systematic, replicate across seeds, and are of the size expected between ENDF/B-VI-era and ENDF/B-VII.1 libraries. The 0.0253 eV deep tail is the most sensitive possible probe of thermal transport (13 cm of pure diffusion) and its 3–4 % diffusion-length difference merits a targeted study; it is not a feature of clinical beams. The photon excess of +1 to +6 % on the epithermal beam (up to +9 % on the monoenergetic neutron beams) points to capture-photon production data (chlorine is a substantial contributor in ICRU 46 tissues) rather than photon transport, since the photon-only beams agree to ~1 %.

The two-evaluation construction was compatible with all recorded non-corner evaluations under the declared uncertainty-tolerant rule. Gaps from the unexpanded edges to sampled extrema were 9–14 % for the tumour and 3–5 % for normal tissue. The finite scan does not identify true extrema or prove coverage of the set. Each corner evaluation (4 × 10^{7} histories, three ROIs, ≤ 2.5 % on every tally) took 11–20 min on eight CPU threads of a partly shared workstation. The construction uses two endpoint runs, but the cost of scoring increasingly detailed spatial models and of establishing coverage has not been measured here.

The divergence window quantifies sensitivity to the uncertainty set and comparator choice. The preset schedule holds T/B at 3.5, whereas DS-2 includes lower and higher ratios from a cohort. Disagreement therefore cannot be attributed solely to the envelope algorithm, nor does a lower edge below threshold demonstrate an actual failing plan. A comparison with a denser scenario schedule varying both blood and tissue ratios is needed to assess computational benefit. Cohort marginal ranges do not establish within-patient uncertainty or a treatment-time joint distribution.

The most useful next steps are a source-resolved photon bound, a separate statistical and model-error budget, and a spatial-dose example using tumour coverage and normal-tissue hotspot endpoints. More successful sampled cases cannot substitute for the missing total-dose argument. These analyses are not reported as completed in this study.

Limitations include a code-to-code comparison, a phantom with a declared tumour, whole-ROI mean doses, fixed biological weights and illustrative thresholds. The composite predeclared validation gate was not met; residual library differences and noise-limited discrepancies remain reported outcomes. The work does not establish clinical plan acceptability or a certified total-dose interval.

## 5. Conclusions

The OpenMC/ENDF-B-VII.1 model reproduces all 184 epithermal reference points within the declared tolerance, while the full validation suite remains incomplete. A two-evaluation total-dose envelope is compatible with the sampled phantom results under the recorded 3σ rule and differs materially from a fixed-ratio preset comparison. Conditional absorption monotonicity supports boron and nitrogen component bounds, but photon-response approximations and unresolved statistical coverage prevent a certified total-dose claim. The evidence supports a feasibility study and targeted methodological follow-up.

## Declaration of competing interest

The author is the founder of Avila Labs LLC, which has a commercial interest in plan-verification software. The author is the named inventor on United States provisional patent application No. 64/132,705, filed August 13, 2026, concerning the method described in this manuscript. The author reports no other competing interests.

## Funding

This research received no external funding.

## Declaration of generative AI and AI-assisted technologies in the writing process

During the preparation of this work the author used Claude (Anthropic) and Codex (OpenAI) to draft and edit the manuscript text from the author's pre-registered study record, to prepare the figures and to format the references. The author is responsible for final verification of the calculations, sources, interpretation and manuscript.

## Data availability

The reproducibility package is available at https://github.com/AvilaLabs/avify-dose-paper/releases/tag/v1.0.0. It includes recorded tally summaries, results, transport harnesses, independent numerical checkers, figure scripts and scientific protocols. The W006 protocol is a labelled public copy with internal patent-strategy discussion omitted; original and public file hashes and checker adaptations are documented in the repository. Full OpenMC statepoint files and transport nuclear-data libraries are not included. The release is publicly accessible with rights reserved; no DOI has been assigned. The reference archive of Goorley et al. (2002) is the publisher's supplementary material, cited here by DOI and SHA-256 (511aa05c4253c77b2778a2bd7aa52334a1fd996cca3d08b7c1899f32db38fd62) and not redistributed.

## Appendix A. Conditional component bounds and unresolved coverage

### A.1 Exact transport statement

Assume a fixed geometry, fixed nonnegative external source, fixed nonnegative scattering kernel and fixed boundary conditions, with no boron-dependent neutron multiplication or secondary neutron source. Let *T*_{c} be the positive free-flight operator whose attenuation includes the added boron absorption, and *S* the fixed scattering operator. Assume the collision expansion converges. If *c*^{−} ≤ *c* ≤ *c*^{+} pointwise, then *T*(*c*^{+}) ≤ *T*(*c*) ≤ *T*(*c*^{−}) on nonnegative sources. Induction preserves this ordering in each term (*T*_{c}*S*)^{n}*T*_{c}*q*; summing gives ψ(*c*^{+}) ≤ ψ(*c*) ≤ ψ(*c*^{−}) in position, direction and energy. The conclusion relies on fixed scattering rather than only a small concentration change.

Any fixed positive nitrogen-dose response preserves that ordering. For a spatially uniform compartment concentration *c*_{r}, write *D*_{B,r}(*c*) = *c*_{r}*K*_{r}[ψ(*c*)], with *K*_{r} a fixed positive response including the energy-dependent reaction weighting. Positivity yields *c*_{r}^{−}*K*_{r}[ψ(*c*^{+})] ≤ *D*_{B,r}(*c*) ≤ *c*_{r}^{+}*K*_{r}[ψ(*c*^{−})]. For a heterogeneous map the corresponding concentration bounds must remain inside the spatial integral. A ratio of total atom inventories cannot generally replace these pointwise factors. These are exact-model component inequalities; they do not establish that every physical approximation or Monte Carlo estimator meets the assumptions.

### A.2 Why the present total-dose envelope remains empirical

A photon-dose proof would need ordered spatial and spectral photon-source bounds and a fixed positive photon-response operator, or explicit bounds on changes in that operator. Fitting two coefficients to global endpoint capture totals does not supply this information and does not isolate other photon-production channels. A fast-component allowance likewise needs a uniform bound or an explicitly limited empirical interpretation. The present manuscript retains the recorded computation for reproducibility rather than asserting these missing results.

### A.3 Statistical interpretation

The recorded 3σ expansion is an approximate uncertainty allowance, not an absolute guarantee or an established simultaneous confidence region. The historical implementation combines component errors in quadrature, uses endpoint photon errors without full propagation through the fitted response coefficients, and incorporates a deterministic fast-component allowance within its standard-error calculation. Correlations, coefficient conditioning and model discrepancies can invalidate a nominal coverage interpretation. The empirical comparison additionally allows 3σ on each direct evaluation. No joint confidence level is assigned to the resulting family of comparisons here.

A future certificate would first bound physical approximation errors separately, then propagate tally uncertainty through the full construction with covariance or a justified conservative alternative, and allocate a declared error probability across all reported inequalities. Those steps require further analysis and potentially additional tally information; relabelling the current expansion does not complete them.

## References

Chadwick, M.B., Herman, M., Obložinský, P., Dunn, M.E., Danon, Y., Kahler, A.C., Smith, D.L., Pritychenko, B., Arbanas, G., Arcilla, R., Brewer, R., Brown, D.A., Capote, R., Carlson, A.D., Cho, Y.S., Derrien, H., Guber, K., Hale, G.M., Hoblit, S., Holloway, S., Johnson, T.D., Kawano, T., Kiedrowski, B.C., Kim, H., Kunieda, S., Larson, N.M., Leal, L., Lestone, J.P., Little, R.C., McCutchan, E.A., MacFarlane, R.E., MacInnes, M., Mattoon, C.M., McKnight, R.D., Mughabghab, S.F., Nobre, G.P.A., Palmiotti, G., Palumbo, A., Pigni, M.T., Pronyaev, V.G., Sayer, R.O., Sonzogni, A.A., Summers, N.C., Talou, P., Thompson, I.J., Trkov, A., Vogt, R.L., van der Marck, S.C., Wallner, A., White, M.C., Wiarda, D., Young, P.G., 2011. ENDF/B-VII.1 nuclear data for science and technology: cross sections, covariances, fission product yields and decay data. Nucl. Data Sheets 112, 2887–2996. https://doi.org/10.1016/j.nds.2011.11.002

Coderre, J.A., Morris, G.M., 1999. The radiation biology of boron neutron capture therapy. Radiat. Res. 151, 1–18. https://doi.org/10.2307/3579742

Goorley, J.T., Kiger, W.S. III, Zamenhof, R.G., 2002. Reference dosimetry calculations for neutron capture therapy with comparison of analytical and voxel models. Med. Phys. 29, 145–156. https://doi.org/10.1118/1.1428758

Hirose, K., Kato, T., Harada, T., Motoyanagi, T., Tanaka, H., Takeuchi, A., Kato, R., Komori, S., Yamazaki, Y., Arai, K., Kadoya, N., Sato, M., Takai, Y., 2022. Determining a methodology of dosimetric quality assurance for commercially available accelerator-based boron neutron capture therapy system. J. Radiat. Res. 63, 620–635. https://doi.org/10.1093/jrr/rrac030

Hu, N., Tanaka, H., Kakino, R., Yoshikawa, S., Miyao, M., Akita, K., Isohashi, K., Aihara, T., Nihei, K., Ono, K., 2021. Evaluation of a treatment planning system developed for clinical boron neutron capture therapy and validation against an independent Monte Carlo dose calculation system. Radiat. Oncol. 16, 243. https://doi.org/10.1186/s13014-021-01968-2

ICRU, 1992. Photon, Electron, Proton and Neutron Interaction Data for Body Tissues. ICRU Report 46. International Commission on Radiation Units and Measurements, Bethesda, MD.

ICRU, 2000. Nuclear Data for Neutron and Proton Radiotherapy and for Radiation Protection. ICRU Report 63. International Commission on Radiation Units and Measurements, Bethesda, MD.

Kim, W., Park, J.-S., Song, J.-H., Yoo, H., Park, K., Kim, H.J., Shin, D.-W., Lee, S.U., Ahn, S., Ha, S., Yi, J., Cho, K., Seo, H.J., Lim, H.-S., Yee, G.-T., 2026. Boron neutron capture therapy in recurrent high-grade gliomas: safety, efficacy, and pharmacokinetics from a multicenter, dose-escalation phase 1 trial. Adv. Radiat. Oncol. 11, 101947. https://doi.org/10.1016/j.adro.2025.101947

Kobayashi, Y., Nakamura, S., Takemori, M., Nakaichi, T., Shuto, Y., Ito, K., Takahashi, K., Kashihara, T., Yonemura, M., Endo, H., Kunito, K., Okamoto, H., Chiba, T., Nakayama, H., Oshika, R., Kishida, H., Itami, J., Kurihara, H., Igaki, H., 2025. Comparison of dose distribution with and without reflecting heterogeneous boron distribution using ^{18}F-BPA positron emission tomography in boron neutron capture therapy. Appl. Radiat. Isot. 219, 111720. https://doi.org/10.1016/j.apradiso.2025.111720

Matsuda, M., Yamamoto, T., Kumada, H., Nakai, K., Shirakawa, M., Tsurubuchi, T., Matsumura, A., 2009. Dose distribution and clinical response of glioblastoma treated with boron neutron capture therapy. Appl. Radiat. Isot. 67, S19–S21. https://doi.org/10.1016/j.apradiso.2009.03.054

Romano, P.K., Horelik, N.E., Herman, B.R., Nelson, A.G., Forget, B., Smith, K., 2015. OpenMC: A state-of-the-art Monte Carlo code for research and development. Ann. Nucl. Energy 82, 90–97. https://doi.org/10.1016/j.anucene.2014.07.048

Stern, R.L., Heaton, R., Fraser, M.W., Goddu, S.M., Kirby, T.H., Lam, K.L., Molineu, A., Zhu, T.C., 2011. Verification of monitor unit calculations for non-IMRT clinical radiotherapy: Report of AAPM Task Group 114. Med. Phys. 38, 504–530. https://doi.org/10.1118/1.3521473

## Table

TABLE1: Table 1. Summary of outcomes (all quantities from the checked record; W006 seeds unless stated).
| Item | Result |
|---|---|
| Epithermal beam vs reference (4 components × 46 depths) | 184/184 within ±10 %; worst 6.0 %; this work ≤ 0.95 % rel. std. error |
| Photon beams 0.2–10 MeV | worst 0.1–1.1 % (both seeds) |
| Neutron beams 100 keV, 1 MeV | worst 6.1–7.9 % and 5.7–7.1 % over the two seeds |
| Neutron beams 1, 2, 10 keV | within 10 % at every point (first seed); five deep fast-neutron points outside with the second seed (ratios 0.73–0.89; both codes at 2–9 % statistics) |
| Neutron beam 0.0253 eV | thermal/boron kerma +10–13 % at 11.4–13.4 cm (both seeds) |
| Depression, 30 µg/g ^{10}B, 2.2–3.8 cm | 7.1–9.9 % (7.8–7.9 % at the 2.6 cm peak); 29 % at 13–14 cm |
| Empirical compatibility (22 non-corner evaluations, 3 ROIs, 2 sets) | 0 failures under the recorded 3σ comparison rule |
| Gap to sampled extrema: tumour lower-bound gap | 8.7 % (DS-1), 14.1 % (DS-2) |
| Gap to sampled extrema: brain / scalp upper-bound gap | 4.2 / 2.5 % (DS-1), 4.9 / 2.8 % (DS-2) |
| Tumour lower envelope edge ÷ preset worst case | 0.685 (DS-1), 0.436 (DS-2) |
| Divergence window (nominal tumour ROI dose) | [25.3, 36.9) Gy-w (DS-1), [24.6, 56.4) Gy-w (DS-2) |
| Corner-evaluation cost | 11–20 min each, 8 threads, 4 × 10^{7} histories |

## Figure captions

![Fig. 1. Left: the declared set DS-2 in (ρ_{T}, *B*), with the two corner evaluations *c*^{−}, *c*^{+}, mixed corners, nominal, and the three preset-comparator points. Right: tumour ROI weighted dose rate from every scored evaluation versus tumour ^{10}B concentration, with the approximate envelope constructed from the two corner runs; the annotated gap is the lower envelope edge relative to the comparator's worst case (0.44).](figures/fig1_concept.png)

![Fig. 2. The analytic modified Snyder head phantom in the *x*–*z* plane, with the beam, the axial tally column and the declared tumour.](figures/fig2_geometry.png)

![Fig. 3. Epithermal-beam depth–kerma-rate profiles (thermal-neutron, boron per ppm, fast-neutron, photon) versus the published reference, and their ratio with the ±10 % band; error bars 1σ combined.](figures/fig3_epithermal.png)

![Fig. 4. Monoenergetic beams: range of this-work/reference over all published depths per beam and component, for two independent seeds.](figures/fig4_mono.png)

![Fig. 5. Thermal-neutron dose depression from 30 µg/g ^{10}B in scalp and brain versus depth, two seeds; the published 8–11 % peak-region range.](figures/fig5_depression.png)

![Fig. 6. Approximate total-dose envelopes and all direct evaluations for both sets and three ROIs, normalised to each set's nominal point; gaps from unexpanded edges to sampled extrema in the panel titles; compatibility permits 3σ on each evaluated dose as well as the envelope expansion.](figures/fig6_bounds.png)

![Fig. 7. Illustrative threshold comparisons using the preset schedule and lower envelope edge as a function of the plan's nominal-assumption tumour ROI dose; the shaded window is where they diverge.](figures/fig7_window.png)
