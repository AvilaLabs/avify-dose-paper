# Avify Dose — paper reproducibility package

Avila Labs · Connor Avila · research@avilalabs.org

Code and recorded simulation results supporting **Two-evaluation dose envelopes for boron-uptake uncertainty in BNCT: a modified Snyder phantom feasibility study and reference-dosimetry comparison**.

## Scope

This is a phantom feasibility study. The reported total-dose envelopes are empirical, not established continuum-wide certificates or clinical plan acceptance tests. Photon-response approximations, fast-component remainders and simultaneous statistical coverage remain unresolved. The historical composite transport validation is **INCOMPLETE-VALIDATION**. The revised manuscript includes a conditional component-bound argument and the limitations. This repository is not a clinical treatment-planning product.

## Reproduce the numerical audit (no OpenMC required)

Python 3.10 or newer; standard library only:

```bash
git clone https://github.com/AvilaLabs/avify-dose-paper.git
cd avify-dose-paper
python3 scripts/verify.py
```

Expected: **93/93 checks passed**, with both historical validation grades still INCOMPLETE-VALIDATION. The audit re-derives quantities from recorded tally means and standard errors. It verifies file hashes and selected manuscript numbers; it does not prove physical accuracy or the total-dose bound. An absent external reference archive is reported as skipped, not verified.

## Rebuild figures and manuscript

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
PAPER_FIG_OUT=paper/manuscript/figures PAPER_FIG_DPI=600 python paper/make_paper_figures.py
python paper/manuscript/build_manuscript.py --figures end --stem manuscript_journal --pdf
python paper/manuscript/build_manuscript.py --figures inline --stem manuscript_preprint --no-line-numbers --spacing 1.15 --pdf
```

PDF conversion additionally requires LibreOffice (`soffice`). Omit `--pdf` to build Word files only. requirements.txt records the versions used for packaging; the standard-library audit does not require them.

## Files

- `controls/w006_results.json`, `w007_results.json`: recorded tally summaries, parameters and results.
- `controls/check_w006.py`, `check_w007.py`: independent numerical re-derivation.
- `controls/w006_run.py`, `w007_run.py`: original OpenMC transport harnesses.
- `controls/w006_reference.json`, `w006_extract_reference.py`: cited reference transcription and extraction code.
- `W006_PROTOCOL.md`, `W007_PROTOCOL.md`: scientific protocols; W006 is a clearly labelled public copy with internal patent discussion omitted.
- `paper/`: revised manuscript source, figure generator, document builder and numerical audit.
- `PROVENANCE.md`, `ledger.md`, `provenance/`: release adaptations and selected historical evidence.
- `external/README.md`: reference-archive hash and external transport-data requirements.

## Optional new transport runs

Full transport requires a matching OpenMC 0.15.3 installation and ENDF/B-VII.1 data; see external/README.md. It can require many hours and substantial storage. Keep new outputs outside the recorded controls directory:

```bash
mkdir -p artifacts/new-w006 artifacts/new-w007
export W006_OUT="$PWD/artifacts/new-w006"
export W006_RUNROOT="$PWD/artifacts/new-w006/runs"
export W006_SCRATCH="$PWD/artifacts/new-w006/scratch"
for stage in c1 c2 configA c3 scored assemble; do
  python controls/w006_run.py "$stage" || break
done
export W007_OUT="$PWD/artifacts/new-w007"
export W007_RUNROOT="$PWD/artifacts/new-w007/runs"
for stage in c2 c3 assemble; do
  python controls/w007_run.py "$stage" || break
done
```

W007 retains the recorded W006 baseline by design. The historical checkers are for the recorded campaign; outputs from new runs may differ and are not automatically a new validated study. Do not replace the recorded JSON to make a new run pass. Full transport was not rerun for this release.

## Citation and access

Use CITATION.cff and cite the versioned release. The GitHub release is a stable version link, not a DOI; no Zenodo deposit is claimed. Reference data remain attributable to their original authors. See RIGHTS.md for release terms.
