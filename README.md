# CHARMM starting-parameter pipeline

Reproducible Tier-1 generation of **starting CHARMM-compatible parameter sets** from MOL2 inputs using the official SwissParam command-line endpoints, followed by technical QC, MATCH-vs-MMFF comparison, manifest/resume logic, and per-molecule ZIP packaging.

> These are starting CHARMM-compatible parameter sets; scientific validation remains separate.

## Scientific boundary

Pipeline statuses are operational/technical only:

- `GENERATED`
- `TECHNICAL_PASS`
- `REVIEW`
- `FAILED`
- `SECOND_LEVEL_RUNNING`
- `SECOND_LEVEL_COMPLETE`

`TECHNICAL_PASS` does **not** mean production-ready or scientifically validated.

The batch also enforces an applicability boundary before SwissParam submission. A MOL2 graph containing a six-coordinate phosphorus atom bonded to exactly four nitrogen and two oxygen atoms (`P-N4-O2`) is routed to `REVIEW` with `PV_PORPHYRIN_CUSTOM_FF_REQUIRED`. This structural signature is an operational detector for the project's specialized P(V)-porphyrin route; it is not, by itself, a general oxidation-state assignment. No SwissParam health, submit, or polling request is made for that route.

The recommended specialized route is a custom CHARMM-compatible P(V)-porphyrin core plus CGenFF phenyl/pyridyl periphery with separate QM/MM validation. Expected molecular charge remains an explicit input-state contract in `config/states.tsv`; the graph detector does not infer charge.

## Repository layout

```text
input/
config/states.tsv
scripts/
results/{raw,match,mmff,review}/
packages/
reports/
manifests/
tests/
.github/workflows/parameterize.yml
```

MOL2 files are ignored by Git by default. This public repository does not automatically publish structures fetched from private Google Drive storage.

## Tier 1: SwissParam

For inputs inside the SwissParam applicability boundary, the client implements the official command-line flow:

1. graph applicability preflight;
2. health check;
3. submit one MOL2 with `approach=both`;
4. parse and persist Session number;
5. bounded status polling;
6. resume an unfinished session from the manifest;
7. retrieve `results.tar.gz`;
8. SHA256 and safe tar extraction;
9. output inventory before branch classification;
10. technical QC;
11. independent MATCH vs MMFF comparison;
12. per-molecule ZIP package.

Terminal SwissParam status payloads are classified before the wait budget expires when they contain known backend failure markers. The manifest preserves a concise summary of the actual backend error lines rather than replacing them with a timeout or generic failure phrase.

MATCH is the primary CHARMM branch. MMFF-based output is retained as an independent comparator. The pipeline never averages or automatically mixes the branches.

## Local setup

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[test]"
pytest -m "not integration"
```

On Windows PowerShell:

```powershell
py -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install -e ".[test]"
pytest -m "not integration"
```

## Environment inventory

Run on the actual workstation/WSL/server:

```bash
python -m scripts.inventory_environment \
  --scope-note "user workstation / WSL"
```

This probes existing Git, Conda/Mamba, Python, VMD, GROMACS and ORCA. It does not reinstall them or modify system Python.

## Input contract

Place MOL2 files in `input/`.

Minimum required sections:

```text
@<TRIPOS>MOLECULE
@<TRIPOS>ATOM
@<TRIPOS>BOND
```

The declared atom and bond counts must match the parsed blocks.

Optional `config/states.tsv`:

```text
molecule	state	expected_charge	headgroup_note
ALC0315_neutral	neutral	0	tertiary amine neutral
ALC0315_protonated	protonated	1	protonated tertiary amine
```

If `states.tsv` has no row for an input, the batch continues with `EXPECTED_CHARGE_UNKNOWN`. Protonation, H atoms, bond orders, and input atom names are never silently edited.

For a specialized `P-N4-O2` route, a configured expected charge is compared with the MOL2 charge sum during preflight. A mismatch is recorded as `INPUT_CHARGE_MISMATCH`; the input is preserved unchanged and remains in `REVIEW`.

## Single explicit integration test

External SwissParam calls are never part of ordinary unit tests.

```bash
export RUN_SWISSPARAM_INTEGRATION=1
export SWISSPARAM_TEST_MOL2=/absolute/path/to/molecule.mol2
pytest -m integration -k swissparam_single_molecule_submission -s
```

## Batch

Only after the single-molecule integration test:

```bash
python -m scripts.batch_parameterize \
  --input-dir input \
  --states config/states.tsv
```

Default polling interval is 15 s. Total wait is bounded and configurable:

```bash
python -m scripts.batch_parameterize \
  --poll-interval 15 \
  --max-total-wait 7200
```

Idempotency key:

```text
input SHA256
+ generator
+ approach
+ pipeline version
```

Completed unchanged inputs are skipped. An unfinished manifest entry with a Session number resumes polling instead of submitting a duplicate job.

## Technical QC

The QC layer checks:

- required MOL2 sections;
- declared vs parsed atom/bond counts;
- duplicate atom names;
- charge readability;
- output inventory;
- MATCH branch presence;
- input/MATCH/MMFF atom counts when output MOL2 is available;
- expected integer charge when configured, tolerance `1e-4`;
- log issue tokens with context.

Duplicate input atom names are preserved in the original MOL2. When they prevent reliable name mapping, the molecule is routed to `REVIEW`.

## Package format

```text
packages/<MOLECULE>_CHARMM_STARTING_PARAMETERS.zip
```

Each package contains the immutable original MOL2, raw SwissParam archive, detected MATCH/MMFF files, comparison report, technical QC, issues table, and molecule-level manifest entry.

## GitHub Actions

The workflow is manual `workflow_dispatch` only. It has optional `input_path` and `force_rerun` inputs and uploads packages, reports, manifest, and raw SwissParam archives as a GitHub artifact.

There is no schedule and no automatic GitHub Release publication.

## Tier 2

Tier 2 is only for `REVIEW` molecules.

The project specification prefers `FFParam-v2`, but the authoritative `FFParam-v2` CLI has not yet been identified unambiguously. `scripts/second_level.py` therefore performs discovery and captures actual `--help` output when a candidate executable exists. It does not invent CLI flags.
