# AGENTS.md

## Quick Start

```bash
python other_helper_files/check_env.py        # validates deps, writes env_report.json
python analyser.py                             # runs analysis (config hardcoded)
python plotter.py                              # generates plots (config hardcoded)
```

Neither `analyser.py` nor `plotter.py` accept CLI arguments — both hardcode their config paths.

## Environment Setup

On CERN machines, source CVMFS ROOT **before** running anything:
```bash
source other_helper_files/make_env.sh          # LCG_108a, ROOT + PyROOT
```

Requires: Python ≥3.9, ROOT ≥6.36, numpy 2.4+, pyyaml 6.0+, IPython 9.7+. `check_env.py` validates all of these.

## Run Order & Data

1. `python other_helper_files/check_env.py` — verify environment
2. `python analyser.py` — produces `./OutHistoFiles/hists_{sample}.root`
3. `python plotter.py` — reads ROOT hists, writes `./plots/` (PNG)
4. `python roc_and_rate/make_roc_and_rate.py` — ROC curves from `.pkl` snapshots
5. `python pypkg/post_analysis_persample.py` runs **inside** `analyser.py` step 2 — exports `.pkl` for ROC

**Output directory is deleted and recreated** each run (`analyser.py:38-41`). All previous hists are lost.

## Configuration

**`analysis_config.yaml`** — controls everything for the analyser:
- `general.input_dir_prefix`: prefix prepended to all sample paths (currently `./`)
- `general.output_dir`: `./OutHistoFiles`
- `general.tree_name`: `Events`
- `samples.{name}.input_file_pattern`: **glob** relative to prefix, e.g. `Phase2EGPuppiIso_151Xv3_Ntuplesv1base/MinBias_PU200/FP/v151Xv1base/perfNano_*.root`
- `samples.{name}.type`: `dytoll` → DY analysis, `qcd` → QCD analysis

Currently active samples: `DY_noPU`, `DY_PU200`, `MinBias`. Commented entries (e.g. `v3` variants) can be enabled.

**`plot_config.yaml`** — controls all plotting:
- `general.hist_folder`: `./OutHistoFiles/`
- `general.samples`: maps plot names to ROOT filenames (must match analyser output names)
- Each top-level key defines a plot group → output directory `./plots/{key}/`
- `hist_tag` uses composite syntax: `col[ele1, ele2]` expands to `col_ele1`, `col_ele2`
- Normalization schemes: `default_no_norm`, `hist_integral`, `summed_components`

## Architecture

**Entry point:** `analyser.py:18` → `analyser()`
- `define_cpp_utils()` on line 21 — loads ALL C++ functions via `ROOT.gInterpreter.Declare()`. **Must run before any RDataFrame code.**
- `SampleRDFManager` per sample (`an_specific_utilities.py:28`) wraps RDataFrame + hists
- Dispatches by type: `dy_to_ll_ana.dy_to_ll_ana_main()` or `qcd_ana.qcd_ana_main()`
- Post-analysis: `post_analysis_persample()` — saves `.pkl` snapshots for ROC curves

**Collection suffixes** (`an_specific_utilities.py:14-16`):
- `sufEl = 'TkEleL2'` — reconstructed track electrons
- `sufGen = 'GenEl'` — generator-level electrons
- `sufPu = 'L1PuppiCands'` — L1 PUPPI candidates

**DY analysis** (`dy_to_ll_ana.py:12`):
1. `GenEl_prompt==2 && |eta|<2.5` → `GenEl_DYP`
2. Filter: `GenEl_DYP_n>0 && TkEleL2_n>0` → `genDYP` dataframe
3. η split: EB (`|η|≤1.5`), EE (`1.5<|η|≤2.5`)
4. Gen-match with dR cuts: EB=0.04, EE=0.05
5. New collection: `*_MCH` (matched), stores `*_recoidx`, `*_genidx` vectors
6. Recalculate puppi isolation on matched electrons
7. **WARNING:** Lines 133-153 in `dy_to_ll_ana.py` are dead code (after `return`). References `dfgenEB` which is undefined.

**QCD analysis** (`qcd_ana.py:13`):
1. `TkEleL2_pt > 0` → `TkEleL2_Pt5`
2. η split: EB, EE
3. Recalculate puppi isolation
4. Line 34 (`anut.make_puppi_by_angdiff_from_tkel`) is a dead TODO comment — not called.

## C++ Functions

`define_cpp_utils.py` is a thin wrapper. Actual C++ code lives in:
- `cpp_utils/cmssw_cpp_utils.py` — base utilities (quantize, etc.), loaded first
- `cpp_utils/calciso_annularcone.py` — `calcisoannularcone()`, `calcisoanncone_singleobj()`
- `cpp_utils/tkelem_puppi_iso_def_legacy2026.py` — legacy 2026 iso: `calciso_legacy26()`, `calciso_legacy26_single()`
- `define_cpp_utils.py` — inline functions: `getminangs()`, `getmatchedidxs()`, `checksorting()`, `getpuppimask_annulardR()`, `getdR()`

**Load order matters** — functions depend on each other. `define_cpp_utils()` respects this order.

## Isolation Variables

`pypkg/calc_puppi_iso.py:recalculate_puppi_iso()` creates columns with **underscored dR values**:
- `reisotot_dRmin0_03_puppiIso` — recalculated total puppi iso (default drmin=0.03, drmax=0.2)
- `reisotot2026_dRmin0_03_puppiIso` — legacy 2026 version
- `reisochg_dRmin0_03_puppiIso`, `reisonut_dRmin0_03_puppiIso` — charged/neutral components
- `reisotot_dRmin0_03_puppiIsoDiff` — difference from Ntuple puppiIso

Dots in dR values become underscores. The regex `:dRmin\d_\d{1,2}` suffix selects them in `add_hists_multiplecolls`.

## Histogram Naming

`varmetadata.py` defines binning for all variables. Histograms follow `{filter_suffix}{collection}_{var}` pattern. Filter suffix comes from RDataFrame `GetFilterNames()[-1]`.

## Post-Processing

`post_analysis_persample.py` runs at the end of `analyser.py` per sample:
- `DY_PU200` → saves `DYPM{EB,EE}` dataframes to `DY_PU200.pkl`
- `MinBias` → saves `{EB,EE}` dataframes to `MinBias.pkl`
- These `.pkl` files feed ROC curve generation (`roc_and_rate/make_roc_and_rate.py`)

## Key Files

| File | Purpose |
|------|---------|
| `analyser.py` | Entry — loads config, C++, dispatches analysis |
| `analysis_config.yaml` | Sample paths, types, I/O settings |
| `plot_config.yaml` | Plot groups, hist tags, normalization |
| `dy_to_ll_ana.py` | DY→ll analysis pipeline |
| `qcd_ana.py` | QCD (MinBias) analysis pipeline |
| `an_specific_utilities.py` | `SampleRDFManager`, gen-match, angdiff |
| `rdf_generic.py` | `define_newcollection`, `add_hists_*`, `save_rdf_snapshot_to_pkl` |
| `varmetadata.py` | Histogram binning + labels for all variables |
| `pypkg/calc_puppi_iso.py` | Puppi isolation recalculation |
| `cpp_utils/*.py` | C++ isolation + utility code |
| `define_cpp_utils.py` | Orchestrator for ALL ROOT.C++ declarations |
| `pypkg/plotBeautifier.py` | Plot formatting (log axes, ranges, labels) |
| `pypkg/post_analysis_persample.py` | Pickle exports for ROC |
| `other_helper_files/make_env.sh` | CVMFS ROOT environment |

## Analysis Step Toggles

Commented blocks marked `STEP_X_X_X:` enable optional analysis stages (gen property plots, selection comparisons, etc.). Uncomment to activate in `dy_to_ll_ana.py` or `qcd_ana.py`.

## Common Pitfalls

- **`check_env.py` is in `other_helper_files/`**, not root — `python check_env.py` will fail
- **No CLI args** for `analyser.py` or `plotter.py` — edit the config files directly
- **C++ order matters** — `cmssw_cpp_utils()` before `calciso_annularcone()` before `tkelem_puppi_iso_def_legacy2026()`
- **`dy_to_ll_ana.py` has dead code** after line 110 — lines referencing `dfgenEB` will fail if uncommented
- **`qcd_ana.py` line 34** is a stale TODO, not functional code
- **Plot sample names** in `plot_config.yaml` must match analyser sample keys (not file paths)
- **`.gitignore`** excludes `*.root`, `*.pkl`, `*.png`, `*.DS_Store` — don't commit these
- **`plots/` and `OutHistoFiles/`** are recreated (deleted first) on each run
