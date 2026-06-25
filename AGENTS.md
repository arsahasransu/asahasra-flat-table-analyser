# AGENTS.md

## Quick Start

```bash
source other_helper_files/make_env.sh          # CERN only: LCG_108a, ROOT + PyROOT
python other_helper_files/check_env.py         # validates deps, writes env_report.json
python analyser.py                             # runs analysis (config hardcoded)
python plotter.py                              # generates plots (config hardcoded)
```

Neither `analyser.py` nor `plotter.py` accept CLI arguments — config file path is hardcoded. Edit `analysis_config.yaml` or `plot_config.yaml` directly.

## Environment Setup

On CERN machines, source CVMFS ROOT **before** running anything:
```bash
source other_helper_files/make_env.sh          # LCG_108a, x86_64-el9-gcc14
```

`check_env.py` validates: `packaging`, `numpy >= 2.1.0`, `pyarrow`, `ROOT`, `uproot`, `pyyaml >= 6.0.0`. Exits non-zero on failure (CI-friendly). Requires: Python ≥ 3.9, ROOT ≥ 6.36.

## Run Order & Data

1. `python other_helper_files/check_env.py` — verify environment
2. `python analyser.py` — produces `./OutHistoFiles/hists_{sample}.root` + `*_snapshot.root` parquet files
3. `python plotter.py` — reads ROOT hists, writes `./plots/` (PNG)
4. `cd roc_and_rate && python make_roc_and_rate.py` — ROC curves from `_snapshot.root` files

**Output directory is deleted and recreated** each run (`analyser.py:38-41`). `plots/` is also recreated by `plotter.py`.

## Configuration

**`analysis_config.yaml`** — controls everything for the analyser:
- `general.input_dir_prefix`: prefix prepended to all sample paths (currently `./`)
- `general.output_dir`: `./OutHistoFiles`
- `general.tree_name`: `Events`
- `samples.{name}.input_file_pattern`: **glob** relative to prefix
- `samples.{name}.type`: `dytoll` → DY analysis, `qcd` → QCD analysis

Active samples: `DY_noPU`, `DY_PU200`, `MinBias`.

**`plot_config.yaml`** — controls all plotting:
- `general.hist_folder`: `./OutHistoFiles/`
- `general.samples`: maps plot names to ROOT filenames (must match analyser sample keys)
- Each top-level key defines a plot group → output directory `./plots/{key}/`
- `hist_tag` uses composite syntax: `col[ele1, ele2]` expands to `col_ele1`, `col_ele2`
- Normalization: `default_no_norm`, `hist_integral`, `summed_components`

## Architecture

**Entry point:** `analyser.py:18` → `analyser()`
- `ROOT.EnableImplicitMT()` on line 20 — multi-threaded RDataFrame
- `define_cpp_utils()` on line 21 — loads ALL C++ via `ROOT.gInterpreter.Declare()`. **Must run before any RDataFrame code.**
- `SampleRDFManager` per sample (`an_specific_utilities.py:28`) wraps RDataFrame + hists
- Dispatches by type: `dy_to_ll_ana.dy_to_ll_ana_main()` or `qcd_ana.qcd_ana_main()`
- Post-analysis: `post_analysis_persample()` — saves `_snapshot.root` + `.parquet` for ROC

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
7. **WARNING:** Lines 133-153 are dead code (after `return` on line 111). References `dfgenEB` which is undefined.

**QCD analysis** (`qcd_ana.py:13`):
1. `TkEleL2_pt > 0` → `TkEleL2_Pt5`
2. η split: EB, EE
3. Recalculate puppi isolation
4. Line 34 (`anut.make_puppi_by_angdiff_from_tkel`) is a dead TODO comment — not called.

## C++ Functions

`define_cpp_utils.py` is a thin wrapper. Actual C++ code lives in:
- `cpp_utils/cmssw_cpp_utils.py` — base utilities, loaded first
- `cpp_utils/calciso_annularcone.py` — `calcisoannularcone()`, `calcisoanncone_singleobj()`
- `cpp_utils/tkelem_puppi_iso_def_legacy2026.py` — legacy 2026 iso
- `define_cpp_utils.py` — inline: `getdR()`, `getminangs()`, `getmatchedidxs()`, `checksorting()`, `getpuppimask_annulardR()`

**Load order matters** — `cmssw_cpp_utils()` → `calciso_annularcone()` → `tkelem_puppi_iso_def_legacy2026()`.

## Isolation Variables

`pypkg/calc_puppi_iso.py:recalculate_puppi_iso()` creates columns with **underscored dR values**:
- `reisotot_dRmin0_03_puppiIso` — recalculated total (default drmin=0.03, drmax=0.2)
- `reisotot2026_dRmin0_03_puppiIso` — legacy 2026 version
- `reisochg_dRmin0_03_puppiIso`, `reisonut_dRmin0_03_puppiIso` — charged/neutral
- Dots become underscores. Regex `:dRmin\d_\d{1,2}` selects them in `add_hists_multiplecolls`.

## Snapshot & ROC

Post-analysis (`pypkg/post_analysis_persample.py`) saves via `rdf_generic.save_rdf_snapshot()`:
- `DY_PU200` → `DY_PU200_{EB,EE}_snapshot.root` + `.parquet`
- `MinBias` → `MinBias_{EB,EE}_snapshot.root` + `.parquet`

`roc_and_rate/make_roc_and_rate.py` reads `_snapshot.root` files via `rdf_generic.load_rdf_snapshot_from_root()`. Requires `matplotlib`, `uproot`. Uses `sys.path` hack to import `rdf_generic` from parent dir. Imports `my_py_ai_utils` from same directory.

## Key Files

| File | Purpose |
|------|---------|
| `analyser.py` | Entry — loads config, C++, dispatches analysis |
| `analysis_config.yaml` | Sample paths, types, I/O settings |
| `plot_config.yaml` | Plot groups, hist tags, normalization |
| `plotter.py` | Plot generator (also has `auto_singlehist_plotter()` — commented) |
| `dy_to_ll_ana.py` | DY→ll analysis pipeline |
| `qcd_ana.py` | QCD (MinBias) analysis pipeline |
| `an_specific_utilities.py` | `SampleRDFManager`, gen-match, angdiff |
| `rdf_generic.py` | `define_newcollection`, `add_hists_*`, snapshots |
| `varmetadata.py` | Histogram binning + labels |
| `define_cpp_utils.py` | Orchestrator for ALL ROOT C++ declarations |
| `pypkg/calc_puppi_iso.py` | Puppi isolation recalculation |
| `pypkg/post_analysis_persample.py` | Snapshot exports for ROC |
| `pypkg/plotBeautifier.py` | Plot formatting |
| `roc_and_rate/make_roc_and_rate.py` | ROC curve generation |
| `roc_and_rate/my_py_ai_utils.py` | ROC helper functions |

## Common Pitfalls

- **`check_env.py` is in `other_helper_files/`** — `python check_env.py` will fail
- **No CLI args** for `analyser.py` or `plotter.py` — edit config YAML files
- **C++ load order matters** — see C++ Functions section
- **`dy_to_ll_ana.py` dead code** after line 111 — references undefined `dfgenEB`
- **`qcd_ana.py` line 34** — stale TODO, not functional
- **`plot_config.yaml` sample names** must match analyser sample keys (not file paths)
- **`.gitignore`** excludes `*.root`, `*.pkl`, `*.png`, `*.DS_Store`, `*.parquet`, `hists_*/`, `puppiIsoData/`, `nb_plots/*`, `env_report.json`
- **`README.md` is stale** — claims CLI args, wrong `check_env.py` location, wrong sample names. Trust this file.
- **`auto_singlehist_plotter()`** in `plotter.py` uncommented creates a plot for every histogram — very I/O heavy
