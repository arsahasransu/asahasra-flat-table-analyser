# Phase2 L1T Electron Isolation Analysis

## Quick Start

**Check environment (Required before running!):**
```bash
python other_helper_files/check_env.py
```

**Run analysis (Config is hardcoded inside script):**
```bash
python analyser.py
```

**Generate plots:**
```bash
python plotter.py
```

**Run ML Inference (Appends isolation scores to snapshots):**
```bash
python train_ml_models/inference.py
```

**Evaluate ROC Curves:**
```bash
cd roc_and_rate
python make_roc_and_rate.py
```

## Project Overview

This is a high-performance HEP analysis framework using ROOT's `RDataFrame` to process flat Ntuples from Phase2 L1T studies. 
The primary goal is electron isolation studies in DY→ll and QCD (MinBias) events with PU200, developing traditional cone-based algorithms as well as ML-based (PyTorch) Soft Isolation networks.

### Recent Updates:
- **Awkward Array Integration:** Fully replaced `numpy` object arrays with native jagged `awkward` arrays for extremely fast event-level filtering, ROC curve generation, and minimum calculations.
- **ML Evaluation Pipeline:** The `SoftIsoSumNetwork` (PyTorch) now natively reads ROOT snapshot trees via `uproot` + `awkward`, evaluates the `weighted_iso_score`, and rapidly updates the ROOT files in place.
- **Upgraded ROC Evaluation:** The evaluation script now compares traditional track/puppi isolations directly against the new ML `weighted_iso_score`.

## Execution Order & Data Flow

1. `python other_helper_files/check_env.py` — Validates environment (creates `env_report.json`). Exits non-zero on failure.
2. `python analyser.py` — Produces main output histograms in `./OutHistoFiles/hists_{sample}.root` and generates sample-level snapshots (`_snapshot.root` + `.parquet` files).
3. `python plotter.py` — Reads `./OutHistoFiles/*.root` and generates PNG plots inside `./plots/`.
4. `python train_ml_models/inference.py` — Loads `_snapshot.root` files, evaluates the trained Soft Isolation model, and inserts the `TkEleL2_*_weighted_iso_score` branches back into the ROOT trees.
5. `cd roc_and_rate && python make_roc_and_rate.py` — High-performance script reading snapshots and generating multi-threshold ROC curves.

## Key Architecture

**Entry point:** `analyser.py` → `analyser()`
- **Multi-threading:** Enables ImplicitMT RDataFrame processing automatically.
- **C++ Interop:** Loads utility C++ scripts (`define_cpp_utils.py`) dynamically into `ROOT.gInterpreter`.
- **Sample Types:** 
  - `dytoll`: e.g. `DY_noPU`, `DY_PU200`
  - `qcd`: e.g. `MinBias`
- Outputs: Everything goes to a freshly generated `./OutHistoFiles/` directory.

**Configuration:**
Neither `analyser.py` nor `plotter.py` accept CLI arguments. The paths are hardcoded to load `analysis_config.yaml` and `plot_config.yaml` directly.

**Main collections:**
- `TkEleL2`: Reconstructed TkElectrons
- `GenEl`: Generator electrons
- `L1PuppiCands`: L1 Puppi candidates by pdgId

## Common Pitfalls

- The output directories (`OutHistoFiles/` and `plots/`) are **deleted and recreated** on each respective run. Do not store manual files inside them.
- All C++ functions must be declared before any `RDataFrame` logic runs. `analyser.py` handles this.
- If editing samples, the sample names in `plot_config.yaml` must exactly match the keys defined in `analysis_config.yaml`.
- The `varmetadata.py` file defines automatic binning rules. If a variable is missing, a default binning is applied.

## Top-Level Files

| File/Folder | Description |
|------|-------------|
| `analyser.py` | Main event loop: loads config, creates SampleRDFManager, dispatches to dy_to_ll_ana or qcd_ana |
| `analysis_config.yaml` | Analyer configuration (Sample definitions, input glob patterns, types) |
| `plot_config.yaml` | Plot configurations (Histogram collections, binning overrides, normalization schemes) |
| `plotter.py` | Plot generator reading from ROOT files |
| `dy_to_ll_ana.py` | DY→ll workflow: gen-matching, η region splitting, puppi isolation |
| `qcd_ana.py` | MinBias (QCD) workflow: TkEle selection, η region splitting, puppi isolation |
| `an_specific_utilities.py` | `SampleRDFManager` class and utility functions |
| `rdf_generic.py` | Generic RDataFrame utilities (define_newcollection, add_hists_*) and `awkward`-based root loaders |
| `varmetadata.py` | Histogram metadata: binning and titles for common variables |
| `define_cpp_utils.py` | Central C++ utility load manager |
| `cpp_utils/` | C++ implementation headers loaded dynamically by ROOT |
| `pypkg/` | Additional Python modules (`calc_puppi_iso.py`, `post_analysis_persample.py`, etc.) |
| `train_ml_models/` | Training and Inference scripts for PyTorch Soft Iso Sum networks |
| `roc_and_rate/` | High-performance ROC curve generation using Awkward Arrays and Matplotlib |
| `other_helper_files/` | Helpers like `make_env.sh` (CVMFS setup) and `check_env.py` (Validation) |
