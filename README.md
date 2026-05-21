# Phase2 L1T Electron Isolation Analysis

## Quick Start

**Run analysis:**
```bash
python analyser.py analysis_config.yaml
```

**Check environment:**
```bash
python check_env.py
```

**Generate plots:**
```bash
python plotter.py
```

## Project Overview

This is a HEP analysis framework using ROOT's RDataFrame to process flat Ntuples from Phase2 L1T studies. Primary analysis: electron isolation studies in DY→ll and QCD events with PU200.

## Key Architecture

**Entry point:** `analyser.py:18` → `analyser()`
- Loads config from `analysis_config.yaml:4-44`
- Creates `SampleRDFManager` per sample (`an_specific_utilities.py:28`)
- Dispatches to `dy_to_ll_ana.py:12` or `qcd_ana.py:13` based on sample type
- Outputs histograms to `./OutHistoFiles/hists_{sample}.root`

**Sample types:**
- `dytoll`: DYToLL_M50_PU0, DYToLL_M50_PU200
- `qcd`: MinBias_PU200

**Main collections:**
- `TkEleL2`: Reconstructed TkElectrons
- `GenEl`: Generator electrons
- `L1PuppiCands`: L1 Puppicandidates by pdgId

## Data Dependencies

**Required dependencies:**
- Python ≥3.9 (tested: 3.13, 3.14)
- ROOT ≥6.36 (tested: 6.38)
- numpy 2.4+, yaml 6.0+, IPython 9.7+

**Data download:**
```bash
scp -r mercury13:/opt/ppd/cms/users/asahasra/Phase2EGPuppiIso_151Xv1_Ntuples ./
```

**Execution order:**
1. `python check_env.py` (creates `env_report.json`)
2. `python analyser.py analysis_config.yaml` → outputs ROOT files
3. `python plotter.py plot_config.yaml` → generates PNG plots

## Output Structure

```
./OutHistoFiles/
├── hists_DY_noPU.root
├── hists_DY_PU200.root
└── hists_MinBias.root
```

Each sample's plots use normalized histogram names: `{rdf_filter}_{collection_var}`

## Common Pitfalls

- Output directory is **deleted and recreated** on each run
- `env_report.json` is written during environment validation
- Custom isolation variables use `dRmin{X_Y}` pattern (dots replaced with underscore)
- C++ functions must be declared before use; `analyser.py:21` runs first

## Top-Level Files

| File | Description |
|------|-------------|
| `analyser.py` | Entry point: loads config, creates SampleRDFManager, dispatches to dy_to_ll_ana or qcd_ana |
| `analysis_config.yaml` | Sample definitions with input patterns and types (dytoll/qcd) |
| `plot_config.yaml` | Plot configurations with histogram collections and normalization schemes |
| `check_env.py` | Environment validation script; writes env_report.json |
| `define_cpp_utils.py` | C++ utility functions (getminangs, getmatchedidxs, calcisoannularcone, quantise_relisoval, checksorting) |
| `dy_to_ll_ana.py` | DY→ll analysis: gen-matching, η region splitting, puppi isolation recalculation |
| `qcd_ana.py` | QCD analysis: TkEle selection, η region splitting, puppi isolation recalculation |
| `plotter.py` | Plot generator: reads ROOT files, applies normalization, calls plotBeautifier |
| `plotBeautifier.py` | Plot formatting utilities: log-axis rebinning, auto x-range, axis labels |
| `calc_puppi_iso.py` | Calculates puppi isolation in annular cones with configurable dRmin/dRmax |
| `an_specific_utilities.py` | SampleRDFManager class and utility functions (gen-matching, angdiff histograms) |
| `rdf_generic.py` | Generic RDataFrame utilities: define_newcollection, add_hists_* |
| `varmetadata.py` | Histogram metadata: binning and titles for common variables (pt, eta, phi, iso, etc.) |
| `my_py_generic_utils.py` | Generic utilities: recreate_dir, time_eval decorator, create_rdf_checkpint |
| `my_py_analysis_utils.py` | Analysis-specific utilities (if any) |
| `my_py_ai_utils.py` | AI/ML utilities (if any) |
| `post_analysis_persample.py` | Post-processing: saves RDataFrame snapshots to pickle for ROC curves |
| `make_roc_and_rate.py` | ROC curve and rate calculation utilities |
| `other_helper_files/make_env.sh` | CVMFS environment setup script (LCG_108a) |
| `other_helper_files/download_data.sh` | Data download command from mercury13 |
