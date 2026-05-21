# AGENTS.md

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
- `TkEleL2`: Reconstructed TkElectrons (defined in `an_specific_utilities.py:14`)
- `GenEl`: Generator electrons (`an_specific_utilities.py:15`)
- `L1PuppiCands`: L1 Puppicandidates by pdgId (`an_specific_utilities.py:16`)

## Critical Paths

**DY analysis flow** (`dy_to_ll_ana.py`):
1. Define collections: `GenEl_prompt==2 && |eta|<2.5` → `GenEl_DYP`
2. Filter: `genDYP_n>0 && TkEleL2_n>0` → `genDYP`
3. Split η regions: EB (|η|≤1.5), EE (1.5<|η|≤2.5)
4. Gen-match with dR cuts: EB=0.04, EE=0.05 (`dy_to_ll_ana.py:56`)
5. Recalculate puppi isolation for gen-matched electrons (`calc_puppi_iso.py:6`)
6. Output: `DYPM{EB,EE}` dataframes with iso variables

**QCD analysis flow** (`qcd_ana.py`):
1. FilterTkEleL2 pT>0 → `TkEleL2_Pt5`
2. Split η regions: EB, EE
3. Recalculate puppi isolation (`qcd_ana.py:47`)
4. Output: `{EB,EE}` dataframes

## C++ Utilities (define_cpp_utils.py)

Auto-loaded via `ROOT.gInterpreter.Declare()` in `analyser.py:21`:
- `getminangs()`: min dη, dφ, dR between two collections
- `getmatchedidxs()`: gen-matching with dR cut
- `calcisoannularcone()`: puppi isolation in annular cones
- `quantise_relisoval()`: custom浮点 quantization (1-8-8 format)
- `checksorting()`: verify descending pT sort

## Configuration

**Input:** `analysis_config.yaml:1-44`
- `input_dir_prefix`: base path for input files (default: `./Phase2EGPuppiIso_151Xv1_Ntuples`)
- `tree_name`: `Events`
- `output_dir`: `./OutHistoFiles`

**Plot config:** `plot_config.yaml:1-201`
- Predefined plot groups in sections (most current are enabled, others commented)
- Key plot: `compare_puppi_iso_consistency_check` (lines 74-84)
- Normalization schemes: `hist_integral`, `default_no_norm`, `summed_components`

## Data Dependencies

**Root level commands to run:**
1. `python check_env.py` (creates `env_report.json`)
2. `python analyser.py analysis_config.yaml` → outputs ROOT files
3. `python plotter.py plot_config.yaml` → generates PNG plots

**Environment:**
- Python ≥3.9 (tested: 3.13, 3.14)
- ROOT ≥6.36 (tested: 6.38)
- Required Python deps: numpy 2.4+, yaml 6.0+, IPython 9.7+

**Data download:**
```bash
scp -r mercury13:/opt/ppd/cms/users/asahasra/Phase2EGPuppiIso_151Xv1_Ntuples ./
```

## Output Structure

```
./OutHistoFiles/
├── hists_DY_noPU.root
├── hists_DY_PU200.root
└── hists_MinBias.root
```

Each sample's plots use normalized histogram names: `{rdf_filter}_{collection_var}`

## Special Conventions

1. **Analysis steps** are marked in comments (`STEP_X_X_X:`) - uncomment to enable
2. Plot config uses composite collection syntax: `col[ele1, ele2, ...]` (see `plotter.py:84-108`)
3. dRmin variants: `reisotot_dRmin0_03` (from `calc_puppi_iso.py:7` default drmin=0.03)
4. Gen-matching: stores `*_recoidx` and `*_genidx` vectors indexing matched pairs
5. η region split uses `define_newcollection()` (`rdf_generic.py:13`)

## Common Pitfalls

- Output directory is **deleted and recreated** on each run (`analyser.py:37-41`)
- `env_report.json` is written by `check_env.py:157` during environment validation
- Custom isolation variables use `dRmin{X_Y}` pattern (dots replaced with underscore)
- C++ functions must be declared before use; `analyser.py:21` must run first
