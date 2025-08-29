# Fast Run-through of Repo

### Dependencies
The repo depends on the following packages, all of which are part of the (LCG package)[https://lcginfo.cern.ch] releases.
 - `python` - Atleast v3.9
 - `pyROOT` - ROOT with python bindings

The repo has been test run successfully with the following releases.
|Package|Version|
|:------|:-----:|
|`python`|3.13|
|`ROOT`|6.36|

### To run the analysis

```shell
# Needs pyROOT, yaml
python analyser.py analysis_config.yaml
```
# Introduction

Note on the files in this directory.

 - `analyser.py` - The starting point, `main()`, of the analysis. Uses the configuration file `analysis_config.yaml`.
 - `analysis_config.yaml` - Analysis configurations for all the data/MC datasets.
 - `code_quality_control_checks.sh` - 
 - `define_cpp_utils.py` - 
 - `download_data.sh` - 
 - `dy_to_ll_ana.py` - 
 - `make_env.sh` - 
 - `plotBeautifier.py` - 
 - `plotter.py` - 
 - `utilities.py` - 

