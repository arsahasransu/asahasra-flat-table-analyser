# import argparse
import os
import pathlib
import shutil
import yaml

import ROOT

import an_specific_utilities as anautil
from define_cpp_utils import define_cpp_utils
import dy_to_ll_ana
import my_py_generic_utils as ut
from post_analysis_persample import post_analysis_persample
import qcd_ana


@ut.time_eval
def analyser():

    ROOT.EnableImplicitMT()
    define_cpp_utils()

    # Get arguments
    # parser = argparse.ArgumentParser(description="Analyse the data")
    # parser.add_argument("config", help="Path to the configuration file")
    # args = parser.parse_args()
    # cfile = args.config
    cfile = "analysis_config.yaml"

    # Open the analysis config file
    with open(cfile, "r") as config_file:
        config = yaml.safe_load(config_file)

    opts = config['general']
    samples = config['samples']

    # Remake the histogram root output directory
    outdir = pathlib.Path(f'{opts['output_dir']}')
    if outdir.exists() and outdir.is_dir():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Loop on simulation samples
    for s_name, s_info in samples.items():

        anamanager = anautil.SampleRDFManager(opts['input_dir_prefix']+'/'+s_info['input_file_pattern'], s_name, tree_name=opts['tree_name'])

        # Run the analysers depending on the samples
        histograms = []
        if s_info['type'] == 'dytoll':
            dy_to_ll_ana.dy_to_ll_ana_main(anamanager)
            histograms = anamanager.get_histograms()
        elif s_info['type'] == 'qcd':
            qcd_ana.qcd_ana_main(anamanager)
            histograms = anamanager.get_histograms()
        else:
            print(f"Unknown sample type {s_info['type']}")
            continue

        # Write histograms to output root file
        anautil.writehists_to_file(f'{opts['output_dir']}/hists_{s_name}.root', histograms)

        # Run post-analysis special code
        post_analysis_persample(anamanager)


# Main function
if __name__ == "__main__":
    analyser()
