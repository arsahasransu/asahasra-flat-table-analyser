# import argparse
import os
import pathlib
import shutil
import yaml

import ROOT

from define_cpp_utils import define_cpp_utils
import dy_to_ll_ana
import my_py_generic_utils as ut
import qcd_ana


@ut.time_eval
def writehists_to_file(ofname, histograms):
    outfile = ROOT.TFile(ofname, 'RECREATE')
    for hist in histograms:
        hist.Write()
    outfile.Close()


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
        try:
            df = ROOT.RDataFrame(opts['tree_name'], opts['input_dir_prefix']+'/'+s_info['input_file_pattern'])
        except FileNotFoundError as err:
            print(f"Could not load the data for sample {s_name}: {err}")
            continue
        print("\n")
        print(f"PROCESSING {df.Count().GetValue()} EVENTS IN SAMPLE {s_name.upper()}")

        # Run the analysers depending on the samples
        histograms = []
        if s_info['type'] == 'dytoll':
            histograms = dy_to_ll_ana.dy_to_ll_ana_main(df)
        elif s_info['type'] == 'qcd':
            histograms = qcd_ana.qcd_ana_main(df)
        else:
            print(f"Unknown sample type {s_info['type']}")
            continue

        # Write histograms to output root file
        writehists_to_file(f'{opts['output_dir']}/hists_{s_name}.root', histograms)


# Main function
if __name__ == "__main__":
    analyser()
