import argparse 
import pathlib
import shutil
import yaml

import ROOT

from define_cpp_utils import define_cpp_utils
import dy_to_ll_ana
import qcd_ana
import utilities as utilities

@utilities.time_eval
def analyser():

    ROOT.EnableImplicitMT()
    define_cpp_utils()

    # Get arguments
    parser = argparse.ArgumentParser(description="Analyse the data")
    parser.add_argument("config", help="Path to the configuration file")
    args = parser.parse_args()

    # Open the analysis config file
    with open(args.config, "r") as config_file:
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
        except:
            print(f"Could not load the data for sample {s_name}")
            continue
        print(f"Processing {df.Count().GetValue()} events in sample {s_name}")

        # Run the analysers depending on the samples
        histograms = []
        hists = []
        if s_info['type'] == 'dytoll':
            hists = dy_to_ll_ana.dy_to_ll_ana_main(df)
        elif s_info['type'] == 'qcd':
            hists = qcd_ana.qcd_ana_main(df)
        else:
            print(f"Unknown sample type {s_info['type']}")
            continue
        histograms.extend(hists)

        # Write histograms to output root file
        outfile = ROOT.TFile(f'{opts['output_dir']}/hists_{s_name}.root', 'RECREATE')
        for hist in histograms:
            hist.Write()
        outfile.Close()


# Main function
if __name__ == "__main__":
    analyser()