import argparse 
import yaml

import ROOT

from define_cpp_utils import define_cpp_utils
from dy_to_ll_ana import dy_to_ll_ana_main
import utilities as utilities

@utilities.time_eval
def analyser():

    ROOT.EnableImplicitMT()
    define_cpp_utils()

    parser = argparse.ArgumentParser(description="Analyse the data")
    parser.add_argument("config", help="Path to the configuration file")
    args = parser.parse_args()

    with open(args.config, "r") as config_file:
        config = yaml.safe_load(config_file)

    opts = config['general']
    samples = config['samples']

    for s_name, s_info in samples.items():
        try:
            df = ROOT.RDataFrame(opts['tree_name'], opts['input_dir_prefix']+'/'+s_info['input_file_pattern'])
        except:
            print(f"Could not load the data for sample {s_name}")
            continue
        print(f"Processing {df.Count().GetValue()} events in sample {s_name}")

        histograms = []
        df,more_hists = dy_to_ll_ana_main(df)
        histograms.extend(more_hists)

        outfile = ROOT.TFile(f'{opts['output_dir']}/hists_{s_name}.root', 'RECREATE')
        for hist in histograms:
            hist.Write()
        outfile.Close()

if __name__ == "__main__":
    analyser()