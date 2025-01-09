import argparse 
import yaml

import ROOT

from dy_to_ll_ana import dy_to_ll_ana
import utilities as utilities

@utilities.time_eval
def analyser():

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

        df = dy_to_ll_ana(df)


if __name__ == "__main__":
    analyser()