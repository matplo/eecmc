#!/usr/bin/env python3

import yaml
import ROOT
import numpy as np
import array
import sys

def logbins(xmin, xmax, nbins):
    xmin = max(xmin, 1e-10)
    xmax = max(xmax, 1e-10)
    lspace = np.logspace(np.log10(xmin), np.log10(xmax), nbins+1)
    arr = array.array('f', lspace)
    return arr

def logbins_from_config(hist_config):
    return logbins(float(hist_config['xrange'][0]), float(hist_config['xrange'][1]), int(hist_config['xnbins']))

# Load the YAML file
yaml_file = 'tdraw_config.yaml'
if len(sys.argv) > 1:
    yaml_file = sys.argv[1]
with open(yaml_file, 'r') as stream:
    try:
        config = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

# Loop over all sections in the config
for section in config:
    # Get the output file and mode from the first item in the current section
    output_file_name = config[section][0]['output']
    output_mode = config[section][0]['mode']

    # Open the output file
    output_file = ROOT.TFile(output_file_name, output_mode)

    # Loop over the rest of the items in the current section
    for item in config[section][1:]:
        # Each item is a dictionary where the key is the histogram name and the value is the histogram configuration
        for hist_name, hist_config in item.items():
            hist_name = hist_name + '_' + section
            # Open the ROOT file and get the tree
            root_file = ROOT.TFile.Open(hist_config['file'])
            tree = root_file.Get(hist_config['tree'])
            print('tree:', type(tree), 'hist_config:', hist_config)
            if isinstance(tree, ROOT.TNtuple) or isinstance(tree, ROOT.TTree) or isinstance(tree, ROOT.TChain):
                print(f'type check for {tree} passed')
                pass
            else:
                print('Could not find tree:', hist_config['tree'])
                continue

            # Create the histogram
            output_file.cd()
            if hist_config['xbins'] == 'log':
                bins = logbins_from_config(hist_config)
                hist = ROOT.TH1F(hist_name, hist_name, len(bins)-1, bins)
            else:
                hist = ROOT.TH1F(hist_name, hist_name, hist_config['xnbins'], *hist_config['xrange'])

            # Draw the variable with the specified condition
            draw_string = "{}>>{}".format(hist_config['var'], hist_name)
            tree.Draw(draw_string, hist_config['cond'])

            # If 'scale' is specified, scale the histogram
            if 'scale' in hist_config:
                scale = eval(hist_config['scale'], {'__builtins__': None}, {'h_jet_pt': ROOT.gDirectory.Get(hist_name)})
                hist.Scale(scale)

            if 'bw' in hist_config:
                hist.Scale(1., 'width')


            print('draw_string:', draw_string, 'entries', hist.GetEntries())
            # Save the histogram to the output file
            hist.Write()

    # Close the output file
    output_file.Close()