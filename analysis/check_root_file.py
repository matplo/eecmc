#!/usr/bin/env python3

import ROOT
import sys
import os

def main():
    if len(sys.argv) < 2:
        print("Usage: {} <file>".format(sys.argv[0]))
        sys.exit(1)
    if not os.path.exists(sys.argv[1]):
        print("File does not exist {}".format(sys.argv[1]))
        sys.exit(1)
    fname = sys.argv[1]
    f = ROOT.TFile(fname)
    if f.IsZombie():
        print("Error opening file {}".format(fname))
        sys.exit(1)
    if f.TestBit(ROOT.TFile.kRecovered):
        print("File {} was not correctly written".format(fname))
        f.Close()
        sys.exit(1)
    keys = f.GetListOfKeys()
    if keys.IsEmpty():
        print("File {} has no keys".format(fname))
        f.Close()
        sys.exit(1)
    else:
        print("File {} was correctly written and contains keys".format(fname))
    f.Close()
  
if __name__ == "__main__":
    main()
    sys.exit(0)