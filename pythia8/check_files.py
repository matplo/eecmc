#!/usr/bin/env python3

import sys
import ROOT
import tqdm

# Suppress the warnings
ROOT.gErrorIgnoreLevel = ROOT.kError

def is_root_file_ok(filename, silent=False):
	try:
			root_file = ROOT.TFile(filename)
	except OSError as e:
		if silent == False:
			print (e)
		return False
	if root_file.IsZombie():
		if silent == False:
			print(f"Error opening file: {filename}")
		return False	
	else:
		if silent == False:
			print(f"File is okay: {filename}")
		return True
			
import os
import fnmatch

def find_root_files(directory):
		root_files = []
		for root, dir, files in os.walk(directory):
				for file in files:
						if fnmatch.fnmatch(file, '*.root'):
								root_files.append(os.path.join(root, file))
		return root_files

# Usage
files_to_remove = []
if len(sys.argv) > 1:
	root_files = find_root_files(sys.argv[1])
	for file in tqdm.tqdm(root_files):
		if not is_root_file_ok(file, silent=True):
			files_to_remove.append(file)
if len(sys.argv) > 2:
  if sys.argv[2] == 'remove':	
    for file in tqdm.tqdm(files_to_remove):
      os.remove(file)
