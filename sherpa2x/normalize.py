#!/usr/bin/env python

from heppyy.util.logger import Logger
log = Logger()
log.set_level('INFO')

import sys
import os

from heec import EEC2file

for arg in sys.argv:
    if '.root' not in arg:
        continue
    if os.path.isfile(arg):
        f = EEC2file(arg)
        print(f.norms())
        f.ls()
        f.normalize()
        f.ls()
        f.SaveAs(arg.replace('.root', '_norm.root'))

def test():
	f = EEC2file('test_eec2.root')

	print(f.norms())

	f.ls()
	f.normalize()
	f.ls()

	f.SaveAs('test_eec2_norm.root')