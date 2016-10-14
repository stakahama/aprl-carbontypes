#!/usr/bin/env python

import os
from glob import glob

header = '''
################################################################################
##
## {scriptname}
## Author: Satoshi Takahama (satoshi.takahama@epfl.ch)
## Oct. 2016
##
## see LICENSE_GPLv3.txt
##
################################################################################
'''

def addheader(filename, header=header):
    filename2 = filename+'~'
    with open(filename, 'r') as finp, open (filename2, 'w') as fout:
        fout.write(header.format(scriptname=filename))
        fout.write('\n')
        for line in finp:
            fout.write(line)
    os.rename(filename2, filename)

addheader('test_file.R')

for rfile in glob('*.R'):
    if 'test_' == rfile[:5]:
        continue
    print rfile
    addheader(rfile)
