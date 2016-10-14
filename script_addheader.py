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
## license: GNU Public License v3.0 (LICENSE_GPLv3.txt)
##
################################################################################
'''

# def addheader(filename, header=header):
#     filename2 = filename+'~'
#     with open(filename, 'r') as finp, open (filename2, 'w') as fout:
#         fout.write(header.format(scriptname=filename))
#         fout.write('\n')
#         for line in finp:
#             fout.write(line)
#     os.rename(filename2, filename)

# addheader('test_file.R')

# for rfile in glob('*.R'):
#     if 'test_' == rfile[:5]:
#         continue
#     print rfile
#     addheader(rfile)



def modifyheader(filename):
    filename2 = filename+'~'
    with open(filename, 'r') as finp, open (filename2, 'w') as fout:
        for line in finp:
            if '## see LICENSE_GPLv3.txt' in line:
                fout.write('## license: GNU Public License v3.0 (LICENSE_GPLv3.txt)'+'\n')
                break
            else:
                fout.write(line)
        for line in finp:
            fout.write(line)
    os.rename(filename2, filename)

# modifyheader('test_file.R')

for rfile in glob('*.R'):
    if 'test_' == rfile[:5]:
        continue
    print rfile
    modifyheader(rfile)
