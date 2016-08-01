#!/usr/bin/env python

## test why pybel does not work in this environment

import os
import subprocess

py = "import sys; print '\\n'.join(sys.path)"
out = subprocess.Popen('python -c "{}"'.format(py), shell=True, stdout=subprocess.PIPE)
print out.communicate()[0]
exit()


out = subprocess.Popen('printenv', stdout=subprocess.PIPE)
print out.communicate()[0]
exit()
