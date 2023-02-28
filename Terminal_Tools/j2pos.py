#!/usr/bin/env python3
"""
Created on Sun Feb  5 10:01:36 2023

@author: coopy
"""

#script to convert .ionpos file in JDFTx to CONTCAR

import subprocess
from ase.io import read, write
import os

path = os.getcwd()
files = os.listdir()
ionpos = False

# this loop checks all the files in the file directory to see if any are named
# .ionpos. If it finds one that doesn't begin with 'init', it stores the name
# that comes before .ionpos as file_name.
for file in files:
    if '.ionpos' in file and 'init' not in file:
        file_name = file.split('.ionpos')[0]
        ionpos = True

if ionpos is False:
    raise Exception('ionpos file not in current directory')
        
# Run the shell command 'createXSF out {file_name} to create an xsf file from the
# converged ionpos and lattice files.
subprocess.run(f'createXSF out {file_name}.xsf', shell=True)

# Use ASE to read in the structure from the .xsf file and write it to a CONTCAR
# file.
structure = read(f'{file_name}.xsf')
write('CONTCAR', structure, format='vasp')
