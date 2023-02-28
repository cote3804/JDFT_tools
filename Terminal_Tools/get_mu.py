#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 07:22:19 2023

@author: coopy
"""

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('bias', type=float, help='The bias in volts')
parser.add_argument('-rhe', type=bool, help='Specify True if you want reversible hydrogen electrode. Requries specifying pH',
                    default=False)
parser.add_argument('-ph', type=float, help='pH for reversible hydrogen electrode',
                    default=0.0)
args = parser.parse_args()


import os
import sys
paths = os.environ['PATH'] #This block looks for the gc_manager path and then 
for path in paths.split(':'): #appends the path to the active python path so
    if 'manager' in path: #that the jdft_helper script functions can be used
        gc_manager_path = path
        break
sys.path.append(gc_manager_path)
from jdft_helper import helper



bias = args.bias
rhe = args.rhe
pH = args.ph

if os.path.exists('./inputs'):
    with open('./inputs', 'r') as f:
        inputs = f.read()
    # print(inputs.split('\n'))
    inputs = helper().read_inputs('.')
    fluid = inputs['fluid'].strip()
    if fluid == 'LinearPCM':  
        pcm_var = inputs['pcm-variant'].strip()
else:
    print('inputs file not found. Solvent model defaulting to CANDLE')
    fluid = 'LinearPCM'
    pcm_var = 'CANDLE'


if fluid == 'LinearPCM' and pcm_var == 'CANDLE':
    Vref = 4.66
elif fluid == 'LinearPCM' and pcm_var == 'GLSSA13':
    Vref = 4.68
elif fluid == 'NonlinearPCM' and (pcm_var == 'GLSSA13'):
    Vref = 4.62
elif fluid == 'SaLSA':
    Vref = 4.54
elif fluid == 'ClassicalDFT':
    Vref = 4.44

rhe_shift = -0.0591 * pH # defaults to no RHE shift
target_mu = -(Vref + bias + rhe_shift)/27.2114
print(f'{target_mu:.4}')