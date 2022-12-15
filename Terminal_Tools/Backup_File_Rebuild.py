#! /usr/bin/env python3
"""
Created on Mon Dec 12 11:20:01 2022

@author: coopy
"""

import json
from os.path import join as opj
import os
from pathlib import Path
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f', help='JSON file name to rebuild from', default='combined_N2R_all_data.json')
args = parser.parse_args()
file=args.f

path = os.getcwd()
file = opj(path,'combined_N2R_all_data.json')

with open(file,'r') as f:
    data = json.load(f)

def build_optlog(optlog_dict):
    text = []
    current_step = 0
    for step_num in optlog_dict.keys():
        text.append(f'Running Convergence Step: {step_num}\n')
        text.append(f'        Step     Time         Energy       fmax\n')
        text.append(f'*Force-consistent energies used in optimization.\n')
        for i, step_dict in enumerate(optlog_dict[step_num]):
            text.append(f'{step_dict["opt_method"]}:    {i:>3} 00:00:00   {step_dict["energy"]:7.6f}* {step_dict["force"]:>9.4f}\n')
        text.append('\n')
    return text

def write_optlog(optlog_text, path_string):
    try:
        with open(opj(path_string,'opt.log'), 'w') as w:
            for line in optlog_text:
                w.write(line)
    except:
        pass
        # print(f'path {path_string} not found')
        
def build_convergence(convergence_dict):
    text = []
    for step in convergence_dict.keys():
        text.append(f'step {step}\n')
        for parameter, parameter_value in convergence_dict[step].items():
            text.append(f'{parameter} {parameter_value}\n')
        text.append('\n')
    return text

def write_convergence(text, path_string):
    try:
        with open(opj(path_string,'convergence'), 'w') as w:
            for line in text:
                w.write(line)
    except:
        pass
        # print(f'path {path_string} not found')
            
def is_bias(string): #function to check if input string is a voltage string
    if string.endswith('V') and '/' not in string:
        return True
    else:
        return False
    
converged_paths= []
for material in data.keys():
    for calc_type in data[material]:
        
        #Adsorbed
        if calc_type == 'adsorbed':
            for molecule in data[material][calc_type].keys():
                for bias in data[material][calc_type][molecule].keys():
                    for site in data[material][calc_type][molecule][bias].keys():
                        if data[material][calc_type][molecule][bias][site]['converged']:
                            path_string = calc_type + '/' + material + '/' + molecule + '/' + bias + '/' + site + '/'
                            path_string = opj(path,path_string)
                            opt_log = build_optlog(data[material][calc_type][molecule][bias][site]['opt'])
                            write_optlog(opt_log, path_string)
                            convergence = build_convergence(data[material][calc_type][molecule][bias][site]['convergence_file'])
                            write_convergence(convergence, path_string)
        #Surfaces
        elif calc_type == 'surf':
            for bias in data[material][calc_type].keys():
                if data[material][calc_type][bias]['converged']:
                    path_string = 'surfs' + '/' + material +  '/' + bias + '/'
                    path_string = opj(path,path_string)
                    opt_log = build_optlog(data[material][calc_type][bias]['opt'])
                    convergence = build_convergence(data[material][calc_type][bias]['convergence_file'])
                    write_optlog(opt_log, path_string)
                    write_convergence(convergence, path_string)
        #Bulks
        elif calc_type == 'bulk':
            if data[material][calc_type]['converged']:
                path_string = 'bulks' + '/' + material +  '/'
                path_string = opj(path,path_string)
                opt_log = build_optlog(data[material][calc_type]['opt'])
                write_optlog(opt_log, path_string)
                convergence = build_convergence(data[material][calc_type]['convergence_file'])
                write_convergence(convergence, path_string)
        
        #Molecules
        elif is_bias(calc_type): #If the items are voltage strings, then the material is a molecule
            bias = calc_type
            if data[material][bias]['converged']:
                path_string = 'molecules' + '/' + material +  '/'
                path_string = opj(path,path_string)
                opt_log = build_optlog(data[material][bias]['opt'])
                write_optlog(opt_log, path_string)
                convergence = build_convergence(data[material][bias]['convergence_file'])
                write_convergence(convergence, path_string)
            
            
            