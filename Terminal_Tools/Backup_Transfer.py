#!/usr/bin/env python
"""
Created on Sat Dec 10 06:23:06 2022

@author: coopy
"""

#This script will copy the tinyouts from a backup directory into a calculation
#directory in gc_manager. Specify the path of the backup folder you want to use

import os
import argparse
from os.path import exists as ope
from  os.path import join as opj
from os.path import isdir as opi
import subprocess
from pathlib import Path 

parser = argparse.ArgumentParser()

parser.add_argument('-p', help='Give path to backup folder',
                    type=str, required=False)

args = parser.parse_args()
path = args.p
managed_path = os.getcwd()

if ope('manager_control.txt'):   
    managed = True #in a gc_manager directory
else:
    raise Exception('Not in a gc_manager directory')
    
if path is None:
    path = opj(os.getcwd(),'backup')

#this function reads the convergence file and returns the largest step number in the
#file. In most cases it will be either 2 or 3
def read_convergence(path):
    with open(opj(path,'convergence'),'r') as f:
        convergence = f.read()
    for line in convergence.split('\n'):
        if line.startswith('step') or line.startswith('Step'):
            step_num = line.split(' ')[1] #store step number in step_num variable
    return step_num

def read_optlog(path):
    with open(opj(path,'opt.log'),'r') as f:
        opt_log = f.read()
    for line in opt_log.split('\n'): #splits at line break
        if line.startswith('Running'):
            step_num = line.split(' ')[-1]
    return step_num

def convergence_compare(path):
    '''
This function takes the path as input and reads the opt.log and convergence files
in that path. If the final steps in each final are the same, the calculation
is converged and it will return a True boolean value. If either file is missing
it will return false.
    '''
    try:
        convergence_step = read_convergence(path)
        opt_log_step = read_optlog(path)
        if convergence_step == opt_log_step or int(opt_log_step) == int(convergence_step) + 1:
            converged = True
        else:
            converged = False
    except: 
        converged = False
    

    return converged


#This part is for handling backup directories that are organized under a directory
#name like 'backup_2'. As long as the directory starts with 'backup', the script
#will recognize that as the backup directory and set it as the name for the 
#delimitter used in the for loop below
backup_path = Path(path).parts
for part in backup_path:
    if part.startswith('backup'):
        delimiter = part
   
    
#loop through directories in the backup folder. For directories that contain a
#tinyout file and that are converged, add their paths to the converged_paths list. 
#Also perform string manipulation on the path to convert the backup folder path
#to the current gc_manager directory path and add it to managed_paths
converged_paths = []
managed_paths = []
for root,folder,file in os.walk(path):
    # if root == '/projects/cote3804/BEAST-data/cooper/backup_test/backup_2/adsorbed/Fe_110/N2/0.00V/01':
        # print(read_convergence(root), read_optlog(root))
        # print(convergence_compare(root))
    if convergence_compare(root) is True or 'molecules' in str(root): # and ope(opj(root,'tinyout')) is True:
        #calc_path converts the backup path to the actively managed gc_manager
        # calcs directory path
        # print(root)
        calc_path = str(managed_path) + '/calcs' + str(root).split(delimiter)[-1]
        converged_paths.append(root)
        managed_paths.append(calc_path)
        

for backup_path, managed_path in zip(converged_paths, managed_paths):
    files_to_transfer = ['tinyout', 'inputs', 'POSCAR', 'CONTCAR', 'opt.log', 'Ecomponents', 'convergence']
    for file in files_to_transfer:
        converged_file = opj(backup_path,file)
        if opi(managed_path):
            cmd = f'cp {converged_file} {managed_path}'
            subprocess.run(cmd, shell=True)
            # if file == 'convergence':
            print(f'copied {file} from {converged_file} to {managed_path}')
        else:
            continue


        
