#!/usr/bin/env python

import os
import sys
import subprocess
from os.path import exists as ope
from os.path import join as opj


def convergence_change(conv_path):
    #load in convergence file
    with open(opj(conv_path,'convergence')) as f:
        data = f.read()
    
    step_dict = {}
    #save convergence information in dictionary.
    #keys correspond to step numbers and each key contains a list of the convergence criteria
    for line in data.split('\n'):
        if line.startswith('step'): #make new empty list at step number read from file
            step_number = int(line.split(' ')[1])
            step_dict[step_number] = []
        else:
            step_dict[step_number].append(line) #append convergence parameters to list stored in associated key
    
    skip = False #False if the rough convergence step is not skipped
    if step_dict[1][0] == 'kpoint-folding 1 1 1':
        skip = True #set true if the rough convergence step is skipped
        step_dict.pop(1) #remove first step from dictionary  

    write_lines = [] #blank list for writing new convergence document
    for step_number in step_dict.keys():
        if skip: #if the first step is skipped, must subtract 1 from step_number
            write_lines.append(f'step {step_number-1}') #write string with step followed by step number
        else:
            write_lines.append(f'step {step_number}')
        for line in step_dict[step_number]:
            write_lines.append(line) #append all lines within dictionary key
    
    with open(opj(conv_path,'convergence'),'w') as r:
        for line in write_lines:
            r.write(line)
            r.write('\n')

            
def file_clean(path):
    file_keep_list = ['convergence', 'POSCAR','inputs']
    for file in os.listdir(path):
        if file in file_keep_list: #don't delete files in the list above
            pass
        else:
            os.remove(opj(path,file)) #delete these files for fresh calculation re-run


if __name__ == '__main__':

    path = os.getcwd()

    #check to see if in gc_manager directory or results sub-directory withing gc_manager
    if ope('manager_control.txt'):
        pass #if in manager directory, initial path variable is correct
    elif os.path.split(path)[1] == 'results':
        path = os.path.split(path)[0] # go to manager directory
    else:
        raise Exception('Not in gc_manager directory')

    results_path = opj(path,'results')
    with open(opj(results_path,'failed.txt'),'r') as f:
        results = f.read()

    paths_list = []
    for line in results.split('\n'):
        paths_list.append(opj(path,line))

    #these lines execute all of the file manipulation 
    print('fixing the following failed calculations:')            
    for path in paths_list:
        convergence_change(path)
        file_clean(path)
        print(path)
