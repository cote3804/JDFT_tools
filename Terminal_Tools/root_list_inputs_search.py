#!/usr/bin/env python3
"""
Created on Wed Jun  7 12:20:32 2023

@author: coopy
"""

import os
import subprocess

root = os.getcwd()

with open("tmp_parallel/root_list.txt") as f:
    root_list = f.read()



for path in root_list.split('\n'):
    target_path = os.path.join(root,path)
    print(target_path)
    # subprocess.run("cd {target_path}", shell=True)
    subprocess.run(f"bands2inputs.py -r {target_path}", shell=True)

    