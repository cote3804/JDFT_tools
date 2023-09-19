# -*- coding: utf-8 -*-
"""
Created on Sat Sep 16 08:42:15 2023

@author: coopy
"""
import sys
import os
script_path = "C:\\Users\\coopy\\OneDrive - UCB-O365\\github\\JDFT_tools\\RxNetwork"
linux_path = "/home/cooper/Research/Binaries"
windows_path = "C:\\Users\\coopy\\OneDrive - UCB-O365\\Research\\Binaries\\all_data"
downloads_path = "C:\\Users\\coopy\\Downloads"
sys.path.append(script_path)

from Analyzer import Analyzer

bulks = ["HfP"]

analysis = Analyzer(data_path=downloads_path,
                    filename="fake_all_data.json", bulks=bulks)
analysis.do_analysis()
data = analysis.get_data()
print(data)
# spans = analysis.get_spans()
# print(data)
print(analysis.get_FED_energies(0))
binding = analysis.get_binding_energies(0)
print(binding)
# analysis.visualize_contcar(0, 'Ni2P_001', 'N2H')