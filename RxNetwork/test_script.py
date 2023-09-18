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
sys.path.append(script_path)

from Analyzer import Analyzer

bulks = ["Ni2P", "MoS2"]

analysis = Analyzer(data_path=linux_path,
                    filename="all_data.json", bulks=bulks)
analysis.do_analysis()
data = analysis.get_data()
# print(data)
print(analysis.get_FED_energies(-0.5))
# analysis.visualize_contcar(0, 'Ni2P_001', 'N2H')