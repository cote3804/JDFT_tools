# -*- coding: utf-8 -*-
"""
Created on Sat Sep 16 08:42:15 2023

@author: coopy
"""
import sys
import os
from matplotlib import pyplot as plt
import matplotlib as mpl
script_path = "C:\\Users\\coopy\\OneDrive - UCB-O365\\github\\JDFT_tools\\RxNetwork"
linux_path = "/home/cooper/Research/Binaries"
windows_path = "C:\\Users\\coopy\\OneDrive - UCB-O365\\Research\\Binaries\\all_data"
downloads_path = "C:\\Users\\coopy\\Downloads"
sys.path.append(script_path)

from Analyzer import Analyzer

# bulks = ["ZrC"]
# bulks = ["HfP"]
bulks = ["HfP"]

analysis = Analyzer(data_path=linux_path,
                    filename="fake_all_data.json", bulks=bulks)
analysis.do_analysis()
# print(data)
# spans = analysis.get_spans()
print(analysis.get_FED_energies(0))
binding = analysis.get_binding_energies(0)
analysis.get_spans("NRR")
FED = analysis.get_FED_energy("HfP_110", "0.00V")
print("FED energies", FED)
# analysis.visualize_contcar(0, 'Ni2P_001', 'N2H')
# analysis.visualize_contcar(0, 'HfP_110', 'N2H')
fig1, ax1 = analysis.plot_FED("NRR", "HfP_110", 0)
fig1, ax1 = analysis.plot_FED("NRR", "HfP_110", -0.5, color="#003EF5", graph_objects=(fig1,ax1))
plt.show()
# plt.savefig(os.path.join(linux_path, "HfP_110_FED.png"))
# fig2, ax2 = analysis.plot_FED("NRR", "Co2Si_100", 'No_bias')




