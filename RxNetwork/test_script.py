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

# bulks = ["Fe2Mo6S8"]
bulks = ['RhB','RhSi','Re3B','TaN','ZrN','MoS2','TaB','CaSi2']
bulks = ['CaSi2']
# bulks = ["CoSe2"]
ben_bulks = ["Co2Si", "HfP", "MoS2", "Nb2C", "TiC", "Zr2Si", "ZrC"]
# bulks = [
#         # "RuP", 
#          "ZrC", 
#         #  "AgO", 
#         #  "WP",
#         "VO2"
#          ]
analysis = Analyzer(data_path=windows_path,
                    filename="combined_4_all_data.json", bulks=bulks)
# print(analysis.get_data("NRR")['TiC_11-1']['NNH2*']['0.00V'])
# spans = analysis.get_spans("NRR", -0.5)
# print("\n spans",spans)
# analysis.get_spans("NRR")
# test = analysis.get_FED_energies(0)
# print(test)
# FED = analysis.get_FED_energy("RhB_2-12", 0, "NRR")
# print(FED)
# charges = analysis.get_pathway_charges(0, "NRR")
# print(charges)
# binding = analysis.get_binding_energies(0, "NRR")
# data = analysis.get_data("NRR")
# print(data)

gmax, span = analysis.get_span("NRR", "CaSi2_001", -0.5)
print(gmax, span)
# band_centers = analysis.get_band_centers()
# print(band_centers)
# min_site = analysis.get_min_site("RhB_2-12", 0, "N2")
struct = analysis.get_structure("FeTe2_100", 0, "N2*", site="min")
struct.visualize()
struct = analysis.get_structure("FeTe2_100", 0, "N*", site="min")
struct.visualize()
# atoms = struct.get_adsorbate_atoms("N2*")
# cnum = analysis.get_coordination_number("ZrN_311", 0, "N2*", site="min")
# print(cnum)
# struct.visualize()

# pzc = analysis.get_surface_pzcs(she=False)
# print(pzc)


# coord_env = struct.coordination_env()
# print(min_site)
# print(analysis.get_spans('HER', 0))
# binding = analysis.get_binding_energies(0, "HER")
# print(binding)
# surface = "Zr2Si_112"
# fig1, ax1 = analysis.plot_FED("NRR", surface, 0)
# fig1, ax1 = analysis.plot_FED("NRR", surface, -0.5, color="#003EF5", graph_objects=(fig1, ax1))
# fig1, ax1 = analysis.plot_FED("HER", surface, -0.5, color="#003EF5", graph_objects=(fig1, ax1))
# plt.show()
# plt.savefig(os.path.join(windows_path, f"C:\\Users\\Cooper\\OneDrive - UCB-O365\\Research\\Binaries\\figures\\{surface}_HER_FED.png"))