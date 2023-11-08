# cluster testing
import sys
import os
from matplotlib import pyplot as plt
import matplotlib as mpl
script_path = "C:\\Users\\Cooper\\OneDrive - UCB-O365\\github\\JDFT_tools\\RxNetwork"
linux_path = "/home/cooper/Research/Binaries"
windows_path = "C:\\Users\\Cooper\\OneDrive - UCB-O365\\Research\\Binaries\\all_data"
figure_path = "C:\\Users\\Cooper\\OneDrive - UCB-O365\\Research\\Binaries\\figures"
downloads_path = "C:\\Users\\coopy\\Downloads"
filename = "fake_all_data (1).json"
sys.path.append(script_path)

from Analyzer import Analyzer
from structure.Site import Site

analysis = Analyzer(data_path=windows_path, filename=filename, bulks=['RhSi'])
surface = "Co2Si_100"
intermediate = "N2*"
reaction = "NRR"    
bias = 0
struct = analysis.get_structure(surface, bias, 'N2H*')
# print(struct.intermediate_neighbors(intermediate))
site = Site.get_site(windows_path, filename=filename, surface=surface, intermediate="NH*", bias=bias)
# print(site.unique_clusters("NH*"), "testing")
# struct.visualize()
# unique_clusters = analysis.count_unique_clusters(surface, intermediate, bias)
# print(unique_clusters)
# FED = analysis.cluster_FED(surface, bias)
# print(FED)
# spans = analysis.cluster_span(surface, bias)
# print(spans)
plot = analysis.plot_cluster_FED(reaction, surface, bias)
fig_str = f"{surface}_{bias}_cluster_FED.png"
plt.savefig(os.path.join(figure_path, fig_str))