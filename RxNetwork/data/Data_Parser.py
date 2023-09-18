from monty.json import MontyDecoder
import os
import json

class Data_Parser:
    def __init__(self, data_path, filename="all_data.json") -> None:
        with open(os.path.join(data_path,filename),'r') as f:
            self.all_data = json.load(f) # load in all_data as dict

    def get_surfaces(self, bulk) -> list:
        surfaces = []
        for mat_key in self.all_data.keys():
            if mat_key.split("_")[0] == bulk and "_" in mat_key:
                surfaces.append(mat_key)
        return surfaces
    
    def surface_data(self, surface) -> dict:
        return self.all_data[surface]
    
    def get_surface_energy(self, surface:str, bias:str) -> float:
        return self.all_data[surface]["surf"][bias]["final_energy"]
    
    def get_adsorbed_energy(self, surface:str, bias:str, intermediate:str, site=None) -> float:
        sites_data = self.all_data[surface]["adsorbed"][intermediate][bias]
        if site != None:
            return sites_data[site]["final_energy"] 
        elif site == None:
            lowest_site, lowest_energy = self.get_lowest_site(surface, intermediate, bias)
            if len(lowest_site) > 0:
                return lowest_energy
            elif len(lowest_site) == 0: # no converged sites
                print(f"{intermediate} on {surface} at {bias} has no converged sites")    
                return 1000
    
    def get_lowest_site(self, surface:str, intermediate:str, bias:str) -> dict:
        site_data = self.all_data[surface]["adsorbed"][intermediate][bias]
        lowest_energy = 1000
        lowest_site = ""
        for site, data in site_data.items():
            if "final_energy" in data.keys() and bool(data["converged"]) != False:
                if data["final_energy"] < lowest_energy:
                    lowest_energy = data["final_energy"]
                    lowest_site = site
            elif "final_energy" not in data.keys():
                continue
        return lowest_site, lowest_energy
    

    def get_sites_data(self, surface:str, bias:str, intermediate:str) -> dict:
        return self.all_data[surface]["adsorbed"][intermediate][bias]
    
    def get_intermediate_biases(self, surface, intermediate):
        biases = []
        for bias in self.all_data[surface]["adsorbed"][intermediate].keys():
            biases.append(bias)
        return biases
    
    def get_intermediates(self, surface):
        return [i for i in self.all_data[surface]["adsorbed"].keys()]
    

    def get_contcar(self, surface, bias, intermediate, site) -> dict:
        return self.all_data[surface]["adsorbed"][intermediate][bias][site]["contcar"]