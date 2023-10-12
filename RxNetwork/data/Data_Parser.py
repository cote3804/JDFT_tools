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
                return None
    
    def get_surface_biases(self, surface:str) -> list:
        # Checks which biases a surface has converged at
        return [bias for bias in self.all_data[surface]["surf"].keys() if bool(self.all_data[surface]["surf"][bias]["converged"]) == True]
    
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
        if "adsorbed" not in self.all_data[surface].keys():
            return []
        return [i for i in self.all_data[surface]["adsorbed"].keys()]
    
    def check_surface_convergence(self, surface:str, bias:str) -> bool:
        # check whether surface calc converged
        # print(surface, self.all_data[surface]["surf"][bias]["converged"])
        if 'surf' not in self.all_data[surface].keys():
            return False
        return bool(self.all_data[surface]["surf"][bias]["converged"])

    def check_any_surface_convergence(self, surface:str) -> bool:
        # check wheter surface converged for any bias
        truth = any([bool(self.all_data[surface]["surf"][bias]["converged"]) for bias in self.all_data[surface]["surf"].keys()])
        return truth
    
    def check_adsorbed_convergence(self, surface:str, bias:str, intermediate:str) -> bool:
        converged = []
        if self.check_if_surface_has_adsorbed(surface) == False:
            return False
        intermediates = self.get_intermediates(surface)
        if intermediate not in intermediates: 
            return False
        if self.check_adsorbed_bias(surface, bias, intermediate) == False:
            return False
        elif self.check_adsorbed_bias(surface, bias, intermediate) == True:
            for site, data in self.get_sites_data(surface, bias, intermediate).items():
                if "converged" in data.keys():
                    converged.append(bool(data["converged"]))
            return any(converged)
    
    def check_adsorbed_bias(self, surface:str, bias:str, intermediate:str) -> bool:
        # checks if a specified adosrbate on a surface has been calculated at a specified bias
        if bias in self.all_data[surface]["adsorbed"][intermediate].keys():
            return True
        else:
            return False
            
    def check_adsorbed_site_convergence(self, surface:str, bias:str, intermediate:str, site:str) -> bool:
        if "converged" in self.all_data[surface]["adsorbed"][intermediate][bias][site].keys():
            if bool(self.all_data[surface]["adsorbed"][intermediate][bias][site]["converged"]) == True:
                return True
            else:
                return False
        else:
            return False

    def get_contcar(self, surface, bias, intermediate, site) -> dict:
        return self.all_data[surface]["adsorbed"][intermediate][bias][site]["contcar"]
    
    def check_if_surface_has_adsorbed(self, surface) -> bool:
        if "adsorbed" in self.all_data[surface].keys():
            return True
        else:
            return False
    
    def get_converged_intermediates(self, surface, bias):
        converged_intermediates = []
        for intermediate in self.get_intermediates(surface):
            if self.check_adsorbed_convergence(surface, bias, intermediate):
                converged_intermediates.append(intermediate)
            else:
                continue
        return converged_intermediates
    
    def get_converged_intermediates_at_all_biases(self, surface):
        pass