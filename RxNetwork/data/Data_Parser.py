from monty.json import MontyDecoder
from structure.Structure import Structure
import os
import json

class Data_Parser:
    def __init__(self, data_path, filename="all_data.json") -> None:
        with open(os.path.join(data_path,filename),'r') as f:
            self.all_data = json.load(f) # load in all_data as dict

    def get_surfaces(self, bulk) -> list:
        surfaces = []
        for mat_key in self.all_data.keys():
            if "_" in mat_key:
                index_key = mat_key.split("_")[1] # need to check that whatever follows the underscore is a surface facet and not a random tag
                characters = [character for character in index_key.replace("-","")] # need to strip negatives from negative surface factes
                is_surface = all([character.isdigit() for character in characters]) # check to see that all are numbers else it's not a surface
            if mat_key.split("_")[0] == bulk and "_" in mat_key and is_surface:
                surfaces.append(mat_key)
        if len(surfaces) == 0:
            raise Exception(f"no surfaces found for {bulk}")
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
        if "adsorbed" not in self.all_data[surface].keys():
            return "", 1000
        if intermediate not in self.all_data[surface]["adsorbed"].keys():
            return "", 1000
        if bias not in self.all_data[surface]["adsorbed"][intermediate].keys():
            return "", 1000
        site_data = self.all_data[surface]["adsorbed"][intermediate][bias]
        lowest_energy = 1000
        lowest_site = ""
        for site, data in site_data.items():
            if "contcar" in data.keys():
                struct_dict = data["contcar"]
                struct = Structure(struct_dict)
                dissociated = struct.check_dissociation(intermediate) # currently doesn't do anything but return False
                if dissociated:
                    continue
            if "final_energy" in data.keys() and bool(data["converged"]) != False:
                if data["final_energy"] < lowest_energy:
                    lowest_energy = data["final_energy"]
                    lowest_site = site
            elif "final_energy" not in data.keys():
                continue
        return lowest_site, lowest_energy

    def get_sites_data(self, surface:str, bias:str, intermediate:str, underscored=False) -> dict:
        # underscored is used to accept/reject underscored sites. Keep it false to reject underscored sites.
        # underscored sites typically notate failed calculations
        site_data = self.all_data[surface]["adsorbed"][intermediate][bias]
        underscored_sites = []
        for site in self.all_data[surface]["adsorbed"][intermediate][bias].keys():
            if site.startswith('_'):
                underscored_sites.append(site)
        if underscored == False:
            for site in underscored_sites:
                site_data.pop(site)
        return site_data
    
    
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
        if bias not in self.all_data[surface]["surf"].keys():
            return False
        return bool(self.all_data[surface]["surf"][bias]["converged"])

    def check_any_surface_convergence(self, surface:str) -> bool:
        # check whether surface converged for any bias
        if 'surf' not in self.all_data[surface].keys():
            return False
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
        if site not in self.all_data[surface]["adsorbed"][intermediate][bias].keys():
            print(site, "not found")
            return False
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

    def check_any_adsorbate_convergence(self, surface, bias):
        if self.check_if_surface_has_adsorbed(surface) == False:
            print(f"{surface} has no adsorbates")
            return False
        truth = any([self.check_adsorbed_convergence(surface, bias, intermediate) for intermediate in self.get_intermediates(surface)])
        return truth
    
    def get_adsorbed_charge(self, surface, bias, intermediate, site=None):
        if site == None:
            site = self.get_lowest_site(surface, intermediate, bias)[0]
        elif site != None:
            pass
        return self.all_data[surface]["adsorbed"][intermediate][bias][site]["nfinal"]
    
    def get_surface_charge(self, surface, bias):
        if "nfinal" in self.all_data[surface]["surf"][bias].keys():
            return self.all_data[surface]["surf"][bias]["nfinal"]
        else:
            print(f"{surface} at {bias} has no surface charge data")
            return None
        
    def minimum_energy_of_sites(self, surface, bias, intermediate, sites:list):
        min_energy = 1000
        sites_data = self.get_sites_data(surface, bias, intermediate)
        min_site = ''
        for site in sites:
            energy = sites_data[site]["final_energy"]
            if energy < min_energy:
                min_energy = energy
                min_site = site

        return min_site, min_energy
    
    def get_band_data(self, surface, bias):
        # might need to convergence check this
        return self.all_data[surface]["surf"][bias]["band_means"]

    def surface_pzc(self, surface):
        if "No_bias" in self.all_data[surface]["surf"].keys() and "eigStats" in self.all_data[surface]["surf"]["No_bias"].keys():
            return self.all_data[surface]["surf"]["No_bias"]["eigStats"]["mu  "]
        else:
            print(f"No pzc data for {surface}")
            return None
    
    def bader_charge(self, surface, bias, adsorbate, atom_index, site="min"):
        if site == "min":
            site = self.get_lowest_site(surface, adsorbate, bias)[0]
        return self.all_data[surface]["adsorbed"][adsorbate][bias][site]["bader"][atom_index]['OXIDATION STATE']