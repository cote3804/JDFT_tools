from monty.json import MontyDecoder
import os
import json

class Data_Parser:
    def __init__(self, data_path, filename="all_data.json") -> None:
        with open(os.path.join(data_path,filename),'r') as f:
            self.all_data = json.load(f, cls=MontyDecoder) # load in all_data as dict

    def get_surfaces(self, bulk) -> list:
        surfaces = []
        for mat_key in self.all_data.keys():
            if mat_key.split("_")[0] == bulk:
                surfaces.append(mat_key)
        return surfaces
    
    def surface_data(self, surface) -> dict:
        return self.all_data[surface]
    
    def get_adsorbate_energy(self, surface:str, bias:str, intermediate:str) -> float:
        return self.all_data[surface]["adsorbed"][bias][intermediate]["final_energy"]