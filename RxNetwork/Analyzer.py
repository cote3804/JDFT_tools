#!/usr/bin/envs python3

from data.Data_Parser import Data_Parser
from data.Materials import Materials
from data.Calculator import Calculator
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from ase.visualize import view


class Analyzer:
    def __init__(self, data_path:str, filename:str, bulks:list) -> None:
        self.bulks = bulks
        self.data = Data_Parser(data_path, filename)

    def calculate_surface_properties(self, materials:Materials, ev=True): # this method needs to return a dictionary with all the necessary values for analysis
        calculator = Calculator(self.data, reaction="NRR")
        # surface_data = materials.get_surface_data()
        calculator.calculate_surface_properties(materials)
        # return materials

    def do_analysis(self) -> None:
        materials = Materials(self.data, self.bulks)
        self.materials = materials
        self.calculate_surface_properties(materials)

    def get_data(self) -> dict:
        data_dict = {}
        for material in self.materials.materials:
            for surface in material.surfaces:
                data_dict.update({surface:material.get_surface_data(surface)})
        return data_dict
    
    def surface_data(self, surface:str) -> dict:
        bulk = surface.split("_")[0]
        material = self.materials.get_material(bulk)
        material_data = material.energies
        surface_data = material_data[surface]
        return surface_data
        
        
    def get_FED_energies(self, bias, referenced="final") -> dict:
        FED_energies = {}
        bias = self.bias_float_to_str(bias)
        for material in self.materials.materials:
            for surface in material.surfaces:
                FED_energies[surface] = material.get_FED_energy(surface, bias, referenced=referenced)
        
        return FED_energies
    
    def bias_float_to_str(self, bias:float) -> str:
        return f"{bias:.2f}" + "V"
                
    def visualize_contcar(self, bias, surface, intermediate, site=None) -> None:
        bias = self.bias_float_to_str(bias)
        if site == None:
            site, energy = self.data.get_lowest_site(surface, intermediate, bias)
        
        struct_dict = self.data.get_contcar(surface, bias, intermediate, site)
        struct = Structure.from_dict(struct_dict)
        atoms = AseAtomsAdaptor.get_atoms(struct)
        view(atoms)