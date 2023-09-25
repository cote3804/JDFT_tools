#!/usr/bin/envs python3

from data.Data_Parser import Data_Parser
from data.Materials import Materials, Material
from data.Calculator import Calculator
from plotting.FEDPlotter import FED_plotter
from network.Network import Network
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from ase.visualize import view
from matplotlib import pyplot as plt


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
            for surface in material.converged_surfaces():
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
            for surface in material.converged_surfaces():
                FED_energies[surface] = material.get_FED_energy(surface, bias, referenced=referenced)
        
        return FED_energies
    
    def get_FED_energy(self, surface:str, bias:str, referenced="final") -> float:
        bulk = surface.split("_")[0]
        material = self.materials.get_material(bulk)
        return material.get_FED_energy(surface, bias, referenced=referenced)
    
    def bias_float_to_str(self, bias:float) -> str:
        return f"{bias:.2f}" + "V"
                
    def visualize_contcar(self, bias, surface, intermediate, site=None) -> None:
        bias = self.bias_float_to_str(bias)
        if site == None:
            site, energy = self.data.get_lowest_site(surface, intermediate, bias)
        
        struct_dict = self.data.get_contcar(surface, bias, intermediate, site)
        print(struct_dict)
        struct = Structure.from_dict(struct_dict)
        atoms = AseAtomsAdaptor.get_atoms(struct)
        view(atoms)

    def get_span(self, reaction:str, material:Material, surface:str, bias:float) -> dict:
        network = Network(reaction)
        network.add_data_to_nodes(material.get_FED_energy(surface, bias))
        intermediates = material.get_converged_intermediates()
        subgraph = network.connected_subgraph(intermediates[surface])
                # network.reconnect(subgraph)
        # subgraph = network.connected_subgraph(intermediates)
        return subgraph
    #TODO implement reaction network graph theory stuff
    
    def get_spans(self, reaction:str) -> dict:
        network = Network(reaction)
        for material in self.materials.materials:
            intermediates = material.get_converged_intermediates()
            for surface in intermediates.keys():
                for bias in material.get_biases(surface):
                    FED_energy = self.get_FED_energy(surface, bias)
                    span = self.get_span(reaction, material, surface, bias)
        pass
        # method to get all spans for all the surfaces and biases

    def get_binding_energies(self, bias, reaction="NRR") -> dict:
        calculator = Calculator(self.data, reaction)
        bias = self.bias_float_to_str(bias)
        binding_energies = calculator.calculate_binding_energies(self.materials, bias)
        return binding_energies
    
    def plot_FED(self, reaction:str, surface:str, bias:float, referenced="final", color="#f00000", graph_objects=None) -> plt.Figure:
        bias = self.bias_float_to_str(bias)
        FED_energy = self.get_FED_energy(surface, bias, referenced=referenced)
        plotter = FED_plotter(reaction, FED_energy)
        fig, ax = plotter.plot(surface, bias, graph_objects=graph_objects, color=color)
        state_length = 1
        connector_length = 1/2
        # for i in range(len(FED_energy)):
        # plt.show()
        return fig, ax

    def save_data(self, path, filename): # save data to json
        pass