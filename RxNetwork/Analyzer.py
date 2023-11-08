#!/usr/bin/envs python3

from data.Data_Parser import Data_Parser
from data.Materials import Materials, Material
from data.Calculator import Calculator, Reaction
from structure.Structure import Structure
from structure.Site import Site, Cluster
from plotting.FEDPlotter import FED_plotter
from network.Network import Network
# from pymatgen.core import Structure
# from pymatgen.io.ase import AseAtomsAdaptor
from ase.visualize import view
from matplotlib import pyplot as plt
from helpers.conversions import bias_float_to_str


class Analyzer:
    def __init__(self, data_path:str, filename:str, bulks:list) -> None:
        self.bulks = bulks
        self.data = Data_Parser(data_path, filename)

    def calculate_surface_properties(self, materials:Materials, reaction:str, ev=True): # this method needs to return a dictionary with all the necessary values for analysis
        calculator = Calculator(self.data, reaction=reaction)
        # surface_data = materials.get_surface_data()
        calculator.calculate_surface_properties(materials)
        # return materials

    def do_analysis(self, reaction="NRR") -> None:
        materials = Materials(self.data, self.bulks)
        self.materials = materials
        self.calculate_surface_properties(materials, reaction)

    def get_data(self, reaction:str) -> dict:
        data_dict = {}
        materials = Materials(self.data, self.bulks)
        calculator = Calculator(self.data, reaction)
        for material in materials.materials:
            material_energy = calculator.calculate_material_energies(material)
            for surface in material.surfaces:
                data_dict.update({surface:material_energy[surface]})
        return data_dict
    
    def surface_data(self, surface:str) -> dict:
        bulk = surface.split("_")[0]
        material = self.materials.get_material(bulk)
        material_data = material.energies
        surface_data = material_data[surface]
        return surface_data
        
        
    def get_FED_energies(self, bias, referenced="final", reaction="NRR") -> dict:
        calculator = Calculator(self.data, reaction="NRR")
        FED_energies = {}
        bias = self.bias_float_to_str(bias)
        materials = Materials(self.data, self.bulks)
        for material in materials.materials:
            for surface in material.converged_surfaces():
                FED_energies[surface] = calculator.get_FED_energy(material, surface, bias, referenced=referenced)
        
        return FED_energies
    
    def get_FED_energy(self, surface:str, bias:float, reaction:str, referenced="final") -> float:
        calculator = Calculator(self.data, reaction)
        bias = self.bias_float_to_str(bias)
        bulk = surface.split("_")[0]
        materials = Materials(self.data, self.bulks)
        material = materials.get_material(bulk)
        return calculator.get_FED_energy(material, surface, bias, referenced=referenced)
    
    def bias_float_to_str(self, bias:float) -> str:
        return f"{bias:.2f}" + "V"

    #TODO get this pymatgen visualization issue working. Pymatgen from_dict() method is failing.   
    # def visualize_contcar(self, bias, surface, intermediate, site=None) -> None:
    #     bias = self.bias_float_to_str(bias)
    #     if site == None:
    #         site, energy = self.data.get_lowest_site(surface, intermediate, bias)
        
    #     struct_dict = self.data.get_contcar(surface, bias, intermediate, site)
    #     print(struct_dict)
    #     struct = Structure.from_dict(struct_dict)
    #     atoms = AseAtomsAdaptor.get_atoms(struct)
    #     view(atoms)

    def get_span(self, reaction:str, surface:str, bias:float) -> dict:
        bias = self.bias_float_to_str(bias)
        calculator = Calculator(self.data, reaction)
        gmax, span = calculator.calculate_span(reaction, surface, bias)
        # network = Network(reaction)
        # network.add_data_to_nodes(material.get_FED_energy(surface, bias))
        # intermediates = material.get_converged_intermediates()
        # subgraph = network.connected_subgraph(intermediates[surface])
                # network.reconnect(subgraph)
        # subgraph = network.connected_subgraph(intermediates)
        # return subgraph
        return gmax, span
    #TODO implement reaction network graph theory stuff
    
    def get_spans(self, reaction:str, bias=0) -> dict:

        # network = Network(reaction)
        materials = Materials(self.data, self.bulks)
        self.materials = materials
        span_dict = {}
        bias_str = self.bias_float_to_str(bias)
        for material in materials.materials:
            for surface in material.surfaces:
                if self.data.check_any_adsorbate_convergence(surface, bias_str) == False:
                    print(f"no converged adsorbates for {surface}")
                    continue
                gmax, span = self.get_span(reaction, surface, bias)
                span_dict.update({surface:(gmax,span)})
        return span_dict
        # method to get all spans for all the surfaces and biases

    def get_binding_energies(self, bias, reaction="NRR") -> dict:
        calculator = Calculator(self.data, reaction)
        reaction = Reaction(reaction)
        adsorbate = reaction.binding_adsorbate()
        bias = self.bias_float_to_str(bias)
        materials = Materials(self.data, self.bulks)
        binding_energies = calculator.calculate_binding_energies(materials, bias, adsorbate=adsorbate.strip("*"))
        return binding_energies
    
    def plot_FED(self, reaction:str, surface:str, bias:float, referenced="final", color="#f00000", graph_objects=None) -> plt.Figure:
        # bias = self.bias_float_to_str(bias)
        FED_energy = self.get_FED_energy(surface, bias, reaction, referenced=referenced)
        if FED_energy == None:
            return None
        plotter = FED_plotter(reaction, FED_energy)
        fig, ax = plotter.plot(surface, bias, graph_objects=graph_objects, color=color)
        state_length = 1
        connector_length = 1/2
        # for i in range(len(FED_energy)):
        # plt.show()
        return fig, ax

    def get_surface_charge(self, surface, bias, reaction):
        calculator = Calculator(self.data, reaction)
        bias = self.bias_float_to_str(bias)
        bulk = surface.split("_")[0]
        materials = Materials(self.data, self.bulks)
        material = materials.get_material(bulk)
        return calculator.get_pathway_charges(material, surface, bias, reference="initial")

    def get_pathway_charges(self, bias:float, reaction:str) -> dict:
        materials = Materials(self.data, self.bulks)
        charge_data = {}
        for material in materials.materials:
            for surface in material.surfaces:
                charge_data[surface] = self.get_surface_charge(surface, bias, reaction)
        return charge_data

    def get_sites_data(self, surface, bias, adsorbate):
        return self.data.get_sites_data(surface, bias, adsorbate)

    def get_min_site(self, surface, bias, adsorbate):
        bias = self.bias_float_to_str(bias)
        return self.data.get_lowest_site(surface, adsorbate, bias)

    def save_data(self, path, filename): # save data to json
        pass

    def get_structure(self, surface:str, bias:float, adsorbate:str, site="min"):
        bias = self.bias_float_to_str(bias)
        intermediate_string = adsorbate.strip('*')
        if site == "min":
            site, energy = self.data.get_lowest_site(surface, intermediate_string, bias)
        else:
            pass
        if site == '':
            print(f"site not found for {surface} with {adsorbate}")
            return None
        sites_data = self.data.get_sites_data(surface, bias, intermediate_string)
        struct_dict = self.data.get_contcar(surface, bias, intermediate_string, site)
        structure = Structure(struct_dict)
        return structure
    
    def get_coordination_number(self, surface:str, bias:float, adsorbate:str, buffer=1, site="min"):
        bias_str = self.bias_float_to_str(bias)
        if self.data.check_adsorbed_convergence(surface, bias_str, adsorbate.strip("*")):
            structure = self.get_structure(surface, bias, adsorbate, site)
            return structure.coordination_number(adsorbate, buffer=buffer)
        else:
            return None

    def count_unique_clusters(self, surface:str, intermediate:str, bias:float) -> int:
        bias_str = self.bias_float_to_str(bias)
        if self.data.check_adsorbed_convergence(surface, bias_str, intermediate.strip("*")):
            site_data = self.data.get_sites_data(surface, bias_str, intermediate.strip("*"))
            site = Site(site_data)
            unique_clusters = site.unique_clusters(intermediate)
            return unique_clusters
            return None

    def cluster_FED(self, surface:str, bias:float) -> dict:
        # method to calculate multiple FEDs using clustering algorithm
        # returns a dictionary with the FEDs for each cluster
        calculator = Calculator(self.data, reaction="NRR")
        bias_str = self.bias_float_to_str(bias)
        data = calculator.aggregate_cluster_FEDs(surface, bias_str, "N2*")
        return data
    
    def cluster_span(self, surface:str, bias:float) -> dict:
        # method to calculate multiple spans using clustering algorithm
        # returns a dictionary with the spans for each cluster
        calculator = Calculator(self.data, reaction="NRR")
        bias_str = self.bias_float_to_str(bias)
        data = calculator.aggregate_cluster_spans(surface, bias_str, "N2*")
        return data
    
    def plot_cluster_FED(self, reaction:str, surface:str, bias:float, referenced="final", graph_objects=None, color="#000000") -> plt.Figure:
        cluster_FEDs = self.cluster_FED(surface, bias)
        colors = ["#33B841","#353FC4","#CF3530","#33B841","#353FC4","#CF3530","#33B841","#353FC4","#CF3530"]
        bias = bias_float_to_str(bias)
        for iclust, cluster in enumerate(cluster_FEDs.keys()):
            surface_str = f"{surface.split('_')[0]} ({surface.split('_')[1]}) {bias} {repr(cluster)}"
            plotter = FED_plotter(reaction, cluster_FEDs[cluster])
            if iclust == 0:
                if graph_objects == None:
                    graph_objects = plotter.plot(surface_str, bias, color=colors[iclust], label_str=surface_str)
                else:
                    graph_objects = plotter.plot(surface_str, bias, graph_objects=graph_objects, color=colors[iclust], label_str=surface_str)
            else:
                graph_objects= plotter.plot(surface_str, bias, graph_objects=graph_objects, color=colors[iclust], label_str=surface_str)
        return graph_objects
    
    def get_band_centers(self, bias=0):
        calculator = Calculator(self.data, reaction="NRR")
        band_centers = {}
        bias = self.bias_float_to_str(bias)
        materials = Materials(self.data, self.bulks)
        for material in materials.materials:
            for surface in material.converged_surfaces():
                band_centers[surface] = calculator.calculate_band_center(surface, bias)
        
        return band_centers
    
    def get_surface_pzcs(self, she=True):
        calculator = Calculator(self.data, reaction="NRR")
        surface_PZCs = {}
        materials = Materials(self.data, self.bulks)
        for material in materials.materials:
            for surface in material.converged_surfaces():
                surface_PZCs[surface] = calculator.get_pzc(surface, she)
        return surface_PZCs
    
    def bias_string(self, bias:float) -> str:
        return bias_float_to_str(bias)
    
    def adsorbate_bader(self, surface, bias, adsorbate):
        # returns oxidation states for adsorbate atoms on specified surface.
        calculator = Calculator(self.data, reaction="NRR")
        structure = self.get_structure(surface, bias, adsorbate)
        if structure == None:
            return None
        atom_indices = structure.get_adsorbate_atoms(adsorbate)
        bader = {}
        bias_str = bias_float_to_str(bias)
        try:
            for atom in atom_indices:
                bader[atom] = self.data.bader_charge(surface, bias_str, adsorbate.strip("*"), str(atom+1))
        except:
            print(f"no bader data for {surface} {bias_str} {adsorbate}")
            return None
        return bader