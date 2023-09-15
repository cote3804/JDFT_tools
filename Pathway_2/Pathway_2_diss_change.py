
"""
Created on Thu Jun  8 17:55:46 2023

@author: coopy
"""
import json
import os
import sys
import numpy as np
import numpy.ma as ma #numpy masked arrays
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from Reaction_Network import Reaction_Network

#TODO currently the energy matrix is stored as an environment variable unreferenced 
# and in Hartrees. I could consider using a method to convert the matrix every time
# I want it referenced and in eV

class Network_Analyzer:
    def __init__(self, path, file, reaction):
        self.reaction = reaction
        self.materials = []
        self.H_to_eV = 27.2114
        network = Reaction_Network()
        network.load_network(reaction)
        self.network_dict = network.get_network_dict()
        with open(os.path.join(path,file)) as f:
            self.all_data = json.load(f)
        self.network = network


    
    def network_to_thermo(self, material, bias, site):
        pathways = self.network_dict['pathways']["labels"]
        network_dimensions = self.rxn_network_dimensions()
        thermo_matrix = np.zeros((len(network_dimensions), max(network_dimensions.values())))
        thermo_mask = np.full((thermo_matrix.shape), False)
        
        for node in self.network_dict["network"]:
            for pathway in node["network_data"]["pathways"]:
                row_index = pathways.index(pathway)
                column_index = node["network_data"]["index"]
                thermo_matrix[row_index, column_index] = self.calculate_intermediate(node["energy"], material, bias, site)
                thermo_mask[row_index, column_index] = True
        return thermo_matrix, thermo_mask
    
    def calculate_intermediate(self, energy_dict: dict, material, bias, site):
        # Calculate one intermediate in the list of intermediates in "network" in
        # the reaction_data dictionary for a given reaction
        E_mol = 0
        E_surf = 0
        
        # loop through molecules and add them to the molecule energy
        for mol, mol_number in energy_dict["molecules"].items():
            E_mol += self.calculate_mol(mol,bias)*mol_number
        
        for ads, ads_number in energy_dict["surfaces"].items():
            E_surf += self.calculate_surface(ads, material, bias)*ads_number #TODO implement site specification
        total_energy = E_mol + E_surf
        return total_energy
        
    def calculate_surface(self, adsorbate, material, bias, site='min'):
        if adsorbate == "*": # loan surface
            energy = self.all_data[material]['surf'][bias]['final_energy']
            return energy
        E_data = self.all_data[material]['adsorbed'][adsorbate.replace('*','')][bias]
        # E_data is the dictionary with the sites keys and all of the calculation
        # data within each site.
        if site == 'min':
            E_old = 0
            for site in E_data.keys():
                try:
                    E_new = E_data[site]["final_energy"]
                except:
                    raise Exception(f"final energy for {adsorbate} on {material} at {bias} on site {site}"
                                    f" not found. Check convergence.")
                if E_new < E_old:
                    E_old = E_new
                    min_site = site #TODO maybe return site number
                else:
                    pass
            energy = E_old
        else:
            try:
                energy = E_data[site]
            except:
                raise Exception(f"site {site} not in all_data for material {material}"
                                f" At {bias} for adsorbate {adsorbate}.")
        return energy
            
    def calculate_mol(self, mol, bias):
        if mol == "H+": #TODO add other proton calculation methods or settle on a good one
            energy = (1/2) * self.all_data["H2"]["0.00V"]["final_energy"]
        else:
            try:
                energy = self.all_data[mol][bias]['final_energy']
            except:
                raise Exception(f"molecule {mol} at bias {bias} not found or not converged. "
                                f"Check all_data.")
        return energy
    
        
    def rxn_network_dimensions(self):
        '''

        Returns
        -------
        pathway_count : dict
            this dictionary specifies the dimensions of the reaction network. Each
            key is another pathway and the value associated with each key is the
            number of intermediate states for that pathway. This includes the initial and final states

        '''
        network_dimensions = {}
        for node in self.network_dict['network']:
            pathways = node["network_data"]["pathways"]
            for pathway in pathways:
                network_dimensions[pathway] = network_dimensions.get(pathway, 0) + 1
        return network_dimensions
    
    def add_material(self, material, bias, site="min"):
        self.materials.append({"material":material, "bias":bias, "site":site})
    
    def calculate_energies(self, ev=True, referenced='final'):
        '''
        Returns
        -------
        energy_array : numpy array
            Contains all energetics data for all materials in the added to the class
            for the specified reaction pathway.
            Index 0: pathway specified in rxn network
            Index 1: states in longest pathway
            Index 2: material in order in which materials were added
        '''
        material_num = len(self.materials)
        matrix_list = []
        mask_list = []
        for material_dict in self.materials:
            thermo_matrix, thermo_mask = self.network_to_thermo(material_dict["material"], material_dict["bias"], material_dict["site"])
            matrix_list.append(thermo_matrix)
            mask_list.append(thermo_mask)
            
        energy_array = np.stack(matrix_list, axis=2)
        energy_mask = np.stack(mask_list, axis=2)
        self.energy_array = energy_array
        self.energy_mask = energy_mask
        
        if referenced == "final":
            subtraction_matrix = energy_array[:,-1:,:]
            energy_array = energy_array - subtraction_matrix 
        if ev == True:
            energy_array = energy_array * self.H_to_eV
        return energy_array, energy_mask
    
    def FED_plot(self):
        #FED stands for free energy diagram
        
        energy_array, energy_mask = self.calculate_energies()
        
        ################# Plot Formatting ######################
        starting_point = 0.5
        number_of_states = len(self.energy_array[0,:,0])
        state_length = 1
        connection_length = 0.75
        
        fig,ax = plt.subplots(dpi=300)
        lines = []
        linestyles = []
        input_colors = ['k','b','r']
        colors = []
        
        for imat, material_dict in enumerate(self.materials):
            material = material_dict["material"]
            mat_matrix = energy_array[:,:,imat]
            mat_mask = energy_mask[:,:,imat]
            for pathway_vector, pathway_mask in zip(mat_matrix, mat_mask):
                for istate, state in enumerate(pathway_vector):
                    # state is an energy value
                    if istate == 0 and pathway_mask[istate]:
                        lines.append([(starting_point, state),
                                      (starting_point+state_length, state)])
                        linestyles.append("solid")
                        colors.append(input_colors[imat])
                    elif istate != 0 and pathway_mask[istate]:
                        #draw connection from previous state
                        if pathway_mask[istate-1]: 
                            lines.append([(starting_point + istate*state_length + 
                                           (istate-1)*connection_length, pathway_vector[istate - 1]),
                                          (starting_point + istate*state_length + (istate)*connection_length, state)])
                            linestyles.append("dashed")
                            colors.append(input_colors[imat])
                        else: #last state was invalid, need to go look for most recent valid state to connect to
                            indices = np.where(pathway_mask)[0]
                            indices = [i for i in indices if i < istate-1]
                            last_valid_index = indices[-1]
                            lines.append([(starting_point + (last_valid_index+1)*state_length + 
                                           (last_valid_index)*connection_length, pathway_vector[last_valid_index]),
                                          (starting_point + istate*state_length + (istate)*connection_length, state)])
                            linestyles.append("dashed")
                            colors.append(input_colors[imat])
                        
                        #now draw current state
                        lines.append([(starting_point + istate*state_length + (istate)*connection_length, state),
                                      (starting_point + (istate+1)*state_length + (istate)*connection_length, state)])
                        linestyles.append("solid")
                        colors.append(input_colors[imat])
        # print(lines)
        line_segments = LineCollection(lines, linewidths = (1), linestyle=linestyles, colors=colors)
        # print(line_segments)
        ax.add_collection(line_segments)
        ax.set_title("plotting test")
        ax.set_xlim(0, mat_matrix.shape[-1]*(state_length+connection_length) + starting_point + 1)
        # ax.set_ylim(-2511.9, -2511.6)
        ax.set_ylim(-4,4)
        plt.show()
                   
    def calculate_energetic_spans(self):

        def energetic_span(energies: np.ndarray, mask: np.ndarray) -> float:
            # Assumes energies are ordered and first energy corresponds to starting state and
            # last energy corresponds to final state such that the reaction energy is
            # E_final - E_initial.
            
            # Energetic span is defined as the largest free energy difference within a catalytic cycle.
            # Read more about it here: https://doi.org/10.1021/ar1000956
            # If the largest free energy difference occurs across the edge of the catalytic cycle,
            # the reaction free energy must be subtracted such that E_span = E_max - E_min - E_rxn.
            energies = ma.masked_array(energies, ~mask)
            # The ~mask flips the boolean values. For some reason they need to be opposite
            imax = np.argmax(energies)
            imin = np.argmin(energies)
            if imax < imin:
                E_span = energies[imax] - energies[imin] - (energies[0] - energies[-1])
            else:
                E_span = energies[imax] - energies[imin]
            return E_span, (imin,imax)
        
        energy_array, energy_mask = self.calculate_energies()
        E_span_dict = {}
        for imat, material_dict in enumerate(self.materials):
            mat_matrix = energy_array[:,:,imat]
            mat_mask = energy_mask[:,:,imat]
            E_span_dict[material_dict["material"]] = {}
            for ipathway, (pathway_vector, pathway_mask) in enumerate(zip(mat_matrix, mat_mask)):
                E_span, span_indices = energetic_span(pathway_vector, pathway_mask)
                pathway_name = self.network.index_to_pathway_name(ipathway)
                min_state_name = self.network.indices_to_state((ipathway, span_indices[0]))
                max_state_name = self.network.indices_to_state((ipathway, span_indices[-1]))
                E_span_dict[material_dict["material"]][pathway_name] = [min_state_name, max_state_name, E_span]
        return E_span_dict

    def get_intermediate_energies(self, intermediate, material):
        if site is None:
            site = "min"

if __name__ == "__main__":
    file = 'test_Combined_all_data.json'
    path = os.path.normpath('C:\\Users\\coopy\\OneDrive - UCB-O365\\Research\\N2R_Scaling\\all_data\\Paper all_data')
    ntwrk = Network_Analyzer(path, file, "NRR")
    ntwrk.add_material("Ru_111", "0.00V", "min")
    ntwrk.add_material("Re_111", "0.00V", "min")
    ntwrk.add_material("Rh_111", "0.00V", "min")
    array, mask = ntwrk.calculate_energies()
    ntwrk.FED_plot()
    print(ntwrk.calculate_energetic_spans())