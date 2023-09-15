# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 19:17:43 2023

@author: coopy
"""

# Pathway script V_2
###############################################################################
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import json
import os
import re
import pprint
H_to_ev = 27.2114
###############################################################################

class Pathway:
    def __init__(self, path=None, filename=None, reaction='N2R'):
        
        if path is None:
            self.path = os.getcwd()
        if filename is None:
            self.filename = 'all_data.json'
            
        # (a) tafel mechanism (b) Heyrovski
        self.rxn_dict = {'state_0': {'a': {'intermediate': '*', 'number': 2, 'molecules': {}, 'protons':2},
                                     'b':{'intermediate': '*', 'number': 1, 'molecules': {}, 'protons':2}},
                        'state_1': {'a': {'intermediate': 'H*','number': 1, 'molecules': {'*':1}, 'protons':1},
                                    'b': {'intermediate': 'H*','number': 1, 'molecules': {}, 'protons':1}},
                        'state_2': {'a': {'intermediate': 'H*', 'number': 2, 'molecules': {}, 'protons':0}},
                        'state_3': {'a': {'intermediate': '*', 'number': 2, 'molecules':{'H2':1}, 'protons':0},
                                    'b': {'intermediate': '*', 'number': 1, 'molecules':{'H2':1}, 'protons':0}}
                   }
        with open(os.path.join(path, filename),'r') as f:
            self.all_data = json.load(f)
            
        self.references = {'H+': (1/2)*-0.8231668041392806, 'H2': -0.8231668041392806}

        
    def get_all_data(self):
        return self.all_data
        
    def add_material(self, material: str):
        try:
            len(self.materials) #TODO better variable check and creation
        except:
            self.materials = []
        self.materials.append(material)
    
    def _calculate_energies(self, bias='0.00V', reference=True, eV = True):
        energies = {}
        all_data = self.all_data
        materials = self.materials
        rxn_dict = self.rxn_dict
        references = self.references
        
        energy_matrices = {} #dict with matrices in each material

        
        for material in materials:
            energy_matrix = np.zeros((len(rxn_dict['state_0']),len(rxn_dict)))
            energies[material] = {}
            for state_index in rxn_dict.keys():
                energies[material][state_index] = {}
                create_dict = True
                for path_branch, state_data in rxn_dict[state_index].items():
                    branch_index = ord(path_branch) - 97 #convert letter to number
                    
                    intermediate = state_data['intermediate']
                    if create_dict is True:
                        energies[material][state_index][intermediate] = {}
                        create_dict = False
                    
                    #----------------------- surfs ---------------------#
                    if intermediate == '*': # clean surface
                        surf_energy =  self.get_surface_energy(material, bias, state_data['number'])
                        proton_energy = state_data['protons'] * self.references['H+']
                        molecule_energy = self.get_molecule_energy(state_data['molecules'], material, bias)
                        energy = surf_energy + proton_energy + molecule_energy
                        
                    #-------------------- adsorbates -------------------#
                    elif intermediate.endswith('*') and intermediate != '*':
                        ads_energy = self.get_adsorbate_energy(material, intermediate, bias, state_data['number'])
                            
                        proton_energy = state_data['protons'] * references['H+']
                        molecule_energy = self.get_molecule_energy(state_data['molecules'], material, bias)
                        energy = ads_energy + molecule_energy + proton_energy
                    
                    energy_matrix[branch_index,int(state_index.split('_')[-1])] = energy
                    energies[material][state_index][intermediate][path_branch] = energy

            energy_matrices[material] = energy_matrix
            if reference:
                energy_matrix = np.subtract(energy_matrix, np.transpose(np.matrix(energy_matrix[:,-1])))
                energy_matrices[material] = np.array(energy_matrix)
            if eV:
                energy_matrix = 27.2114 * energy_matrix
                energy_matrices[material] = energy_matrix
                
        return energies, energy_matrices
    def get_molecule_energy(self, mol_dict, material, bias):
        #takes molecules dictionary from rxn_dict as input
        #returns molecule energy
        molecule_energy = 0
        if len(mol_dict) > 0:
            for molecule, number in mol_dict.items():
                if molecule == '*': #reference surface
                    molecule_energy += self.get_surface_energy(material, bias, number)
                elif molecule.endswith('*') and molecule != '*': # reference adsorbate
                    molecule_energy += self.get_adsorbate_energy(material, molecule, bias, number)
                else:
                    molecule_energy += self.references[molecule] * number
        else: molecule_energy = 0
        return molecule_energy
    
    def get_surface_energy(self, material, bias, number):
        # input material string, bias, and state_data dict from rxn_dict
        # output surface energy
        try:
            surf_energy = number * self.all_data[material]['surf'][bias]['final_energy']
        except:
            print(f"final energy not found for {material} surface at {bias}")
        

        return surf_energy
        
    def get_adsorbate_energy(self, material, intermediate, bias, number):
        intermediate = intermediate.replace('*','')
        
        try:
            ads_energy = 0
            for site_index, sites in self.all_data[material]['adsorbed'][intermediate][bias].items():
                if sites['final_energy'] < ads_energy: # find minimum site energy
                    ads_energy = sites['final_energy']
                    lowest_site = site_index #TODO could pull out lowest energy sites
        except:
            print(f"final energy for adsorbate {intermediate} on"
                  f" {material} at {bias} not found")
            
        ads_energy = number * ads_energy
        return ads_energy
    
        
        

if __name__ == "__main__":
    file = 'combined_N2R_all_data.json'
    path = os.path.normpath('C:\\Users\\coopy\\OneDrive - UCB-O365\\Research\\N2R_Scaling')
    path = Pathway(path, file)
    materials = ["Pt_111","Au_111","Ru_111","Re_111","Ag_111","Ir_111",
                 "Rh_111","Co_111_mag","Ni_111"]
    for material in materials:
        path.add_material(material)
    energies, energy_matrices = path._calculate_energies(reference=True)
    
    import matplotlib.pyplot as plt
    for material in materials:
        print(type(energy_matrices[material]))
        plt.plot([0,1,2,3])
    






# def read_rxn(file_name):
#     with open(os.path.join("C:\\Users\\coopy\\OneDrive\\Desktop\\", "rxn_data.txt"), 'r') as f:
#         text = f.read()
#     rxn_dict = {}
#     states = []
#     state_index = 0
#     branch_index = 0
#     intermediates = []
#     for line in text.split("\n"):
#         for istate,state in enumerate(line.split("->")): # split at reaction arrow
#             # rxn_dict[f"state_{state_index}"]
#             for species in state.split(" + "):
#                 if istate == 0: # left side of reaction arrow 
#                     if species.strip().endswith("*") and species.strip() not in intermediates: #if adsorbed intermediate
#                         intermediates.append(species.strip())
#                         rxn_dict[f"state_{state_index}"] = {}
#                         rxn_dict[f"state_{state_index}"][species.strip()] = branch_index
#                         state_index += 1
#                 elif istate == 1: # right side of reaction arrow
#                     print(species.strip())
                    
        
#     return rxn_dict