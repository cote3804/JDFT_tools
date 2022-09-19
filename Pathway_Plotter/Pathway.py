# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 12:50:59 2022

@author: coopy
"""
#This script provides the functionality to parse through the all_data.json 
#file generated by gc_manager and plot desired energy pathways using 
#matplotlib

###############################################################################
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import json
import os
import re
H_to_ev = 27.2114

###############################################################################
# with open(os.path.join(path,json_file),'r') as f:
#     all_data = json.load(f)

class Pathway:
    def __init__(self, path=None, filename=None, reaction='N2R'):
        if path is None:
            self.path = os.getcwd()
        else:
            self.path = path
            
        if filename is None:
            self.filename = 'all_data.json'
        else:
            self.filename = filename
        self.pathways = {}
        self.count = 0
        self.reaction = reaction
        self.reaction_data = {
                            'N2R': 
                          {'pathway' : ['N2','N2H','NNH2','NNH3','N','NH','NH2','NH3'],
                          'protons' : [6, 5, 4, 3, 3, 2, 1, 0],
                          'molecules' : ['','','','','NH3','NH3','NH3','NH3'],
                          'entropy' : {'N2': 0.021720172, 'NH3': 0.021830347, 
                                       'H3O': 0.021892558, 'H2O': 0.021404843},
                          'ZPE' : {'N2': 0.005677, 'NH3': 0.034492, 
                                       'H3O': 0.035438, 'H2O': 0.021586 }
                          },
            
                            'HER':
                          {'pathway' : ['H'],
                           'protons' : [1],
                           'molecules' : [''],
                           'entropy' : {'H3O': 0.035438, 'H2O': 0.021586,
                                        'H2': 0.014789719},
                           'ZPE' : {'H3O': 0.035438, 'H2O': 0.021586,
                                    'H2': 0.010088}
                           }
                             
            
            
            }
        self.potentials = {'0.00': -0.1713, '-0.25': -0.1621, '-0.50': -0.1529}
        
    def get_json(self):
        '''
        This method will return the json as a dictionary to be investigated in 
        the IDE variable pane
        '''
        with open(os.path.join(self.path,self.filename),'r') as f:
            all_data = json.load(f)
        return all_data
    
    def get_energy(self, material, adsorbates, bias, site, ev=True, free_energy=True, proton='smart'):
        with open(os.path.join(self.path,self.filename),'r') as f:
            all_data = json.load(f)
        if material in all_data.keys():
            pass
        else:
            raise Exception('Material must match the materials listed in the all_data.json file. '
                            f'{material} was entered which does not match the json')
        
        energy= []
        step = []
        for intermediate in adsorbates: 
            if intermediate in all_data[material]['adsorbed'].keys() and all_data[material]['adsorbed'][intermediate][f'{bias}V'][f'0{site}']['converged']:
                energy.append(all_data[material]['adsorbed'][intermediate][f'{bias}V'][f'0{site}']['final_energy'])
                step.append(f'{intermediate}*')
            else:
                raise Exception(f'intermediate {intermediate}* at bias {bias}V ' 
                                f'and site 0{site} on {material} not in all_data.json. '
                                f'Check to see if it\'s converged')
        energy, step = self.reaction_references(all_data, energy, step, material, 
                                                adsorbates, bias, site, free_energy, proton)        
        if ev:
            energy = list(np.array(energy)*27.2114)
        #Set reference point to final state
        energy = list(np.array(energy)-energy[-1])
        return energy, step
    
    def reaction_references(self, all_data, energies, steps, material, adsorbates, 
                            bias, site, free_energy=True, proton='smart'):
        #start by building initial states for pathway with clean surface and 
        #infinite separation of molecules and protons
        
        #calculated from H3O+ - H2O. Calculation can be found in Cooper's projects/
        #directory /projects/cote3804/proton_energy
        if proton == 'smart': #use solution proton reference
            E_proton = -17.7048336646746769 - (-17.2809934718248783) 
        elif proton == 'stupid':
            E_proton = (1/2) * all_data['H2'][f'{bias}V']['final_energy']
        S_proton = 0.000487715
        ZPE_proton = 0.013852
        
        if self.reaction == 'N2R':
            initial_mol = 'N2'
            final_mol = 'NH3'
            #set up initial and final steps in pathway which for N2R is surface + N2
            #and surface + NH3
        # TODO: add other reaction pathways into reference generation logic
        elif self.reaction == 'CO2R':
            pass  
        elif self.reaction == 'HER':
            initial_mol = None
            final_mol = 'H2'
        elif self.reaction == 'OER':
            pass
        elif self.reaction == 'ORR':
            pass
        
        clean_surf_energy = all_data[material]['surf'][f'{bias}V']['final_energy']
        #HER pathway does not have a starting reference molecule so the logic must
        #be able to handle an initial_mol value of None
        if self.reaction == 'HER':
            initial_mol_energy = 2*E_proton
        else:
            initial_mol_energy = all_data[initial_mol][f'{bias}V']['final_energy']
        final_mol_energy = all_data[final_mol][f'{bias}V']['final_energy']
        
        #Build list of added molecules/proton energies to be added to list of 
        #adsorption energies
        #===================================================================#
        full_pathway = self.reaction_data[self.reaction]['pathway']
        proton_number = self.reaction_data[self.reaction]['protons']
        #molecules includes the list of molecules necessary to include to make
        #sure the number of atoms is constant across the pathway
        molecules = self.reaction_data[self.reaction]['molecules']
        entropy = self.reaction_data[self.reaction]['entropy']
        ZPE = self.reaction_data[self.reaction]['ZPE']
        i = 0
        for adsorbate, energy in zip(steps, energies):
            try:
                path_index = full_pathway.index(adsorbate.split('*')[0])
            except ValueError:
                print('References could not be added because the provided reaction '
                      'intermediates are not in agreement with the currently '
                      'supported intermediates')
            adsorbate_energy = energy
            proton_energy = proton_number[path_index] * E_proton
            if free_energy:
                proton_energy = proton_energy + proton_number[path_index]*(ZPE_proton - S_proton)
                #ZPE and Entropy 
                # proton_energy = proton_energy + proton_number[path_index]*(-S_proton)
            molecule = molecules[path_index]
            if molecule == '':
                molecule_energy = 0
            else:
                molecule_energy = all_data[molecule][f'{bias}V']['final_energy']
            if free_energy and molecule in entropy.keys():
                # molecule_energy = molecule_energy + ZPE[molecule] - entropy[molecule]
                #ZPE and Entropy
                molecule_energy = molecule_energy - entropy[molecule]
                #Entropy alone
                
            energies[i] = adsorbate_energy + proton_energy + molecule_energy
            i+=1
            
        #Add initial and final reference states to steps and energies lists 
        #=================================================================# 
        if self.reaction == 'N2R':
            initial_total_energy = clean_surf_energy + initial_mol_energy + proton_number[0]*E_proton
        elif self.reaction == 'HER':
            initial_total_energy = clean_surf_energy + initial_mol_energy
        if self.reaction == 'N2R':
            final_total_energy = clean_surf_energy + 2*final_mol_energy + proton_number[-1]*E_proton
        elif self.reaction == 'HER':
            final_total_energy = clean_surf_energy + final_mol_energy
        if free_energy:
            
            
            # final_total_energy = final_total_energy + 2*(ZPE['NH3'] - entropy['NH3']) + proton_number[-1]*(ZPE_proton - S_proton)
            # initial_total_energy = initial_total_energy + 2*(ZPE['N2'] - entropy['N2']) + proton_number[0]*(ZPE_proton - S_proton)
            #ZPE and Entropy
            
            
            if self.reaction == 'N2R':
                final_total_energy = final_total_energy + 2*(-entropy['NH3']) + proton_number[-1]*(-S_proton)
                initial_total_energy = initial_total_energy + 2*(-entropy['N2']) + proton_number[0]*(-S_proton)
            elif self.reaction == 'HER':
                #TODO fix the entropies for other reactions
                final_total_energy = final_total_energy + (-entropy['H2']) + ZPE['H2']
                initial_total_energy = initial_total_energy + (-entropy['H2']) -2 *S_proton +2*ZPE_proton
            #Entropy only
        
        #TODO: Add logic to specify 2 final NH3 molecules instead of hard coding above
        #On second thought, hard coding each pathway might be better
        steps.insert(0, initial_mol)
        if self.reaction == 'N2R':
            steps.append(f'2{final_mol}')
        else:
            steps.append(final_mol)
        energies.insert(0, initial_total_energy)
        energies.append(final_total_energy)
        
        return energies, steps
        
    def add_pathway(self, material, adsorbates, bias, site, free_energy=True, proton='smart'):
        self.pathways[self.count] = [material, adsorbates, bias, site, free_energy, proton]
        self.count += 1
    
    def pathway_plot(self, ylim=None, legend=True, alpha=1, annotation=True, custom_colors=None):
        plt.figure(dpi=300)
        plt.ylabel('$\Delta \Phi (eV)$')
        plt.xticks([])
        plt.xlabel('Reaction Coordinate')
        # colors = ['#CE3B7E','#E6552E','#9E50F4','#3AC2BC','#ECB42E','#D1DB49']
        if custom_colors is None:
            colors = ['#390099','#9e0059','#ff0054','#ff5400','#ffbd00','#D1DB49']
        else:
            colors = custom_colors
        labels = []
        legend_entries = []

        for ipath in range(len(self.pathways.keys())):
            path_list = self.pathways[ipath]
            color = colors[ipath]
            labels.append(f'{path_list[2]}V')
            energies, steps = self.get_energy(path_list[0],path_list[1],path_list[2],
                                              path_list[3],free_energy=path_list[4],proton=path_list[5])
            for i,energy in enumerate(energies):
                plt.hlines(energy, i, i+0.5, color=color, label=labels[ipath],
                           alpha=alpha)
                if i > 0:
                    plt.plot([i,i-0.5], [energies[i],energies[i-1]], linestyle='--', color=color, label=labels[ipath],
                             alpha=alpha)
                if ipath == len(self.pathways.keys())-1 and annotation: #only annotate last pathway
                    plt.annotate(steps[i], ((i+0),energy+0.1), color=color)
            #Create legend entry objects 
            free_energy_term = 'free energy' if path_list[4] else 'enthalpy'
            legend_entries.append(Line2D([0], [0], 
                                         color=colors[ipath], lw=2,
                                         label=f'{path_list[0]} {labels[ipath]} {free_energy_term}'))
        if legend:
            plt.legend(handles=legend_entries, fontsize=6)
        else:
            pass
        
        if ylim is None:
            pass
        else:
            plt.ylim(ylim)
        plt.savefig(os.path.join(self.path,'N2R_Pathway_0V'))
        
    def energies(self,ev=True):
        energies_dict = {}
        for ipath in range(len(self.pathways.keys())):
            path_list = self.pathways[ipath]
            energies, steps = self.get_energy(path_list[0],path_list[1],path_list[2],path_list[3],
                                              free_energy=path_list[4],ev=ev,proton=path_list[5])
            energy_dict = {}
            for i,energy in enumerate(energies):
                energy_dict[steps[i]] = energy
            energies_dict[path_list[0]] = energy_dict
        return energies_dict
    
    def adatom_energy(self, all_data, material, site, bias='0.00', ev=True):
        H_to_ev = 27.2114
        E_proton = -17.7048336646746769 - (-17.2809934718248783)
        if self.reaction == 'N2R':
            adatom= 'N'
        else:
            pass
        #TODO add more reaction support
        molecular_energy = (all_data['N2'][f'{bias}V']['final_energy'] #N2 molecule
                            +3*E_proton) #3 protons
        adsorbate_energy = (all_data[material]['adsorbed'][adatom][f'{bias}V'][f'0{site}']['final_energy'] #adsorbate complex
                            + all_data['NH3'][f'{bias}V']['final_energy'])
        surface_energy = all_data[material]['surf'][f'{bias}V']['final_energy']
        adatom_energy = adsorbate_energy - (surface_energy + molecular_energy)
        
        if ev:
            adatom_energy = adatom_energy * H_to_ev
        return adatom_energy
    
    
if __name__ == '__main__':
    
    json_file = 'all_data.json'
    path = os.path.normpath('C:\\Users\\coopy\\OneDrive - UCB-O365\\Research\\N2R_Scaling')
    pathway = Pathway(path=path)
    all_data = pathway.get_json()
    pathway.add_pathway('Re_111',['N2','N2H','NH3'],'0.00',1,free_energy=False)
    pathway.add_pathway('Ru_111',['N2','N2H','NH3'],'0.00',1,free_energy=False)
    pathway.add_pathway('Ir_111',['N2','N2H','NH3'],'0.00',1,free_energy=False)
    pathway.add_pathway('Co_111_mag',['N2','N2H','NH3'],'0.00',1,free_energy=False)
    # pathway.add_pathway('F_110',['N2','N2H','NNH2','NNH3','N','NH','NH2','NH3'],'0.00',1,free_energy=False)
    
    energies = pathway.energies()
    pathway.pathway_plot()
