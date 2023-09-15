# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 17:55:46 2023

@author: coopy
"""
import json
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

#TODO currently the energy matrix is stored as an environment variable unreferenced 
# and in Hartrees. I could consider using a method to convert the matrix every time
# I want it referenced and in eV

#TODO get dissociative working

class reaction_network:
    def __init__(self, path, file, reaction):
        self.rxn_data = {"NRR": {
                            "initial state": {"molecules": {"H+": 6, "N2": 1},
                                              "surfaces": {"*"},
                                              "label": "N2"},
                            "final state": {"molecules": {"NH3": 2},
                                            "surfaces": {"*"},
                                            "label": "2NH3"},
                            "pathways": {"labels": ["A", "B"],
                                         "Names": ["assoc. dist.", "assoc. alt."]},
                            "network": [
                                {
                                "energy": {"molecules": {"H+": 6},
                                          "surfaces": {"N2*": 1}
                                              },
                                "network_data": {"pathways": ["A","B"],
                                            "index": 1,
                                            "label": "N2*"
                                            }
                                },
                                {
                                "energy": {"molecules": {"H+":5},
                                          "surfaces": {"N2H*": 1}
                                              },
                                "network_data": {"pathways": ["A","B"],
                                            "index": 2,
                                            "label": "N2H*"
                                            }
                                },
                                {
                                "energy": {"molecules": {"H+":4},
                                          "surfaces": {"NNH2*": 1}
                                              },
                                "network_data": {"pathways": ["A"],
                                            "index": 3,
                                            "label": "NNH2*"
                                            }
                                },
                                {
                                "energy": {"molecules": {"H+":4},
                                          "surfaces": {"NHNH*": 1}
                                              },
                                "network_data": {"pathways": ["B"],
                                            "index": 3,
                                            "label": "NHNH*"
                                            }
                                },
                                {
                                "energy": {"molecules": {"H+":3, "NH3":1},
                                          "surfaces": {"N*": 1}
                                              },
                                "network_data": {"pathways": ["A"],
                                            "index": 4,
                                            "label": "N*"
                                            }
                                },
                                {
                                "energy": {"molecules": {"H+":3},
                                          "surfaces": {"NHNH2*": 1}
                                              },
                                "network_data": {"pathways": ["B"],
                                            "index": 4,
                                            "label": "NHNH2*"
                                            }
                                },
                                {
                                "energy": {"molecules": {"H+":2, "NH3":1},
                                          "surfaces": {"NH*": 1}
                                              },
                                "network_data": {"pathways": ["A"],
                                            "index": 5,
                                            "label": "NH*"
                                            }
                                },
                                {
                                   "energy": {"molecules": {"H+":2},
                                             "surfaces": {"NH2NH2*": 1}
                                                 },
                                   "network_data": {"pathways": ["B"],
                                               "index": 5,
                                               "label": "NH2NH2*"
                                               } 
                                },
                                {
                                    "energy": {"molecules": {"H+":1, "NH3":1},
                                              "surfaces": {"NH2*": 1}
                                                  },
                                    "network_data": {"pathways": ["A","B"],
                                                "index": 6,
                                                "label": "NH2*"
                                                } 
                                },
                                {
                                    "energy": {"molecules": {"H+":0, "NH3":1},
                                              "surfaces": {"NH3*": 1}
                                                  },
                                    "network_data": {"pathways": ["A","B"],
                                                "index": 7,
                                                "label": "NH3*"
                                                } 
                                }
                                        ]
            },
            'NRR_test': {
                                "initial state": {"molecules": {"H+": 6, "N2": 1},
                                                  "surfaces": {"*"},
                                                  "label": "N2"},
                                "final state": {"molecules": {"NH3": 2},
                                                "surfaces": {"*"},
                                                "label": "2NH3"},
                                "pathways": {"labels": ["A", "B", "C"],
                                             "Names": ["assoc. dist.", "assoc. alt.", "diss"]},
                                "network": [
                                    {
                                    "energy": {"molecules": {"H+": 6},
                                              "surfaces": {"N2*": 1}
                                                  },
                                    "network_data": {"pathways": ["A","B"],
                                                "index": 1,
                                                "label": "N2*"
                                                }
                                    },
                                    {
                                    "energy": {"molecules": {"H+":5},
                                              "surfaces": {"N2H*": 1}
                                                  },
                                    "network_data": {"pathways": ["A","B"],
                                                "index": 2,
                                                "label": "N2H*"
                                                }
                                    },
                                    {
                                    "energy": {"molecules": {"H+":4},
                                              "surfaces": {"NNH2*": 1}
                                                  },
                                    "network_data": {"pathways": ["A"],
                                                "index": 3,
                                                "label": "NNH2*"
                                                }
                                    },
                                    {
                                    "energy": {"molecules": {"H+":4},
                                              "surfaces": {"NHNH*": 1}
                                                  },
                                    "network_data": {"pathways": ["B"],
                                                "index": 3,
                                                "label": "NHNH*"
                                                }
                                    },
                                    {
                                    "energy": {"molecules": {"H+":3, "NH3":1},
                                              "surfaces": {"N*": 1}
                                                  },
                                    "network_data": {"pathways": ["A"],
                                                "index": 4,
                                                "label": "N*"
                                                }
                                    },
                                    {
                                    "energy": {"molecules": {"H+":3},
                                              "surfaces": {"NHNH2*": 1}
                                                  },
                                    "network_data": {"pathways": ["B"],
                                                "index": 4,
                                                "label": "NHNH2*"
                                                }
                                    },
                                    {
                                    "energy": {"molecules": {"H+":2, "NH3":1},
                                              "surfaces": {"NH*": 1}
                                                  },
                                    "network_data": {"pathways": ["A"],
                                                "index": 5,
                                                "label": "NH*"
                                                }
                                    },
                                    {
                                       "energy": {"molecules": {"H+":2},
                                                 "surfaces": {"NH2NH2*": 1}
                                                     },
                                       "network_data": {"pathways": ["B"],
                                                   "index": 5,
                                                   "label": "NH2NH2*"
                                                   } 
                                    },
                                    {
                                        "energy": {"molecules": {"H+":1, "NH3":1},
                                                  "surfaces": {"NH2*": 1}
                                                      },
                                        "network_data": {"pathways": ["A","B"],
                                                    "index": 6,
                                                    "label": "NH2*"
                                                    } 
                                    },
                                    {
                                        "energy": {"molecules": {"H+":0, "NH3":1},
                                                  "surfaces": {"NH3*": 1}
                                                      },
                                        "network_data": {"pathways": ["A","B"],
                                                    "index": 7,
                                                    "label": "NH3*"
                                                    } 
                                    }
                                            ]
                }
            }
        
        self.reaction = reaction
        self.materials = []
        self.H_to_eV = 27.2114
        
        with open(os.path.join(path,file)) as f:
            self.all_data = json.load(f)
    
    def network_to_thermo(self, material, bias, site):
        rxn_data = self.rxn_data
        reaction = self.reaction
        pathways = self.rxn_data[reaction]['pathways']["labels"]
        network_dimensions = self.rxn_network_dimensions()
        thermo_matrix = np.zeros((len(network_dimensions), max(network_dimensions.values())))
        
        thermo_matrix[:,0] = self.calculate_end_state('initial', material, bias, site)
        for node in rxn_data[reaction]["network"]:
            for pathway in node["network_data"]["pathways"]:
                row_index = pathways.index(pathway)
                column_index = node["network_data"]["index"]
                thermo_matrix[row_index, column_index] = self.calculate_intermediate(node["energy"], material, bias, site)
        thermo_matrix[:,-1] = self.calculate_end_state('final', material, bias, site)
        return thermo_matrix

    def calculate_end_state(self, end, material, bias, site):
        state_string = f"{end} state"
        state_data = self.rxn_data[self.reaction][state_string]
        E_surf = 0
        E_mol = 0
        for surface in state_data['surfaces']:
            E_surf += self.calculate_surface(surface, material, bias, site)
        for mol in state_data['molecules']:
            E_mol += self.calculate_mol(mol, bias) * state_data['molecules'][mol]
        
        return E_mol + E_surf
    
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
                    min_site = site #TODO maybe return site number purposes
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
        for node in self.rxn_data[self.reaction]['network']:
            pathways = node["network_data"]["pathways"]
            for pathway in pathways:
                network_dimensions[pathway] = network_dimensions.get(pathway, 2) + 1
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
        for material_dict in self.materials:
            thermo_matrix = self.network_to_thermo(material_dict["material"], material_dict["bias"], material_dict["site"])
            matrix_list.append(thermo_matrix)
        energy_array = np.stack(matrix_list, axis=2)
        self.energy_array = energy_array
        
        if referenced == "final":
            subtraction_matrix = energy_array[:,-1:,:]
            energy_array = energy_array - subtraction_matrix 
        if ev == True:
            energy_array = energy_array * self.H_to_eV
        return energy_array
    
    def FED_plot(self):
        #FED stands for free energy diagram
        
        energy_array = self.calculate_energies()
        
        ################# Plot Formatting ######################
        starting_point = 0.5
        number_of_states = len(self.energy_array[0,:,0])
        state_length = 1
        connection_length = 0.25
        
        fig,ax = plt.subplots(dpi=300)
        lines = []
        linestyles = []
        colors = ['k','b','r']
        
        for imat, material_dict in enumerate(self.materials):
            material = material_dict["material"]
            mat_matrix = energy_array[:,:,imat]
            for pathway_vector in mat_matrix:
                for istate, state in enumerate(pathway_vector):
                    # state is an energy value
                    if istate == 0:
                        lines.append([(starting_point, state),
                                      (starting_point+state_length, state)])
                        linestyles.append("solid")
                    else:
                        #need to draw connection from previous state and new state horizontal line
                        #draw connection first
                        lines.append([(starting_point + istate*state_length + (istate-1)*connection_length, pathway_vector[istate - 1]),
                                      (starting_point + istate*state_length + (istate)*connection_length, state)])
                        linestyles.append("dashed")
                        
                        #now draw current state
                        lines.append([(starting_point + istate*state_length + (istate)*connection_length, state),
                                      (starting_point + (istate+1)*state_length + (istate)*connection_length, state)])
                        linestyles.append("solid")
        # print(lines)
        line_segments = LineCollection(lines, linewidths = (1), linestyle=linestyles, colors=colors)
        # print(line_segments)
        ax.add_collection(line_segments)
        ax.set_title("plotting test")
        ax.set_xlim(0, mat_matrix.shape[-1]*(state_length+connection_length) + starting_point + 1)
        # ax.set_ylim(-2511.9, -2511.6)
        ax.set_ylim(-4,4)
        plt.show()
                   
    def calculate_energetic_span():
        pass
                        
        
    
if __name__ == "__main__":
    file = 'test_Combined_all_data.json'
    path = os.path.normpath('C:\\Users\\coopy\\OneDrive - UCB-O365\\Research\\N2R_Scaling\\all_data\\Paper all_data')
    ntwrk = reaction_network(path, file, "NRR")
    # ntwrk.add_material("Ru_111", "0.00V", "min")
    # ntwrk.add_material("Re_111", "0.00V", "min")
    ntwrk.add_material("Rh_111", "0.00V", "min")
    # array = ntwrk.calculate_energies()
    ntwrk.FED_plot()