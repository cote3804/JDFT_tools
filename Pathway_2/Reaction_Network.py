#%% Main Loop
import os
import json
import numpy as np
import tkinter as tk

# TODO make a bunch of methods to handle getting and setting of rxn_data parameters
class Reaction_Network:
    def __init__(self):
        with open('./network_defaults.json', 'r') as f:
            self.network_data = json.load(f)
        
    def load_network(self, reaction:str):
        self.network_dict = self.network_data[reaction]
        self.reaction = reaction
        self.dimensions = self.get_dimensions()
        self.reaction_mask = self.get_reaction_mask()
        self.lookup_matrix = self.get_lookup_matrix() #indexed to 1 so that zeros represent invalid states

    def gui_network(self):
        pass # TODO need to make this whole gui network addition work

    def get_dimensions(self) -> tuple:
        """
        returns a tuple of the network dimensions (rows, columns)
        """
        network_dimensions = {}
        for node in self.network_data[self.reaction]['network']:
            pathways = node["network_data"]["pathways"]
            for pathway in pathways:
                # dictionary method .get(foo, val) checks if the key foo exitsts, and if it doesn't,
                # it creates the key foo and stores val in it. This is a way to guarantee that by stumbling
                # upon a new key, you won't throw an error but instead add the key with a default value.
                # In this case, val is 1 because it has counted one occurence of a new reaction pathway.
                network_dimensions[pathway] = network_dimensions.get(pathway, 0) + 1
        network_dimensions = (len(network_dimensions),max(network_dimensions.values()))
        return network_dimensions
    
    def get_reaction_mask(self):
        pathways = self.network_dict["pathways"]["labels"]
        reaction_mask = np.full((self.dimensions), False)
        for node in self.network_dict["network"]:
            for pathway in node["network_data"]["pathways"]:
                row_index = pathways.index(pathway)
                column_index = node["network_data"]["index"]
                reaction_mask[row_index, column_index] = True
        return  reaction_mask
    
    def get_lookup_matrix(self) -> dict:
        """
        This function is used to create an array whose indices correspond to the indices specified in network_dict.
        The column index is specified in network_dict directly with each intermediate and the row index corresponds to the
        pathways specified in network_dict. The ordering in which the pathways were specified is preserved here.
        Example: a reaction network with 2 pathways each with 4 intermediates would produce a 2x4 array, where each entry in
        the array references the index of the corresponding intermediate dictionary in the list of dictionaries in 
        network_dict["network"].
        Also it's indexed to 1 so that zeros can signify an invalid intermediate that only exists for some pathways but not all
        Why did I structure things this way? Because I'm bad at coding. My failures have now become your problem :) -Coop
        """
        pathways = self.network_dict["pathways"]["labels"]
        lookup_matrix = np.zeros(self.dimensions)
        for inode, node in enumerate(self.network_dict["network"]):
            for pathway in node["network_data"]["pathways"]:
                row_index = pathways.index(pathway)
                column_index = node["network_data"]["index"]
                lookup_matrix[row_index, column_index] = inode + 1
        return lookup_matrix

    def indices_to_state(self, indices:tuple):
        if self.lookup_matrix[indices] == 0:
            raise Exception(f"Invalid state with indices {indices}")
        state_data = self.network_dict["network"][int(self.lookup_matrix[indices]-1)]
        return state_data["network_data"]["label"]
    
    def index_to_pathway_name(self, index):
        pathway_name = self.network_dict["pathways"]["names"][index]
        return pathway_name

    def get_network_dict(self):
        return self.network_dict
    
    def add_network(self, network_dict):
        self.network_data.update(network_dict)

if __name__ == "__main__":
    test = Reaction_Network()
    test.load_network("NRR")
    state = test.indices_to_state((0,7))
    network_dict = test.get_network_dict()

'''
self.rxn_data = {"NRR": {
                            "pathways": {"labels": ["A", "B", "C"],
                                         "names": ["assoc. dist.", "assoc. alt.", "diss"]},
                            "network": [
                                {
                                "energy": {"molecules": {"H+": 6, "N2": 1},
                                          "surfaces": {"*": 1}
                                              },
                                "network_data": {"pathways": ["A","B"],
                                            "index": 0,
                                            "label": "N2"
                                            }
                                    
                                },
                                {
                                "energy": {"molecules": {"H+": 6, "N2": 1},
                                          "surfaces": {"*": 2}
                                              },
                                "network_data": {"pathways": ["C"],
                                            "index": 0,
                                            "label": "N2"
                                            }
                                },
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
                                "energy": {"molecules": {"H+": 6},
                                          "surfaces": {"N*": 2}
                                              },
                                "network_data": {"pathways": ["C"],
                                            "index": 1,
                                            "label": "2N*"
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
                                "energy": {"molecules": {"H+": 4},
                                          "surfaces": {"NH*": 2}
                                              },
                                "network_data": {"pathways": ["C"],
                                            "index": 3,
                                            "label": "2NH*"
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
                                "energy": {"molecules": {"H+": 2},
                                          "surfaces": {"NH2*": 2}
                                              },
                                "network_data": {"pathways": ["C"],
                                            "index": 5,
                                            "label": "2NH2*"
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
                                },
                                {
                                "energy": {"molecules": {"H+":0},
                                          "surfaces": {"NH3*":2}
                                              },
                                "network_data": {"pathways": ["C"],
                                            "index": 7,
                                            "label": "2NH3*"
                                            }
                                },
                                {
                                    "energy": {"molecules": {"H+":0, "NH3":2},
                                              "surfaces": {"*": 1}
                                                  },
                                    "network_data": {"pathways": ["A","B"],
                                                "index": 8,
                                                "label": "2NH3"
                                                } 
                                },
                                {
                                "energy": {"molecules": {"H+":0, "NH3":2},
                                          "surfaces": {"*":2}
                                              },
                                "network_data": {"pathways": ["C"],
                                            "index": 8,
                                            "label": "2NH3*"
                                            }
                                }
                                        ]
            }
            }

'''