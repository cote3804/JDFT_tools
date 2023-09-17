from data.Data_Parser import Data_Parser
from data.Materials import Material, Materials

class Calculator:
    def __init__(self, data:Data_Parser, reaction:str) -> None:
        self.data = data
        self.reaction = Reaction(reaction)
    
    def calculate_surface_properties(self, materials:Materials) -> None:
        surface_properties = {}
        for material in materials.materials:
            self.set_material_energies(material)
            
    def intermediate_energy(self, surface:str, intermediate:str, bias:str, site:str) -> float:
        adsorbed_energy = self.data.get_adsorbed_energy(surface, bias, intermediate)
        reference_energy = self.reaction.reference_energy(intermediate, bias)
        calculated_energy = adsorbed_energy + reference_energy
        return calculated_energy
    
    def set_material_energies(self, material: Material) -> None:
        energies = {}
        surfaces = material.surfaces
        for surface in surfaces:
            energies[surface] = {}
            for intermediate in self.reaction.find_reaction_intermediates(self.data, surface):
                bias_energies = {}
                for bias in self.data.get_intermediate_biases(surface, intermediate):
                    bias_energies[bias] = {}
                    for site in self.data.get_sites_data(surface, bias, intermediate).keys():
                        calculated_energy = self.intermediate_energy(surface, intermediate, bias, site)
                        bias_energies[bias].update({site: calculated_energy})
                    energies[surface].update({intermediate: bias_energies})
        material.energies = energies

#TODO Don't need to define list of intermediates in both the References and Reaction
# objects

class Reaction:
    def __init__(self, reaction:str) -> None:
        intermediate_dict = {"NRR": ["N2", "N2H", "NNH2", "N", "NH", "NH2", "NH3"],
                             "HER": ["H"]}
        self.references = References(reaction)
        self.intermediates = intermediate_dict[reaction]
    
    def reference_energy(self, intermediate, bias):
        reference_energy = self.references.get_intermediate_references(intermediate, bias)
        return reference_energy
    
    def find_reaction_intermediates(self, data:Data_Parser, surface) -> list:
        '''
        Compares reaction intermediates in the all_data to reaction intermediates
        defined in the reaction
        '''
        
        data_intermediates = data.get_intermediates(surface)
        return [i for i in data_intermediates if i in self.intermediates]

class References:
    def __init__(self, reaction:str) -> None:
        self.reference_energies = {"N2":{"G": 10, "nelec":14}, "NH3":{"G":8, "nelec":10}, 
                                    "H+":{"G":1, "nelec":0}}
        self.reference_table = {"N2": {"H+":6}, "N2H": {"H+":5}, "NNH2": {"H+":4},
                           "N": {"NH3":1, "H+":3}, "NH": {"NH3":1, "H+":2},
                           "NH2": {"NH3":1, "H+":1}, "NH3":{"NH3":1}}
        self.reaction_terminations = {"NRR":{"initial":{"N2":1, "H+":6}, "final":{"NH3":2}}
            }
        self.reaction = reaction
        
    def get_references(self, intermediate):
        pass
        
    
    def get_intermediate_references(self, intermediate, bias):
        reference_energy = 0
        for molecule, quantity in self.reference_table[intermediate].items():
            # print(self.reference_energies[molecule]["G"], self.reference_energies[molecule]["nelec"])
            reference_energy += ((self.reference_energies[molecule]["G"]
                                + self.reference_energies[molecule]["nelec"] * self.bias_string_to_float(bias))
                                 ) * quantity
        return reference_energy
        
    def bias_string_to_float(self, bias):
        return float(bias.split('V')[0])
    
