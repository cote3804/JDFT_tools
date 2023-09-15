from Data_Parser import Data_Parser

class Calculator:
    def __init__(self, data:Data_Parser, reaction:str) -> None:
        self.data = data
        self.reaction = Reaction(reaction)
    
    def calculate_surface_properties(self, materials) -> None:
        surface_properties = {}
        for surface in materials.get_surfaces():

        surfaces_data = materials.get_surface_data()
        for surface_name, surface_dict in surfaces_data.items():
            for bias_str, bias_name in surface_dict["adsorbed"].items
                for intermediate in self.reaction.intermediates:
                    surface_properties.update

class Reaction:
    def __init__(self, reaction:str) -> None:
        intermediate_dict = {
            "NRR": ["N2", "N2H", "NNH2", "N", "NH", "NH2", "NH3"],
            "HER": ["H"]
        },
        self.references = References(reaction)
        self.intermediates = intermediate_dict[reaction]

class References:
    def __init__(self, reaction:str) -> None:
        self.reaction = reaction
    
    def get_intermediate_references(self, intermediate):

class SurfaceProperties:
    def __init__(self):
        self.surface_properties = {}
    
    def build_empty_struct(self, materials):
        
    
