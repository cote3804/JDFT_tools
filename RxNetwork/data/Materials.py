
from data.Data_Parser import Data_Parser

# Material class to store data with each material

class Material:
    def __init__(self, bulk_name, data: Data_Parser):
        self.bulk_name = bulk_name
        self.surfaces = data.get_surfaces(self.bulk_name)
        self.energies = {}

    def get_FED_energy(self, surface:str, bias:str, referenced:str):
        if referenced == "final":
            reference_energy = self.energies[surface]["final"][bias]
        elif referenced == "initial":
            reference_energy = self.energies[surface]["initial"][bias]
        FED_energies = {}
        surface_energy = self.energies[surface]
        for intermediate in surface_energy.keys():
            if intermediate not in ["initial", "final"]:
                FED_energy = min([E for (s,E) in surface_energy[intermediate][bias].items()])
                FED_energies[intermediate] = FED_energy - reference_energy
            elif intermediate in ["initial", "final"]:
                continue
        FED_energies["initial"] = surface_energy["initial"][bias] - reference_energy
        FED_energies["final"] = surface_energy["final"][bias] - reference_energy
        return FED_energies
    # @property
    # def surfaces(self) -> list:
    #     return self.surfaces
    
    # @surfaces.setter
    # def surfaces(self, surfaces:list) -> None:
    #     self.surfaces = surfaces

    @property
    def surface_data(self) -> dict:
        return self.surface_data
    
    @surface_data.setter
    def surface_data(self, surface_data:dict) -> None:
        self.surface_data = surface_data
        
    def get_surface_data(self, surface):
        return self.energies[surface]


class Materials:
    def __init__(self, data:Data_Parser, bulks:list) -> None:
        self.bulks = bulks
        self.materials = [Material(bulk, data) for bulk in bulks] # creates list of Material objects

    # @property
    # def bulks(self) -> list:
    #     return self.bulks
    
    # @bulks.setter
    # def bulks(self, bulks:list) -> None:
    #     self.bulks = bulks

    def get_surface_data(self, bulk:str) -> dict:
        surface_data = {}
        for material in self.materials:
            surface_data[material.bulk_name] = material.surface_data
        return surface_data
    
    def get_surfaces(self) -> list:
        surfaces = []
        for material in self.materials:
            surfaces.extend(material.surfaces)
        return surfaces
    
    def get_material(self, bulk:str) -> Material:
        for material in self.materials:
            if material.bulk_name == bulk:
                return material
            else:
                raise Exception(f"material {bulk} not previously initialized")

    def get_bias_energy(self, surface:str, bulk:str):
        bulk_name = surface.split("_")[0]
        material = self.get_material(bulk_name)
        


    

        
        