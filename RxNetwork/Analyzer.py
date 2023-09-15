#!/usr/bin/envs python3

from data.Data_Parser import Data_Parser
from data.Materials import Materials
from data.Calculator import Calculator

class Analyzer:
    def __init__(self, data_path:str, filename:str, bulks:list) -> None:
        self.bulks = bulks
        self.data = Data_Parser(data_path, filename)

    def calculate_surface_properties(self, materials:Materials) -> None: # this method needs to return a dictionary with all the necessary values for analysis
        calculator = Calculator(self.data)
        # surface_data = materials.get_surface_data()
        calculator.calculate_surface_properties(materials)

    def do_analysis(self) -> None:
        materials = Materials(self.data, self.bulks)

