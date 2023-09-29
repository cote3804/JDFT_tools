from data.Data_Parser import Data_Parser
from data.Materials import Materials, Material

class Calculator:
    def __init__(self, data:Data_Parser, reaction:str, ev=True) -> None:
        self.data = data
        self.reaction = Reaction(reaction)
        if ev == True:
            self.ev = 27.2114 #H to eV conversion
        elif ev == False:
            self.ev = 1
    
    def calculate_surface_properties(self, materials) -> None:
        surface_properties = {}
        for material in materials.materials:
            self.set_material_energies(material)
            
    def intermediate_energy(self, surface:str, intermediate:str, bias:str, site:str) -> float:
        converged = self.data.check_adsorbed_site_convergence(surface, bias, intermediate, site)
        if converged == False:
            return None
        elif converged == True:
            adsorbed_energy = self.data.get_adsorbed_energy(surface, bias, intermediate)
            reference_energy = self.reaction.reference_energy(intermediate, bias, self.ev)
            calculated_energy = adsorbed_energy + reference_energy
            # print(calculated_energy, adsorbed_energy, reference_energy, surface, intermediate, bias, site)
            return calculated_energy * self.ev
    
    def terminal_energies(self, surface:str, bias:str) -> (float, float):
        surface_energy = self.data.get_surface_energy(surface, bias) * self.ev
        initial_energy = self.reaction.reference_energy("initial", bias, self.ev) * self.ev
        # print(bias, initial_energy)
        final_energy = self.reaction.reference_energy("final", bias, self.ev) * self.ev
        # print(final_energy)
        return (surface_energy + initial_energy), (surface_energy + final_energy)

    def set_material_energies(self, material) -> None:
        energies = {}
        surfaces = material.surfaces
        for surface in surfaces:
            if self.data.check_any_surface_convergence(surface) == False:
                print(f"{surface} has no converged sites")
                continue
            elif self.data.check_any_surface_convergence(surface) == True:
                energies[surface] = {}
                energies[surface]["initial"] = {}
                energies[surface]["final"] = {}
                for intermediate in self.reaction.find_reaction_intermediates(self.data, surface):
                    bias_energies = {}
                    for bias in self.data.get_intermediate_biases(surface, intermediate):
                        bias_energies[bias] = {}
                        for site in self.data.get_sites_data(surface, bias, intermediate).keys():
                            calculated_energy = self.intermediate_energy(surface, intermediate, bias, site)
                            print(intermediate, bias, site, calculated_energy)
                            bias_energies[bias].update({site: calculated_energy})
                        energies[surface].update({intermediate: bias_energies})
                for bias in self.data.get_surface_biases(surface):
                    energies[surface]["initial"].update({bias: self.terminal_energies(surface, bias)[0]})
                    energies[surface]["final"].update({bias: self.terminal_energies(surface, bias)[1]})
        material.energies = energies

    def calculate_binding_energies(self, materials:Materials, bias='0.00V') -> dict:
        binding_energies = {}
        for material in materials.materials:
            for surface in material.converged_surfaces():
                reference_molecule, reference_energy = self.reaction.binding_reference(self.ev)
                if self.data.check_adsorbed_convergence(surface, bias, "N") == False:
                    continue
                surface_energy = self.data.get_surface_energy(surface, bias)
                bound_energy = self.data.get_adsorbed_energy(surface, bias, reference_molecule)
                binding_energy = (bound_energy - (reference_energy + surface_energy)) * self.ev
                binding_energies[surface] = binding_energy
        return binding_energies
#TODO Don't need to define list of intermediates in both the References and Reaction
# objects

class Reaction:
    def __init__(self, reaction:str) -> None:
        intermediate_dict = {"NRR": ["N2", "N2H", "NNH2", "N", "NH", "NH2", "NH3"],
                             "HER": ["H"]}
        self.references = References(reaction)
        self.intermediates = intermediate_dict[reaction]
    
    def reference_energy(self, intermediate, bias, ev):
        if intermediate in ["initial", "final"]:
            reference_energy = self.references.get_termination_reference(intermediate, bias, ev)
        else:
            reference_energy = self.references.get_intermediate_references(intermediate, bias, ev)
        return reference_energy
    
    def find_reaction_intermediates(self, data:Data_Parser, surface) -> list:
        '''
        Compares reaction intermediates in the all_data to reaction intermediates
        defined in the reaction
        '''
        
        data_intermediates = data.get_intermediates(surface)
        return [i for i in data_intermediates if i in self.intermediates]
    
    def binding_reference(self, ev):
        return self.references.get_binder(ev)

class References:
    def __init__(self, reaction:str) -> None:
        self.reference_energies = {"N2":{"G": -18.2700562, "nelec":10}, "NH3":{"G":-10.4283238, "nelec":8}, 
                                    "H+":{"G":-0.419025298, "nelec":0}}
        self.reference_table = {"N2": {"H+":6}, "N2H": {"H+":5}, "NNH2": {"H+":4},
                           "N": {"NH3":1, "H+":3}, "NH": {"NH3":1, "H+":2},
                           "NH2": {"NH3":1, "H+":1}, "NH3":{"NH3":1}, "H":{"H+":1}
        self.reaction_terminations = {"NRR":{"initial":{"N2":1, "H+":6}, "final":{"NH3":2}}
            }
        self.reaction = reaction
        self.binders = {"NRR": {"N": {"N2": 0.5}}}
        
    def get_references(self, intermediate):
        pass
        
    
    def get_intermediate_references(self, intermediate, bias, ev):
        if bias == 'No_bias':
            bias = '0.00V' # convert bias to zero so that it isn't used in the reference calculation
        reference_energy = 0
        bias = self.bias_string_to_float(bias)
        mu = self.bias_to_mu(bias, ev)
        for molecule, quantity in self.reference_table[intermediate].items():
            # print(self.reference_energies[molecule]["G"], self.reference_energies[molecule]["nelec"])
            reference_energy += ((self.reference_energies[molecule]["G"]
                                - self.reference_energies[molecule]["nelec"] * mu)
                                 ) * quantity
            # print(molecule, quantity)
            # print("references for ", intermediate, ": ", reference_energy, " ", molecule, " ", quantity)
            # print( "MuN: " , mu * self.reference_energies[molecule]['nelec'])
        # print(reference_energy)
        return reference_energy
    
    #TODO: these should be the same method

    def get_termination_reference(self, termination, bias, ev):
        if bias == 'No_bias':
            bias = '0.00V' # convert bias to zero so that it isn't used in the reference calculation
        reference_energy = 0
        bias = self.bias_string_to_float(bias)
        mu = self.bias_to_mu(bias, ev)
        for molecule, quantity in self.reaction_terminations[self.reaction][termination].items():
            reference_energy += ((self.reference_energies[molecule]["G"]
                                - self.reference_energies[molecule]["nelec"] * mu)
                                 ) * quantity
            # print("references for ", termination, ": ", reference_energy, " ", molecule, " ", quantity)
            # print( "MuN: " , mu * self.reference_energies[molecule]['nelec'], mu)
        return reference_energy
    
    def bias_string_to_float(self, bias):
        return float(bias.split('V')[0])
    
    def bias_to_mu(self, bias:int, ev:float) -> float:
        return (-bias) / ev
    
    def get_binder(self, ev) -> (str,float):
        molecule = list(self.binders[self.reaction].keys())[0]
        reference_energy = 0
        for ref_molecule, quantity in self.binders[self.reaction][molecule].items():
            reference_energy += self.reference_energies[ref_molecule]["G"] * quantity
        return molecule, reference_energy
    
