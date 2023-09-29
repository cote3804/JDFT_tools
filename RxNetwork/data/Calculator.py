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
        for material in materials.materials:
            energies = self.calculate_material_energies(material)
            for surface in material.surfaces:
                material.energies.update({surface:energies[surface]})

    def calculate_reference_energy(self, surface:str, bias:str, reference="final") -> float:
        if reference == "final":
            molecule_reference_energy = self.reaction.reference_energy("final", bias, self.ev) * self.ev
            surface_reference_energy = self.data.get_surface_energy(surface, bias) * self.ev
            reference_energy = molecule_reference_energy + surface_reference_energy
        elif reference == "initial":
            molecule_reference_energy = self.reaction.reference_energy("initial", bias, self.ev) * self.ev
            surface_reference_energy = self.data.get_surface_energy(surface, bias) * self.ev
            reference_energy = molecule_reference_energy + surface_reference_energy
        return reference_energy
            
    def intermediate_energy(self, surface:str, intermediate:str, bias:str, site:str) -> float:
        allowed_bias = self.data.check_adsorbed_bias(surface, bias, intermediate)
        if allowed_bias:
            converged = self.data.check_adsorbed_site_convergence(surface, bias, intermediate, site)
            if converged == False:
                return None
            elif converged == True:
                adsorbed_energy = self.data.get_adsorbed_energy(surface, bias, intermediate)
                reference_energy = self.reaction.reference_energy(intermediate, bias, self.ev)
                calculated_energy = adsorbed_energy + reference_energy
                return calculated_energy * self.ev
        elif allowed_bias == False:
            return None

    
    def terminal_energies(self, surface:str, bias:str) -> (float, float):
        surface_energy = self.data.get_surface_energy(surface, bias) * self.ev
        initial_energy = self.reaction.reference_energy("initial", bias, self.ev) * self.ev
        final_energy = self.reaction.reference_energy("final", bias, self.ev) * self.ev
        return (surface_energy + initial_energy), (surface_energy + final_energy)

    def calculate_material_energies(self, material) -> None:
        energies = {}
        surfaces = material.surfaces
        for surface in surfaces:
            if self.data.check_any_surface_convergence(surface) == False:
                print(f"{surface} has no converged sites")
                continue
            elif self.data.check_any_surface_convergence(surface) == True:
                energies[surface] = {}
                initial_state, final_state = self.reaction.terminal_to_states()
                energies[surface][initial_state] = {}
                energies[surface][final_state] = {}
                for intermediate in self.reaction.find_reaction_intermediates(self.data, surface):
                    bias_energies = {}
                    for bias in self.data.get_intermediate_biases(surface, intermediate):
                        bias_energies[bias] = {}
                        for site in self.data.get_sites_data(surface, bias, intermediate).keys():
                            calculated_energy = self.intermediate_energy(surface, intermediate, bias, site)
                            bias_energies[bias].update({site: calculated_energy})
                        energies[surface].update({f"{intermediate}*": bias_energies}) # adding * here for adsorbed states
                for bias in self.data.get_surface_biases(surface):
                    energies[surface][initial_state].update({bias: self.terminal_energies(surface, bias)[0]})
                    energies[surface][final_state].update({bias: self.terminal_energies(surface, bias)[1]})
        return energies

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
    
    def calculate_span(self, reaction, material, surface, bias):
        reaction = Reaction(reaction)
        intermediates = self.data.get_converged_intermediates(surface, bias)
        energies = []
        E_init, E_final = self.terminal_energies(surface, bias)
        E_rxn = E_final - E_init
        energies.extend([E_init, E_final])
        index_list = [reaction.intermediate_index("initial") , reaction.intermediate_index("final")]
        for intermediate in intermediates:
            min_site, min_energy = self.data.get_lowest_site(surface, intermediate, bias)
            energies.append(self.intermediate_energy(surface, intermediate, bias, min_site))
            index_list.append(reaction.intermediate_index(intermediate))
        E_span = -0.0001
        span_indices = []
        for E_i, index_i in zip(energies, index_list):
            for E_j, index_j in zip(energies, index_list):
                if index_i < index_j:
                    if E_j - E_i > E_span:
                        E_span = E_j - E_i
                        span_indices = [index_i, index_j]
                elif index_i > index_j:
                    if E_i - E_j - E_rxn > E_span:
                        E_span = E_i - E_j - E_rxn
                        span_indices = [index_i, index_j]
        span_intermediates = (reaction.index_to_state(span_indices[0]), reaction.index_to_state(span_indices[1]))
        return E_span, span_intermediates
    
    def get_FED_energy(self, material:Material, surface:str, bias:str, referenced="final"):
        if self.data.check_surface_convergence(surface, bias) == False:
                print(f"no converged surfaces for {surface} at {bias}")
                return None
        elif self.data.check_surface_convergence(surface, bias) == True:
            reference_energy = self.calculate_reference_energy(surface, bias, reference=referenced)
        FED_energies = {}
        initial_state, final_state = self.reaction.terminal_to_states()
        surface_energy = self.calculate_material_energies(material)[surface]
        for intermediate in surface_energy.keys():
            if intermediate not in [initial_state, final_state]:
                if self.data.check_adsorbed_convergence(surface, bias, intermediate.strip('*')) == False:
                    continue
                float_energies = [i for i in surface_energy[intermediate][bias].values() if type(i) == float]
                # float_energies is needed to filter out the None values from unconverged calculations
                FED_energy = min(float_energies)
                FED_energies[intermediate] = FED_energy - reference_energy
            elif intermediate in [initial_state, final_state]:
                continue
        FED_energies[initial_state] = surface_energy[initial_state][bias] - reference_energy
        FED_energies[final_state] = surface_energy[final_state][bias] - reference_energy
        return FED_energies


#TODO Don't need to define list of intermediates in both the References and Reaction
# objects

class Reaction:
    def __init__(self, reaction:str) -> None:
        intermediate_dict = {"NRR": ["N2*", "N2H*", "NNH2*", "N*", "NH*", "NH2*", "NH3*"],
                             "HER": ["H*"]}
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
        # data_intermediates = [f"{i}*" for i in data_intermediates]
        return [i for i in data_intermediates if f"{i}*" in self.intermediates]
    
    def binding_reference(self, ev):
        return self.references.get_binder(ev)

    def intermediate_index(self, intermediate:str) -> int:
        index_data = self.references.get_state_indices()
        return index_data[intermediate][0]
    
    def index_to_state(self, index:int) -> str:
        for intermediate, index_tuple in self.references.get_state_indices().items():
            if index_tuple[0] == index:
                return intermediate
            
    def terminal_to_states(self) -> tuple:
        # returns the state string for the initial and final states
        return (self.references.get_initial_state(), self.references.get_final_state())




class References:
    def __init__(self, reaction:str) -> None:
        self.reference_energies = {"N2":{"G": -18.2700562, "nelec":10}, "NH3":{"G":-10.4283238, "nelec":8}, 
                                    "H+":{"G":-0.419025298, "nelec":0}, "H2":{"G":-0.838050596, "nelec":2},
                                }
        self.reference_table = {"N2*": {"H+":6}, "N2H*": {"H+":5}, "NNH2*": {"H+":4},
                           "N*": {"NH3":1, "H+":3}, "NH*": {"NH3":1, "H+":2},
                           "NH2*": {"NH3":1, "H+":1}, "NH3*":{"NH3":1}, "H*":{"H+":1}
                            }
        self.reaction_terminations = {"NRR":{"initial":{"N2":1, "H+":6}, "final":{"NH3":2}},
                                      "HER":{"initial":{"H+":2}, "final":{"H2":1}}
                                    }
        self.termination_labels = {"NRR":{"initial":"N2", "final":"2NH3"},
                                   "HER":{"initial":"2H+", "final":"H2"}
                                }
        self.state_indices = {"NRR": {"initial":(0,0),
                                      "N2*": (1,0),
                                        "N2H*": (2,0),
                                        "NNH2*": (3,0),
                                        "N*": (4,0),
                                        "NH*": (5,0),
                                        "NH2*": (6,0),
                                        "NH3*": (7,0),
                                        "final": (8,0)
                                        },
                                "HER": {"initial":(0,0),
                                        "H*": (1,0),
                                        "final": (2,0)
                                        }
                                }
        self.reaction = reaction
        self.binders = {"NRR": {"N*": {"N2": 0.5}},
                        "HER": {"H*": {"H+": 1}}
                    }
        
    def get_references(self, intermediate):
        pass
        
    
    def get_intermediate_references(self, intermediate, bias, ev):
        if bias == 'No_bias':
            bias = '0.00V' # convert bias to zero so that it isn't used in the reference calculation
        reference_energy = 0
        bias = self.bias_string_to_float(bias)
        mu = self.bias_to_mu(bias, ev)
        for molecule, quantity in self.reference_table[f"{intermediate}*"].items():
            reference_energy += ((self.reference_energies[molecule]["G"]
                                - self.reference_energies[molecule]["nelec"] * mu)
                                 ) * quantity
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
    
    def get_state_indices(self):
        return self.state_indices[self.reaction]
    
    def get_final_state(self) -> str:
        return self.termination_labels[self.reaction]["final"]
    
    def get_initial_state(self) -> str:
        return self.termination_labels[self.reaction]["initial"]
    
