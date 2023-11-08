from data.Data_Parser import Data_Parser
from data.Materials import Materials, Material
from structure.Structure import Structure
from structure.Site import Site, Cluster
from helpers.conversions import bias_float_to_str, mu_to_she, H_to_eV
# Adsorbed intermediates in the all_data are stored without the * symbol, but
# intermdiate keys in the data structures in this script for calculating
# reaction references are labeled with the * symbol. Whenever a reaction 
# reference is being calculated, the * must be added to the intermediate string
# from the all_data.

#TODO move this explanation to a more logical spot. Maybe in the doc string for the Calculator class.


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
                adsorbed_energy = self.data.get_adsorbed_energy(surface, bias, intermediate, site)
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

    def calculate_surface_energies(self, surface:str) -> dict:
        energies = {}
        if self.data.check_any_surface_convergence(surface) == False:
            print(f"{surface} has no converged sites")
        elif self.data.check_any_surface_convergence(surface) == True:
            initial_state, final_state = self.reaction.terminal_to_states()
            energies[initial_state] = {}
            energies[final_state] = {}
            for intermediate in self.reaction.find_reaction_intermediates(self.data, surface):
                bias_energies = {}
                for bias in self.data.get_intermediate_biases(surface, intermediate):
                    bias_energies[bias] = {}
                    for site in self.data.get_sites_data(surface, bias, intermediate).keys():
                        calculated_energy = self.intermediate_energy(surface, intermediate, bias, site)
                        bias_energies[bias].update({site: calculated_energy})
                    energies.update({f"{intermediate}*": bias_energies}) # adding * here for adsorbed states
            for bias in self.data.get_surface_biases(surface):
                energies[initial_state].update({bias: self.terminal_energies(surface, bias)[0]})
                energies[final_state].update({bias: self.terminal_energies(surface, bias)[1]})
        return energies

    def calculate_binding_energies(self, materials:Materials, bias='0.00V', adsorbate="N") -> dict:
        binding_energies = {}
        for material in materials.materials:
            for surface in material.converged_surfaces():
                reference_molecule, reference_energy = self.reaction.binding_reference(self.ev)
                if self.data.check_adsorbed_convergence(surface, bias, adsorbate) == False:
                    continue
                surface_energy = self.data.get_surface_energy(surface, bias)
                bound_energy = self.data.get_adsorbed_energy(surface, bias, reference_molecule.strip('*'))
                binding_energy = (bound_energy - (reference_energy + surface_energy)) * self.ev
                binding_energies[surface] = binding_energy
        return binding_energies
    
    def calculate_span(self, reaction, surface, bias):
        reaction = Reaction(reaction)
        intermediates = reaction.find_reaction_intermediates(self.data, surface)
        intermediates = [f"{i}*" for i in intermediates] # convert intermediates to the format necessary for calculating references
        energies = []
        E_init, E_final = self.terminal_energies(surface, bias)
        E_rxn = E_final - E_init
        energies.extend([E_init, E_final]) # first two energies are initial and final energies
        initial_state, final_state = reaction.terminal_to_states() # get the strings for the initial and final states
        index_list = [reaction.intermediate_index(initial_state) , reaction.intermediate_index(final_state)] # set first two indices to initial and final
        for intermediate in intermediates:
            if self.data.check_adsorbed_convergence(surface, bias, intermediate.strip('*')) == False:
                continue
            min_site, min_energy = self.data.get_lowest_site(surface, intermediate.strip('*'), bias)
            energies.append(self.intermediate_energy(surface, intermediate.strip('*'), bias, min_site))
            index_list.append(reaction.intermediate_index(intermediate))
        E_span = -0.001
        span_indices = []
        for E_i, index_i in zip(energies, index_list):
            for E_j, index_j in zip(energies, index_list):
                if index_i < index_j:
                    if E_j - E_i > E_span:
                        E_span = E_j - E_i
                        span_indices = [index_i, index_j]
                elif index_i > index_j:
                    if E_j - E_i + E_rxn > E_span:
                        E_span = E_j - E_i + E_rxn
                        span_indices = [index_i, index_j]
        if len(span_indices) == 0:
            E_span = 0
            span_intermediates = (initial_state, final_state)
        else:
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

    def get_pathway_charges(self, material:Material, surface:str, bias:str, reference="final") -> dict:
        FED_energies = self.get_FED_energy(material, surface, bias)
        pathway_charges = {}
        initial_state, final_state = self.reaction.terminal_to_states()  
        initial_charge = self.calculate_terminal_charge(surface, bias, initial_state)
        final_charge = self.calculate_terminal_charge(surface, bias, final_state)
        if initial_charge == None or final_charge == None:
            return None
        if reference == "final":
            reference_charge = final_charge
        elif reference == "initial":
            reference_charge = initial_charge
        elif reference == "absolute":
            reference_charge = 0

        pathway_charges[initial_state] = initial_charge - reference_charge
        pathway_charges[final_state] = final_charge - reference_charge
        for intermediate, energy in FED_energies.items():
            if intermediate not in [initial_state, final_state]:
                pathway_charges[intermediate] = self.calculate_intermediate_charge(surface, bias, intermediate) - reference_charge
        return pathway_charges

    def calculate_intermediate_charge(self, surface:str, bias:str, intermediate:str):
        reference_charge = self.reaction.charge(intermediate)
        adsorbed_charge = self.data.get_adsorbed_charge(surface, bias, intermediate.strip('*'))
        total_charge = adsorbed_charge + reference_charge
        return total_charge
    
    def calculate_terminal_charge(self, surface:str, bias:str, terminal:str):
        reference_charge = self.reaction.charge(terminal)
        surface_charge = self.data.get_surface_charge(surface, bias)
        if surface_charge == None:
            return None
        total_charge = surface_charge + reference_charge
        return total_charge

    def set_common_clusters(self, surface, bias, intermediate) -> list:
        # method that gets all the cluster instances for the binding configurations of the specified intermediates
        # returns a list of cluster instances
        site = Site(self.data.get_sites_data(surface, bias, intermediate.strip('*')))
        unique_clusters = site.unique_clusters(intermediate)
        return unique_clusters

    def aggregate_cluster_FEDs(self, surface, bias, intermediate) -> dict:
        # method to calculate the FED of multiple clusters on a surface for a given reaction.
        # if multiple sites are within the same cluster, the ground state is chosen.
        reference_clusters = self.set_common_clusters(surface, bias, intermediate)
        cluster_FEDs = {}
        print("need to check these clusters:", reference_clusters)
        for cluster in reference_clusters:
            cluster_str = repr(cluster)
            print("checking cluster ", cluster_str, "\n")
            cluster_FED = self.calculate_cluster_FED(cluster, surface, bias)
            cluster_FEDs[cluster_str] = cluster_FED
        
        return cluster_FEDs
    
    def aggregate_cluster_spans(self, surface, bias, intermediate) -> dict:
        # method to calculate the span of multiple clusters on a surface for a given reaction.
        # if multiple sites are within the same cluster, the ground state is chosen.
        reference_clusters = self.set_common_clusters(surface, bias, intermediate)
        cluster_spans = {}
        for cluster in reference_clusters:
            cluster_span = self.calculate_cluster_span(cluster, surface, bias)
            cluster_str = repr(cluster)
            cluster_spans[cluster_str] = cluster_span
        return cluster_spans

    def calculate_cluster_FED(self, reference_cluster:Cluster, surface:str, bias:str, referenced="final"):
        # method to calculate the FED of a single cluster on a surface for a given reaction.
        # if multiple sites are within the same cluster, the ground state is chosen.
        if self.data.check_surface_convergence(surface, bias) == False:
                print(f"no converged surfaces for {surface} at {bias}")
                return None
        elif self.data.check_surface_convergence(surface, bias) == True:
            reference_energy = self.calculate_reference_energy(surface, bias, reference=referenced)
        surface_energies = self.calculate_surface_energies(surface)
        initial_state, final_state = self.reaction.terminal_to_states()
        FED_energy_dict = {}
        for intermediate in surface_energies.keys():
            if intermediate in [initial_state, final_state]:
                continue
            elif intermediate not in [initial_state, final_state]:
                site = Site(self.data.get_sites_data(surface, bias, intermediate.strip('*')))
                ground_sites = reference_cluster.matching_ground_site(site, intermediate) # list of matching sites
                min_site, min_energy = self.data.minimum_energy_of_sites(surface, bias, intermediate.strip('*'), ground_sites)
                if min_site == None:
                    site_energy = None
                else:
                    site_energy = surface_energies[intermediate][bias][min_site] - reference_energy
                FED_energy_dict[intermediate] = site_energy
        FED_energy_dict[initial_state] = surface_energies[initial_state][bias] - reference_energy
        FED_energy_dict[final_state] = surface_energies[final_state][bias] - reference_energy
        return FED_energy_dict

    def calculate_cluster_span(self, reference_cluster:Cluster, surface:str, bias:str):
        reaction = self.reaction
        intermediates = reaction.find_reaction_intermediates(self.data, surface)
        intermediates = [f"{i}*" for i in intermediates] # convert intermediates to the format necessary for calculating references
        energies = []
        E_init, E_final = self.terminal_energies(surface, bias)
        E_rxn = E_final - E_init
        energies.extend([E_init, E_final]) # first two energies are initial and final energies
        initial_state, final_state = reaction.terminal_to_states() # get the strings for the initial and final states
        index_list = [reaction.intermediate_index(initial_state) , reaction.intermediate_index(final_state)] # set first two indices to initial and final
        for intermediate in intermediates:
            if self.data.check_adsorbed_convergence(surface, bias, intermediate.strip('*')) == False:
                continue
            site = Site(self.data.get_sites_data(surface, bias, intermediate.strip('*')))
            ground_site = reference_cluster.matching_ground_site(site, intermediate)
            energies.append(self.intermediate_energy(surface, intermediate.strip('*'), bias, ground_site))
            index_list.append(reaction.intermediate_index(intermediate))
        E_span = -0.001
        span_indices = []
        for E_i, index_i in zip(energies, index_list):
            for E_j, index_j in zip(energies, index_list):
                if index_i < index_j:
                    if E_j - E_i > E_span:
                        E_span = E_j - E_i
                        span_indices = [index_i, index_j]
                elif index_i > index_j:
                    if E_i - E_j + E_rxn > E_span:
                        E_span = E_j - E_i + E_rxn
                        span_indices = [index_i, index_j]
        if len(span_indices) == 0:
            E_span = 0
            span_intermediates = (initial_state, final_state)
        else:
            span_intermediates = (reaction.index_to_state(span_indices[0]), reaction.index_to_state(span_indices[1]))
        return E_span, span_intermediates
            
    def calculate_band_center(self, surface:str, bias:str) -> float:
        # method to calculate the band center of a surface
        transition_metals = [
            'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
            'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
            'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg'
        ]
        band_centers = {}
        s_total = 0
        p_total = 0
        d_total = 0
        s_count = 0
        p_count = 0
        d_count = 0
        if self.data.check_surface_convergence(surface, bias) == False:
                print(f"no converged surfaces for {surface} at {bias}")
                return None
        else:
            band_data = self.data.get_band_data(surface, bias)
            for spin in band_data.keys():
                for element in band_data[spin].keys():
                    if element in transition_metals:
                        for atom_number in band_data[spin][element].keys():
                            if "s" in band_data[spin][element][atom_number].keys():
                                s_total += band_data[spin][element][atom_number]["s"]
                                s_count += 1
                            if "p" in band_data[spin][element][atom_number].keys():
                                p_total += band_data[spin][element][atom_number]["p"]
                                p_count += 1
                            if "d" in band_data[spin][element][atom_number].keys():
                                d_total += band_data[spin][element][atom_number]["d"]
                                d_count += 1
        
        if s_count > 0:
            band_centers["s"] = s_total / s_count * 27.2114
        else:
            band_centers["s"] = None
        if p_count > 0:
            band_centers["p"] = p_total / p_count * 27.2114
        else:
            band_centers["p"] = None
        if d_count > 0:
            band_centers["d"] = d_total / d_count * 27.2114
        else:
            band_centers["d"] = None
        return band_centers

    def get_pzc(self, surface:str, she=True) -> float:
        mu = self.data.surface_pzc(surface)
        if mu == None:
            return None
        else:
            if she:
                return mu_to_she(mu)
            else:
                return mu
            

#TODO Don't need to define list of intermediates in both the References and Reaction
# objects

class Reaction:
    def __init__(self, reaction:str) -> None:
        intermediate_dict = {"NRR": ["N2*", "N2H*", "NNH2*", "N*", "NH*", "NH2*", "NH3*"],
                             "HER": ["H*"]}
        self.references = References(reaction)
        self.intermediates = intermediate_dict[reaction] # these intermediates don't include the terminal states
    
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
        return self.references.get_binder()

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
    
    def charge(self, intermediate:str) -> int:
        return self.references.get_charge(intermediate)
    
    def binding_adsorbate(self) -> str:
        # returns adsorbate string
        return self.references.get_binder()[0]
    




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
        self.state_indices = {"NRR": {"N2":(0,0),
                                      "N2*": (1,0),
                                        "N2H*": (2,0),
                                        "NNH2*": (3,0),
                                        "N*": (4,0),
                                        "NH*": (5,0),
                                        "NH2*": (6,0),
                                        "NH3*": (7,0),
                                        "2NH3": (8,0)
                                        },
                                "HER": {"2H+":(0,0),
                                        "H*": (1,0),
                                        "H2": (2,0)
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
        elif bias in ["No_bias_freeze_none","No_bias_nofreeze"]:
            bias = '0.00V'
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
    
    def get_binder(self) -> (str,float):
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
    
    def get_charge(self, intermediate:str) -> int:
        charge = 0
        # calculate the charge of any state's references
        if intermediate == self.get_initial_state():
            for molecule in self.reaction_terminations[self.reaction]["initial"].keys():
                charge += self.reference_energies[molecule]["nelec"] * self.reaction_terminations[self.reaction]["initial"][molecule]
            return charge
        elif intermediate == self.get_final_state():
            for molecule in self.reaction_terminations[self.reaction]["final"].keys():
                charge += self.reference_energies[molecule]["nelec"] * self.reaction_terminations[self.reaction]["final"][molecule]
            return charge
        else:
            for ref_molecule in self.reference_table[intermediate].keys():
                charge += self.reference_energies[ref_molecule]["nelec"] * self.reference_table[intermediate][ref_molecule]
            return charge
    
