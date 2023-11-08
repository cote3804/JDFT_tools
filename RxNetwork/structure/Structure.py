import numpy as np
import json
from pymatgen.core import Structure as PyMatStruct
from pymatgen.io.ase import AseAtomsAdaptor
from ase.geometry.analysis import Analysis as aseAnalysis
from ase.visualize import view
from ase import Atoms
from ase.data import covalent_radii, vdw_radii
from ase.neighborlist import NeighborList, natural_cutoffs

class Structure:
    def __init__(self, struct_dict):
        if type(struct_dict) == str:
            struct_dict = json.loads(struct_dict)
        self.pymatstruct = PyMatStruct.from_dict(struct_dict)
        symbols = [str(site.specie.symbol) for site in self.pymatstruct]
        positions = [site.coords for site in self.pymatstruct]
        if hasattr(self.pymatstruct, "lattice"):
            pbc = True
            cell = self.pymatstruct.lattice.matrix
        else:  # Molecule without lattice
            pbc = False
            cell = None
        self.atoms = Atoms(symbols=symbols, positions=positions, pbc=pbc, cell=cell)
    
    # def __repr__(self): # formats what gets printed when calling print() on an instance of this class
    #     return repr(self.atoms)
    
    def visualize(self):
        view(self.atoms)
    
    def count_neighbor_bonds(self, atom_index, buffer=0.1):
        neighbor_indices = self.get_neighbors(atom_index, buffer=buffer)
        return len(neighbor_indices)
    
    def get_neighbors(self, atom_index, buffer):
        cutoffs = natural_cutoffs(self.atoms, mult=buffer) # list of covalent radii
        neighbor_list = NeighborList(cutoffs, self_interaction=False, bothways=True)
        neighbor_list.update(self.atoms)
        neighbors = neighbor_list.get_neighbors(atom_index) # list of neighbor indices within vocalent radii times a buffer value
        return neighbors[0]
        
    def list_of_atomic_symbols(self, intermediate) -> list:
        atoms_data = {"N2*": {"N": 2}, "N2H*": {"N": 2, "H": 1},
                      "N*": {"N": 1}, "NH*": {"N": 1, "H": 1}, 
                      "NH2*": {"N": 1, "H": 2}, "NH3*": {"N": 1, "H": 3},
                      "NNH2*": {"N": 2, "H": 2}} # defines number of atoms in adsorbate
        atomic_symbols = []
        for atom, number in atoms_data[intermediate].items():
            atomic_symbols.extend([atom]*number)
        return atomic_symbols
    
    def get_adsorbate_atoms(self, intermediate):
        atomic_symbols = self.list_of_atomic_symbols(intermediate)
        # assuming that the adsorbate atoms are highest on the surface.
        # first convert to cubic to make sure we're finding highest atoms relative to the surface plane
        cubic_atoms = self.cubic_cell()
        structure_symbols = cubic_atoms.get_chemical_symbols()
        cubic_positions = cubic_atoms.get_positions()
        # now loop though the atoms defined in atomic_symbols, find the highest atoms, and return their indices
        highest_atom_indices = []
        for ads_atom in atomic_symbols:
            for iatom, atom in enumerate(structure_symbols):
                z_position = -1000 # initialize to a very low value
                if atom in atomic_symbols and cubic_positions[iatom][2] > z_position:
                    index = iatom
                    z_position = cubic_positions[iatom][2]
            structure_symbols.pop(index) # get rid of value we just found
            highest_atom_indices.append(index)
        
        return highest_atom_indices

    def bond_distance(self, atom_indices):
        if len(atom_indices) != 2:
            raise Exception("bond distance requires two atoms")
        distance = self.atoms.get_distance(atom_indices[0], atom_indices[1])
        return distance
    
    def cubic_cell(self):
        # returns the cubic cell of the structure
        atomic_numbers = self.atoms.get_atomic_numbers()
        S = np.eye(3) # cubic lattice 
        L = np.array(self.atoms.get_cell()) # current lattice
        np.fill_diagonal(S, np.linalg.norm(L, axis=1)) # resize cubic lattice to match current lattice vector lengths
        positions = self.atoms.get_positions()
        new_positions = []
        for i, vec in enumerate(positions):
            new_positions.append(np.linalg.inv(L.T).dot(vec)) # convert to fractional coordinates in original lattice
        new_positions = np.array(new_positions) 
        new_positions = S.dot(new_positions.T).T # convert to coordinates in new lattice basis
        cubic_atoms = Atoms(atomic_numbers, positions=new_positions, cell=S, pbc=True)
        return cubic_atoms
    
    def coordination_number(self, intermediate, buffer=1.1):
        unique_bonds = self.intermediate_neighbors(intermediate, buffer=buffer)
        return len(unique_bonds)

    def intermediate_neighbors(self, intermediate:str, buffer=1.05) -> list:
        # function to get the neighbors of an adsorbate on a surface given only the adsorbate string
        adsorbate_indices = self.get_adsorbate_atoms(intermediate)
        neighbor_indices = []
        for atom_index in adsorbate_indices:
            neighbor_indices.extend(self.get_neighbors(atom_index, buffer=buffer))
        
        neighbor_indices = [iatom for iatom in neighbor_indices if iatom not in adsorbate_indices] # remove duplicates  
        unique_bonds = list(set(neighbor_indices))
        return unique_bonds
    
    def check_dissociation(self, intermediate:str, buffer=1.05) -> bool:
        # function to check if an adsorbate is dissociated on a surface
        # returns True if dissociated, False if not
        # currently nothing is actually implemented here
        return False

    @property
    def positions(self):
        positions = self.atoms.get_positions()
        return positions
    
    @positions.setter
    def positions(self, positions):
        self.atoms.set_positions(positions)

    @property
    def cell(self):
        cell = self.atoms.get_cell()
        return cell
    
    @cell.setter
    def cell(self, cell):
        self.atoms.set_cell(cell)