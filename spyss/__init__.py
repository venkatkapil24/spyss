import numpy as np
from ase.io import read, write
from ase.constraints import FixAtoms


class spyss:
    """
    A stochastic pythonic structure search class to generate a metastable
    structure of compound. A sypss object generates a set of metastable 
    strcuture by
    - creating a set of lattice vectors (predefined or random using a seed)
    - adding atoms (deterministic or stochastic using a seed)
    - optimizing the structure (at fixed or variable lattice vectors)
    """

    def __init__(self, cell, atoms, seed, total_structures, optimizer):
        """
        ! Initializes a spyss object.
        @param cell The name of the spyss cell object.
        @param atoms The name of the spyss atoms object.
        @param seed The random number seed to reproduce structures.
        @param total_structures The total number of structures to be generated.
        @param optimizer The name of the spyss optimizer object.
        """

        self.cell = cell
        self.atoms = atoms
        self.seed = seed
        self.total_structures = total_structures
        self.optimizer = optimizer

        self.count = 0 
        self.stochastic_cells = []
        self.stochastic_structures = []
        self.stochastic_metastable_structures = []


    def update_count(self):
        """
        ! Keeps track of the number of generated structures.
        """

        self.count = len(self.stochastic_cells)


    def generate_stochastic_cell(self):
        """
        ! Generates a cell with fixed or (stochastic) variable shape.
        """

        if self.cell.mode == "fixed-cell":
            tmp_atoms = self.cell.structure.copy()

        if self.cell.mode == "variable-cell":
            raise NotImplementedError("Variable cell stochastic cell generation is not yet implemented.") 
       
        self.stochastic_cells.append(tmp_atoms)
        self.update_count()
        del tmp_atoms


    def generate_stochastic_structure(self):
        """
        ! Generates a structure with a stochasic cell and with  initial and 
        stochastic atoms.
        """

        self.stochastic_structures = self.stochastic_structures[0:self.count-1]
        
        tmp_atoms = self.stochastic_cells[self.count - 1].copy()
        
        # Checks if initial_atoms exist
        # If yes, adds them to the structure
        # and fixes all the initial_atoms
        
        if self.atoms.initial_atoms is not None: 
            tmp_atoms += self.atoms.get_initial_atoms()
            c = FixAtoms(indices=[atom.index for atom in tmp_atoms])
            tmp_atoms.set_constraint(c)
    
        # Checks if stochastic_atoms exist
        # If yes, adds them to the structure
        
        if len(self.atoms.stochastic_atoms_list) !=0:
            tmp_atoms += self.atoms.get_stochastic_atoms()
            self.stochastic_structures.append(tmp_atoms)
        
        del tmp_atoms

        
    def optimize_stochastic_structure(self):
        """
        ! Generates a structure with a stochasic cell and with  initial and 
        stochastic atoms.
        """

        self.stochastic_metastable_structures = self.stochastic_metastable_structures[0:self.count-1]
        tmp_atoms = self.optimizer.optimize(self.stochastic_structures[self.count - 1].copy())
        self.stochastic_metastable_structures.append(tmp_atoms)
        del tmp_atoms


    def generate_metastable_structure(self):
        """
        ! Generates a stochastic metastable structure for a given random seed.
        """

        self.generate_stochastic_cell()
        self.generate_stochastic_structure()
        self.optimize_stochastic_structure()
        
        
    def generate_all_metastable_structures(self):
        """
        ! Generates all stochastic metastable structures.
        """

        for index in range(self.total_structures):
            self.generate_metastable_structure()
            
            
    def sort_metastable_structures(self):
        """
        ! Sorts all stochastic metastable structures by energy.
        """

        tmparray = np.argsort([atoms.info['energy'] for atoms in self.stochastic_metastable_structures])
        self.sorted_metastable_energies = [self.stochastic_metastable_structures[index].info['energy'] for index in tmparray]
        self.sorted_metastable_structures = [self.stochastic_metastable_structures[index].copy() for index in tmparray]
        
