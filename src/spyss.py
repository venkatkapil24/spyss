import numpy as np
from ase.io import read, write


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
        Initializes a spyss object.
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
        Keeps track of the number of generated structures.
        """

        self.count = len(self.stochastic_cells)


    def generate_stochastic_cell(self):
        """
        Generates a cell with fixed or (stochastic) variable shape.
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
        Generates a structure with a stochasic cell and with  initial and 
        stochastic atoms.
        """

        tmp_atoms = self.stochastic_cells[self.count - 1]
        tmp_atoms += self.atoms.get_initial_atoms()
        tmp_atoms += self.atoms.get_stochastic_atoms()
        self.stochastic_structures.append(tmp_atoms)
        del tmp_atoms


    def generate_metastable_structure():
        """
        Generates a stochastic metastable structure for a given random seed.
        """

        self.generate_stochastic_cell()
        self.generate_stochastic_structure()
        self.optimize_stochastoc_structure()
