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

    def __init__(self, cell, atoms, seed, count, optimizer):
        """
        Initializes a spyss object.
        """

        self.cell = cell
        self.atoms = atoms
        self.seed = seed
        self.count = count
        self.optimizer = optimizer

        self.stochastic_cells = []
        self.stochastic_structures = []
        self.stochastic_metastable_structures = []

    def generate_stochastic_cell(self):
        """
        Generates a stochastic cell using the spyss cell object.
        """

        if self.cell.mode == "fixed-cell":
            tmp_atoms = self.cell.structure.copy()

        if self.cell.mode == "variable-cell":
            raise NotImplementedError("Variable cell stochastic cell generation is not yet implemented.") 
       
        self.stochastic_cells.append(tmp_atoms)
        del tmp_atoms

    def generate_metastable_structure():
        """
        Generates a stochastic metastable structure for a given random seed.
        """

        self.generate_stochastic_cell()
        self.generate_stochastic_structure()
        self.optimize_stochastoc_structure()
