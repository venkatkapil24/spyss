import ase
from ase.io import read
from ase.build import make_supercell

class cell:
    """
    A sypyss cell object is described by an ASE atoms object 
    of a structure used to initialize the cell, a cell mode 
    used to contruct lattice vectors, and a pseudo random 
    number seed.
    """

    def __init__(self, ase_atoms, P, mode):
        """
        Initializes a spyss cell.
        """

        self.ase_atoms = ase_atoms
        self.P = P
        self.mode = mode
        
        self.volume = self.ase_atoms.get_volume()

