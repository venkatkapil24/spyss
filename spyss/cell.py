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

    def __init__(self, structure_file, P, mode):
        """
        Initializes a spyss cell.
        """

        self.structure_file = structure_file
        self.P = P
        self.mode = mode

        with open(structure_file, 'r') as tmp_file:
            tmp_atoms = read(tmp_file)
            self.structure = ase.Atoms(cell=tmp_atoms.cell, pbc=tmp_atoms.pbc)
            self.structure = make_supercell(self.structure, P)
            del tmp_atoms

        self.volume = self.structure.get_volume()

