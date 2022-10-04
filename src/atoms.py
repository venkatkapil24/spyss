import ase
import numpy as np

from ase.io import read
from ase import neighborlist

def check_for_overlap(atoms, cutoff_dict=None):
    """
    Reads an ASE-atoms object and determines if there are atoms that overlap.
    """
    
    if cutoff_dict is None:
        cutoff_dict = np.asarray(neighborlist.natural_cutoffs(atoms)) * 0.5
    
    number_of_clashes = len(neighborlist.neighbor_list('i', atoms, cutoff=cutoff_dict))

    if number_of_clashes > 0:
        return True
    else:
        return False


class stochastic_region:
    """
    Determines geometric space and associated functions where stochastic 
    atoms are to be added.
    """

    def __init__(self, mode, space_dict):
        """
        Defines the mode (the nature of the region) and a dictionary that
        with paramters that aid transformation of a pesudo random numbers
        in [0,1]^3 to the described space.
        """

        self.mode = mode
        self.space_dict = space_dict

        self.read_space_dict()


    def _ufunc_transformation_cuboidal(self, r3):
        """"
        Transforms a randfom number in [0,1]^3 to a cuboidal region.
        """

        tmparray = np.zeros(3)

        for i, il  in zip(range(3), ['A', 'B', 'C']):

            l = self.space_dict['cell_matrix'][i]
            lim = self.space_dict[il + '_limits']
            tmparray += (lim[0] + r3[i] * (lim[1] - lim[0])) * l
            del l, lim

        return tmparray


    def read_space_dict(self):
        """
        Checks if mode and space_dict are commensurate
        """

        if self.mode == "full":

            self.transformation_function = self._ufunc_transformation_cuboidal

            try:
                with open(self.space_dict['cell_file'], 'r') as tmpfile:
                    self.space_dict['cell_matrix'] = np.asarray(read(tmpfile).cell)
            except:
                raise ValueError('key: cell_file not found in space_dict.')

            self.space_dict['A_limits'] = [0, 1]
            self.space_dict['B_limits'] = [0, 1]
            self.space_dict['C_limits'] = [0, 1]
                
        elif self.mode == 'slab':
            
            self.transformation_function = self._ufunc_transformation_cuboidal

            # Reads the cell matrix
            try:
                with open(self.space_dict['cell_file'], 'r') as tmpfile:
                    self.space_dict['cell_matrix'] = np.asarray(read(tmpfile).cell)
            except:
                raise ValueError('key: cell_file not found in space_dict.')

            # Checks if the bounds for the A,B or C cell vector are valid
            for l in ['A', 'B', 'C']:
                try:

                    limits = self.space_dict[l + '_limits']

                    if limits == 'full':
                        print (l + ' full')
                        self.space_dict[l + '_limits'] = [0,1]

                    elif len(limits) == 2:

                        assert (type(limits[0]) is int or type(limits[0]) is float)
                        assert (limits[0] > 0 and limits[0] <1)

                        assert (type(limits[1]) is int or type(limits[1]) is float)
                        assert (limits[1] > 0 and limits[1] <1)

                        assert (limits[1] > limits[0])

                    else:
                        raise ValueError(l + '_limits' + ' should either be a string "full" or a list of two numbers indicating fractional bounds.')
                        
                    del limits
            
                except KeyError:
                    raise ValueError('key: A_limits not found in space_dict.')


class atoms:
    """
    A sypyss atoms object is described by an ASE atoms object 
    of initial atoms, a list of ASE atoms objects of stochastic 
    atoms, and a list of the number of individual stochastic atoms, 
    a random number seed, and the region where stochastic atoms 
    are to be added. 
    """

    def __init__(self, initial_atoms_file, stochastic_atoms_files, number_of_stochastic_atoms, stochastic_region, seed, cutoff_dict):
        """
        Initializes a spyss cell.
        """

        self.initial_atoms_file = initial_atoms_file
        self.stochastic_atoms_files = stochastic_atoms_files
        self.number_of_stochastic_atoms = number_of_stochastic_atoms
        self.stochastic_region = stochastic_region
        self.seed = seed
        self._max_iterations = 100000
        self.cutoff_dict = cutoff_dict

        with open(self.initial_atoms_file) as tmp_file:
            self.initial_atoms = read(tmp_file)

        self.stochastic_atoms_list = []
        for stochastic_atoms_file in stochastic_atoms_files:
            with open(stochastic_atoms_file) as tmp_file:
                self.stochastic_atoms_list.append(read(tmp_file))

    def get_initial_atoms(self):
        """
        Returns the ASE atoms object of the initial atoms.
        """

        return self.initial_atoms.copy()


    def get_stochastic_atoms(self):
        """
        Returns a single ASE atoms object with the stochastic atoms.
        """

        i = 0
        for x in range(self._max_iterations):
            
            self.stochastic_atoms = ase.Atoms()

            for index in range(len(self.stochastic_atoms_list)):

                for index_stochastic_atoms in range(self.number_of_stochastic_atoms[index]):

                    stochastic_atoms = self.stochastic_atoms_list[index].copy()
                    # translates the centre of maass of the atoms
                    random_displacement = np.asarray([np.random.rand(), np.random.rand(), np.random.rand()])
                    random_displacement = self.stochastic_region.transformation_function(random_displacement)
                    stochastic_atoms.translate(random_displacement)
                    del random_displacement

                    # rotates the set of atoms about their centre-of-mass
                    random_angle = np.random.rand() * 180.0
                    random_direction = [np.random.rand(), np.random.rand(), np.random.rand()]
                    stochastic_atoms.rotate(random_angle, random_direction, center='COM')
                    del random_angle
                    del random_direction

                    self.stochastic_atoms += stochastic_atoms.copy()
                    i += 1

            if check_for_overlap(self.stochastic_atoms, self.cutoff_dict) is True:
                continue
            else:
                print ('Generated structure(s) after', i, ' tries.')
                break

        return self.stochastic_atoms.copy()


