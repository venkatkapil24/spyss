class atoms:
    """
    A sypyss atoms object is described by an ASE atoms object 
    of initial atoms, a list of ASE atoms objects of stochastic 
    atoms, and a list of the number of individual stochastic atoms, 
    a random number seed, and the region where stochastic atoms 
    are to be added. 
    """

    def __init__(self, initial_atoms, stochastic_atoms, number_of_stochastic_atoms, stochastic_region, seed):
        """
        Initializes a spyss cell.
        """

        self.initial_atoms = initial_atoms
        self.stochastic_atoms = stochastic_atoms
        self.number_of_stochastic_atoms = number_of_stochastic_atoms
        self.stochastic_region = stochastic_region
        self.seed = seed
