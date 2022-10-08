import ase
from ase.io import read
from ase.optimize import BFGS
#from ase.constraints import FixAtoms

class optimizer:
    """
    A sypyss optimizer object is described by an ASE calculator
    and an optimization algorithm. 
    """

    def __init__(self, ASE_calculator, optimization_algorithm, optimization_parameter_dict):
        """
        Initializes a spyss optimizer.
        """

        self.ASE_calculator = ASE_calculator
        self.optimization_algorithm = optimization_algorithm
        self.optimization_parameter_dict = optimization_parameter_dict
        
        self._implemented_optimization_algorithms = {'BFGS' : BFGS}
        
        if self.optimization_algorithm not in self._implemented_optimization_algorithms:
            raise NotImplementedError('Currently only allowing these algorithms: ', self._implemented_optimization_algorithms.keys())
              
    def optimize(self, structure):
        """
        Optimizes a structure.
        """
        
        tmpatoms = structure.copy()
        tmpatoms.set_calculator(self.ASE_calculator)
        
#         c = FixAtoms(indices=[atom.index for atom in tmpatoms])
#         tmpatoms.set_constraint(c)
        
        tmpobject = self._implemented_optimization_algorithms[self.optimization_algorithm](tmpatoms)
        tmpobject.run(**self.optimization_parameter_dict)       
        tmpatoms.info['energy'] = tmpatoms.get_potential_energy()
        del tmpobject
        
        return tmpatoms.copy()
        
        
        
        
