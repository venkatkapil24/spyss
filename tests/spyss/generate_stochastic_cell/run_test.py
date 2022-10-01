from ase.io import read
from cell import cell
from spyss import spyss

spyss_cell = cell(structure_file='square.pdb', mode='fixed-cell')
spyss_object = spyss(cell=spyss_cell, atoms=None, seed=31514, count=None, optimizer=None)
spyss_object.generate_stochastic_cell()

assert spyss_object.stochastic_cells[0] == read('square.pdb')
print ('test passed.')
