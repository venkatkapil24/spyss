from ase.io import read
from cell import cell

spyss_cell = cell(structure_file='square.pdb', mode='fixed-cell')

assert spyss_cell.structure_file == 'square.pdb'
assert spyss_cell.mode == 'fixed-cell'
assert spyss_cell.volume == read('square.pdb').get_volume()
print ('test passed.')
