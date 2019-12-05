from ase.build import mx2
from ase import Atoms
from ase.io import read, write
from ase.calculators.espresso import Espresso
import math

def monolayer_MX2(formula, a, thickness,vacuum):
    crystal=mx2(formula=formula, kind='2H', a=a, thickness=thickness, size=(1,1,1), vacuum=None)
    #print(crystal.positions)
    #convert to all PBC cell
    slab=Atoms(crystal)
    slab.set_cell([crystal.cell[0],crystal.cell[1],[0.0,0.0,vacuum]], scale_atoms=False)
    slab.set_pbc([True, True, True])
    return slab

def monolayer_Xene(formula, a, buckling, vacuum):
    if vacuum is not None:
        buckling=buckling/vacuum

    positions=[[2/3, 1/3,buckling/2.0],[1/3,2/3,-buckling/2.0]]
    cell=[[a, 0, 0], [-a/2, a * 3**0.5 / 2, 0], [0, 0, 0]]
    atoms = Atoms(formula, positions=positions, cell=cell, pbc=(1, 1, 0))
    atoms.set_scaled_positions(positions)
    if vacuum is not None:
        atoms.center(vacuum, axis=2)
    return atoms




