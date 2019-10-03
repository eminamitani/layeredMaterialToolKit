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
    positions=[[math.sqrt(3.0)/3.0,0.0,buckling/2.0],[math.sqrt(3.0)/3.0*2.0,0.0,-buckling/2.0]]



