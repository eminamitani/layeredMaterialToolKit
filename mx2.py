from ase.build import mx2
from ase import Atoms
from ase.io import read, write
from ase.calculators.espresso import Espresso

def monolayer_MX2(formula, a, thickness,vacuum):
    crystal=mx2(formula=formula, kind='2H', a=a, thickness=thickness, size=(1,1,1), vacuum=None)
    print(crystal.positions)
    #convert to all PBC cell
    slab=Atoms(crystal)
    slab.set_cell([crystal.cell[0],crystal.cell[1],[0.0,0.0,vacuum]], scale_atoms=False)
    slab.set_pbc([True, True, True])
    return slab

if __name__ == '__main__':
    test=monolayer_MX2('MoS2', 3.3,3.3,20.0)
    write("test.vasp", test, format='vasp', vasp5=True)
    input_data = {
        'system': {
            'ecutwfc': 64,
            'ecutrho': 576},
        'disk_io': 'low'}

    calc=Espresso(input_data=input_data)
    calc.write_input(test)

