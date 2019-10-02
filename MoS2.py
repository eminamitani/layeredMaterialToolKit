from ase.build import mx2
from ase import Atoms
from ase.io import read, write
from ase.calculators.espresso import Espresso
from mx2 import monolayer_MX2
import os


formula = 'MoS2'
MoS2 = monolayer_MX2(formula, 3.3, 3.3, 20.0)
write("test.vasp", MoS2, format='vasp', vasp5=True)

input_data = {
    'control': {
        'calculation': 'relax',
        'outdir': './',
        'pseudo_dir': './',
        'prefix': formula,
        'wf_collect': True
    },
    'system': {
        'ecutwfc': 70,
        'ecutrho': 560},
}
psuedopotentials = {'Mo': 'Mo.upf',
                    'S': 'S.upf'}

#seems strange, but the psuedopoteitial data is not update
calc = Espresso(input_data=input_data, psuedopotentials=psuedopotentials, kpts=(12, 12, 1), label=formula)

calc.write_input(MoS2)

with open(formula+".pwi", "r") as input:
    lines=input.readlines()

datas=[l.rstrip('\n') for l in lines]

for pot in psuedopotentials.keys():
    print(pot)
    print(psuedopotentials[pot])
    for l in range(len(datas)):
        if (datas[l].find(pot+"_dummy.UPF")) > 0:
            print('find')
            ll=datas[l].replace(pot+"_dummy.UPF", psuedopotentials[pot])
            print(ll)
            datas[l]=ll

with open(formula+".pwi.rev","w") as input:
    for l in datas:
        input.write(l+'\n')


#MoS2.set_calculator(calc)
#MoS2.get_potential_energy()
