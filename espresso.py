from ase.io import read, write

class espresso:

    #self.system: ASE Atoms
    def __init__(self, templete):
        self.templete=read(templete,format='espresso.in')

    def optimize_lattice_constant(self):
        return



