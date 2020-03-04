from ast import literal_eval
import glob
import re
import os
import shutil
import numpy as np
from layeredMaterialToolKit.monolayer import *
import datetime
import spglib
from ase.calculators.vasp import Vasp

class vasp:

    #input_data: text file include dictionary type data for INCAR

    def __init__(self, input_data, config, PBS_header):
        with open(input_data, 'r') as inf:
            self.input_data=literal_eval(inf.read())

        print('VASP base setup')
        print(self.input_data)

        with open(config, 'r') as inf:
            self.config=literal_eval(inf.read())

        with open(PBS_header, 'r') as header:
            self.PBS_header=header.read()

        print('PBS header')
        print(self.PBS_header)

    def generate_different_lattice_constant_inputs(self, initial_lattice_constant, range, samples, kpts, vacuum):
        '''

        :param initial_lattice_constant: center value to search the lattice constant
        :param range: percentage to search the lattice constant
        :param samples: number of samples
        :param kpts: k-point sampling like (12,12,1)
        :param vacuum: unit cell size to z direction
        :return:
        '''


        #in MX2, I set initial guess of thickness of the layer = initial_lattice_constant
        #this is Ok for typical MX2 like MoS2, WSe2, but be careful


        #set range to seardh
        start=initial_lattice_constant-initial_lattice_constant*range/100.0
        last=initial_lattice_constant+initial_lattice_constant*range/100.0

        latta = np.linspace(start, last, samples)
        dirs=[]
        for a in latta:
            dirname=self.config['formula']+"-a-"+str(a)
            dirs.append(dirname)
            os.makedirs(dirname, exist_ok=True)
            atoms=None

            if(self.config['type']=='MX2'):
                atoms=monolayer_MX2(self.config['formula'], a, initial_lattice_constant, vacuum)
            elif(self.config['type'] == 'Xene'):
                atoms = monolayer_Xene(self.config['formula'], a, self.config['buckling'], vacuum)
            #write('CONTCAR',atoms,format='vasp', vasp5=True)
            write(str(dirname) + '/POSCAR', atoms, format='vasp', vasp5=True)
            if(self.input_data['luse_vdw']):
                calc=Vasp(xc=self.input_data['xc'],encut=self.input_data['encut'],
                          kpts=kpts,ismear=self.input_data['ismear'],
                          luse_vdw=self.input_data['luse_vdw'],
                          ibrion=self.input_data['ibrion'],
                          potim=self.input_data['potim'],
                          nsw=self.input_data['nsw'],
                          gga=self.input_data['gga'],
                          param1=self.input_data['param1'],
                          param2=self.input_data['param2'],
                          aggac=self.input_data['aggac'],
                          zab_vdw=self.input_data['zab_vdw'],
                          lasph=self.input_data['lasph'],
                          sigma=self.input_data['sigma'],setups='recommended')
                shutil.copy('vdw_kernel.bindat',str(dirname)+'/vdw_kernel.bindat')
            else:
                calc = Vasp(xc=self.input_data['xc'], encut=self.input_data['encut'], kpts=kpts,
                            ismear=self.input_data['ismear'],
                            ibrion=self.input_data['ibrion'],
                            potim=self.input_data['potim'],
                            nsw=self.input_data['nsw'],
                            sigma=self.input_data['sigma'], setups='recommended')
            atoms.set_calculator(calc)
            calc.initialize(atoms)
            calc.write_kpoints(directory=str(dirname))
            calc.write_potcar(directory=str(dirname))
            calc.write_incar(atoms,directory=str(dirname))


            #shutil.copyfile(INCAR, os.path.join(dirname,INCAR))

        #generate run script
        with open("run.sh" ,"w") as f:
            f.write(self.PBS_header+'\n')
            directory=self.config['formula']+"-a-"
            f.write("for dir in ./"+directory+"*; do \n")
            f.write("echo $dir \n")
            f.write("cd $dir \n")
            f.write(self.config['mpicommand'] +" "+self.config['path']+"/vasp_std " + " > scf.out \n")
            f.write("cd ../ \n")
            f.write("done \n")

        #create log file
        logfilename="log_lattice_optim-"+str(datetime.datetime.now())
        with open(logfilename,"w") as log:
            for d in dirs:
                log.write(d+"\n")

        return

    #TODO not updated for VASP
    def gather_optimize_results(self,logfile):
        with open(logfile,'r') as f:
            dirs=f.readlines()

        results=[]

        for d in dirs:
            #parse directory name
            #directory name is prefix-a-latticeConstant
            latticeConstant=float(d.split('-')[-1].rstrip('\n'))

            #in generate_different_lattice_constant_inputs, standard output is set to scf.out
            scfout=os.path.join(d.rstrip('\n'),'scf.out')
            with open(scfout,'r') as f:
                lines=f.readlines()

            final = [s for s in lines if re.match('\!.+total.*', s)]
            tf = final[-1].strip("\n").rstrip("Ry").split("=")
            fene = float(tf[1])

            results.append([latticeConstant,fene])

        sortedData = sorted(results, key=lambda x: (x[0]))

        fdata = open("results.txt", "w")
        fdata.write("# a   energy(Ry) \n")
        for data in sortedData:
            fdata.write(str(data[0]) + "  " + str(data[1]) +  " \n")
        fdata.close()



    def create_band_input(self,scf_file, band_config_file):
        '''
        :param scf_file: pw.x input file after relaxation
        :param band_config: setup for band calculation
        :return:
        '''

        band_config={}
        with open(band_config_file, 'r') as inf:
            band_config=literal_eval(inf.read())

        scf=read(scf_file,format='espresso-in')
        lat=scf.cell.get_bravais_lattice()
        #get special point list
        sps=lat.get_special_points()

        #remove point in z-direction (because 2D system)
        keys=sps.keys()
        sps_copy=sps.copy()
        for key in keys:
            if sps[key][2] > 0.0:
                del sps_copy[key]


        sp_inplane=''.join(list(sps_copy))
        #better to finish by Gamma-point
        sp_inplane=sp_inplane+'G'
        print("band-path:"+sp_inplane)
        path=scf.cell.bandpath(sp_inplane, npoints=band_config['nbandkpts'])

        self.input_data['control'].update({'calculation': 'scf'})
        calc = Espresso(input_data=self.input_data, psuedopotentials=self.pseudopotentials,
                        kpts=band_config['nscfkpts'],label=self.config['formula']+'-scf')
        calc.write_input(scf)
        scf_file_name = self.config['formula'] + "-scf.pwi"
        self.rename_psuedo(scf_file_name)



        self.input_data['control'].update({'calculation':'bands',
                                           'restart_mode':'restart',
                                           'verbosity':'high'})
        calc = Espresso(input_data=self.input_data, psuedopotentials=self.pseudopotentials,
                        kpts=path,label=self.config['formula']+'-band')
        calc.write_input(scf)
        band_file_name = self.config['formula'] + "-band.pwi"
        self.rename_psuedo(band_file_name)


        #input file for bands.x
        with open('bandsx.in','w') as bf:
            bf.write("&bands \n")
            bf.write("outdir=" + "\'"+self.input_data['control']['outdir'] +"\' \n")
            bf.write("filband=" + "\'"+ self.config['formula']+".band\' \n")
            bf.write("lsym=" + "."+str(band_config['lsym']).lower()+". \n")
            bf.write("no_overlap=" + "." + str(band_config['no_overlap']).lower() + ". \n")
            bf.write("/")



        #generate run script
        with open("run.sh" ,"w") as f:
            f.write(self.PBS_header+'\n')
            #scf
            f.write(self.config['mpicommand'] +" "+self.config['path']+"/pw.x " +self.config['qeoption']+ " -input " + scf_file_name + " > scf.out \n")
            #then band
            f.write(self.config['mpicommand'] + " " + self.config['path'] + "/pw.x " + self.config[
                'qeoption'] + " -input " + band_file_name + " > band.out \n")
            f.write(self.config['mpicommand'] + " " + self.config['path'] + "/bands.x " + self.config[
                'qeoption'] + " -input " + "bandsx.in" + " > bandsx.out \n")
