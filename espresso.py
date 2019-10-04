from ase.io import read, write
from ast import literal_eval
import glob
import re
import os
import shutil
import numpy as np
from ase.calculators.espresso import Espresso
from monolayer import *
import datetime

class espresso:

    #input_data: text file include dictionary type data for espresso input

    def __init__(self, input_data, config, pseudo_dir, PBS_header):

        with open(input_data, 'r') as inf:
            self.input_data=literal_eval(inf.read())

        print('espresso base setup')
        print(self.input_data)

        #constracting psuedopotentials dictionary from psuedo_dir
        #assume elemental name_option.upf(# UPF)
        self.pseudopotentials={}
        files=glob.glob(pseudo_dir+"/*")
        for f in files:
            sp=re.split('[._]',os.path.basename(f))
            self.pseudopotentials[sp[0]]=os.path.basename(f)

        print('avairable pseudopotentials')
        print(self.pseudopotentials)

        with open(config, 'r') as inf:
            self.config=literal_eval(inf.read())

        with open(PBS_header, 'r') as header:
            self.PBS_header=header.read()

        print('PBS header')
        print(self.PBS_header)


    #minor patch program to fix pseudo potential name in pw-input
    def rename_psuedo(self, input_file):

        with open(input_file, "r") as input:
            lines = input.readlines()

        datas = [l.rstrip('\n') for l in lines]

        for pot in self.pseudopotentials.keys():
            for l in range(len(datas)):
                if (datas[l].find(pot + "_dummy.UPF")) > 0:
                    #print('find')
                    ll = datas[l].replace(pot + "_dummy.UPF", self.pseudopotentials[pot])
                    datas[l] = ll

        shutil.copyfile(input_file, input_file+'.org')

        with open(input_file, "w") as input:
            for l in datas:
                input.write(l + '\n')



    def optimize_lattice_constant(self, initial_lattice_constant, range, samples, kpts, vacuum):
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

        #add relax tag
        self.input_data['control']['calculation']='relax'

        calc = Espresso(input_data=self.input_data, psuedopotentials=self.pseudopotentials, kpts=kpts,
                        label=self.config['formula'])

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

            calc.write_input(atoms)
            input_file_name=self.config['formula']+".pwi"
            self.rename_psuedo(input_file_name)
            shutil.copyfile(input_file_name , os.path.join(dirname,self.config['formula']+".pwi"))

        input_file_name=self.config['formula']+".pwi"
        #generate run script
        with open("run.sh" ,"w") as f:
            f.write(self.PBS_header+'\n')
            directory=self.config['formula']+"-a-"
            f.write("for dir in ./"+directory+"*; do \n")
            f.write("echo $dir \n")
            f.write("cd $dir \n")
            f.write(self.config['mpicommand'] +" "+self.config['path']+"/pw.x " +self.config['qeoption']+ " -input " + input_file_name + " > scf.out \n")
            f.write("cd ../ \n")
            f.write("done \n")

        #create log file
        logfilename="log_lattice_optim-"+str(datetime.datetime.now())
        with open(logfilename,"w") as log:
            for d in dirs:
                log.write(d+"\n")

        return

    def gather_optimize_results(self,logfile):
        with open(logfile,'r') as f:
            dirs=f.readlines()

        results=[]

        for d in dirs:
            #parse directory name
            #directory name is prefix-a-latticeConstant
            latticeConstant=float(d.split('-')[-1].rstrip('\n'))

            #in optimize_lattice_constant, standard output is set to scf.out
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

        #create scf calculation input from most stable structure
        energySort=sorted(results,key=lambda x:(x[1]))
        targetdir=self.config['formula']+"-a-"+str(energySort[0][0])
        mostStable=os.path.join(targetdir,'scf.out')
        print(mostStable)

        relaxresult=read(mostStable,format='espresso-out')
        #change tag to scf
        self.input_data['control']['calculation']='scf'
        #create dummy input
        calc = Espresso(input_data=self.input_data, psuedopotentials=self.pseudopotentials,
                        kpts=(16,16,1),label=self.config['formula']+'-optimized')
        #kpoints & psuedo is dummy!!
        calc.write_input(relaxresult)

    def create_band_input(self,scf_file, band_config):
        '''
        :param scf_file: pw.x input file after relaxation
        :param band_config: setup for band calculation
        :return:
        '''

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
        path=scf.cell.bandpath(sp_inplane, npoints=100)

        self.input_data['control'].update({'calculation': 'scf'})
        calc = Espresso(input_data=self.input_data, psuedopotentials=self.pseudopotentials,
                        kpts=(24,24,1),label=self.config['formula']+'-scf')
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

        #generate run script
        with open("run.sh" ,"w") as f:
            f.write(self.PBS_header+'\n')
            #scf
            f.write(self.config['mpicommand'] +" "+self.config['path']+"/pw.x " +self.config['qeoption']+ " -input " + scf_file_name + " > scf.out \n")
            #then band
            f.write(self.config['mpicommand'] + " " + self.config['path'] + "/pw.x " + self.config[
                'qeoption'] + " -input " + band_file_name + " > band.out \n")

















