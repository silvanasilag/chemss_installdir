#!/usr/local/anaconda3/bin/python
# -*- coding: utf-8 -*-
#programa que extraeel xyz optimizado de un out de gausian


class AtomNMR:
    def __init__(self, atomic_symbol,xx,yy,zz):
        self.s=atomic_symbol
        self.x=xx
        self.y=yy
        self.z=zz
#------------------------------------------------------------------------------------------
class MoleculeNMR:
    def __init__(self,name_molecule):
        self.atoms = []
        self.n = 0
        self.i_m = name_molecule
    def add_atom(self, atom):
        self.atoms.append(atom)
        natoms=len(self.atoms)
        self.n = natoms
    def __repr__(self):
        self=sorted(self, key=lambda atomm: atomm.s) 
#------------------------------------------------------------------------------------------
def get_geo(out):
    fg=0
    fg_nl=[]
    fg_n=0
    with open(out, 'r') as f2:    
        for line in f2:
            line = line.strip()
            lin = line.split()
            if len(lin)>0:
                if fg_n==1:fg_nl.append(str(line))
                if "******************************************" in line: fg_n=1
                if "Symbolic Z-matrix:" in line: 
                    fg_n=0
                    # print(fg_nl[-3])
                    name=fg_nl[-3]
                
                if len(lin)==6 and fg==3:
                    # print(lin)
                    s=str(lin[1])
                    x=float(lin[3])
                    y=float(lin[4])
                    z=float(lin[5])
                    ai=AtomNMR(s,x,y,z)
                    mol0.add_atom(ai)
                if "---------------------------------------------------" in line and fg==3:
                    fg=0
                if "---------------------------------------------------" in line and fg==2:fg=3 
                if "Standard orientation:" in line: 
                    fg=1 #primero en aparecer
                    mol0=MoleculeNMR(name)
                if "---------------------------------------------------" in line and fg==1:fg=2
    return mol0

#------------------------------------------------------------------------------------------
def write_xyz(xyz,mol):
    with open(xyz,'w') as f2:
        # print(os.getcwd())
        f2.write(mol.i_m)
        f2.write("\n  \n")
        f2.write("0   1\n")
        for i, iatom in enumerate(mol.atoms):
             f2.write('%s  %16f  %11f %11f \n' % ( iatom.s,iatom.x, iatom.y, iatom.z))
