import os.path
import numpy as np
from glomos_utils.atomic import get_chemical_symbol,get_covalent_radius

#------------------------------------------------------------------------------------------
hartree2kcalmol=float(627.509391)
#------------------------------------------------------------------------------------------
class Atom:
    def __init__(self, symbol, xcartesian, ycartesian, zcartesian, xfix='T', yfix='T', zfix='T'):
        self.s=symbol
        self.xc=xcartesian
        self.yc=ycartesian
        self.zc=zcartesian
        self.xf=xfix
        self.yf=yfix
        self.zf=zfix
    def vec(self):
        vector=np.array([self.xc,self.yc,self.zc])
        return vector
#------------------------------------------------------------------------------------------
class Molecule:
    def __init__(self, name, energy, matrix=[], comments=[]):
        self.atoms = []
        self.i = name
        self.e = energy
        self.m = matrix
        self.c = comments
        self.n = 0
    def add_atom(self, atom):
        self.atoms.append(atom)
        natoms=len(self.atoms)
        self.n = natoms
    def clear(self):
        delattr(self, 'i')
        delattr(self, 'e')
        delattr(self, 'm')
        delattr(self, 'c')
        delattr(self, 'n')
        delattr(self, 'atoms')
def neighbor_finder(adj_mtx):
    '''
    This function builds a dictionary that contains the neighbors of all atoms

    in: adj_mtx (list); a list of list elements that represent the adjacency matrix
        of the molecule as a graph.
    out: dict_neig (dict); a dictionary containing as keys the position of the atom
        in the molecule.atoms list, and as values a list of all its neighbors.
        ({0:[1,3,7,11], 1:[0,7],...})
    '''
    dict_neig = {}
    conta = 0
    for i in adj_mtx:
        neig = []
        contb = 0
        for j in i:
            if j == 1:
                neig.append(contb+1)
            contb = contb + 1
        dict_neig[conta+1] = neig
        # print(conta,neig)
        conta = conta + 1
    return dict_neig

def conectmx(moleculein, factor=1.2):
    natoms=moleculein.n
    matrixc=np.zeros(shape=(natoms,natoms),dtype=np.int)
    for iatom in range(natoms):
        si=moleculein.atoms[iatom].s
        ri=get_covalent_radius(si)
        xi=moleculein.atoms[iatom].xc
        yi=moleculein.atoms[iatom].yc
        zi=moleculein.atoms[iatom].zc
        for jatom in range(iatom+1,natoms):
            sj=moleculein.atoms[jatom].s
            rj=get_covalent_radius(sj)
            xj=moleculein.atoms[jatom].xc
            yj=moleculein.atoms[jatom].yc
            zj=moleculein.atoms[jatom].zc
            xx,yy,zz=xj-xi, yj-yi, zj-zi
            distance = np.sqrt(xx*xx+yy*yy+zz*zz)
            dist = distance/(ri + rj)
            if dist <= factor :
                matrixc[iatom][jatom] = int(1)
                matrixc[jatom][iatom] = int(1)
    return matrixc

def get_geometry_and_energy_gaussian( filename, zpe=1, ccsd='F'):
    # Open and read the entire file into memory
    with open(filename, 'r') as gaufile:
        lines = gaufile.readlines()
    namein = filename.split('.')[0]
    moleculeout = None
    energy = None
    enehartree=float(0.0)
    enehartree_correction=float(0.0)
    for line in lines:
        # Retrieve energy information
        if "SCF Done" in line and ccsd == 'F':
            scf = line.split()
            enehartree = float(scf[4])
        if "CCSD(T)= " in line and ccsd == 'T':
            zerosplit = line.split()
            firstsplit = zerosplit[1].split('D')
            ener = float(firstsplit[0])
            exp = float(firstsplit[1])
            enehartree = float(ener * pow(10, exp))
        if "Zero-point correction=" in line and zpe == 1:
            zpc = line.split()
            enehartree_correction = float(zpc[2])
        if "Thermal correction to Energy=" in line and zpe == 2:
            zpc = line.split()
            enehartree_correction = float(zpc[4])
        if "Thermal correction to Enthalpy=" in line and zpe == 3:
            zpc = line.split()
            enehartree_correction = float(zpc[4])
        if "Thermal correction to Gibbs Free Energy=" in line and zpe == 4:
            zpc = line.split()
            enehartree_correction = float(zpc[6])

        # Retrieve molecular geometry
        if line.strip() in ("Input orientation:", "Standard orientation:"):
            energy = enehartree + enehartree_correction  # Calculate total energy
            moleculeout = Molecule(namein, energy)
            start_index = lines.index(line) + 5  # Skip header lines
            for geom_line in lines[start_index:]:
                if geom_line.startswith(" --------"):
                    break
                ls = geom_line.split()
                if len(ls) == 6 and ls[0].isdigit() and ls[1].isdigit() and ls[2].isdigit():
                    ss = get_chemical_symbol(int(ls[1]))
                    xc, yc, zc = float(ls[3]), float(ls[4]), float(ls[5])
                    ai = Atom(ss, xc, yc, zc)
                    moleculeout.add_atom(ai)
    # Convert energy to kcal/mol
    if energy is not None:
        energy_kcalmol = energy * hartree2kcalmol
        moleculeout.energy = energy_kcalmol  # Update molecule energy

    return moleculeout