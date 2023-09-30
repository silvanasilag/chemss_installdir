import sys
import os
from os import path
import glob
import datetime


class AtomNMR:
    def __init__(self,num_in_mol,atomic_symbol,signal_exp, signal_theo, residuo, correccion):
        self.nz=num_in_mol
        self.s=atomic_symbol
        self.e=signal_exp
        self.t=signal_theo
        self.r=residuo
        self.c=correccion
#------------------------------------------------------------------------------------------
class MoleculeNMR:
    def __init__(self, name,cpu_time,elapse_time,rmsd,mol_id):
        self.atoms = []
        self.i = name
        self.n = 0
        self.ct=cpu_time
        self.et=elapse_time
        self.rmsd=rmsd
        self.im=mol_id
    def add_atom(self, atom):
        self.atoms.append(atom)
        natoms=len(self.atoms)
        self.n = natoms
    def __repr__(self):
        self=sorted(self, key=lambda atomm: atomm.s) 
#------------------------------------------------------------------------------------------
def flotante(variable):
    try:
        float(variable)
        return True
    except:
        return False

def from_reader(txt, out, mol_id):
    oslist, otlist=[], []
    sum = datetime.timedelta()
    time=0
    (d1,h1,m1,s1)=(0,0,0,0)
    (d2,h2,m2,s2)=(0,0,0,0)
    fl=0
    with open(out, 'r') as f2:    #abre el out
        for line in f2:
            line = line.strip()
            lin = line.split()
            os = 'X'
            if "Isotropic" in line and len(lin)==8:
                ot=float(lin[4])
                os=str(lin[1])
                oslist.append(os) #appendiza simbolo del atomo
                otlist.append(ot) #appendiza las valor isotropico
            if "Job cpu time:" in line:
                (d,h,m,s)=(int(lin[3]), int(lin[5]), int(lin[7]), float(lin[9]))
                (d1,h1,m1,s1)=(d+d1,h+h1,m+m1,s+s1)
            if "Elapsed time:" in line:
                (d,h,m,s)=(int(lin[2]), int(lin[4]), int(lin[6]),float(lin[8]))
                (d2,h2,m2,s2)=(d+d2,h+h2,m+m2,s+s2)         
        timecpu = datetime.timedelta(days=int(d1),hours=int(h1), minutes=int(m1), seconds=round(float(s)))
        timelaps = datetime.timedelta(days=int(d2),hours=int(h2), minutes=int(m2), seconds=round(float(s2)))
    nz=0   
    with open(txt, 'r') as f:  #txt
        for line in f:
            line = line.strip()
            lin = line.split()
            if fl == 0:
                name=line
                mol0=MoleculeNMR(name,timecpu,timelaps,0.0,mol_id)
                fl=1
            if len(lin)>1:
                if flotante(lin[1])==True or "=" in line or "Nan" in line:
                    nz=nz+1
                    s=str(lin[0])
                    ex=str(lin[1])
                    #if s =="H" or s=="C":
                    ai=AtomNMR(nz,s,ex,0.0,0.0,0.0)
                    mol0.add_atom(ai)

    for i, iatom in enumerate(mol0.atoms):   
        if oslist[i] != iatom.s:
            print('-------------------------ERROR!!',mol0.im)
            print(oslist[i],iatom.s)
        else:
            iatom.t = otlist[i]
        #print(iatom.s,iatom.e,iatom.t,iatom.r,iatom.c)
    
    mol1 = filter_data(mol0)
    for i, iatom in enumerate(mol1.atoms):
        if iatom.s == "H" and iatom.e > 15:
            sys.exit("Error")
    return {'mol':mol1}
#----------------------------------------------------------------------------------------------
def filter_data(mol0): 
    av,ex,iso,s,nz,nn=[],[],[],[],[],[] #l
    is_v=-200
    name=""
    mol1=MoleculeNMR(mol0.i,mol0.ct,mol0.et,mol0.rmsd,mol0.im)
    for i,iatom in enumerate(mol0.atoms):
        ex.append(iatom.e)
        iso.append(iatom.t)
        s.append(iatom.s)
        nn.append(iatom.nz)
    for j in nn:
        i=j-1
        if ex[i]=="Nan":pass
        #---------prueva
        elif "=" in ex[i]:
            va=int(ex[i].split("=")[1])
            valor=ex[va-1]
            if flotante(valor)==False:
                sys.exit("Equivocacion en el =")
        #---------------
        else:
            #print(nn[i],"exp:",ex[i],"Iso:",iso[i],s[i])
            av.append(iso[i])
            nz.append(nn[i])
            n1=nn
            for jj in n1:
                ii=jj-1
                if "=" in ex[ii]:
                    exj=ex[ii].split("=")[1]
                    if int(exj) == int(j):
                        av.append(iso[ii])
                        nz.append(nn[ii])
                        if s[i]!= s[ii]:
                            sys.exit("Error")
            is_v=sum(av)/len(av)
            is_v=round(is_v, 4)
            for ii in range(len(nz)):
                    name=name+str(nz[ii])+","
            a=AtomNMR(name[:-1],s[i],float(ex[i]),is_v,0.0,0.0)  
            av,nz=[],[]
            name=""
            mol1.add_atom(a)
    return mol1

#------------------------------------------------MAIN------------MAIN------MAIN-----------
def molecules_data(path1,path2):
    filesout = []
    filestxt = []
    namemol =[]
    moleculeout=[] 
    os.chdir(path1)
    for iifile in glob.glob("*.out"): 
        mol = os.path.splitext(os.path.basename(iifile))[0]
        path_out = path1 + "/" +iifile
        filesout.append(path_out)
        namemol.append(mol)
    os.chdir(path2)
    if len(glob.glob("*.dnmr2"))==1: 
        dnmr=[str(glob.glob("*.dnmr2")[0]) for i in range(len(filesout))]
    else: 
        dnmr=glob.glob("*.dnmr2")
    for ifile in dnmr:
        # mol = os.path.splitext(os.path.basename(ifile))[0]
        path_txt = path2 + "/" +ifile
        filestxt.append(path_txt)
        # namemol.append(mol)
    filesout.sort()
    filestxt.sort()
    namemol.sort()
    if len(filesout) != len(filestxt):
        print("Missing files")
        print(filesout,len(filesout))
        print(filestxt,len(filestxt))
        sys.exit(1)
    for i, j, k in zip(filestxt,filesout,namemol):
        reader=from_reader(i, j, k)
        molx = reader['mol']
        moleculeout.extend([molx])    
    return moleculeout
#------------------------------------------------FILTRAR
