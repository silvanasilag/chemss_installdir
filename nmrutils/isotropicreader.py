import sys
import os
from os import path
import glob
import datetime
from nmrutils.terminacion     import terminacion,key_compare,key_compare_all
"""
class Level_theory:
    def __init__(self,opt_fun,opt_base,nmr_fun,nmr_base, dispersion, solvent, solv_method):
        self.opf = opt_fun
        self.opb = opt_base
        self.nmrf= nmr_fun
        self.nmrb = nmr_base
        self.g0 = dispersion
        self.sol= solvent
        self.sol_m= solv_method
"""

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
    def __init__(self, name,cpu_time,cpu_time_nmr,rmsd,mol_id,chk):
        self.atoms = []
        self.i = name
        self.n = 0
        self.ct=cpu_time
        self.ct_nmr=cpu_time_nmr
        self.rmsd=rmsd
        self.im=mol_id
        self.chk=chk
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

def extrac_time(lin):
    (d,h,m,s)=(int(lin[3]), int(lin[5]), int(lin[7]), float(lin[9]))
    return (d,h,m,s)

def from_reader(txt, out, mol_id,filesout,key_opt,key_nmr):
    oslist, otlist=[],[]
    time=0
    timecpu=(0,0,0,0)
    timecpu_nmr=(0,0,0,0)
    fl=0
    chk=0
    with open(out, 'r') as f2:    #abre el out
        nmr,opt="-","-"
        for line in f2:
            line = line.strip()
            lin = line.split()
            os = 'X'
            if str("#")  and str("OPT") in line and opt=="-":
                opt =line
                if  "Geom=Check" in line:
                    chk=1
                    opt=opt.strip(' Guess=Read')
                    opt=opt.strip(' Geom=Check')
            if str("#") and str("NMR")in line and nmr=="-":
                nmr=line
                nmr=nmr.strip(' Guess=Read')
                nmr=nmr.strip(' Geom=Check')
            if str("#") and str("IOP")in line.upper() : 
                if "(3/76=1000007400,3/77=0999900001,3/78=0000109999)" in line: w="WC04"
                elif "(3/76=1000001189,3/77=0961409999,3/78=0000109999)" in line: w="WP04"
                else: w="Iop"
                if nmr == "-": opt = opt=opt.replace("BLYP",w) 
                if opt != "-" and nmr != "-": nmr =nmr.replace("BLYP",w) 
            if "Isotropic" in line and len(lin)==8:
                fl=1
                ot=float(lin[4])
                os=str(lin[1])
                oslist.append(os) #appendiza simbolo del atomo
                otlist.append(ot) #appendiza las valor isotropico
            if "Job cpu time:" in line and fl==0:
                timecpu=extrac_time(lin)
            if "Job cpu time:" in line and fl==1: 
                timecpu_nmr=extrac_time(lin)

    key_opt,key_nmr=key_compare(key_opt,key_nmr,nmr,opt,out)
    if fl==0:terminacion(out,filesout)
    nz,fl=0,0
    with open(txt, 'r') as f:  #txt
        for line in f:
            line = line.strip()
            lin = line.split()
            if fl == 0:
                name=line
                mol0=MoleculeNMR(name,timecpu,timecpu_nmr,0.0,mol_id,chk)
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
    return {'mol':mol1,'kys':[opt,nmr],'keys':[key_opt,key_nmr]}
#----------------------------------------------------------------------------------------------
def filter_data(mol0): 
    av,ex,iso,s,nz,nn=[],[],[],[],[],[] #l
    is_v=-200
    name=""
    mol1=MoleculeNMR(mol0.i,mol0.ct,mol0.ct_nmr,mol0.rmsd,mol0.im,mol0.chk)
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
    key_opt,key_nmr=[],[]
    key_opt_f,key_nmr_f=[],[]
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
        reader=from_reader(i, j, k,filesout,key_opt,key_nmr)
        keys=reader['kys']
        molx = reader['mol']
        key_opt,key_nmr=reader['keys'][0],reader['keys'][1]
        moleculeout.extend([molx])
        key_opt_f.append((keys[0],k))
        key_nmr_f.append((keys[1],k))
    if len(key_opt)!=1 : key_compare_all(key_opt,key_opt_f)
    if len(key_nmr)!=1 :key_compare_all(key_nmr,key_nmr_f)
    return moleculeout,keys
#------------------------------------------------FILTRAR
