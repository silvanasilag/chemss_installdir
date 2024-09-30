import sys
import os
import glob
import datetime
from os import path
from glomos_utils.glomos_utils import get_geometry_and_energy_gaussian,neighbor_finder,conectmx
from nmrutils.terminacion     import terminacion,key_compare,key_compare_all
from runrun.qpbs import printf

class AtomNMR:
    def __init__(self,num_in_mol,atomic_symbol,signal_exp, signal_theo, residuo, correccion,neighbors):
        self.nz=num_in_mol
        self.s=atomic_symbol
        self.e=signal_exp
        self.t=signal_theo
        self.r=residuo
        self.c=correccion
        self.nb=neighbors
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
    return d,h,m,s

def from_reader(dnmr, out, mol_id,filesout,key_opt,key_nmr):
    # Function that opens the Gaussian output (.out) and NMR experimental data (.dnmr) file simultaneously
    oslist, otlist=[],[]
    timecpu,timecpu_nmr=(0,0,0,0),(0,0,0,0)
    chk,fl,fl2=0,0,0
    with open(out, 'r') as f2:    #opens goussian out
        nmr,opt="-","-"
        full_line = ""
        for line in f2:
            line = line.strip()
            lin = line.split()
            if full_line != "":
                fl2=4
                full_line += line
                line = full_line
                full_line=''
            if "iop(" in line.lower() and not line.endswith(")") and nmr=="-":
                full_line += line
                fl2=2
            if "#" in line:
                if "OPT" in line and opt=="-":
                    opt =line.strip(' Guess=Read')
                    opt =opt.strip(' Geom=Check')
                    if  "Geom=Check" in line:
                        chk=1
                if "NMR" in line and nmr=="-":
                    if fl2 == 2:
                        pass
                    else:
                        nmr = line.strip(' Guess=Read')
                        nmr = nmr.strip(' Geom=Check')
                if "IOP" in line.upper():
                    if fl2 == 2:
                        pass
                    else:
                        if "(3/76=1000007400,3/77=0999900001,3/78=0000109999)" in line:
                            opt = opt.replace("BLYP", "WC04")
                            nmr = nmr.replace("BLYP", "WC04") if nmr != "-" else nmr
                        elif "(3/76=1000001189,3/77=0961409999,3/78=0000109999)" in line:
                            opt = opt.replace("BLYP", "WP04")
                            nmr = nmr.replace("BLYP", "WP04") if nmr != "-" else nmr
                        full_line = ""
                        fl2=0
            if "Isotropic" in line and len(lin)==8:
                fl=1
                ot=float(lin[4])
                os=str(lin[1])
                oslist.append(os) #atom symbol
                otlist.append(ot) #isotropic value
            if "Job cpu time:" in line and fl==0:
                timecpu=extrac_time(lin)
            if "Job cpu time:" in line and fl==1:
                timecpu_nmr=extrac_time(lin)
    key_opt,key_nmr=key_compare(key_opt,key_nmr,nmr,opt,out)
    if fl==0:terminacion(out,filesout)
    nz,fl,simb=0,0,[]
    # -----VECINOSSS
    x=get_geometry_and_energy_gaussian(out,1,'F')
    xmtx = conectmx(x)
    nbs = neighbor_finder(xmtx)
    # -----VECINOSSS
    with open(dnmr, 'r') as f:  #opens dnmr file
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
                    ai=AtomNMR(nz,s,ex,0.0,0.0,0.0,nbs[nz])
                    mol0.add_atom(ai)
                    simb.append(s)
    for i, iatom in enumerate(mol0.atoms):
        if oslist[i] != iatom.s:
            print('-------------------------ERROR!!',mol0.im)
            print(oslist[i],iatom.s)
        else:
            iatom.t = otlist[i]
        #print(iatom.s,iatom.e,iatom.t,iatom.r,iatom.c)
    mol1 = filter_data(mol0)
    for i, iatom in enumerate(mol1.atoms):
        nnb = [simb[k - 1] for k in iatom.nb]
        iatom.nb = nnb
        if iatom.s == "H" and iatom.e > 15:
            sys.exit("Error")
    return {'mol':mol1,'kys':[opt,nmr],'keys':[key_opt,key_nmr]}
#----------------------------------------------------------------------------------------------
def filter_data(mol0):
    av,ex,iso,s,nz,nn,nbh=[],[],[],[],[],[],[] #l
    is_v=-200
    name=""
    mol1=MoleculeNMR(mol0.i,mol0.ct,mol0.ct_nmr,mol0.rmsd,mol0.im,mol0.chk)
    for i,iatom in enumerate(mol0.atoms):
        ex.append(iatom.e)
        iso.append(iatom.t)
        s.append(iatom.s)
        nn.append(iatom.nz)
        nbh.append(iatom.nb)
    for j in nn:
        i=j-1
        if ex[i]=="Nan":pass
        #---------prueva
        elif "=" in ex[i]:
            va=int(ex[i].split("=")[1])
            valor=ex[va-1]
            if flotante(valor)==False:
                print()
                sys.exit("Potential error in chemical equivalence '=' missplace")
        #---------------
        else:
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
            a=AtomNMR(name[:-1],s[i],float(ex[i]),is_v,0.0,0.0,nbh[i])
            av,nz=[],[]
            name=""
            mol1.add_atom(a)
    return mol1

#------------------------------------------------MAIN------------MAIN------MAIN-----------
def molecules_data(path1,path2):
    filesout = []
    filesdnmr = []
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
        path_dnmr = path2 + "/" +ifile
        filesdnmr.append(path_dnmr)
    filesout.sort()
    filesdnmr.sort()
    namemol.sort()
    if len(filesout) != len(filesdnmr):
        print("Missing files")
        print(filesout,len(filesout))
        print(filesdnmr,len(filesdnmr))
        sys.exit(1)
    for i, j, k in zip(filesdnmr,filesout,namemol):
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
