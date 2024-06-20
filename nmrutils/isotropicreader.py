import sys
import os
from os import path
from glomos_utils.get_geometry  import get_geometry_gaussian
from glomos_utils.connectivity import conectmx
import glob

# import datetime
class Level_ot: #Level of theory
    def __init__(self, opt_fun,opt_base, nmr_fun, nmr_base, dispersion, solv_model, solvent):
        self.of=opt_fun
        self.ob=opt_base
        self.nf=nmr_fun
        self.nb= nmr_base
        self.d=dispersion
        self.sm=solv_model
        self.sl=solvent

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
    def __init__(self, name,cpu_time,elapse_time,cpu_time_nmr,elapse_time_nmr,rmsd,mol_id,chk):
        self.atoms = []
        self.i = name
        self.n = 0
        self.ct=cpu_time
        self.et=elapse_time
        self.ct_nmr=cpu_time_nmr
        self.et_nmr=elapse_time_nmr
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

def flotante(variable):
    try:
        float(variable)
        return True
    except:
        return False

def extrac_time(lin):
    (d,h,m,s)=(int(lin[3]), int(lin[5]), int(lin[7]), float(lin[9]))
    return (d,h,m,s)

def from_reader(txt, out, mol_id):
    oslist, otlist=[], []
    time=0
    timecpu=(0,0,0,0)
    timelaps=(0,0,0,0)
    timecpu_nmr=(0,0,0,0)
    timelaps_nmr=(0,0,0,0)
    fl=0
    chk=0
    with open(out, 'r') as f2:    #abre el out
        for line in f2:
            line = line.strip()
            lin = line.split()
            os = 'X'
            if "#" and "OPT" and "Geom=Check" in line: 
                chk=1
            if "Isotropic" in line and len(lin)==8:
                fl=1
                ot=float(lin[4])
                os=str(lin[1])
                oslist.append(os) #appendiza simbolo del atomo
                otlist.append(ot) #appendiza las valor isotropico
            if "Job cpu time:" in line and fl==0: timecpu=extrac_time(lin)
            if "Elapsed time:" in line and fl==0: 
                lin.insert(0,"-")
                timelaps=extrac_time(lin)
            if "Job cpu time:" in line and fl==1: timecpu_nmr=extrac_time(lin)
            if "Elapsed time:" in line and fl==1: 
                lin.insert(0,"-")
                timelaps_nmr=extrac_time(lin)
    #TERMINACIONNNNN
    if fl==0:
        Terminacin ()
    #TERMINACIONNN

    # -----VECINOSSS
    x = get_geometry_gaussian("",out,0)
    xmtx = conectmx(x)
    nbs = neighbor_finder(xmtx)
    print(nbs)
    # -----VECINOSSS

    nz,fl,simb=0,0,[]
    with open(txt, 'r') as f:  #txt
        print(txt)
        for line in f:
            line = line.strip()
            lin = line.split()
            if fl == 0:
                name=line
                mol0=MoleculeNMR(name,timecpu,timelaps,timecpu_nmr,timelaps_nmr,0.0,mol_id,chk)
                fl=1
            if len(lin)>1 and lin != 0:
                print(line)
                if flotante(lin[1])==True or "=" in line or "Nan" in line:
                    nz=nz+1
                    s=str(lin[0])
                    ex=str(lin[1])
                    #if s =="H" or s=="C":
                    print(nz,s,ex,nbs[nz])
                    ai=AtomNMR(nz,s,ex,0.0,0.0,0.0,nbs[nz])
                    mol0.add_atom(ai)
                    simb.append(s)
                print(line)
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
        print(iatom.nz,"",iatom.s,":",iatom.nb)
    return {'mol':mol1}
#----------------------------------------------------------------------------------------------
def filter_data(mol0): 
    av,ex,iso,s,nz,nn,nbh=[],[],[],[],[],[],[] #l
    is_v=-200
    name=""
    mol1=MoleculeNMR(mol0.i,mol0.ct,mol0.et,mol0.ct_nmr,mol0.et_nmr,mol0.rmsd,mol0.im,mol0.chk)
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
            a=AtomNMR(name[:-1],s[i],float(ex[i]),is_v,0.0,0.0,nbh[i])  
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
    return moleculeout,"borrar"
#------------------------------------------------FILTRAR
