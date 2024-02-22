import os.path
from glomos_utils.atomic import get_chemical_symbol
from glomos_utils.libmoleculas import Molecule, Atom
# from nmrutils.getbilparam  import get_a_str
#------------------------------------------------------------------------------------------
log_file =''
hartree2kcalmol=float(627.509391)
#------------------------------------------------------------------------------------------
""" THIS ROUTINE RETURNS 1 (OR 2) IF THE TERMINATION IS NORMAL AND 0 OTHERWISE.
    HOW TO USE:
    >>> from gaussian.get_geometry import get_normaltermination_gaussian
    >>> n=get_normaltermination_gaussian('opt0001.out')
    >>> print (n)
    """

def get_normaltermination_gaussian(filename):
    if os.path.isfile(filename):
        gaufile=open(filename,'r')
        normal=0
        for line in gaufile:
            if "Normal termination" in line: normal=normal+1
        gaufile.close()
        return normal
    else:
        return False
#------------------------------------------------------------------------------------------
""" THIS ROUTINE GETS THE ELECTRONIC ENERGY (HARTREES) FROM GAUSSIAN OUTPUT FILES (*.out)
    IF zpe = 1 THE ZPE CORRECTION IS INCLUDED. 
    HOW TO USE:
    >>> from gaussian.get_geometry import get_energy_gaussian
    >>> ene1=get_energy_gaussian('opt0001.out',0)
    >>> ene2=get_energy_gaussian('opt0001.out',1)
    >>> print (ene1, ene2)
    """

def get_energy_gaussian(filename, correction=1, ccsd='F'):
    enehartree=float(0.0)
    enehartree_correction=float(0.0)
    gaufile=open(filename,'r')
    for line in gaufile:
        if "SCF Done" in line and ccsd=='F':
            scf=line.split()
            enehartree=float(scf[4])
        if "CCSD(T)= " in line and ccsd=='T':
            zerosplit=line.split()
            firstsplit=zerosplit[1].split('D')
            ener=float(firstsplit[0])
            exp=float(firstsplit[1])
            enehartree=float(ener*pow(10,exp))
        if "Zero-point correction=" in line: 
            zpc=line.split()
            if correction == 1: enehartree_correction=float(zpc[2])
        if "Thermal correction to Energy=" in line: 
            zpc=line.split()
            if correction == 2: enehartree_correction=float(zpc[4])
        if "Thermal correction to Enthalpy=" in line: 
            zpc=line.split()
            if correction == 3: enehartree_correction=float(zpc[4])
        if "Thermal correction to Gibbs Free Energy=" in line: 
            zpc=line.split()
            if correction == 4: enehartree_correction=float(zpc[6])
    gaufile.close()
    enehartree=enehartree + enehartree_correction
    enekcalmol=enehartree*hartree2kcalmol
    return enekcalmol
#------------------------------------------------------------------------------------------
def get_ecorrection_gaussian(filename, correction=1):
    enehartree_correction=float(0.0)
    gaufile=open(filename,'r')
    for line in gaufile:
        if "Zero-point correction=" in line:
             zpc=line.split()
             if correction == 1: enehartree_correction=float(zpc[2])
        if "Thermal correction to Energy=" in line:
             zpc=line.split()
             if correction == 2: enehartree_correction=float(zpc[4])
        if "Thermal correction to Enthalpy=" in line:
             zpc=line.split()
             if correction == 3: enehartree_correction=float(zpc[4])
        if "Thermal correction to Gibbs Free Energy=" in line:
             zpc=line.split()
             if correction == 4: enehartree_correction=float(zpc[6])
    gaufile.close()
    enekcalmol=enehartree_correction*hartree2kcalmol
    return enekcalmol
#------------------------------------------------------------------------------------------
""" THIS ROUTINE GETS THE GEOMETRY (X, Y, Z COORDINATES) FROM GAUSSIAN OUTPUT FILES (*.out).
    HOW TO USE:
    >>> from gaussian.get_geometry import get_geometry_gaussian
    >>> geo1=get_geometry_gaussian('path', 'filename', 0, 0)
    >>> print (geo1)
    """

def get_geometry_gaussian(path, filename, zpe=1):
    gaufile=open(path+filename,'r')
    namein=filename.split('.')[0]
    energy=get_energy_gaussian(path+filename,zpe)
    for line in gaufile:
        if line.strip() in ("Input orientation:", "Standard orientation:"):
            moleculeout=Molecule(namein,energy)
            for ii in range(4): line=gaufile.readline()
            line=gaufile.readline()
            while not line.startswith(" --------"):
                ls = line.split()
                if (len(ls) == 6 and ls[0].isdigit() and ls[1].isdigit() and ls[2].isdigit()):
                    ss=get_chemical_symbol(int(ls[1]))
                    xc,yc,zc = float(ls[3]), float(ls[4]), float(ls[5])
                    ai=Atom(ss,xc,yc,zc)
                    moleculeout.add_atom(ai)
                else:
                    break
                line=gaufile.readline()
                ls = line.split()
    gaufile.close()
    return moleculeout
#------------------------------------------------------------------------------------------
""" HOW TO USE:
    >>> from gaussian.calculator_all import get_nt_geometry_gaussian
    >>> from inout.xyz import writexyz
    >>> mol=get_nt_geometry_gaussian('./','opt0001.out',0,0)
    >>> writexyz(mol,"output.xyz")
    """

def get_nt_geometry_gaussian(path, filename, zpe=1):
    r=get_normaltermination_gaussian(path+filename)
    if r is False:
        fopen = open(log_file,'a')
        print("  %s do NOT EXIST!!" %(filename), file=fopen)
        fopen.close()
        return False
    elif (r == 0):
        fopen = open(log_file,'a')
        print("ABNORMAL Termination (0 with zpe=%d) in %s" %(zpe, filename), file=fopen)
        fopen.close()
        return False
    #elif (r==1 and zpe==1):
    #    fopen = open(log_file,'a')
    #    print("ABNORMAL Termination (1 with zpe=%d) in %s" %(zpe, filename), file=fopen)
    #    fopen.close()
    #    return False
    else:
        fopen = open(log_file,'a')
        print("Normal Termination in %s" %(filename), file=fopen)
        fopen.close()
        moleculeout=get_geometry_gaussian(path, filename, zpe)
        moleculeout.c=[r]
        return moleculeout
#------------------------------------------------------------------------------------------
def get_at_geometry_gaussian(path, filename, zpe=1):
    r=get_normaltermination_gaussian(path+filename)
    if r is False: return False
    ##LOGIC FOR FREQUENCY CALCULATIONS HAVE (ZPE=1) r=2 NECESSARILY
    elif (r == 0 and zpe==0) or (r == 1 and zpe==1):
        if get_energy_gaussian(path+filename)==float(0.0): return False
        else:
            moleculeout=get_geometry_gaussian(path,filename, zpe)
            moleculeout.c=[r]
            return moleculeout
    else: return False
#------------------------------------------------------------------------------------------
def get_xt_geometry_gaussian(path, filename, zpe=1):
    r=get_normaltermination_gaussian(path+filename)
    if r is False: return False
    elif r == 0:
        if get_energy_gaussian(path+filename)==float(0.0): return False
        else:
            moleculeout=get_geometry_gaussian(path, filename, zpe)
            moleculeout.c=[r]
            return moleculeout
    else:
        moleculeout=get_geometry_gaussian(path, filename, zpe)
        moleculeout.c=[r]
        return moleculeout
#------------------------------------------------------------------------------------------
def get_all_nt_geometry_gaussian(path, moleculelist, zpe=1):
    moleculeout=[]
    for imol in moleculelist:
        basename=imol.i
        mol01=get_nt_geometry_gaussian(path, basename+'.out', zpe)
        if mol01 is not False:
            moleculeout.extend([mol01])
    return moleculeout
#------------------------------------------------------------------------------------------
def get_all_at_geometry_gaussian(path, moleculelist, zpe=1, erase=0):
    moleculeout=[]
    count=1
    for imol in moleculelist:
        basename=imol.i
        #print(basename)
        mol01=get_at_geometry_gaussian(path, basename+'.out', zpe)
        if mol01 is not False:
            moleculeout.extend([mol01])
            if (erase==1):
                scount=str(count).zfill(4)
                fopen = open(log_file,'a')
                print("%s Remove the files %s.inp and %s.out" %(scount, basename, basename), file=fopen)
                fopen.close()
                count=count+1
                os.system('rm -f '+path+basename+'.inp')
                os.system('rm -f '+path+basename+'.out')
    return moleculeout
#------------------------------------------------------------------------------------------
def get_all_xt_geometry_gaussian(path, moleculelist, zpe=1):
    moleculeout,count,fall=[],0,0
    for imol in moleculelist:
        basename=imol.i
        #print(basename)
        mol01=get_xt_geometry_gaussian(path, basename+'.out', zpe)
        if mol01 is not False:
            moleculeout.extend([mol01])
            count=count+1
        else:
            fall=fall+1
    fopen = open(log_file,'a')
    print("Number of NON-recoverable calculations = %d" %(fall), file=fopen)
    print("Recoverable + SUCCESSFUL calculations  = %d" %(count), file=fopen)
    fopen.close()
    return moleculeout
#------------------------------------------------------------------------------------------
def get_freqneg_gaussian(filename):
    gaufile=open(filename,'r')
    freq_negative=0
    freq_neg_list=[]
    for line in gaufile:
        if "Frequencies" in line:
            freq=line.split()
            freq.pop(0)
            freq.pop(0)
            for ifreq in range(len(freq)):
                if float(freq[ifreq]) < 0.0:
                    freq_negative=freq_negative+1
                    freq_neg_list.append(float(freq[ifreq]))
    gaufile.close()
    freq_neg_list.sort()
    freq_sample= freq_neg_list[0] if (freq_negative > 0) else 0.0
    return freq_negative, freq_sample
#------------------------------------------------------------------------------------------
def cutter_negfreq_gaussian(moleculein,path):
    moleculeout, count=[], 0
    for imol in moleculein:
        basename=imol.i
        freq,fsam=get_freqneg_gaussian(path+basename+'.out')
        if ( freq == 0 ):
            moleculeout.extend([imol])
        else:
            count=count+1
            fopen = open(log_file,'a')
            print("%s ... DISCRIMINATED by Negative Frequency (%f)" %(basename, fsam), file=fopen)
            fopen.close()
    if count==0:
        fopen = open(log_file,'a')
        print("No elements (or new elements) discriminated by Negative Frequency", file=fopen)
        fopen.close()
    return  moleculeout
#------------------------------------------------------------------------------------------
def cutter_normterm_gaussian(moleculein, path, nt=1):
    if nt==0: return moleculein
    moleculeoutnt1, moleculeoutnt2, count=[], [], 0
    for imol in moleculein:
        basename=imol.i
        norm=get_normaltermination_gaussian(path+basename+'.out')
        if ( norm == 0 ):
            count=count+1
            if nt==1 or nt==2:
                fopen = open(log_file,'a')
                print("%3d) %s ... DISCRIMINATED by Normal Termination (0)" %(count,basename), file=fopen)
                fopen.close()
        if ( norm >= 1 ):
            if nt==2 and norm==1:
                count=count+1
                fopen = open(log_file,'a')
                print("%3d) %s ... DISCRIMINATED by Normal Termination (1)" %(count,basename), file=fopen)
                fopen.close()
            moleculeoutnt1.extend([imol])
        if ( norm == 2 ):
            moleculeoutnt2.extend([imol])
    if nt==1: return moleculeoutnt1
    if nt==2: return moleculeoutnt2
#------------------------------------------------------------------------------------------
def get_full_point_group(filename):
    gaufile=open(filename,'r')
    ans='X'
    for line in gaufile:
        if "Full point group" in line:
            fpg=line.split()
            ans=str(fpg[3])
    gaufile.close()
    return ans
#------------------------------------------------------------------------------------------
def display_info_gaussian(moleculein, path, logfile):
    fopen = open(log_file,'a')
    print("%s %18s %18s %9s %7s %6s" % ("#","File"," Delta E ", "Norm","Point","Img "), file=fopen)
    print("%s %18s %18s %9s %7s %6s" % ("N","Name","(kcal/mol)","Term","Group","Freq"), file=fopen)
    fopen.close()
    count=1
    ene0=moleculein[0].e
    for imol in moleculein:
        istr=str(count).zfill(4)
        filename=imol.i+'.out'
        deltaene=imol.e -  ene0
        norm=get_normaltermination_gaussian(path+filename)
        freq,fsam=get_freqneg_gaussian(path+filename)
        fpg=get_full_point_group(path+filename)
        if freq==0:
            fopen = open(log_file,'a')
            print("%s %20s %16.9f %5d %6s %6d        " % (istr,filename,deltaene,norm,fpg,freq), file=fopen)
            fopen.close()
        else:
            fopen = open(log_file,'a')
            print("%s %20s %16.9f %5d %6s %6d (%6.1f)" % (istr,filename,deltaene,norm,fpg,freq,fsam), file=fopen)
            fopen.close()
        count=count+1
#------------------------------------------------------------------------------------------
def get_traj(filename):
    start, end, ene = [], [], []
    openold = open(filename,"r")
    rline = openold.readlines()
    for i in range(len(rline)):
        if "Standard orientation:" in rline[i]:
            start.append(i)
            for m in range (start[-1] + 5, len(rline)):
                if "---" in rline[m]:
                    end.append(m)
                    break
        if "SCF Done" in rline[i]:
            eneline = rline[i].split()
            ene.append(eneline[4])
    moleculeout=[]
    for i,iStart in enumerate(start):
        ei=str(ene[i-1])
        enekcalmol=float(ei)*hartree2kcalmol
        name='traj'+str(i).zfill(4)
        molx=Molecule(name, enekcalmol)
        for line in rline[start[i]+5 : end[i]] :
            words = line.split()
            ss = get_chemical_symbol(int(words[1]))
            xc,yc,zc = float(words[3]), float(words[4]), float(words[5])
            ai=Atom(ss,xc,yc,zc)
            molx.add_atom(ai)
        moleculeout.extend([molx])
    openold.close()
    return moleculeout
#------------------------------------------------------------------------------------------
def run_sample():
    list=['stage100001', 'stage100002', 'stage100003']
    mol1=get_all_nt_geometry_gaussian('./', list, 0)
    display_info_gaussian(mol1, './')
#------------------------------------------------------------------------------------------
