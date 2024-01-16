# -*- coding: utf-8 -*-

import sys
import os
import glob
class Level_ot: #Level of theory
    def __init__(self, opt_fun,opt_base, nmr_fun, nmr_base, dispersion, solv_model, solvent):
        self.of=opt_fun
        self.ob=opt_base
        self.nf=nmr_fun
        self.nb= nmr_base
        self.d=dispersion
        self.sm=solv_model
        self.sl=solvent

#--------------------------DATA READER-----------------------------------------
def from_reader(inp,key_opt,key_nmr):
    ter =  "not normal"
    linea ='normal'
    fl=0
    with open(inp, 'r') as f:
        print("---------",inp)
        for line in f:
            line = line.strip()
            if len(line) != 0:
                if str("#")  and str("OPT") in line and key_opt=="-":
                    key_opt =line
                    key_opt=key_opt.strip(' Guess=Read')
                    key_opt=key_opt.strip(' Geom=Check')
                    fl=1
                if str("#") and str("NMR")in line and key_nmr=="-":
                    key_nmr=line
                    key_nmr=key_nmr.strip(' Guess=Read')
                    key_nmr=key_nmr.strip(' Geom=Check')
                    fl=1
                if str("#") and str("IOP")in line.upper() and fl==2:
                    print(line)
                    if "(3/76=1000007400,3/77=0999900001,3/78=0000109999)" in line: w="WC04"
                    elif "(3/76=1000001189,3/77=0961409999,3/78=0000109999)" in line: w="WP04"
                    else: w="Iop"
                    print(inp)
                    print("line 37:",key_nmr,key_opt)
                    if key_nmr == "-": key_opt = key_opt=key_opt.replace("BLYP",w) 
                    if key_opt != "-" and key_nmr != "-": key_nmr =key_nmr.replace("BLYP",w) 
                if "Magnetic shielding" in line:
                    ter = "Normal"  
                if "Error" in line: 
                    linea = line
                fl=fl+1
    keys.of=key_opt
    keys.nf=key_nmr
    return ter,linea,keys

def key_compare(key_opt,key_nmr,nmr,opt):
    opt=opt.strip("=CARTESIAN")
    f=[0,0,0]
    if nmr.upper() not in key_nmr: 
        f[1]=1
        key_nmr.append(nmr.upper())
    if opt.upper() != "-" and opt.upper() not in key_opt: 
        f[0]=1
        key_opt.append(opt.upper())
    if opt =="-": f[2]=1
    return key_opt,key_nmr,f

#----------------------------------------------------------------------
def terminacion(out):
    path = out
    os.chdir(path)
    error = [ ]
    filesinp = [ ]
    key_opt=[]
    key_nmr=[]
    key_opt_f=[]
    key_nmr_f=[]
    for iifile in glob.glob("*.out"):  
        filesinp.append(iifile)
        filesinp.sort()
    nofalla = [ ]
    falla =[ ]
    print(out)
    for x in range(len(filesinp)):
        innp = filesinp[x]
        ter,linea,keys = from_reader(innp,"-","-")
        nmr=keys[1]
        opt=keys[0]
        key_opt,key_nmr,f=key_compare(key_opt,key_nmr,nmr,opt)
        if  ter == "not normal":
            falla.append(innp)
            error.append(str(linea))
        if  ter == "normal":
            nofalla.append(innp)
        if f[0]==1 and x != 0 :print("Not the same OPT key words",filesinp[x],x)
        if f[1]==1 and x != 0:print("Not the same NMR key words",filesinp[x],x)
        if f[2]==1: print("Not OPT key check it out",filesinp[x],x)
        key_opt_f.append((opt,filesinp[x]))
        key_nmr_f.append((nmr,filesinp[x]))
    if len(key_opt)!=1 :
        print(key_opt) 
        for i in key_opt_f:print(i)
    if len(key_nmr)!=1 :
        print(key_nmr) 
        for ii in key_nmr_f:print(ii)

    os.chdir("..")

    if len(falla) != 0:
        out = str(path)+ "mmmmistake" + ".txt"       #Name of the files with error
        with open(out,'a') as out:
            out.write("OH NO!\n\n")
            out.write("Calculation data\n")
            out.write( "%s \n %s \n"%(opt,nmr))
            out.write("Termination with error:%i \n" %(len(error)))
            for i in range(len(falla)):
                out.write("|%s || %-20s|\n " %(falla[i],error[i]))
        print("Termination with erro:%i \n" %(len(error)))
        for i in range(len(falla)):
            print("|%s || %-20s|\n " %(falla[i],error[i]))
        sys.exit(1)

    keys = Level_ot("na", "na", "", "", "", "", "")
    return falla,error,keys
    
