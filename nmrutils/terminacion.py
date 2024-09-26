# -*- coding: utf-8 -*-

import sys
import os
import glob

#--------------------------DATA READER-----------------------------------------   
def from_reader(inp,key_opt,key_nmr):
    ter =  "not normal"
    linea ='normal'
    with open(inp, 'r') as f:
        for line in f:
            line = line.strip()
            if len(line) != 0:
                if str("#")  and str("OPT") in line and key_opt=="-":
                    key_opt =line
                    key_opt=key_opt.strip(' Guess=Read')
                    key_opt=key_opt.strip(' Geom=Check')
                if str("#") and str("NMR")in line and key_nmr=="-":
                    key_nmr=line
                    key_nmr=key_nmr.strip(' Guess=Read')
                    key_nmr=key_nmr.strip(' Geom=Check')
                if str("#") and str("IOP")in line.upper() : 
                    if "(3/76=1000007400,3/77=0999900001,3/78=0000109999)" in line: w="WC04"
                    elif "(3/76=1000001189,3/77=0961409999,3/78=0000109999)" in line: w="WP04"
                    else: w="Iop"
                    if key_nmr == "-": key_opt = key_opt=key_opt.replace("BLYP",w) 
                    if key_opt != "-" and key_nmr != "-": key_nmr =key_nmr.replace("BLYP",w) 
                if "Magnetic shielding" in line:
                    ter = "Normal"  
                if "Error" in line: 
                    linea = line
    keys=[key_opt,key_nmr]
    return ter,linea,keys

def key_compare(key_opt,key_nmr,nmr,opt,glist):   
    glist=glist.split("/")[-1]
    opt=opt.strip("=CARTESIAN")
    if nmr.upper() not in key_nmr : 
        key_nmr.append(nmr.upper())
        if len(key_nmr)>1: print("--*--Not the same NMR key words:",glist)
    if opt.upper() != "-" and opt.upper() not in key_opt : 
        key_opt.append(opt.upper())
        if len(key_opt)>1: print("--*--Not the same OPT key words:",glist)
    if opt =="-": print("Not OPT key check it out",glist)
    

    return key_opt,key_nmr

def key_compare_all(key,key_f):
    print("--*--The diferent key words--*--")
    for ii in key:print(ii)
    print("--*--Key words per file--*--")
    for i in key_f:print(i)
#----------------------------------------------------------------------
def terminacion(f_err,glist):
    error = []
    nofalla = []
    falla =[]
    key_opt,key_nmr=[],[]
    key_opt_f,key_nmr_f=[],[]
    run_index=glist.index(f_err)
    for x in range(len(glist))[run_index:]:
        innp = glist[x]
        ter,linea,keys = from_reader(innp,"-","-")
        nmr=keys[1]
        opt=keys[0]
        if  ter == "not normal":
            falla.append(innp)
            error.append(str(linea))
        if  ter == "normal":
            nofalla.append(innp)
        key_opt,key_nmr=key_compare(key_opt,key_nmr,nmr,opt,glist[x])
        key_opt_f.append(opt,glist[x])
        key_nmr_f.append(nmr,glist[x])
    if len(key_opt)!=1 : key_compare_all(key_opt,key_opt_f)
    if len(key_nmr)!=1 :key_compare_all(key_nmr,key_nmr_f)

    path = os.path.dirname(f_err)  # Obtener el directorio de la ruta
    os.chdir(path)
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

    
