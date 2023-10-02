# -*- coding: utf-8 -*-

#Silvana Silva Aguirre, CINVESTAV Unidad Mérida,2019 
#El programa genera archivos inp con base a las cordenadas del archivo xyz, y archivos pbs
#PARAMETROS:Funcional OPT,  funcional RMN, carpeta, datos para los pbs

##modulos

import os
from os import path 
import glob
import shutil
from runrun.make_pbs import make_apbs_gaussian
from nmrutils.getbilparam     import  get_a_str, read_block_of_inp, get_a_int, get_a_float,key_norm
from nmrutils.get_geometry    import get_geo, write_xyz
#-------------------------------------------------------------------------data extraction
# nmr = read_block_of_inp('gaussian nmr')
# opt = read_block_of_inp('gaussian opt')
cola = get_a_str('queue','qintel')
nproc = get_a_int('nodes',4)
time_sleep = get_a_float('timesleep',1.0)
ram      = get_a_int('memory_in_gb',16)
walltime = get_a_str('walltime','08:00:00')

#-------------------------------------------------------------------------

def from_reader(xyz): #function that extracts the coordinates and name of xyz files
    nombre  = [ ] 
    atomos  = [ ]
    xx = [ ]
    yy = [ ]
    zz = [ ]
    with open(xyz, 'r') as f2:
        for line in f2:
            line = line.strip()
            lin = line.split()
            if len(lin) > 0:
                if len(lin) != 4:
                    nombre.append(str(lin[0]))
                if len(lin)==4 :
                    atomos.append(str(lin[0]))
                    xx.append(float(lin[1]))
                    yy.append(float(lin[2]))
                    zz.append(float(lin[3]))                                                  
    from_dict = { }
    from_dict = {'xx':xx,'yy':yy,'zz':zz,'atomos':atomos,'nombre':nombre }
    return from_dict

#-------------------------------------------------------------------------MAIN
def inp_pbs_writer(path,xyzpath,chkp,nmr,opt):
    print(nmr,opt)
    filesout = [ ]
    optk=key_norm(opt[0],"OPT")
    optk_l=optk.split("_")
    if len(opt)==1 and os.path.exists(chkp+"/"+str(optk)) == True : 
        f=2
        chks= chkp+"/"+str(optk)
        print("1")
    elif len(opt)==1 and len(optk_l)==3 and os.path.exists(chkp+"/"+str(optk_l[0])+"_"+str(optk_l[1])) == True :
        f=0
        chks= chkp+"/"+str(optk_l[0])+"_"+str(optk_l[1])
        print("2")
    elif os.path.exists(chkp+"/"+"B3LYP_6-31G(d,p)")==True:
        f=0
        chks= chkp+"/"+"B3LYP_6-31G(d,p)"
        print("3")
    else: 
        f=1
    os.chdir(xyzpath)
    for iifile in glob.glob("*.xyz"):  
        filesout.append(xyzpath + "/" + iifile)
    npro = str("%NProc=") + str(nproc)
    mem = str("%MEM=")+str(ram)+str("GB") #----gaussian memory  *****
    chk = str("%Chk=")
    for ii in filesout:           
        base = os.path.basename(ii)
        name = os.path.splitext(base)[0]
        reader = from_reader(ii) #---gets xyz
        nombre = reader['nombre']
        Z = reader['atomos']
        xx= reader['xx']
        yy= reader['yy']
        zz= reader['zz']
        #--------------------------------------------------------------------------------------------
        os.chdir(path)     
        # make_apbs_gaussian(nproc, 16, cola, walltime, name)
        out = str(name) + ".inp"       
        with open(out,'w') as out:  #--------------writes inp in the working file
            out.write('%s \n'%npro)
            out.write('%s \n'%mem) 
            out.write('%s' %chk)
            out.write('%s \n' %name)
            if f==0:  
                for ix in range(len(opt)):
                    if ix == 0:out.write('%s  Geom=Check Guess=Read \n' % opt[ix])
                    else: out.write('%s \n' % opt[ix])            
                out.write("       \n")
                out.write(str(nombre[0]) + " \n")            
                out.write("       \n")
                out.write("0      1\n")
                out.write("       \n")
                out.write('--Link1-- \n')
                out.write('%s \n'%npro)
                out.write('%s' %chk)
                out.write('%s \n' %name)
                for ix in range(len(nmr)):
                    if ix == 0:out.write('%s  Geom=Check Guess=Read \n' % nmr[ix])
                    else: out.write('%s \n' % nmr[ix])
                out.write('       \n')
                out.write(' TMS NMR        \n')
                out.write('       \n')
                out.write('0 1\n')
                out.write("       \n")
                original = str(chks+"/"+name+".chk")
                target =str(path+"/"+name+".chk")
                shutil.copyfile(original, target)
            if f==2:  
                for ix in range(len(nmr)):
                    if ix == 0:out.write('%s  Geom=Check Guess=Read \n' % nmr[ix])
                    else: out.write('%s \n' % nmr[ix])
                out.write('       \n')
                out.write(' TMS NMR        \n')
                out.write('       \n')
                out.write('0 1\n')
                out.write("       \n")
                original = str(chks+"/"+name+".chk")
                target =str(path+"/"+name+".chk")
                shutil.copyfile(original, target)
       
            if f==1:
                 for ix in range(len(opt)):
                     if ix == 0:out.write('%s  \n' % opt[ix])
                     else: out.write('%s \n' % opt[ix])
                 out.write("       \n")
                 out.write(str(nombre[0]) + " \n")
                 out.write("       \n")
                 out.write("0      1\n")
                 for ix in range(len(Z)):
                     out.write( '%s  %16f  %11f %11f \n' % ( Z[ix],xx[ix], yy[ix], zz[ix]))                              
                 out.write("       \n")
                 out.write('--Link1-- \n')
                 out.write('%s \n'%npro)
                 out.write('%s' %chk)
                 out.write('%s \n' %name)
                 for ix in range(len(nmr)):
                     if ix == 0:out.write('%s  Geom=Check Guess=Read \n' % nmr[ix])
                     else: out.write('%s \n' % nmr[ix])
                 out.write('       \n')
                 out.write(' TMS NMR        \n')
                 out.write('       \n')
                 out.write('0 1\n')
                 out.write("       \n")
        os.chdir(xyzpath)
