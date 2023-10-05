#!/usr/local/anaconda3/bin/python
# -*- coding: utf-8 -*-

#Silvana Silva Aguirre, CINVESTAV Unidad Mérida,2019

##modulos
import sys
import os
from os import path 
import glob
import time
import shutil

from nmrutils.getbilparam     import get_a_str, get_a_int, get_a_float
# from nmrutils.regression      import regression

from runrun.qpbs              import send_pbs_files_to_queue
from nmrutils.terminacion     import terminacion
from nmrutils.isotropicreader import molecules_data
from nmrutils.get_geometry    import get_geo, write_xyz
from nmrutils.getbilparam     import read_block_of_inp
from runrun.chemssinput_writer import chemss_inp,chemssinp_edit

chms_path= "/Users/silvana/installdir/chemss_installdir/dataset" #path de la base de datos del xyz el dnmr
work_dir=str(os.getcwdb()) #path de la carpeta desde donde se ejecuta
work_dir=work_dir.replace("b'","'")
work_dir=work_dir.strip("'")

dnmr=str(chms_path +"/dnmr") 
xyz=str(chms_path+"/xyz")
start_time = time.time()



def exist(dirct,n): # hacesr estas banderas más a pueba de tontos
    m=0
    if n == True : warning  =str( "\t Existing file ")
    if n == False : warning  =str( "\t Not existing file") # mensaje de inicio 
    if n== "start" : 
        m=1
        n=False
    if n == "cyc":
        m=2 
        n=False
    if os.path.exists(dirct) == n:
        if m==1:
            print("CHEMSS 1.0")
            print("Please fill in the INTPUT.txt file ")
            chemss_inp()
            if os.path.exists("nw_ds")==False:
                os.mkdir("nw_ds")
        elif m==2:
            f = open("sequence.txt", "a")
            f.write("dir_name/ dir_name_new/ OPT Funcional and basis / NMR funcional and basis")
            f.close()
            print("Please fill in the sequence.txt file ")
        else:
            print(warning)
        sys.exit(1)

def chk_files(key_opt,path,nn):
    outdata=[]
    key_opt=key_opt.replace("/","_")
    os.chdir(path)
    print("1")
    if len(glob.glob("*.chk")) != 0:
        print("2")
        if nn==1: chk_p = chms_path
        else : chk_p = str(work_dir+"/nw_ds")
        os.chdir(chk_p)
        print(os.getcwd())
        if os.path.exists(key_opt) == False:
            os.mkdir(key_opt)
            os.chdir(path) 
            print(os.getcwd())
            for ifile in glob.glob("*.chk"): outdata.append(ifile)
            outdata.sort()
            for iifile in outdata:
                original = str(path+"/"+iifile)
                target =str(chk_p+"/"+key_opt+"/"+iifile)
                shutil.copyfile(original, target)
            os.system("rm *.chk")


def CHESMS():
    r=0
    fname = get_a_str('foldername_base','test')
    fname2 = get_a_str('foldername_new','test2')
    njobs = get_a_int('njobs',5)
    stat = get_a_str('statistics','YES')
    calc = get_a_str('calculus','YES')
    time_sleep = get_a_float('timesleep',1.0)
    new_mol = get_a_str('new_mol_statistics','NO')
    new_molg = get_a_str('new_mol_calculus','NO')
    tbl_comp=get_a_str('comparative','-') #bandera oculta
    nmr = read_block_of_inp('gaussian nmr')
    opt = read_block_of_inp('gaussian opt')
    path = str(work_dir) + "/"+fname
    path2 = str(work_dir) + "/"+fname2
    print(new_mol,fname2)
#---------------------------------------------------
    if calc =='YES':
        r=1
        exist(path,True)
        os.mkdir(path)
        chk=chms_path
        inp_pbs_writer(path,xyz,chk,nmr,opt) 
       	os.chdir(path)
        # send_pbs_files_to_queue(njobs, time_sleep)
        #----------------------------------------------
        print('Calculation complete')
    if new_molg == "YES":
        r=2
        exist(path2,True)
        os.mkdir(path2)
        chk=work_dir+"/nw_ds"
        inp_pbs_writer(path2,chk,chk,nmr,opt) #*******
        os.chdir(path2)
        send_pbs_files_to_queue(njobs, time_sleep)
        #----------------------------------------------
    if stat == 'YES' and new_mol =="NO":
        r=3
        exist(path,False)
        falla,error,opt,nmr = terminacion(path)
        data = molecules_data(path,dnmr)
        data2  = []
        key_opt=out_w(path,data,data2,tbl_comp)
        chk_files(key_opt,path,1)
    # ·····························································
    if new_mol =="YES":
        r=4
        exist(path,False)
        exist(path2,False)
        falla,error,opt,nmr = terminacion(path)
        falla,error,opt,nmr = terminacion(path2)
        data = molecules_data(path,dnmr)
        data2 = molecules_data(path2,str(work_dir+"/nw_ds"))
        key_opt=out_w(path,data,[],"-")
        key_opt_n=out_w(path2,data,data2,tbl_comp)
        """
        if key_opt != key_opt_n : 
            print("Not the same level of theory")
            print(path+"    :   "+key_opt)
            print(path2+"   :  "+key_opt_n)
            sys.exit("")
        """
        chk_files(key_opt,path,1)
        chk_files(key_opt_n,path2,0)
        
    # ·····························································
    if  stat == 'NO' and  calc == 'NO'and new_mol == 'NO'and new_molg =='NO':
        print("Nothing required")

    return r

    # os.system("rm "+str(path)+"/*.chk")  #scrach de gaussian
    # os.system("rm ../*.ER")               #scrach de gaussian            
    # os.system("rm ../*.OU")               #scrach de gaussian

exist("INPUTNMR.txt","start")

from runrun.inpwriter         import inp_pbs_writer
from nmrutils.out_writer      import out_w

scuenc = get_a_str('sequential','NO')
if scuenc=="YES":
    exist("sequence.txt","cyc")
    fil=[]
    with open("sequence.txt", 'r') as f:
        for line in f:
            line = line.strip()
            lin=line.split("/")
            if len(lin)==3:
                fil.append(lin) 
            elif len(lin)!=0:
                fil.append(lin) 
        if len(fil)==1 and "all" in fil[0][0] :
            print("filtro1")
            fil_temp=[]
            dirs = os.listdir(work_dir)
            dir_list = [i for i in dirs if os.path.isdir(os.path.join(work_dir, i))]
            dirr_list=dir_list.sort()
            print(dirr_list)
            for i in dir_list:
                if fil[0][0] == "all_base":
                    ld=[i,'-','-','-']
                    fil_temp.append(ld)
                if fil[0][0] == "all_new":
                    ld=['-',i,'-','-']
                    fil_temp.append(ld)
            fil=fil_temp
    for idata in fil:
        os.chdir(work_dir)
        chemssinp_edit(idata)
        time.sleep(3)
        r=CHESMS()
else:
    r=CHESMS()

# if r==3 or r==4:
    # e = int(time.time() - start_time)
    # file1 = open(str(path)+"/sclsf_" + str(path)+".txt", "a" )
    # file1.write('{:02d}:{:02d}:{:02d}'.format(e // 3600, (e % 3600 // 60), e % 60))
    # file1.close()


