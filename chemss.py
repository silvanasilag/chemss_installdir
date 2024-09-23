#!/usr/local/anaconda3/bin/python
# -*- coding: utf-8 -*-
#Silvana Silva Aguirre, CINVESTAV Unidad Mérida

##modules
import sys
import os
import glob
import time
import shutil

from nmrutils.getbilparam     import get_a_str, get_a_int, get_a_float
from runrun.qpbs import send_pbs_files_to_queue, printf
from nmrutils.terminacion     import terminacion
from nmrutils.isotropicreader import molecules_data
from nmrutils.getbilparam     import read_block_of_inp
from runrun.chemssinput_writer import chemss_inp,chemssinp_edit

chms_path= "/Users/silvana/installdir/chemss_installdir/dataset" #xyz and dnmr database path
work_dir=str(os.getcwdb()) #working path
work_dir=work_dir.replace("b'","'")
work_dir=work_dir.strip("'")

dnmr=str(chms_path +"/dnmr") 
xyz=str(chms_path+"/xyz")
silvana=93
start_time = time.time()

def exist(dirct,n): # haces estas banderas más a pueba de tontos
    messages = {
        True: "\t Existing directory",
        False: "\t Not existing directory",
        "start": "CHEMSS 1.0\nPlease fill in the INTPUT.txt file",
        "cyc": "Please fill in the sequence.txt file"
    }
    m=0
    if n != False and n!= True:
        m=n
        n = False
    if os.path.exists(dirct) == n:
        if m=="start":
            print(messages[m])
            chemss_inp()
            if os.path.exists("nw_ds")==False:
                os.mkdir("nw_ds")
        elif m=="cyc":
            f = open("sequence.txt", "a")
            f.write("dir_name/ dir_name_new/ OPT Funcional and basis / NMR funcional and basis")
            f.close()
            print(messages[m])
        else:
            print(messages[n],dirct.split("/")[-1])
        sys.exit(1)

def chk_files(key_opt,path,nn):
    outdata=[]
    key_opt=key_opt.replace("/","_")
    os.chdir(path)
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

def CHEMSS():
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
#---------------------------------------------------
    if calc =='YES':
        exist(path,True)
        os.mkdir(path)
        inp_pbs_writer(path,xyz,chms_path,nmr,opt)
        os.chdir(path)
        # send_pbs_files_to_queue(njobs, time_sleep)
        #----------------------------------------------
        print('Calculation complete')
    if new_molg == "YES":
        exist(path2,True)
        os.mkdir(path2)
        chk=work_dir+"/nw_ds"
        inp_pbs_writer(path2,chk,chk,nmr,opt) #*******
        os.chdir(path2)
        send_pbs_files_to_queue(njobs, time_sleep)
        #----------------------------------------------
    if stat == 'YES' and new_mol =="NO":
        exist(path,False)
        #falla,error,keys = terminacion(path)
        data,keys = molecules_data(path,dnmr)
        key_opt=out_w(path,data,[],tbl_comp,keys)
        chk_files(key_opt,path,1)
    # ·····························································
    if new_mol =="YES":
        # exist(path,False)
        exist(path2,False)
        # falla,error,opt,nmr = terminacion(path)
        #falla,error,keys= terminacion(path2)
        keys=["-","-"]
        data,keys= molecules_data(path,dnmr)
        data2,keys = molecules_data(path2,str(work_dir+"/nw_ds"))
        # key_opt_n=out_w(path2,data2,[],tbl_comp,keys)   #scale factor perform with just the new data
        key_opt_n=out_w(path2,data,data2,tbl_comp,keys) #scale factor perform with just the base dataset
        # key_opt_n=out_w(path2,data2,data,tbl_comp) #scale factor aplaided to base dataset performed with just the new data 
        # data3=data+data2
        # key_opt_n=out_w(path2,data3,data,tbl_comp) #scale factor perform with the base dataset + the new data
        if stat == 'YES':  out_w(path,data,[],tbl_comp,keys)
        """
        if key_opt != key_opt_n : 
            print("Not the same level of theory")
            print(path+"    :   "+key_opt)
            print(path2+"   :  "+key_opt_n)
            sys.exit("")
        """
        #chk_files(key_opt,path,1)
        chk_files(key_opt_n,path2,0)
    # ·····························································
    if  stat == 'NO' and  calc == 'NO'and new_mol == 'NO'and new_molg =='NO':
        print("Nothing required")
    return r

    # os.system("rm "+str(path)+"/*.chk")  #scrach de gaussian
    # os.system("rm ../*.ER")               #scrach de gaussian            
    # os.system("rm ../*.OU")               #scrach de gaussian


if __name__ == '__main__':
    exist("INPUTNMR.txt","start")
    from runrun.inpwriter import inp_pbs_writer
    from nmrutils.out_writer import out_w
    scuenc = get_a_str('sequential','NO')
    if scuenc=="YES":
        exist("sequence.txt","cyc")
        fil=[]
        with open("sequence.txt", 'r') as f:
            for line in f:
                line = line.strip()
                lin=line.split("/")
                if len(lin)==3: fil.append(lin)
                elif len(lin)!=0: fil.append(lin)
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
            r=CHEMSS()
    else:
        r=CHEMSS()

# if r==3 or r==4:
    # e = int(time.time() - start_time)
    # file1 = open(str(path)+"/sclsf_" + str(path)+".txt", "a" )
    # file1.write('{:02d}:{:02d}:{:02d}'.format(e // 3600, (e % 3600 // 60), e % 60))
    # file1.close()


