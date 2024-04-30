# -*- coding: utf-8 -*-
#codigo con el calculo del python

##modulos
import os
import numpy as np
import matplotlib.pyplot as plt
import datetime
from datetime import datetime as dt
import random
import glob
import pandas as pd
import re

from scipy import stats
from nmrutils.terminacion import from_reader
from nmrutils.regression import scale,stat,splot,xy
from sklearn.metrics import mean_squared_error
from nmrutils.isotropicreader import molecules_data
from nmrutils.getbilparam import get_a_str, read_block_of_inp, get_a_int, get_a_float,key_norm
from nmrutils.table_ngbh import tabla_alldataset


def fechayhora():
    months=["Jan","Feb","March","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
    now = dt.now()
    mon=months[now.month-1]
    noww=str(now.day)+"-"+mon+"-"+str(now.year)+"   "+str(now.hour)+":"+str(now.minute)+":"+str(now.second)
    
    return noww
#--------------------------------------------------- 
class DataNMR:
    def __init__(self, atom, slope, intercept, r_squere, rmsd):
        self.s=atom
        self.m=slope
        self.b=intercept
        self.r2=r_squere
        self.rmsd=rmsd

nmr = read_block_of_inp('gaussian nmr') #
opt = read_block_of_inp('gaussian opt') #
cola = get_a_str('queue','qintel') #
nproc = get_a_int('nodes',4) #

def kys_changes(lev):
    if "PBE1PBE" in lev: lev=lev.replace("PBE1PBE","PBE0") 
    if "TPSSTPSS" in lev: lev=lev.replace("TPSSTPSS","TPSS")
    if "GD3BJ" in lev: lev=lev.replace("_GD3BJ","")
    if "GD3" in lev: lev=lev.replace("_GD3","")
    if "GD2" in lev: lev=lev.replace("_GD2","")
    if "Def2SVP" in lev: lev=lev.replace("Def2SVP","def2-SVP") 
    lev=lev.replace("_","/")
    return lev


#--------------------------------------------------- 

def comp_table(tbl_comp,data_c,data_h,keys,path):
    os.chdir("..")
    # print(os.getcwd())
    kop=str(keys[0]).split("_")
    knmr=re.sub(r'_', '/', (keys[1]), 1)
    knmr=str(knmr).split("_")
    if len(kop)==3 :
        gd=str(kop[2])
        if "GD3BJ" in gd: gd=gd.replace("GD3BJ","D3BJ")
        if "GD3" in gd: gd=gd.replace("GD3","D3")
        if "GD2" in gd: gd=gd.replace("GD2","D2")
    else: 
        gd="No dispersion added"
    if len(knmr)>1:
        nmr=knmr[0]
        solvent=knmr[1]
        if len(knmr)>2:
            solv_method=knmr[2]
        else:
            solv_method="PCM"
    else:
        nmr=kys_changes(str(keys[1]))
        solvent="Gas phase"
        solv_method=""
    opt=kys_changes(str(keys[0]))
    tbl_compf=str(tbl_comp+".csv")
    df1 = pd.DataFrame({" ":[path],
                        "Geometric OPT": opt,
                        "NMR": nmr,
                        "Dispersion":gd,
                        "Solvent": solvent,
                        "Solvent model": solv_method,
                        "slope H1": [round(data_h.m,4)],
                        "intercept H1":[round(data_h.b,4)],
                        "RMSD H1" :[round(data_h.rmsd,4)],
                        "R^2 H1" :[round(data_h.r2,4)],
                        "slope C13": [round(data_c.m,4)],
                        "intercept C13":[round(data_c.b,4)],
                        "RMSD C13" :[round(data_c.rmsd,4)],
                        "R^2 C13" :[round(data_c.r2,4)],
                        })
    if os.path.isfile(tbl_compf)== True:
        df = pd.read_csv (tbl_compf)
        dff=pd.concat([df,df1])
        dff=dff.drop_duplicates()
        dff.to_csv(tbl_compf,index=False)
    if os.path.isfile(tbl_compf)== False: 
        df1.to_csv(tbl_compf,index=False)

#---------------------------------------------------    
def out_w(path,data,new,tbl_comp,keys):
    cputime=datetime.timedelta(days=int(0),hours=int(0), minutes=int(0), seconds=round(float(0),2))
    elepsetime=datetime.timedelta(days=int(0),hours=int(0), minutes=int(0), seconds=round(float(0),2))
    cputime_nmr=datetime.timedelta(days=int(0),hours=int(0), minutes=int(0), seconds=round(float(0),2))
    elepsetime_nmr=datetime.timedelta(days=int(0),hours=int(0), minutes=int(0), seconds=round(float(0),2))
    xhh,xcc,yhh,ycc,cputime,elepsetime,cputime_nmr,elepsetime_nmr = xy(data,cputime,elepsetime,cputime_nmr,elepsetime_nmr)
    #----------------------------------------------------
    data_h =DataNMR("H", 0, 0, 0, 0)
    data_c =DataNMR("C", 0, 0, 0, 0) 
    print ("\n----Data for H1-----")
    stat(xhh,yhh,data_h)
    print ("----Data for C13----")
    stat(xcc,ycc,data_c)
    scale(data,data_h,data_c)
    #---------------------------------------------------------
    xn_h,xn_c,yn_h,yn_c=[],[],[],[]
    if len(new) != 0:
        xn_h,xn_c,yn_h,yn_c,cputime,elepsetime,cputime_nmr,elepsetime_nmr =xy(new,cputime,elepsetime,cputime_nmr,elepsetime_nmr)
        scale(new, data_h, data_c)
    #---------------------------------------------------Escribe los  archivos out finales     
    os.chdir(path)
    out=(glob.glob("*.out"))
    key_opt=keys[0]
    key_nmr=keys[1]
    keys=[key_norm(key_opt,"OPT"),key_norm(key_nmr,"NMR")]
    if key_opt == "-": key_opt=opt
    plot_name = str(path) +"/scaled_H1" +  ".png"     
    plotc_name =  str(path)+"/scaled_C13" +".png" 
    splot(xhh,yhh,plot_name,data_h.m,data_h.b,data_h.r2,"Hidrogeno-1","steelblue")
    if len(xn_h) != 0:splot(xn_h,yn_h,plot_name,data_h.m,data_h.b,data_h.r2,"Hidrogeno-1","salmon")
    plt.clf()
    splot(xcc,ycc,plotc_name,data_c.m,data_c.b,data_c.r2,"Carbono-13","steelblue")
    if len(xn_h) != 0:splot(xn_c,yn_c,plotc_name,data_c.m,data_c.b,data_c.r2,"Carbono-13","salmon")
    #---------------------------------------------------    
    # writer table for xy all data. 
    #---------------------------------------------------
    out =str(path)+"/scaled_hormons.txt"
    out=open(out,'a')
    out.write("CHEMical Shift Scaler\n\n")
    out.write("This software is provided by TheoChem Merida. \n")
    out.write("Code by Silvana Silva-Aguirre, Filiberto Ortiz-Chi and Gabriel Merino \n")
    out.write("Username: "+os.path.split(os.path.expanduser('~'))[1]+"\n\n")
    out.write("Date: "+fechayhora()+"\n\n")
    out.write("Total Job cpu time: %s\n\n"%(str(cputime+cputime_nmr)))
    out.write("Total Job cpu opt: %s\n\n"%(str(cputime)))
    out.write("Total Job cpu nmr: %s\n\n"%(str(cputime_nmr)))
    out.write("------------------------------------------------------------------------------------------------\n")
    out.write("The Geometry optimization and the isotropic shielding constants were computed using Gaussian G16\n")
    out.write("------------------------------------------------------------------------------------------------\n")
    out.write("{:<24}{:>47}\n".format("Geometry OPT","NMR"))                  #Datos a extraer del out
    out.write("------------------------------------------------------------------------------------------------\n")
    out.write("{:<47s}{:>24s}\n".format(str(key_opt),str(key_nmr)))
    out.write("\n\n||*********************Data for H1******************||\n\n")
    out.write("Scaling Factors H1           \n") #Datos a extraer del out
    out.write("------------------\n") 
    out.write('{0:5} {1:10} \n'.format("slope","     intercept"))
    out.write("------------------\n")
    out.write('{0:5.4f} {1:-10.4f}\n\n'.format(data_h.m, data_h.b))
    out.write("Performance     H1                \n") #Datos a extraer del out
    out.write("------------------\n")
    out.write('{0:5} {1:10} \n'.format("RMSD", "     R^2"))
    out.write("----------------- \n") 
    out.write('{0:5.4f} {1:-10.4f} \n'.format(data_h.rmsd,data_h.r2))    
    out.write("\n\n||*********************Data for C13******************|| \n\n")
    out.write("Scaling Factors C13\n") #Datos a extraer del out
    out.write("-----------------\n") 
    out.write('{0:5} {1:10} \n'.format("slope","     intercept"))
    out.write("------------------ \n")
    out.write('{0:5.4f} {1:-10.4f} \n\n'.format(data_c.m, data_c.b))
    out.write("Performance   C13\n") #Datos a extraer del out
    out.write("-----------------\n")
    out.write('{0:5} {1:10} \n'.format("RMSD", "     R^2"))
    out.write("----------------- \n") 
    out.write('{0:5.4f} {1:-10.4f} \n'.format(data_c.rmsd,data_c.r2))
    out.write("\n\n|******************************************|\n")
    out.write("|                Test Set Results          |\n") #Datos a extraer del out
    out.write("|******************************************|\n")   
    out.write("Name                                         Isotropic Values    Exp.Values    Shift Value(scaled)    Residual^2 \n")
    if len(new)==0:
        for imol in data:
            out.write("------------------------------------\n")
            out.write("%s\n"%(imol.im))
            out.write("%s\n"%(imol.i))
            ph=[]
            pc=[]
            for iatom in imol.atoms:
                symbol= iatom.s
                if symbol == "H":
                    out.write("%-5s %-30s" %(iatom.s,iatom.nz))
                    out.write('{0:-18.4} {1:16.2f} {2:-17.2f} {3:-18.4f} \n'.format(iatom.t, iatom.e, iatom.c, iatom.r))
                    ph.append((iatom.r))
                if symbol == "C":            
                    out.write("%-5s %-30s" %(iatom.s,iatom.nz))
                    out.write('{0:-18.4} {1:16.2f} {2:-17.2f} {3:-18.4f} \n'.format(iatom.t, iatom.e, iatom.c, iatom.r))
                    pc.append((iatom.r))
            if len(ph)!= 0 :
                rmsdh = (sum(ph)/len(ph))**0.5
                out.write("rmsd of H: %25.4F \n"%(rmsdh))
            if len(pc)!= 0 :   
                rmsdc = (sum(pc)/len(pc))**0.5
                out.write("rmsd of C: %25.4f \n"%(rmsdc))
            timecpu=datetime.timedelta(days=int(imol.ct[0]),hours=int(imol.ct[1]), minutes=int(imol.ct[2]), seconds=round(float(imol.ct[3])))
            timecpu_nmr=datetime.timedelta(days=int(imol.ct_nmr[0]),hours=int(imol.ct_nmr[1]), minutes=int(imol.ct_nmr[2]), seconds=round(float(imol.ct_nmr[3])))
            out.write("Job cpu time for optimization: %s\n"%(timecpu))
            if imol.chk==1: out.write("Chk: YES \n")
            out.write("Job cpu time for NMR: %s\n"%(timecpu_nmr))
        tabla_alldataset(path, data)

    else:
        pph=[]
        ppc=[]
        for imol in new:
            print("prev:",imol.i)
            out.write("------------------------------------\n")
            out.write("%s\n"%(imol.im))
            out.write("%s\n"%(imol.i))
            ph=[]
            pc=[]
            for iatom in imol.atoms:
                symbol= iatom.s
                if symbol == "H":
                    out.write("%-5s %-30s" %(iatom.s,iatom.nz))
                    out.write('{0:-18.4} {1:16.2f} {2:-17.2f} {3:-18.4f} \n'.format(iatom.t, iatom.e, iatom.c, iatom.r))
                    ph.append((iatom.r))
                    pph.append((iatom.r))
                if symbol == "C":            
                    out.write("%-5s %-30s" %(iatom.s,iatom.nz))
                    out.write('{0:-18.4} {1:16.2f} {2:-17.2f} {3:-18.4f} \n'.format(iatom.t, iatom.e, iatom.c, iatom.r))
                    pc.append((iatom.r))
                    ppc.append((iatom.r))
            if len(ph)!= 0 :
                rmsdh = (sum(ph)/len(ph))**0.5
                out.write("rmsd of H: %25.4F \n"%(rmsdh))
            if len(pc)!= 0 :   
                rmsdc = (sum(pc)/len(pc))**0.5
                out.write("rmsd of C: %25.4f \n"%(rmsdc))
        if len(pph)!= 0 :
            rmsdh = (sum(pph)/len(pph))**0.5
            out.write("\n\ngeneral rmsd of H: %25.4F \n"%(rmsdh))
        if len(ppc)!= 0 :   
            rmsdc = (sum(ppc)/len(ppc))**0.5
            out.write("general rmsd of C: %25.4f \n"%(rmsdc))
        timecpu=datetime.timedelta(days=int(imol.ct[0]),hours=int(imol.ct[1]), minutes=int(imol.ct[2]), seconds=round(float(imol.ct[3])))
        timecpu_nmr=datetime.timedelta(days=int(imol.ct_nmr[0]),hours=int(imol.ct_nmr[1]), minutes=int(imol.ct_nmr[2]), seconds=round(float(imol.ct_nmr[3])))
        out.write("Job cpu time for optimization: %s\n"%(timecpu))
        if imol.chk==1: out.write("Chk: YES \n")
        out.write("Job cpu time for NMR: %s\n"%(timecpu_nmr))
        print("new")
        tabla_alldataset(path, new)


    out.write("\n\nThat might sound boring, but I think the boring stuff is the stuff I remember the most\n\n")
    out.write("**** Don't get amxiaty, there was no problem whatsoever ****\n")
    out.close
    if tbl_comp != "-": comp_table(tbl_comp,data_c,data_h,keys,str(path.split('/')[-1]))
    print("outt")
    return keys[0]
