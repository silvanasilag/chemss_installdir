# -*- coding: utf-8 -*-
#codigo con el calculo del python

##modulos
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
import os
import sys

from scipy import stats
from sklearn.metrics import mean_squared_error
from datetime import datetime, timedelta

#--------------------------------------------------- 
class DataNMR:
    def __init__(self, atom, slope, intercept, r_squere, rmsd):
        self.s=atom
        self.m=slope
        self.b=intercept
        self.r2=r_squere
        self.rmsd=rmsd

#---------------------------------------------------
def bias(element,all_df,data):
    z = str(element + '-')
    df = all_df[all_df[z] > 0]
    x = list(df['Chemical Shift'])
    y = list(df['Isotropic Values'])
    print(z,"---------leng data---------",len(y))
    xscal = []
    n, p = 0, 0
    b = data.b
    for i in range(len(y)):
        xs = (y[i] - b) / data.m
        xscal.append(xs)
        p, n = (p + 1, n) if xs > x[i] else (p, n + 1)
    rmsd = (mean_squared_error(x, xscal)) ** 0.5
    print("Scale RMSD", z, "=", rmsd)
    c = -0.05 if p > n else 0.05
    while True:
        bi=b
        xscal = []
        b = b + c
        for i in range(len(y)):
            xs = (y[i] - b) / data.m
            xscal.append(xs)
        rmsd_n = (mean_squared_error(x, xscal)) ** 0.5
        #print("Scale RMSD", z, "=", rmsd_n,b)
        if rmsd_n > rmsd:
            print("The minimals is Scale RMSD", z, "=", rmsd,bi)
            break
        else:
            rmsd = rmsd_n
    return bi

def xy(df,n,f):
    dfn = df[df['Z'] == n]
    x = list(dfn['Chemical Shift'])
    if f=="y":
        y = list(dfn['Isotropic Values'])
    if f=="yc":
        y = list(dfn['Scale_Shift'])
    return x,y
def df_dataset(data,cputime,cputime_nmr):
    atoms_vecinos = []
    for imol in data:
        timecpu = timedelta(days=int(imol.ct[0]), hours=int(imol.ct[1]), minutes=int(imol.ct[2]),seconds=round(float(imol.ct[3])))
        timecpu_nmr = timedelta(days=int(imol.ct_nmr[0]), hours=int(imol.ct_nmr[1]), minutes=int(imol.ct_nmr[2]),seconds=round(float(imol.ct_nmr[3])))
        cputime = timecpu + cputime
        cputime_nmr = timecpu_nmr + cputime_nmr
        for iatom in imol.atoms:
            nbs=iatom.nb
            for n in nbs:
                atoms_vecinos.append(n)
    nb=list(set(atoms_vecinos))
    print("Atomos vecinoss",nb)
    dic_df = {
        "Z": [],
        "Isotropic Values": [],
        "Chemical Shift": [],
        "Single": [],
        "Dobles": [],
        "Triple": []
    }
    for elemento in atoms_vecinos:
        dic_df[f"{elemento}-"] = []
    df1 = pd.DataFrame(dic_df)
    for imol in data:
        for iatom in imol.atoms:
            if not (iatom.s == 'H' or iatom.s == 'C'):
                print(imol.im,imol.i,iatom.s,iatom.t,iatom.e)
                print("PANICO SATANICOOO",iatom.s )
            nbs=iatom.nb
            frec={} # Diccionario para contar la frecuencias
            for i in nbs:
                if i in frec:
                    frec[i]+= int(1)
                else:
                    frec[i]= int(1)
            if iatom.s=='H': #hay que borrar esto probablemente
                if iatom.t  == 681.3693:print("AHHHHHHH!!!!---------------------------------11-")
                s=1
                frec["single"]=int(1)
                if len(nbs)!=1:
                    sys.exit(nbs)
            if iatom.s=='C':
                if iatom.t  == 681.3693:print("AHHHHHHH!!!!---------------------------------2-")
                s=6
                if len(nbs)==4:
                    frec["single"]=int(4)
                elif len(nbs)==3:
                    frec["doble"]=int(1)
                    frec["single"]=int(2)
                elif len(nbs)==2:
                    frec["triple"]=int(1)
                    frec["single"]=int(1)
                else:
                    sys.exit(nbs)
            #out.write('\n'+iatom.nz+" "+iatom.s+" : "+str(iatom.nb)+"  s:"+str(frec.get('single', 0))+"  d:"+str(frec.get('doble', 0))+"  t:"+str(frec.get('triple', 0)))
            df_1 = {
                "Z": s,
                "Isotropic Values": iatom.t,
                "Chemical Shift": iatom.e,
                "Single": frec.get('single', 0),
                "Dobles": frec.get('doble', 0),
                "Triple": frec.get('triple', 0)
            }
            for elemento in atoms_vecinos:
                df_1[f"{elemento}-"] = frec.get(elemento, 0)
            df_2 = pd.DataFrame([df_1]) # Crear un nuevo DataFrame con los nuevos datos
            df1 = pd.concat([df1, df_2], ignore_index=True) # Concatenar el nuevo DataFrame con el DataFrame existente

    return df1,cputime,cputime_nmr
#--------------------------------------------------- 
def stat(x,y,data):
    yscal = []
    print("largo de los datos",len(x))
    data.m, data.b, r_value, p_value, std_err = stats.linregress(x, y)
    data.r2 = r_value**2
    #------------------------
    for i in range(len(y)):
        ys = (data.b - y[i])/-data.m
        yscal.append(ys)
    data.rmsd = (mean_squared_error(x,yscal))**0.5
    print("slope: %.4f  \nintercept:%.4f   \nrmsd: %.4f    \nR**2:%.4f \n" % (data.m, data.b,data.rmsd,data.r2))
    return data
#---------------------------------------------------
def scale(data,data_h,data_c,df,dir_bias):
    resl=[]
    '''

    :param data: Molecules out
    :param data_h: stadistic  informacion from the regression of the 1H chimical shift vs isotropic value
    :param data_c:  stadistic  informacion from the regression of the 13C chimical shift vs isotropic value
    :param df:
    :return:
    '''
    for imol in data:
        for iatom in imol.atoms:
            symbol= iatom.s
            if symbol == "H":
                ys = (data_h.b - iatom.t)/-data_h.m
                res = (abs(iatom.e-ys))**2
                iatom.r =res
                iatom.c =ys
            if symbol == "C":
                b=data_c.b
                for element in iatom.nb:
                    if element not in ["H", "C"]:
                        if element not in dir_bias:
                            dir_bias[element]=bias(element,df,data_c)
                        b = dir_bias[element]
                yc = (b - iatom.t)/-data_c.m
                res = (abs(iatom.e-yc))**2
                iatom.r =res
                iatom.c =yc
            resl.append(iatom.c)
    df.loc[:, "Scale_Shift"]=resl
    return dir_bias,df
#--------------------------------------------------- 
def splot(x, y,fname,slope,inter,r2,name,clr):
    c = []
    xx = sorted(x)  
    l = np.linspace(xx[0],xx[-1], 100)
    for ix in l:
        cc = inter + slope * ix
        c.append(cc)   
    plt.clf()
    plt.plot(x, y, 'o',color=clr)
    plt.plot(l, c,'-k', lw=1.2)
    plt.xlabel('Experimental Chimical Shifts')
    plt.ylabel('Computed Isotropic Values')
    plt.title(name)
    plt.savefig(fname, dpi=700, facecolor='w', edgecolor='w', format= 'png')
#---------------------------------------------------


