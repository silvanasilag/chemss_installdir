# -*- coding: utf-8 -*-
#codigo con el calculo del python

##modulos
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
import os

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
def xy(df,n):
    dfn = df[df['Z'] == n]
    x = list(dfn['Chemical Shift'])
    y = list(dfn['Isotropic Values'])
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
            if iatom.s=='H':
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
    #-----------------------------------create cvs 
    datos = list(zip(x, y))
    ruta_archivo = "datos_xy.csv"
    with open(ruta_archivo, 'w', newline='') as archivo_csv:
    # Crear un objeto escritor CSV
        escritor_csv = csv.writer(archivo_csv)
    # Escribir los datos en el archivo CSV
        escritor_csv.writerow(["Chemical Shift", "Isotropic Values"])  # Escribir encabezados
        escritor_csv.writerows(datos)
    #-----------------------------------create cvs 
    yscal = []
    
    print("largo de los datos",len(x))
    data.m, data.b, r_value, p_value, std_err = stats.linregress(x, y)
    data.r2 = r_value**2
    for i in range(len(y)):
        ys = (data.b - y[i])/-data.m
        yscal.append(ys)
    data.rmsd = (mean_squared_error(x,yscal))**0.5
    print("slope: %.4f  \nintercept:%.4f   \nrmsd: %.4f    \nR**2:%.4f \n" % (data.m, data.b,data.rmsd,data.r2))
    return data
#--------------------------------------------------- 
def scale(data, data_h,data_c):
    for imol in data:
        for iatom in imol.atoms:
            symbol= iatom.s
            if symbol == "H":
                ys = (data_h.b - iatom.t)/-data_h.m
                res = (abs(iatom.e-ys))**2
                iatom.r =res
                iatom.c =ys
            if symbol == "C":
                yc = (data_c.b - iatom.t)/-data_c.m
                res = (abs(iatom.e-yc))**2
                iatom.r =res
                iatom.c =yc
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


