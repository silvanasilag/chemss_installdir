# -*- coding: utf-8 -*-
#codigo con el calculo del python

##modulos
import numpy as np
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
def xy(data,cputime,elepsetime,cputime_nmr,elepsetime_nmr):
    xn_h,xn_c,yn_h,yn_c=[],[],[],[]
    for imol in data:
        timecpu=timedelta(days=int(imol.ct[0]),hours=int(imol.ct[1]), minutes=int(imol.ct[2]), seconds=round(float(imol.ct[3])))
        timelaps=timedelta(days=int(imol.et[0]),hours=int(imol.et[1]), minutes=int(imol.et[2]), seconds=round(float(imol.et[3])))
        timecpu_nmr=timedelta(days=int(imol.ct_nmr[0]),hours=int(imol.ct_nmr[1]), minutes=int(imol.ct_nmr[2]), seconds=round(float(imol.ct_nmr[3])))
        timelaps_nmr=timedelta(days=int(imol.et_nmr[0]),hours=int(imol.et_nmr[1]), minutes=int(imol.et_nmr[2]), seconds=round(float(imol.et_nmr[3])))
        cputime= timecpu+cputime
        cputime_nmr= timecpu_nmr+cputime_nmr
        elepsetime= timelaps + elepsetime
        elepsetime_nmr= timelaps_nmr + elepsetime_nmr
        for iatom in imol.atoms:
            symbol= iatom.s
            if symbol == "H":
                xn_h.append(iatom.e)
                yn_h.append(iatom.t)
            if symbol == "C":
                xn_c.append(iatom.e)
                yn_c.append(iatom.t)
    return xn_h,xn_c,yn_h,yn_c,cputime,elepsetime,cputime_nmr,elepsetime_nmr
#--------------------------------------------------- 
def stat(x,y,data):
    #-----------------------------------create cvs 
    datos = list(zip(x, y))
    print(os.getcwd())
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
    """
    data.m, data.b, r_value, p_value, std_err = stats.linregress(y,x)
    data.r2 = r_value**2
    for i in range(len(y)):
        ys = (data.b +(y[i])*data.m)
        yscal.append(ys)
    data.rmsd = (mean_squared_error(x,yscal))**0.5
    print("slope: %.4f  \nintercept:%.4f   \nrmsd: %.4f    \nR**2:%.4f \n" % (data.m, data.b,data.rmsd,data.r2))    
    """
    return data
#--------------------------------------------------- 
def scale(data, data_h,data_c):
    print("data")
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
