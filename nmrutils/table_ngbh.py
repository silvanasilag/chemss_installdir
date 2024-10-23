import pandas as pd
import sys
import os
from datetime import datetime, timedelta

def tabla_alldataset(data,cputime,cputime_nmr):
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
#    out =name+"_neighbors_analisis.txt"
#    out=open(out,'a')
#    out.write("CHEMical Shift Scaler neighbors analisis \n\n")
#    out.write("Atom   Nuclei   neighbors   Single bond  Double bond   Triple bond")
    for imol in data:
#        out.write("\n------------------------------------\n")
#       out.write("%s\n"%(imol.im))
#        out.write("%s\n"%(imol.i))
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
    """tbl_compf = str(name + ".csv")
    print(tbl_compf)
    print(os.getcwd())
    if os.path.isfile(tbl_compf):
        df = pd.read_csv(tbl_compf)
        dff = pd.concat([df, df1])
        dff = dff.drop_duplicates()
        dff.to_csv(tbl_compf, index=False)
    else:
        df1.to_csv(tbl_compf, ind