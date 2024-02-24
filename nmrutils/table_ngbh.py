import pandas as pd
import sys
import os

def tabla_alldataset(name, data):
    df1 = pd.DataFrame({"Z":[],
                        "Isotropic Values":[],
                        "Chemecal shift":[],
                        "-H": [],
                        "-C": [],
                        "-Cl": [],
                        "-O":[],
                        "-N" :[],
                        "-F" :[],
                        "Single":[],
                        "Dobles" :[],
                        "Triple":[]
                        })
    org_e = ['C','H','Cl','O','N','F']
    for imol in data:
        for iatom in imol.atoms:
            nbs=iatom.nb
            frec={} # Diccionario para contar la frecuencias
            for i in nbs:
                if i in frec:
                    frec[i]+= 1
                else:
                    frec[i]= 1

            if iatom.s=="H":
                s=1
                frec["single"]=1
                if len(nbs)!=1:
                    print(imol.i,iatom.nz)
                    sys.exit(nbs)

            if iatom.s=="C":
                s=6
                if len(nbs)==4:
                    frec["single"]=4
                elif len(nbs)==3:
                    frec["doble"]=1
                    frec["single"]=2
                elif len(nbs)==2:
                    frec["triple"]=1
                    frec["single"]=1
                else: 
                    print(imol.i,iatom.nz)
                    sys.exit(nbs)

            otros=([i for i in nbs if i not in org_e])
            df_1 = {
                    "Z": s,
                    "Isotropic Values":iatom.t,
                    "Chemecal shift":iatom.e,
                    "-H": frec.get('H', 0),
                    "-C": frec.get('C', 0),
                    "-Cl": frec.get('Cl', 0),
                    "-O":frec.get('O', 0),
                    "-N" :frec.get('N', 0),
                    "-F" :frec.get('F', 0),
                    "Single":frec.get('doble', 0),
                    "Dobles" :frec.get('doble', 0),
                    "Triple":frec.get('triple', 0),
                    }
        
            df_2 = pd.DataFrame([df_1]) # Crear un nuevo DataFrame con los nuevos datos
            df1 = pd.concat([df1, df_2], ignore_index=True) # Concatenar el nuevo DataFrame con el DataFrame existente

    tbl_compf = str(name + ".csv")
    print(tbl_compf)
    if os.path.isfile(tbl_compf):
        df = pd.read_csv(tbl_compf)
        dff = pd.concat([df, df1])
        dff = dff.drop_duplicates()
        dff.to_csv(tbl_compf, index=False)
    else:
        df1.to_csv(tbl_compf, index=False)