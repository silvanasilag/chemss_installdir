import os
import numpy as np
from utils.libmoleculas import copymol, translate_to_cm, molecular_radius
from utils.atomic import get_covalent_radius, get_chemical_color
#------------------------------------------------------------------------------------------
### DEFINICION DE LOS ENLACES VALIDOS
iones=['Ti','Al', 'Zr']
fd=2.8
bwd1='0.04'                   # 0.04 GRUESO DE LOS ENLACES EN iones
bwd2='0.10'                   # 0.10 GRUESO DE LOS ENLACES
sphere_factor=float(0.55)     # 0.5 (g 1.5)
projection='orthographic'     # 'perspective'
reflection='reflection 0.0'
reflection_model='phong 1.0'  # 'specular'
radio_factor=float(1.243)     # 1.3, 1.4 1.6 #ALCANCE DE LOS ENLACES
tmit = 0.35                   # transmit (transparencia) 0.35 0.75
#-----------------------------------------------------------------------------------------
cadena1="cylinder {<%9.6f, %9.6f, %9.6f> <%9.6f, %9.6f, %9.6f> %s pigment {color rgb %s transmit %5.2f} }"
cadena2="cylinder {<%9.6f, %9.6f, %9.6f> <%9.6f, %9.6f, %9.6f> %s pigment {color rgb %s } finish {%s %s}}" 
cadena3="sphere {<%9.6f, %9.6f, %9.6f>, %9.6f pigment {color rgb %s } finish {%s %s}}" 
cadena4="sphere {<%9.6f, %9.6f, %9.6f>, %9.6f pigment {color rgb %s transmit %5.2f}}" 
def write_povray(poscarin, basename, f=fd, quality='md'):
    poscarxx=copymol(poscarin)
    translate_to_cm(poscarxx)
    a,b,c=molecular_radius(poscarxx)
    factor=f*c
    background1,background2,background3='2.0','2.0','2.0'
    light11,light12,light13='1','1','1'
    light21,light22,light23='1','1','1'
    camara_rotate='y*0.0'
    camara_loca1,camara_loca2,camara_loca3='0','0',str(factor)
    camara_look1,camara_look2,camara_look3='0','0','0'
    name_pov=basename+'.pov' 
    name_png=basename+'.png'
    opnew = open(name_pov,"w")
    print("#version 3.7;\n global_settings {\n assumed_gamma 1.0\n }\n", file=opnew)
    print("# include \"colors.inc\"", file=opnew)
    print("# include \"textures.inc\"\n", file=opnew)
    print("background{color rgb<2.0, 2.0, 2.0>}\n", file=opnew)
    print("light_source {< 10, -8, -8> color rgb <1, 1, 1>}", file=opnew)
    print("light_source {< -8,  8,  8> color rgb <1, 1, 1>}\n", file=opnew)
    print("camera {%s location <0, 0, %4.2f> look_at <0, 0, 0> rotate y*0.0}\n" %(projection, factor), file=opnew)
    for iatom in range(poscarxx.n):
        symi=poscarxx.atoms[iatom].s
        radii_i= get_covalent_radius(symi)
        xxi=poscarxx.atoms[iatom].xc
        yyi=poscarxx.atoms[iatom].yc
        zzi=poscarxx.atoms[iatom].zc
        izf=poscarxx.atoms[iatom].zf
        for jatom in range(iatom+1,poscarxx.n):
            symj=poscarxx.atoms[jatom].s
            radii_j= get_covalent_radius(symj)
            xxj=poscarxx.atoms[jatom].xc
            yyj=poscarxx.atoms[jatom].yc
            zzj=poscarxx.atoms[jatom].zc
            jzf=poscarxx.atoms[jatom].zf
            rr=np.sqrt((xxi-xxj)**2+(yyi-yyj)**2+(zzi-zzj)**2)
            uijx=(xxj-xxi)/rr
            uijy=(yyj-yyi)/rr
            uijz=(zzj-zzi)/rr
            rt=float(rr)/(radii_i+radii_j)
            if rt < radio_factor:
                bwd=bwd1 if ((symi in iones) and (symj in iones)) else bwd2
                kolori=get_chemical_color(symi)
                kolorj=get_chemical_color(symj)
                otf=float(0.85)
                xxic=xxi+sphere_factor*otf*radii_i*uijx
                yyic=yyi+sphere_factor*otf*radii_i*uijy
                zzic=zzi+sphere_factor*otf*radii_i*uijz
                xxjc=xxj-sphere_factor*otf*radii_j*uijx
                yyjc=yyj-sphere_factor*otf*radii_j*uijy
                zzjc=zzj-sphere_factor*otf*radii_j*uijz
                xpm=(xxic+xxjc)/2.0
                ypm=(yyic+yyjc)/2.0
                zpm=(zzic+zzjc)/2.0
                if (izf == 'F'):
                    lista1=(-xxic,yyic,zzic,-xpm, ypm, zpm, bwd,kolori,tmit)
                    print(cadena1 %lista1, file=opnew)
                else:
                    lista2=(-xxic,yyic,zzic,-xpm, ypm, zpm, bwd,kolori,reflection_model,reflection)
                    print(cadena2 %lista2, file=opnew)
                if (jzf == 'F'):
                    lista3=(-xpm, ypm, zpm, -xxjc,yyjc,zzjc,bwd,kolorj,tmit)
                    print(cadena1 %lista3, file=opnew)
                else:
                    lista4=(-xpm, ypm, zpm, -xxjc,yyjc,zzjc,bwd,kolorj,reflection_model,reflection)
                    print(cadena2 %lista4, file=opnew)
    for iatom in poscarxx.atoms:
        sym=iatom.s
        radii_div2=sphere_factor*(get_covalent_radius(sym))
        kolor= get_chemical_color(sym)
        xx, yy, zz=iatom.xc, iatom.yc, iatom.zc
        s0="pigment {color rgb"
        if iatom.zf == 'T':
            lista5=(-xx,yy,zz,radii_div2,kolor,reflection_model,reflection)
            print(cadena3 %lista5, file=opnew)
        if iatom.zf == 'F':
            lista6=(-xx,yy,zz,radii_div2,kolor,tmit)
            print(cadena4 %lista6, file=opnew)
    opnew.close()
    if quality == 'hd': width,height=6000,4500
    if quality == 'md': width,height=1600,1200
    if quality == 'ld': width,height=1200,900
    comp="povray +A Display=Off Output_File_Type=N All_Console=Off Width="+str(width)+" Height="+str(height)+" "+name_pov+" > /dev/null 2>&1"
    os.system(comp)
    print("Output= %s %s" %(name_pov, name_png))
    doma="rm -f "+str(name_pov)
    os.system(doma)
#-----------------------------------------------------------------------------------------
