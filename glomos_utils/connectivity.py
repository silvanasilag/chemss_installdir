import numpy as np
from glomos_utils.atomic         import get_covalent_radius
#------------------------------------------------------------------------------------------
""" >>> from discriminate.connectivity import conectmx
    >>> conectmx(mol,0)
"""
def conectmx(moleculein, factor=1.2):
    natoms=moleculein.n
    matrixc=np.zeros(shape=(natoms,natoms),dtype=np.int)
    for iatom in range(natoms):
        si=moleculein.atoms[iatom].s
        ri=get_covalent_radius(si)
        xi=moleculein.atoms[iatom].xc
        yi=moleculein.atoms[iatom].yc
        zi=moleculein.atoms[iatom].zc
        for jatom in range(iatom+1,natoms):
            sj=moleculein.atoms[jatom].s
            rj=get_covalent_radius(sj)
            xj=moleculein.atoms[jatom].xc
            yj=moleculein.atoms[jatom].yc
            zj=moleculein.atoms[jatom].zc
            xx,yy,zz=xj-xi, yj-yi, zj-zi
            distance = np.sqrt(xx*xx+yy*yy+zz*zz)
            dist = distance/(ri + rj)
            if ( dist <= factor ):
                matrixc[iatom][jatom] = int(1)
                matrixc[jatom][iatom] = int(1)
    return matrixc
#------------------------------------------------------------------------------------------
def distance_matrix(moleculein):
    natoms=moleculein.n
    matrixd=np.zeros(shape=(natoms,natoms),dtype=np.float)
    for iatom in range(natoms):
        si=moleculein.atoms[iatom].s
        ri=get_covalent_radius(si)
        xi=moleculein.atoms[iatom].xc
        yi=moleculein.atoms[iatom].yc
        zi=moleculein.atoms[iatom].zc
        for jatom in range(iatom+1,natoms):
            sj=moleculein.atoms[jatom].s
            rj=get_covalent_radius(sj)
            xj=moleculein.atoms[jatom].xc
            yj=moleculein.atoms[jatom].yc
            zj=moleculein.atoms[jatom].zc
            xx,yy,zz=xj-xi, yj-yi, zj-zi
            distance = np.sqrt(xx*xx+yy*yy+zz*zz)
            dist = distance/(ri + rj)
            matrixd[iatom][jatom] = float(dist)
            matrixd[jatom][iatom] = float(dist)
    matrixd=np.array(matrixd)
    return matrixd
#------------------------------------------------------------------------------------------
""" HOW TO USE:
    >>> from inout.xyz import readxyz, writexyz
    >>> from discriminate.connectivity import conectmx, degree_matrix
    >>> mol00=readxyz('mol.xyz')
    >>> mx=conectmx(mol00,0)
    >>> my=degree_matrix(mx)
"""
def degree_matrix(adjmatrix):
    matrixadj=np.array(adjmatrix)
    vectordeg=matrixadj.sum(axis=1)
    n=len(matrixadj)
    matrixdeg=np.zeros(shape=(n, n),dtype=np.int)
    for iatoms in range(n): matrixdeg[iatoms,iatoms]=vectordeg[iatoms]
    return matrixdeg
#------------------------------------------------------------------------------------------
""" HOW TO USE:
    >>> from discriminate.connectivity import check_connectivity
    >>> print (check_connectivity(molb01,0,adjmatrix))
"""
def check_connectivity(moleculeina, adjmatrixref):
    mxx=conectmx(moleculeina)
    return ( np.array_equiv(mxx,adjmatrixref) or np.array_equal(mxx,adjmatrixref) )
#------------------------------------------------------------------------------------------
""" THIS ROUTINE RETURNS A LIST OF MOLECULES WITH THE SAME CONNECTIVITY 
    THAT THE REFERENCE MOLECULE.
    HOW TO USE:
    >>> from inout.xyz import readxyz, writexyz
    >>> from discriminate.connectivity import cutter_connectivity
    >>> molecule=readxyz('all_optim.xyz')
    >>> reference=readxyz('ethanol.xyz')
    >>> mol=cutter_connectivity(molecule,reference)
    >>> writexyz(mol,"output.xyz")
"""
def cutter_connectivity(moleculein, adjmatrixref):
    moleculeout, count=[],0
    for imol in moleculein:
        if ( check_connectivity(imol, adjmatrixref) ):
            moleculeout.extend([imol])
        else:
            count=count+1
            print("%4d %s ... DISCRIMINATED by Connectivity" %(count, imol.i))
    if count==0: print("No elements (or new elements) discriminated by Connectivity")
    return moleculeout
#------------------------------------------------------------------------------------------
""" THIS ROUTINE RETURNS 1 IF ALL THE NODES IN A MOLECULAR GRAPH ARE CONNECTED (FORM A
    SINGLE CONNECTED GRAPH), AND 0 OTHERWISE. 
    HOW TO USE:
    >>> from discriminate.connectivity import connectedmoleculargraph
"""
def connectedmoleculargraph(moleculein, factor=1.2):
    matrixadj=conectmx(moleculein,factor)
    matrixadjp=matrixadj.copy()
    nd=len(matrixadjp)
    vectord=np.zeros(shape=(nd),dtype=np.int)
    vectordp=np.array(vectord)
    vectord[0]=int(1)
    sumd=1
    while sumd != 0:
        vectord = np.dot(matrixadjp, vectord)
        vectord = vectord + vectordp
        for i, element in enumerate(vectord):
            if element > 1:
                vectord[i] = int(1)
        sumd=sum(vectord - vectordp)
        vectordp = vectord
    scg = 1 if ( nd==sum(vectord) ) else 0
    return scg
#------------------------------------------------------------------------------------------
""" THIS ROUTINE RETURNS A LIST OF MOLECULES CONNECTED SIMPLY.
    HOW TO USE:
    >>> from inout.xyz import readxyz, writexyz
    >>> from discriminate.connectivity import cutter_nonconnectedmoleculargraph
    >>> mol00=readxyz('all_optim.xyz')
    >>> mol01=cutter_nonconnectedmoleculargraph(mol00)
    >>> writexyz(mol01,'output.xyz')
"""
def cutter_nonconnectedmoleculargraph(moleculein, silence=1):
    from inout.getbilparam    import get_a_str
    log_file=get_a_str('output_file','glomos_out.txt')
    if silence==0:
        fopen = open(log_file,'a')
        print("disc_unconnected_opt = on", file=fopen)
        fopen.close()
    moleculeout, count=[], 0
    for imol in moleculein:
        if ( connectedmoleculargraph(imol)==1 ):
            moleculeout.extend([imol])
        else:
            count=count+1
            jj=str(count).zfill(5)
            fopen = open(log_file,'a')
            print ("%s %15s ... DISCRIMINATED: NON-SINGLE Connection" %(jj, imol.i), file=fopen)
            fopen.close()
    if count==0:
        fopen = open(log_file,'a')
        print ("ZERO elements discriminated by NON-SINGLE connection in molecular graph", file=fopen)
        fopen.close()
    elif moleculeout==[]:
        fopen = open(log_file,'a')
        print ("All the elements were discriminated by NON-SINGLE connection in molecular graph", file=fopen)
        print ("(-) You can turn the option disc_nonconnected to off", file=fopen)
        print ("(-) Or increase the initial_population to have a greater diversity", file=fopen)
        fopen.close()
        exit()
    return moleculeout
#------------------------------------------------------------------------------------------
