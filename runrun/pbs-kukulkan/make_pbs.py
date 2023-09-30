import os.path
from conf.conf import config_file
#-----------------------------------------PBS KUKULKAN-------------------------------------------------
def make_apbs_gaussian(ppj, ram, queue_name, walltime, basename):
    filepbs=basename+'.pbs'
    filein =basename+'.inp'
    fileout=basename+'.out'
    #-------------------------------------------------------
    conf=config_file('gaussian')
    #-------------------------------------------------------
    if not os.path.isfile(fileout):
        fh=open(filepbs,'w')
        print("#!/bin/bash", file=fh)
        print("#PBS -N %s" %(basename), file=fh)
        print("#PBS -l nodes=1:ppn=%d" %(ppj), file=fh)
        print("#PBS -l mem=%dgb" %(ram), file=fh)
        print("#PBS -q %s" %(queue_name), file=fh)
        print("#PBS -l walltime=%s\n" %(walltime), file=fh)
        print("cd $PBS_O_WORKDIR\n", file=fh)
        #-------------------------------------------------------
        for ii in range(len(conf)): fh.write(conf[ii])
        #-------------------------------------------------------
        print("mkdir -p $GAUSS_SCRDIR\n", file=fh)
        print("g16 < %s > %s\n" %(filein, fileout), file=fh)
        print ("rm -rf $GAUSS_SCRDIR", file=fh)
        fh.close()
    else:
        os.system("rm -f -v "+filepbs)
#------------------------------------------------------------------------------------------
