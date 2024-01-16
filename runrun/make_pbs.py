import os.path
from conf.conf import config_file
#------------------------------------PBS IXCEL------------------------------------------------------
def make_apbs_gaussian(nproc, ram, queue, walltime, basename):
    filepbs=basename+'.pbs'
    filein =basename+'.inp'
    fileout=basename+'.out'
    gaussian=16
    #-------------------------------------------------------
    conf=config_file('gaussian')
    #-------------------------------------------------------
    if not os.path.isfile(fileout):
        with open(filepbs,'w') as pbs:
            pbs.write("#!/bin/bash                                                                      \n")
            pbs.write("#PBS -N %s                                                            \n" % basename)
            pbs.write("#PBS -l nodes=1:ppn=%s                                                   \n" % nproc)
            pbs.write("#PBS -q %s                                                               \n" % queue)
            pbs.write("#PBS -l walltime=%s                                                   \n" % walltime)
            pbs.write("                                                                                 \n")
            pbs.write("cd $PBS_O_WORKDIR                                                                \n")
            pbs.write("                                                                                 \n")
            #pbs.write("env > entorno-$PBS_JOBID.txt                                                     \n")
            pbs.write("export g%sroot=/usr/local/g16_AVX/                                    \n" % gaussian)
            pbs.write("source $g%sroot/g%s/bsd/g%s.profile                            \n" % ((gaussian,)))
            #pbs.write("export GAUSS_SCRDIR=$PWD/$USER-scratch/$PBS_JOBNAME-$PBS_JOBID                   \n")
            pbs.write("export GAUSS_SCRDIR=/scratch/$USER/$PBS_JOBNAME-$PBS_JOBID                       \n")
            pbs.write("mkdir -p $GAUSS_SCRDIR                                                           \n")
            pbs.write("sed -i '1i %%NProcShared=%s' %s                               \n" % (nproc, filein))
            pbs.write("                                                                                 \n")
            pbs.write("g%s < %s > %s                                     \n" % (16, filein, fileout))
            pbs.write("                                                                                 \n")
            pbs.write("mv $GAUSS_SCRDIR $PBS_O_WORKDIR                                                  \n")
            pbs.write("rm -rf $PBS_JOBNAME-$PBS_JOBID                                                   \n")
            #pbs.write("rm entorno-$PBS_JOBID.txt                                                        \n")
            pbs.write("rm $PBS_JOBNAME.pbs                                                              \n")
    else:
        os.system("rm -f -v "+filepbs)
#------------------------------------------------------------------------------------------
