import subprocess
from sys import stdout
import time 
import glob
import os.path
#------------------------------------------------------------------------------------------
printf = stdout.write 
#------------------------------------------------------------------------------------------
def queued_processes():
    workid=[]
    qs = os.popen('qstat')
    lines = qs.readlines()
    for line in lines[2:]:
        pid,job,user,time,status,queue = line.split()
        workid.append(pid)
    return workid
#------------------------------------------------------------------------------------------
def qsub(args,iii):
    ii=str(iii+1).zfill(4)
    mytime=time.strftime("%c")
    base_cmd = "qsub"
    cmd = [base_cmd] + args
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, errmsg = p.communicate()
    if p.returncode != 0: raise OSError(p.returncode, "command %r failed: %s" % (" ".join(cmd), errmsg))
    anspbs=output.strip()
    anspbs=anspbs.decode('ascii')
    print("## %s Date %s --> %s JobID: %s" %(ii, mytime, args, anspbs))
    return anspbs
#------------------------------------------------------------------------------------------
def clean_pbs():
    print("Cleaning pbs files ...")
    for file in sorted(glob.glob("*.pbs")):
        base=file.split('.')
        namein=base[0]
        if os.path.isfile(file): os.remove(file)
        if os.path.isfile(namein+".e[0-9][0-9][0-9][0-9]*"): os.remove(namein+".e[0-9][0-9][0-9][0-9]*")
        if os.path.isfile(namein+".o[0-9][0-9][0-9][0-9]*"): os.remove(namein+".o[0-9][0-9][0-9][0-9]*")
#------------------------------------------------------------------------------------------
def send_pbs_files_to_queue(njobs, time_sleep):#----------------------Este es el main que manda los pbs
    jobslist, pids, jobindex, inqueue = [], [], 0, 1
    for file in sorted(glob.glob("*.pbs")): jobslist.append(file)
    if len(jobslist) == 0:
        print("There are not pbs files type or they already have been done")
        return 0 
    totaljobs=len(jobslist)
    print("Total jobs found = %d" %(totaljobs))
    timein=time.strftime("%c")
    print("Enter to the pool at : Date %s" %(timein))
    inqueue=totaljobs
    while inqueue >= 1:
        time.sleep(float(time_sleep))
        workid=queued_processes() 
        work_queue=[]
        for elem in pids:
            if elem in workid: work_queue.append(elem)
        inqueue=len(work_queue)
        if (njobs > inqueue) and (jobindex < totaljobs): 
            time.sleep(1.0)
            inqueue=1 ## at least one in queue
            spearhead=[jobslist[jobindex]]
            jobid=qsub(spearhead, jobindex)
            pids.append(jobid)
            jobindex=jobindex+1
    clean_pbs()
    timeout=time.strftime("%c")
    #print("Out of the pool at : Date %s" %(timeout))
    return 0
#------------------------------------------------------------------------------------------
