# -*- coding: utf-8 -*-

def chemss_inp():
    out ="INPUTNMR.txt"
    out=open(out,'a')
    out.write("memory_in_gb\t\t1\n")
    out.write("queue\t\t\tqfast\n")
    out.write("nodes\t\t\t5\n")
    out.write("walltime\t\t00:59:00\n")
    out.write("njobs\t\t\t10\n")
    out.write("timesleep\t\t2.0\n\n\n")

    out.write("foldername_base\t\t \n")
    out.write("foldername_new\t\t \n\n\n")
    
    out.write("calculus\t\tNO\n")
    out.write("statistics\t\tNO\n")
    out.write("new_mol_calculus\tNO\n")
    out.write("new_mol_statistics\tNO\n\n")
    out.write("---GAUSSIAN OPT---\n")
    out.write("# \n")
    out.write("---GAUSSIAN OPT---\n\n\n")
    out.write("---GAUSSIAN NMR---\n")
    out.write("# \n")
    out.write("---GAUSSIAN NMR---\n\n\n")

def chemssinp_edit(ii):
    fil=[]
    flag=1
    f=0
    replacement=''
    changes=''
    with open('INPUTNMR.txt','r') as f2:
        for line in f2:
            line = line.strip()
            lin=line.split()
            if "foldername_base" in line:
                if len(lin)==1:
                    changes=line + "\t\t" + str(ii[0])
                    replacement = replacement + changes + "\n"
                else:
                    r=lin[1]
                    changes=line.replace(str(r),str(ii[0]))
                    replacement = replacement + changes + "\n"
            elif "foldername_new" in line:
                if len(lin)==1:
                    changes=line + "\t\t" + str(ii[1])
                    replacement = replacement + changes + "\n" 
                else:
                    r=lin[1]
                    changes=line.replace(str(r),str(ii[1]))
                    replacement = replacement + changes + "\n"
            elif f==1 and "#" in line:
                changes= "#" +str(ii[2])
                replacement= replacement + changes+ "\n"
            elif f==2 and "#" in line:
                changes= "#" +str(ii[3])
                replacement= replacement + changes+ "\n"
            elif "---GAUSSIAN OPT---" in line and f==0 and flag==1:
                replacement= replacement + line+ "\n"
                f=1
            elif "---GAUSSIAN NMR---" in line and f==1 and flag==1:
                replacement= replacement + line+ "\n"
                f=2
            else: replacement= replacement + line+ "\n"
    fout = open("INPUTNMR.txt", "w")
    fout.write(replacement)
    fout.close()