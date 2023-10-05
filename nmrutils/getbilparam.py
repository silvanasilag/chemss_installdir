import re
#------------------------------------------------------------------------------------------
def get_a_float(strchain, defaultvalue):
    bilatu_file='INPUTNMR.txt'
    bilfile=open(bilatu_file,"r")
    finalvalue=float(defaultvalue) #DEFAULT VALUE
    for line in bilfile:
        line=line.strip(' \t\n\r')
        if len(line.strip()) != 0 :
           li = line.lstrip()
           if not li.startswith("#"):
              readline=line.split()
              if len(readline) == 2:
                 data0=readline[0].strip('\t\n\r') 
                 data1=readline[1].strip('\t\n\r')
                 if data0.lower() == str(strchain): finalvalue=float(data1)
    bilfile.close()
    return finalvalue
#------------------------------------------------------------------------------------------
def get_a_int(strchain, defaultvalue):
    bilatu_file='INPUTNMR.txt'
    bilfile=open(bilatu_file,"r")
    finalvalue=int(defaultvalue) #DEFAULT VALUE
    for line in bilfile:
        line=line.strip(' \t\n\r')
        if len(line.strip()) != 0 :
           li = line.lstrip()
           if not li.startswith("#"):
              readline=line.split()
              if len(readline) == 2:
                 data0=readline[0].strip('\t\n\r') 
                 data1=readline[1].strip('\t\n\r')
                 if data0.lower() == str(strchain): finalvalue=int(data1)
    bilfile.close()
    return finalvalue
#------------------------------------------------------------------------------------------
def get_a_str(strchain, defaultvalue):
    bilatu_file='INPUTNMR.txt'
    bilfile=open(bilatu_file,"r")
    finalvalue=str(defaultvalue)
    for line in bilfile:
        line=line.strip(' \t\n\r')
        if len(line.strip()) != 0 :
           li = line.lstrip()
           if not li.startswith("#"):
              readline=line.split()
              if len(readline) == 2:
                 data0=readline[0].strip('\t\n\r') 
                 data1=readline[1].strip('\t\n\r')
                 if data0.lower() == str(strchain): finalvalue=str(data1)
    bilfile.close()
    return finalvalue
#------------------------------------------------------------------------------------------
def read_block_of_inp(id):
    bilatu_file='INPUTNMR.txt'
    bilfile=open(bilatu_file,"r")
    chainchar='---'+id.upper()+'---'
    printer=0
    data_block=[]
    for line in bilfile:
         line=line.strip(' \t\n\r')
         if len(line.strip()) != 0 :
            lin = line.lstrip()
            if lin.startswith(chainchar): printer=1+printer
            if printer == 1 and not lin.startswith(chainchar): data_block.append(line)
    bilfile.close()
    return data_block 
#------------------------------------------------------------------------------------------
"""
def get_key(i):
  d="-"
  l=0
  i=i.replace(" =","=")
  i=i.replace("= ","=")
  i=i.strip(" ")
  ii=i.strip("#")
  ii=ii.split()
  if "dispersion" in i.lower():
    ld=i.split("=")
    for j in ld :
      if l==1:
        d=str(j)
        break
      if "dispersion" in j.lower():l=1
  if "/" in ii[0]: k=str(ii[0])
  else: k=str(ii[0]+"/"+ii[1])
  if d != "-": k= k+"_"+d
  
  k=k.replace("/","_")
  return k
"""
#------------------------------------------------------------------------------------------

def key_norm(key,fun):
    key=key.strip("#")
    fun_n=re.escape(fun)
    pt=re.compile(f'(?i){fun_n}')
    w_key=pt.split(key)[0]
    w_key=w_key.strip()
    w_key=w_key.replace(" ","_")
    if "empiricaldispersion" in key.lower() :
      key=key.replace(" =","=")
      key=key.replace("= ","=")
      l=key.split()
      for i in l:     
        i=i.strip()
        if "empiricaldispersion" in i.lower():
          ii=i.split("=")[1]
          w_key=w_key+"_"+ ii
    return w_key

