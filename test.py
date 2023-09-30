#!/usr/local/anaconda3/bin/python
# -*- coding: utf-8 -*-

#Silvana Silva Aguirre, CINVESTAV Unidad Mérida,2019

##modulos
import sys
import os



from nmrutils.getbilparam     import read_block_of_inp,get_key
from nmrutils.out_writer 		import key_norm


key="# M062X 6-311+G(2d,p)   opt=CARTESIAN EmpiricalDispersion = GD3 #int=ultrafin"
kenmr="# M062X 6-311+G(2d,p)   NMR=CARTESIAN EmpiricalDispersion = GD3 #int=ultrafin"

print(get_key(kenmr))
print(key_norm(kenmr,"NMR"))
