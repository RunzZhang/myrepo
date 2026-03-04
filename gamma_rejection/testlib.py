# routine to test librefprop.so

#import os
#os.chdir("D:\\Pico\\PICOcode\\REFPROP")
#exec(open("testlib.py").read())

import os

import ctREFPROP.ctREFPROP as ct

RPPath = "D:/Pico/PICOcode/REFPROP"
RP = ct.REFPROPFunctionLibrary(RPPath)
print(RP.RPVersion())
RP.SETPATHdll(RPPath)

ierr, herr = RP.SETUPdll(1, "D:/Pico/PICOcode/REFPROP/FLUIDS/R218.FLD", "", "DEF")
print(ierr)
print(herr)

p, rhol, rhov, xliq, xvap, ierr, herr = RP.SATTdll(293, [1], 1)
print(p, rhol, rhov, xliq, xvap, ierr, herr)
