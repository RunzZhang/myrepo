# routine to test librefprop.dylib
#
#   Built on macOS 10.14 using https://github.com/usnistgov/REFPROP-cmake
#   with homebrew-provided gcc/gfortran (8.3.0), cmake (3.14.2)
#
#   prereq: pip install ctREFPROP

import os
import ctREFPROP.ctREFPROP as ctRP

#RPpath = os.path.join(os.path.expanduser("~"), "github/picoexperiment/PICOcode/REFPROP")
RPpath = os.path.dirname(os.path.abspath(__file__))
lib_path = os.path.join(RPpath, "librefprop.dylib")

RP = ctRP.REFPROPFunctionLibrary(lib_path)
print("Loaded REFPROP", RP.RPVersion())

fluid_path = os.path.join(RPpath, "FLUIDS/R218.FLD")#"FLUIDS/C3F8.FLD")
ierr, herr = RP.SETUPdll(1, fluid_path, "", "DEF")

fluid_info = RP.INFOdll(1)
# fluid_info_dict = fluid_info._asdict()

print("INFOdlloutput (iterable namedtuple) for C3F8")
for prop in fluid_info._fields:
    print(prop, ":", getattr(fluid_info, prop))
