import pandas as pd
import matplotlib.pyplot as plt
import csv
import numpy as np
import re

before_KE_path = "./Cf252_spontanous.txt"
after_KE_path = "/data/runzezhang/Geant4Simulaions/g411_TN/CF252Sap7.5Poly6Pad10KE.dat"
plot_path= "/data/runzezhang/result/TN_sims_D/plot/"

def dat_to_list(path):
    energy_list = [] # in MeV
    intensity_list = [] # I don't know unit
    with open(path, "r") as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) == 3:
                energy_list.append(float(parts[1]))
                intensity_list.append(float(parts[2]))

    return energy_list, intensity_list

def raw_to_list(path):
    energy = []
    intensity = []

    def parse_endf_float(s):
        """Convert ENDF-style float (e.g., '2.000000-5') to Python float."""
        s = s.strip()
        if not s:
            return None
        if '+' in s[1:]:
            s = s.replace('+', 'E+', 1)
        elif '-' in s[1:]:
            s = s.replace('-', 'E-', 1)
        try:
            return float(s)
        except ValueError:
            return None

    with open(path, "r") as f:
        lines = f.readlines()

    # Skip first 8 lines (metadata headers)
    for line in lines[8:]:
        if re.match(r'\s*\d+\.\d{6}[+-]\d\s+\d+\.\d{6}[+-]\d', line):
            for i in range(0, 66, 22):
                s1 = line[i:i + 11]
                s2 = line[i + 11:i + 22]
                val1 = parse_endf_float(s1)
                val2 = parse_endf_float(s2)
                if val1 is not None and val2 is not None:
                    energy.append(val1)
                    intensity.append(val2)

    energy = np.array(energy)
    intensity = np.array(intensity)

    nonzero = intensity > 0
    energy = energy[nonzero]
    intensity = intensity[nonzero]

    intensity /= np.trapz(intensity, energy)

    return energy, intensity

#fast neutron fs thermal neutron tn
(ene_fn, intens_fn)=raw_to_list(before_KE_path)

(ene_tn,intens_tn) = dat_to_list(after_KE_path)
print((ene_fn,intens_fn))

plt.plot(ene_fn,intens_fn,label='Fast Neutron Spectrum')
# plt.plot(ene_tn,intens_tn,label='Thermal Neutron Spectrum')
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Energy/eV")
plt.ylabel("Intensity")
plt.xlim()
# plt.ylim()
plt.legend()
plt.savefig(plot_path+"Cf_filter.png")
