import pandas as pd
import matplotlib.pyplot as plt
import csv
import numpy as np
import re

after_KE_path = "./Cf252_spontanous.txt"
before_KE_path = "/data/runzezhang/Geant4Simulaions/g411_TN/CF252Sap7.5Poly6Pad10KE.dat"
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

    with open(path, "r") as f:
        for line in f:
            if re.match(r'\s*\d+\.\d{6}[+-]\d\s+\d+\.\d{6}[+-]\d', line):
                # Each line has 3 (energy, intensity) pairs
                numbers = [float(line[i:i + 11]) * 10 ** float(line[i + 11:i + 13])
                           for i in range(0, 66, 22)]
                for i in range(0, len(numbers), 2):
                    energy.append(numbers[i])
                    intensity.append(numbers[i + 1])

    # Convert to NumPy arrays
    energy = np.array(energy)
    intensity = np.array(intensity)

    # Optional: remove trailing zeros or normalize
    nonzero = intensity > 0
    energy = energy[nonzero]
    intensity = intensity[nonzero]

    # Normalize if needed
    intensity /= np.trapz(intensity, energy)
    return energy,intensity
#fast neutron fs thermal neutron tn
(ene_fn, intens_fn)=dat_to_list(before_KE_path)
(ene_tn,intens_tn) = raw_to_list(after_KE_path)

plt.plot(ene_fn,intens_fn,label='Fast Neutron Spectrum')
plt.plot(ene_tn,intens_tn,label='Thermal Neutron Spectrum')
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Energy/eV")
plt.ylabel("Intensity")
plt.xlim()
# plt.ylim()
plt.legend()
plt.savefig(plot_path+"Cf_filter.png")
