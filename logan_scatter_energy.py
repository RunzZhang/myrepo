import matplotlib.pyplot as plt
import numpy as np
import os, pickle
import sys
import pandas as pd
import uproot,csv
import math
#from matplotlib.ticker import (LogLocator, MultipleLocator, AutoMinorLocator)
import matplotlib.ticker as ticker
filename = "/data/runzezhang/result/CF252BigEl6,10.root" #EdepSap5cmx7.5cm3cmPad10cm25cmAmLi.root, BigEl.root
file = uproot.open(filename)["tree1"]
print(file.keys())
df = file.arrays(["EDep", "EventNo"], library="pd")
EDep = df['EDep'].tolist()
EDep = [i*10**6 for i in EDep]
#print(EDep)
EventNo = df['EventNo'].tolist()
threshold = 1000
highethreshold = 16700
#print(EventNo)
OverThresholdScatters = []
HighEnergyScatters = []
ElasticTot = []
n=0
for i in EDep: #Creates a 2d list with event no. and EDep only if EDep of the scatter > 1keV
    if i > threshold:
        OverThresholdScatters.append([EventNo[n],i])
    if i > highethreshold:
        HighEnergyScatters.append([EventNo[n],i])
    n+=1
#for i in range(0, len(EDep)): #Creates a 2d list with event no. and EDep only if EDep of the scatter > 1keV
#    if EDep[i] > 1000:
#        OverThresholdScatters.append([EventNo[i],i])
#print(OverThresholdScatters)
#print(HighEnergyScatters)
print(len(OverThresholdScatters))
print(len(HighEnergyScatters))
OverThresholdSingleScatters = []
#169, 192, 280, 626, 970, 1870, 1911, 1936, 2043, 2075, 2210, 2368, 2511, 2558,
#2695!!!!!, 2744!!!!!!, 2753, 3131!!!!, 3225, 3599
for i in range(0, len(OverThresholdScatters)): #Counts the amount of single over threshold scatters
    if i == 0:
        if OverThresholdScatters[i][0] != OverThresholdScatters[i+1][0]:
            OverThresholdSingleScatters.append(OverThresholdScatters[i][0])
    elif i == len(OverThresholdScatters)-1:
        if OverThresholdScatters[i][0] != OverThresholdScatters[i-1][0]:
            OverThresholdSingleScatters.append(OverThresholdScatters[i][0])
    else:
        if (OverThresholdScatters[i][0] != OverThresholdScatters[i-1][0]) and (OverThresholdScatters[i][0] != OverThresholdScatters[i+1][0]):
            OverThresholdSingleScatters.append(OverThresholdScatters[i][0])
HighEnergySingleScatters = []
for i in range(0, len(HighEnergyScatters)): #Counts the amount of single high energy (>10keV) scatters
    if i == 0:
        if HighEnergyScatters[i][0] != HighEnergyScatters[i+1][0]:
            HighEnergySingleScatters.append([HighEnergyScatters[i][0], HighEnergyScatters[i][1]])
    elif i == len(HighEnergyScatters)-1:
        if HighEnergyScatters[i][0] != HighEnergyScatters[i-1][0]:
            HighEnergySingleScatters.append([HighEnergyScatters[i][0], HighEnergyScatters[i][1]])
    else:
        if (HighEnergyScatters[i][0] != HighEnergyScatters[i-1][0]) and (HighEnergyScatters[i][0] != HighEnergyScatters[i+1][0]):
            HighEnergySingleScatters.append([HighEnergyScatters[i][0], HighEnergyScatters[i][1]])
print(len(OverThresholdSingleScatters))
print(len(HighEnergySingleScatters))
first20over = []
first20high = []
first20 = []
for i in range(0,20):
    first20over.append(OverThresholdSingleScatters[i])
    first20high.append(HighEnergySingleScatters[i])
    first20.append(OverThresholdScatters[i])
print(first20over)
print(first20high)
print(first20)
# Counts the number of multiple over threshold scatters
OverThresholdScattersEventNoNoDuplicates = []
for i in OverThresholdScatters:
    OverThresholdScattersEventNoNoDuplicates.append(i[0])
from collections import Counter
TotalOverThresholdEvents = 0
c = Counter(OverThresholdScattersEventNoNoDuplicates)
for i in c:
    TotalOverThresholdEvents += 1
OverThresholdScattersEventNoNoDuplicates = list(set(OverThresholdScattersEventNoNoDuplicates))
MultipleOverThresholdScattersCount = len(OverThresholdScattersEventNoNoDuplicates)-len(OverThresholdSingleScatters)
#print(OverThresholdSingleScatters)
#print(EventNo)
print(len(EventNo))
print(len(EDep))
print(len(OverThresholdScatters))
#print(OverThresholdScatters)
#print(OverThresholdSingleScatters)
#print(EventNo)
#print("# of captures is", len(NCaptureTot))
print("# of scatterers over the 1keV threshold is", len(OverThresholdScatters))
print("# of events with only one 1keV scatter is", len(OverThresholdSingleScatters))
print('# of events with more than one 1keV scatter is', MultipleOverThresholdScattersCount)
print("# of single scatters >", highethreshold/1000, "keV (these can produce scintillation light) is", len(HighEnergySingleScatters))
print(TotalOverThresholdEvents)
bottombin = 1e-4
binnumke = 100
bins = np.logspace(np.log10(bottombin),np.log10(max(EDep)), binnumke)
bins1 = np.linspace(0, 1000000, 1000)
fig1, ax1 = plt.subplots()
#ax1.hist(KEe, bins1, histtype = "step", label = "Escape")
ax1.hist(EDep, bins, histtype = "step", label = "Escape")
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel('AmLi Neutron Deposition Energy (eV)', fontsize = 14)
ax1.set_ylabel('Number of Events', fontsize = 14)
ax1.set_title('Energy Depositions in LAr from Elastic Scatterers', fontsize = 16)
#ax1.xaxis.set_ticks(np.logspace(np.log10(bottombin), 1, num=1-int(np.log10(bottombin))+1))
ax1.xaxis.set_major_locator(ticker.LogLocator(numticks=999))
ax1.xaxis.set_minor_locator(ticker.LogLocator(numticks=999, subs="auto"))
plt.axvline(x = 1000, color = 'red', linestyle = '--', alpha = 0.5, label = "Elastic Max")
plt.show()
#Plots the >16.7keV Recoils
bottombinOver = 16700
binnumkeOver = 100
binsOver = np.logspace(np.log10(bottombinOver),np.log10(max(EDep)), binnumkeOver)
#binsOver = np.linspace(bottombinOver, max(EDep), binnumkeOver)
SingleScatterList = [HighEnergySingleScatters[i][1] for i in range(0, len(HighEnergySingleScatters))]
print(min(SingleScatterList))
print(max(SingleScatterList))
fig2, ax2 = plt.subplots()
#ax1.hist(KEe, bins1, histtype = "step", label = "Escape")
ax2.hist(SingleScatterList, binsOver, histtype = "step", label = "Single Scatters >16.7keV")
with open("/data/runzezhang/scatter_spectrum.csv", 'w', newline='') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(SingleScatterList)
SingleScatterList.to_csv("/data/runzezhang/scatter_spectrum.csv",index=False)
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel('Deposition Energy (eV)', fontsize = 16)
ax2.set_ylabel('Number of Events', fontsize = 16)
ax2.set_title('Distribution of Single Scatter Energies >16.7keV from CF252 Thermal Neutron Source', fontsize = 18)
#ax1.xaxis.set_ticks(np.logspace(np.log10(bottombin), 1, num=1-int(np.log10(bottombin))+1))
ax2.xaxis.set_major_locator(ticker.LogLocator(numticks=999))
ax2.xaxis.set_minor_locator(ticker.LogLocator(numticks=999, subs="auto"))
plt.axvline(x = 16700, color = 'red', linestyle = '--', alpha = 0.5, label = "16.7keV")
plt.legend()
plt.show()