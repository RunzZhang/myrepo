import os

import matplotlib.pyplot as plt
import numpy as np

import pandas as pd

plot_path = '/data/runzezhang/result/TN_sims_D/plot/'
df_list = ["Cold-Cs-11_17-18_exposures","Cold-Cs-12_01_exposures","Cold-Cs-12_10-11_exposures","60Co-12_15-16_exposures"]
fmt_list = ['o','-s','--^',':d']
colors = [
    ['#8B0000', '#FF0000', '#FF7F7F'],         #'red'
     ['#006400', '#228B22', '#90EE90'],        #'green'
     ['#00008B', '#0000FF', '#ADD8E6'],        #'blue'
     ['#3E2723', '#795548', '#D7CCC8']]        #'brown'
fig, ax = plt.subplots(1,3)
# 3 graphs, one for original plot(Gray) about rate, 2 and 3 for rejections
for i in range(len(df_list)):
    path = os.path.join(plot_path, df_list[i] + "_output.txt")
    df = pd.read_csv(path,index_col=0)
    doc_label = df_list[i].rstrip("_exposures")
    # signal
    
    ax[0].errorbar(df['Updated Setiz [keV]'],df["Exp Rate [mHz]"],
                   yerr = df["Exp Sigma [mHz]"],label=doc_label+"signal",fmt=fmt_list[i],color = colors[i][0])
    ax[0].errorbar(df[ 'Updated Setiz [keV]'], df[ "Bkg Rate [mHz]"],
                   yerr=df[ "Bkg Sigma [mHz]"], label=doc_label + "bkg",fmt=fmt_list[i],color = colors[i][1])
    ax[0].errorbar(df[ 'Updated Setiz [keV]'], df[ "Clean Rate [mHz]"],
                   yerr=df[ "Clean Sigma [mHz]"], label=doc_label + "clean", fmt=fmt_list[i],color = colors[i][2])

    ax[1].errorbar(df[ 'Updated Setiz [keV]'], df[ "Rejection Rate Scattering[mHz]"],
                   yerr=df[ "Rejection Sigma Scattering[mHz]"], label=doc_label,fmt=fmt_list[i])

    ax[2].errorbar(df[ 'Updated Setiz [keV]'], df[ "Rejection Rate KeV[mHz]"],
                   yerr=df[ "Rejection Sigma KeV[mHz]"], label=doc_label,fmt=fmt_list[i])

ax[0].set_xlabel("Setiz [keV]")
ax[0].set_ylabel("Rate [mHz]")
ax[0].set_title("Signal/BKG Rates ")

ax[1].set_xlabel("Setiz [keV]")
ax[1].set_ylabel("Gamma Rejection Per Scattering []")
ax[1].set_title("Gamma Rejection Per Scattering ")


ax[1].set_xlabel("Setiz [keV]")
ax[1].set_ylabel("Gamma Rejection Per keV [/keV]")
ax[1].set_title("Gamma Rejection Per keV ")

plt.legend()
plt.savefig(plot_path + "gamma_rejection.pdf")