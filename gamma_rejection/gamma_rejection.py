import os

import matplotlib.pyplot as plt
import numpy as np

import pandas as pd

plot_path = '/data/runzezhang/result/TN_sims_D/plot/'
df_co_list = ["60Co-12_15-16_exposures"]
# df_cs_list = ["Cold-Cs-11_17-18_exposures_mix","Cold-Cs-12_01_exposures_mix","Cold-Cs-12_10-11_exposures_mix"]
df_cs_list = ["Cold-Cs-1_20-21_exposures",
                         "Cold-Cs-2_2-3_exposures"]
df_full_list = df_co_list+df_cs_list
# df_list = df_co_list
# df_list = df_cs_list
df_list = df_full_list
fmt_list = ['o','s','^','d']
colors = [
    ['#8B0000', '#FF0000', '#FF7F7F'],         #'red'
     ['#006400', '#228B22', '#90EE90'],        #'green'
     ['#00008B', '#0000FF', '#ADD8E6'],        #'blue'
     ['#3E2723', '#795548', '#D7CCC8']]        #'brown'
fig, ax = plt.subplots(1,3, figsize=(22, 4))
# 3 graphs, one for original plot(Gray) about rate, 2 and 3 for rejections

for i in range(len(df_list)):
    path = os.path.join(plot_path, df_list[i] + "_output.txt")
    df = pd.read_csv(path)
    # print(df.columns)
    doc_label = df_list[i].rstrip("_exposures")
    # signal
    # drop 2.75,3.25, 3.75 bara pressure
    pressure_drop_list = [2.75,3.25,3.75]
    df = df[~df['Pressure [bara]'].isin(pressure_drop_list )]
    # for j in pressure_drop_list:
    #     df.drop(df[df['Pressure [bara]'] == j].index)
    print(df['Pressure [bara]'])

    # ax[0].errorbar(df['Pressure [bara]'],df["Exp Rate [mHz]"],
    #                yerr = df["Exp Sigma [mHz]"],label=doc_label+"signal",fmt=fmt_list[i],color = colors[i][0])
    # ax[0].errorbar(df['Pressure [bara]'], df[ "Bkg Rate [mHz]"],
    #                yerr=df[ "Bkg Sigma [mHz]"], label=doc_label + "bkg",fmt=fmt_list[i],color = colors[i][1])
    ax[0].errorbar(df[ 'Pressure [bara]'], df[ "Clean Rate [mHz]"],
                   yerr=df[ "Clean Sigma [mHz]"], label=doc_label + "clean", fmt=fmt_list[i],color = colors[i][2])

    ax[1].errorbar(df[ 'Updated Setiz [keV]'], df[ "Rejection Rate Scattering[]"],
                   yerr=df[ "Rejection Sigma Scattering[]"], label=doc_label,fmt=fmt_list[i])

    ax[2].errorbar(df[ 'Eion_rl-1_rhol-1 [10GeVcm**2 g-1]'], df[ "Rejection Rate KeV[/keV]"],
                   yerr=df[ "Rejection Sigma KeV[/keV]"], label=doc_label,fmt=fmt_list[i])

    print("count", doc_label, df[ "Rejection Rate Scattering[]"])
    print("count*energy", doc_label, df[ "Rejection Rate KeV[/keV]"])
    # ax[2].errorbar(df['Eion_rl-1_rhol-1 [10GeVcm**2 g-1]'], df["Rejection Rate Scattering[mHz]"],
    #                yerr=df["Rejection Sigma Scattering[mHz]"], label=doc_label, fmt=fmt_list[i])



print()
ax[0].set_xlabel("Pressure [bara]")
ax[0].set_ylabel("Rate [mHz]")
ax[0].set_title("Signal/BKG Rates ")
ax[0].legend()

ax[1].set_xlabel("Setiz [keV]")
ax[1].set_ylabel("Gamma Rejection Per Scattering []")
ax[1].set_title("Gamma Rejection Per Scattering ")
ax[1].set_ylim(1.0e-12,1.0e-2)
ax[1].set_xlim(0,6)
ax[1].set_yscale("log")

ax[1].legend()


ax[2].set_xlabel("Eion_rl-1_rhol-1 [10GeVcm**2 g-1]")
ax[2].set_ylabel("Gamma Rejection Per keV [/keV]")
ax[2].set_title("Gamma Rejection Per keV ")
ax[2].set_ylim(1.0e-14,1.0e-4)
ax[2].set_xlim(0.08,0.15)
ax[2].set_yscale("log")
ax[2].legend()


plt.savefig(plot_path + "gamma_rejection_2026_Cs.pdf")