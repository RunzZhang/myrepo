import os

import matplotlib.pyplot as plt
import numpy as np

import pandas as pd

plot_path = '/data/runzezhang/result/TN_sims_D/plot/'
df_co_list = ["60Co-12_15-16_exposures_zoom"]
df_cs_list = ["Cold-Cs-11_17-18_exposures_mix","Cold-Cs-12_01_exposures_mix","Cold-Cs-12_10-11_exposures_mix"]
df_full_list = df_co_list+df_cs_list
fmt_list = ['o','s','^','d']
colors = [
    ['#8B0000', '#FF0000', '#FF7F7F'],         #'red'
     ['#006400', '#228B22', '#90EE90'],        #'green'
     ['#00008B', '#0000FF', '#ADD8E6'],        #'blue'
     ['#3E2723', '#795548', '#D7CCC8']]        #'brown'
fig, ax = plt.subplots()
# 3 graphs, one for original plot(Gray) about rate, 2 and 3 for rejections
df_list = df_co_list
# df_list = df_cs_list
# df_list = df_full_list
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

    ax.errorbar(df['Updated Setiz [keV]'], df["Rejection Rate Scattering[mHz]"],
                   yerr=df["Rejection Sigma Scattering[mHz]"], label=doc_label+ "per scattering",fmt=fmt_list[i])
    ax2 = ax.twinx()
    ax.errorbar(df['Updated Setiz [keV]'], df["Rejection Rate KeV[mHz]"],
                yerr=df["Rejection Sigma KeV[mHz]"], label=doc_label+"per ER",fmt=fmt_list[i])
    # ax.errorbar(df['Pressure [bara]'], df[ "Bkg Rate [mHz]"],
    #                yerr=df[ "Bkg Sigma [mHz]"], label=doc_label + "bkg",fmt=fmt_list[i],color = colors[i][1])
    # ax.errorbar(df[ 'Pressure [bara]'], df[ "Clean Rate [mHz]"],
    #                yerr=df[ "Clean Sigma [mHz]"], label=doc_label + "clean", fmt=fmt_list[i],color = colors[i][2])


ax.set_xlabel("Setiz [keV]")
ax.set_ylabel("Ratio [1]")
ax.set_title("Rejection compare")

ax2.set_ylabel("Ratio [1/keV]")
ax.legend()



plt.savefig(plot_path + "gamma_rejection_Co.pdf")