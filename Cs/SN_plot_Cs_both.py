import pandas as pd
import matplotlib.pyplot as plt
import csv
import numpy as np
import os
class SN():
    def __init__(self):
        # after generate new files, you need to select the capture ratio/source for different configs in read_files function.
        # then choose the correct signal/noise of with clause in read files.
        # at last change the self.name and plot_name in plot function
        # v2: change back to 2 backgrounds but with finer definitions
        # v4 kill duplicated NRERs
        self.base_path = "/data/runzezhang/result/TN_sims_D/chunked_root_files_Cs_1E7/"
        self.base_path2 = "/data/runzezhang/result/TN_sims_D/chunked_root_files_Cs_1E7/" # for gamma path

        # self.base_path = "/data/runzezhang/result/TN_sims_D/chunked_root_files_pn_outside_1E7/"
        # self.base_path2 = "/data/runzezhang/result/TN_sims_D/chunked_root_files_pn_1E7_outside_gamma/"  # without lead

        self.plot_path = '/data/runzezhang/result/TN_sims_D/plot/'
        self.false_1 = "PN_false1.csv"
        self.false_2 = "PN_false2.csv"
        self.signal = "PN_sig.csv"
        self.name1 = "Correlated ER Background"
        self.name2 = "Hard Scatter Background"
        self.name = "Backgrounds"
        self.plot_name = self.name+"PN_1E6_wt_gamma_outside.pdf"
        self.gamma_bkg_path = self.plot_path+"Cs_gamma position_distribution.pdf"
        self.pho_threshold = 100

        cols = ["Event","name", "R/mm", "Z/mm", "Volume", "Process", "ER_near/eV", "Multiplicity"]

        self.df_list =[]



        #982 statics false 1
        for i in range(1,51):
        # for i in range(1, 11):
            self.main_body(i)
        # self.main_body(1)
        self.combine_df()
        self.data_analysis()


    def main_body(self,i):
        print(i)
        self.false_1 = f"PN_1E7_false1_part{i}.csv"
        self.false_2 = f"PN_1E7_false2_part{i}.csv"
        self.false_3 = f"PN_1E7_false3_part{i}.csv"
        self.false_gamma_1 = f"PN_gamma_1E7_false1_part{i}.csv"
        self.signal = f"PN_1E7_sig_part{i}.csv"

        self.info_path = self.base_path + f"Cs_gamma_1E6_info_scube_part{i}.csv"

        self.false_1_path = self.base_path + self.false_1
        self.false_2_path = self.base_path + self.false_2
        self.false_3_path = self.base_path + self.false_3
        self.false_gamma_1_path = self.base_path2 + self.false_gamma_1
        self.signal_path = self.base_path + self.signal

        self.ini_path = self.base_path+ f"PN_1E7_ini_part{i}.csv"
        self.ar_ke_path = self.base_path + f"PN_1E7_ke_part{i}.csv"
        self.ar_ke_alter_path = self.base_path + f"PN_1E7_ke_part{i}.csv"


        self.read_files()


        # self.read_files_s_to_N1()

# main funtion we use
    def read_files(self):
        self.original_Activity = 5 # original activity in the paper
        self.Activity = 5  # source practical activity in mivro curie for 50 bubbles/hour

        self.gamma_rate = 2.44e6 # /s


        self.G4_events_gamma =  5E6 # only 50 chunks
        self.ambient_bubble = 5 # /h


        temp_df = pd.read_csv(self.info_path)

        self.df_list.append(temp_df)


    def combine_df(self):
        self.merged_df = pd.concat(self.df_list, ignore_index=True)

    def data_analysis(self):
        #position distributions histogram, dependisng on step number
        # self.read_positions()
        # self.read_positions_2d_hist()
        # self.read_positions_zslice()
        #mulitipliciy distribtuion depending on events
        # self.read_multiplicity()
        # self.read_Ar_multiplicity()
        # ER distribution per row
        # self.read_ER_Ar_CF()
        # self.read_ER_Ar_CF_per_deposit_rate()
        # self.read_ER_Ar_CF_per_deposit_rate_cumulative()
        # self.read_ER_CF_per_deposit_rate_cumulative()
        # self.read_ER_Ar_CF_1d_sum()
        # self.read_ER_Ar_CF_2d_sum()
        # self.read_ER_Ar_CF_1d_sum_rate()
        # self.read_ER_Ar_CF_1d_sum_rate_cummulative()

        self.gamma_rejection_rate_per_keV_vs_Setiz()



    def read_positions(self):
        max_multi_num = self.merged_df["Multiplicity"].max()
        print(max_multi_num)
        self.mutiplicity_list =[]
        for i in range(1,max_multi_num+1,1):
            temp_df = self.merged_df[self.merged_df["Multiplicity"]==i][["R/mm","Z/mm"]]
            temp_list = temp_df.to_numpy()
            self.mutiplicity_list.append(temp_list)

        # self.cmap =  plt.cm.plasma
        fig, ax = plt.subplots()

        sc=ax.scatter(self.merged_df["R/mm"],self.merged_df["Z/mm"],
        c=self.merged_df["Multiplicity"],   # color comes from data
        cmap="plasma",
        s=5,
        alpha=0.7)

        ax.plot([189.95,189.95, 0.8485], [0,663.22, 714.03], color="red")
        ax.plot([114.98,114.98, 0.75575], [0,587.01, 617.78], color="blue")
        ax.plot([99.01,99.01, 4.34], [0,366.49, 399.82], color="red")

        ax.set_xlabel("R [mm]")
        ax.set_ylabel("Z [mm]")
        cbar = plt.colorbar(sc, ax=ax)
        cbar.set_label("Scatter Order")
        plt.savefig(self.plot_path+"Cs_1E5_position.pdf")

    def read_positions_zslice(self):
        z_slice = self.merged_df[(self.merged_df["Z/mm"]>=450)&(self.merged_df["Z/mm"]<=550)]

        fig, ax = plt.subplots()
        # self.cmap =  plt.cm.plasma
        sc = ax.hist2d(z_slice["X/mm"], z_slice["Y/mm"], bins=50,
                       cmap="plasma", norm="log", alpha=0.7)



        ax.set_xlabel("X [mm]")
        ax.set_ylabel("Y [mm]")

        cbar = plt.colorbar(sc[3], ax=ax)
        cbar.set_label("Counts(log)")
        plt.savefig(self.plot_path + "Cs_1E5_zslice_position_density.pdf")

    def read_positions_2d_hist(self):
        max_multi_num = self.merged_df["Multiplicity"].max()
        print(max_multi_num)
        self.mutiplicity_list =[]
        for i in range(1,max_multi_num+1,1):
            temp_df = self.merged_df[self.merged_df["Multiplicity"]==i][["R/mm","Z/mm"]]
            temp_list = temp_df.to_numpy()
            self.mutiplicity_list.append(temp_list)

        # self.cmap =  plt.cm.plasma
        fig, ax = plt.subplots()

        sc=ax.hist2d(self.merged_df["R/mm"],self.merged_df["Z/mm"],bins=50,
        cmap="plasma",norm="log",alpha=0.7)

        ax.plot([189.95,189.95, 0.8485], [0,663.22, 714.03], color="red")
        ax.plot([114.98,114.98, 0.75575], [0,587.01, 617.78], color="blue")
        ax.plot([99.01,99.01, 4.34], [0,366.49, 399.82], color="red")

        ax.set_xlabel("R [mm]")
        ax.set_ylabel("Z [mm]")
        ax.set_xlim(0,200)
        ax.set_ylim(-100,800)
        cbar = plt.colorbar(sc[3], ax=ax)
        cbar.set_label("Counts(log)")
        plt.savefig(self.plot_path+"Cs_1E5_position_density.pdf")


    def read_multiplicity(self):
        multiplicity= self.merged_df.groupby("Event")["Multiplicity"].max()
        max_m = self.merged_df["Multiplicity"].max()
        print(max_m)
        bins = np.arange(1, max_m + 2)
        fig,ax = plt.subplots()
        ax.hist(multiplicity, bins= bins, align="left", rwidth=0.9, density=True)

        ax.set_xlabel("Multiplicity")
        ax.set_ylabel("Probability")
        plt.savefig(self.plot_path+"Cs_1E7_multi.pdf")

    def read_Ar_multiplicity(self):
        multiplicity = self.merged_df[self.merged_df["Volume"]=="LAr_phys"].groupby("Event")["Multiplicity"].max()
        max_m = self.merged_df["Multiplicity"].max()
        print(max_m)
        bins = np.arange(1, max_m + 2)
        fig, ax = plt.subplots()
        ax.hist(multiplicity, bins=bins, align="left", rwidth=0.9, density=True)

        ax.set_xlabel("Argon Multiplicity")
        ax.set_ylabel("Probability")
        plt.savefig(self.plot_path + "Cs_1E7_argon_multi.pdf")
    def read_ER_Ar_CF(self):
        # per energy deposition
        # the sum is per event
        ER_Ar = self.merged_df[self.merged_df["Volume"]=="LAr_phys"]["ER_near/eV"]/1000
        ER_CF4 = self.merged_df[self.merged_df["Volume"] == "hydraulic_fluid_phys"]["ER_near/eV"] / 1000
        ER_sum = self.merged_df.groupby("Event")["ER_near/eV"].sum()/1000

        ER_huge = self.merged_df.groupby("Event")["ER_near/eV"].sum()
        events_keep = ER_huge[ER_huge>9e5].index
        filtered_df = self.merged_df[self.merged_df["Event"].isin(events_keep)]

        print(filtered_df.head(20))

        fig, ax = plt.subplots(1,3, figsize=(14, 4))
        ax[0].hist(ER_Ar, bins=40,align="left")
        ax[0].set_xlabel("ER/keV per scattering LAr")
        ax[0].set_ylabel("Counts")

        ax[1].hist(ER_CF4, bins=40,align="left")
        ax[1].set_xlabel("ER/keV per scattering CF4")
        ax[1].set_ylabel("Counts")

        ax[2].hist(ER_sum, bins=40, align="left")
        ax[2].set_xlabel("ER/keV per event")
        ax[2].set_ylabel("Counts")

        plt.savefig(self.plot_path + "Cs_1E5_ER_all.pdf")

    def read_ER_Ar_CF_per_deposit_rate(self):
        # per energy deposition and total
        # MHz
        Rate_factor = self.gamma_rate*1000 / (self.G4_events_gamma)
        ER_Ar = self.merged_df[self.merged_df["Volume"]=="LAr_phys"]["ER_near/eV"]/1000
        ER_CF4 = self.merged_df[self.merged_df["Volume"] == "hydraulic_fluid_phys"]["ER_near/eV"] / 1000
        print("ER CF4 counts", len(ER_CF4))
        ER_sum = self.merged_df["ER_near/eV"]/1000

        hist_array = [None] * 3

        hist_array[0] = np.histogram(ER_Ar, bins=100,range=(0, 660))
        hist_array[1] = np.histogram(ER_CF4, bins=100,range=(0, 660))
        hist_array[2] = np.histogram(ER_sum, bins=100,range=(0, 660))

        # find if compton edge exist in CF4 cumulative spectrum
        for j in range(len(hist_array[1][0])):
            if hist_array[1][1][j] > 500:
                print(hist_array[1][1][j],"keV edge", hist_array[1][0][j])


        fig, ax = plt.subplots(1, 3, figsize=(16, 4))
        ax[0].bar(hist_array[0][1][:-1], Rate_factor * hist_array[0][0],width=np.diff(hist_array[0][1]),
        align="edge",
        edgecolor="black")
        bin0_len = int(hist_array[0][1][1] - hist_array[0][1][0])
        ax[0].set_xlabel("ER/keV per deposition in LAr")
        ax[0].set_ylabel(" Rate mHz/(bin[" + str(bin0_len) + " keV])")
        # ax[0].ticklabel_format(axis="y",style="sci", scilimits=(0, 0) )
        ax[0].set_yscale("log")
        ax[0].minorticks_on()
        # ax[0].grid(which="major", linestyle="-", linewidth=0.8, alpha=0.7)
        # ax[0].grid(which="minor", linestyle=":", linewidth=0.5, alpha=0.4)

        ax[1].bar(hist_array[1][1][:-1], Rate_factor * hist_array[1][0],width=np.diff(hist_array[1][1]),
        align="edge",
        edgecolor="black")
        ax[1].set_xlabel("ER/keV per deposition in CF4")
        bin1_len = int(hist_array[1][1][1] - hist_array[1][1][0])
        ax[1].set_ylabel("Rate mHz/([" + str(bin1_len) + " keV])")
        ax[1].set_yscale("log")
        # ax[1].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        # ax[1].grid(True)
        ax[1].minorticks_on()
        # ax[1].grid(which="major", linestyle="-", linewidth=0.8, alpha=0.7)
        # ax[1].grid(which="minor", linestyle=":", linewidth=0.5, alpha=0.4)

        ax[2].bar(hist_array[2][1][:-1], Rate_factor * hist_array[2][0],width=np.diff(hist_array[2][1]),
        align="edge",
        edgecolor="black")
        ax[2].set_xlabel("ER/keV per deposition ")
        bin2_len = int(hist_array[2][1][1] - hist_array[2][1][0])
        ax[2].set_ylabel("Rate mHz/([" + str(bin2_len) + " keV])")
        # ax[2].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        ax[2].set_yscale("log")
        # ax[2].grid(True)
        ax[2].minorticks_on()
        # ax[2].grid(which="major", linestyle="-", linewidth=0.8, alpha=0.7)
        # ax[2].grid(which="minor", linestyle=":", linewidth=0.5, alpha=0.4)

        plt.savefig(self.plot_path + "Cs_1E7_ER_perdepostion_coldrate.pdf")


    def read_ER_Ar_CF_per_deposit_rate_cumulative(self):
        # rate factor in mHz
        Rate_factor = self.gamma_rate*1000 / (self.G4_events_gamma)
        ER_Ar = self.merged_df[self.merged_df["Volume"] == "LAr_phys"]["ER_near/eV"] / 1000
        ER_CF4 = self.merged_df[self.merged_df["Volume"] == "hydraulic_fluid_phys"]["ER_near/eV"] / 1000
        ER_sum = self.merged_df["ER_near/eV"] / 1000

        hist_array = [None] * 3

        hist_array[0] = np.histogram(ER_Ar, bins=100,range=(0, 660))
        hist_array[1] = np.histogram(ER_CF4, bins=100,range=(0, 660))
        hist_array[2] = np.histogram(ER_sum, bins=100,range=(0, 660))
        cumulative_threshold_array = [None]*3

        cumulative_threshold_array[0] = np.array([sum(hist_array[0][0][i:]) for i in range(len(hist_array[0][0]))])
        cumulative_threshold_array[1] = np.array([sum(hist_array[1][0][i:]) for i in range(len(hist_array[1][0]))])
        cumulative_threshold_array[2] = np.array([sum(hist_array[2][0][i:]) for i in range(len(hist_array[2][0]))])
        # find first 2 bins and rate for argon
        print("argon bin", hist_array[0][1][:4])
        print("argon rate", Rate_factor * cumulative_threshold_array[0][:3])
        # find if compton edge exist in LAr cumulative spectrum
        for j in range(len(cumulative_threshold_array[0])):
            if hist_array[0][1][j]>500:
                print("500 keV edge Ar",cumulative_threshold_array[0][j])
                break
        #CF4
        for j in range(len(cumulative_threshold_array[1])):
            if hist_array[1][1][j]>500:
                print("500 keV edge CF4",cumulative_threshold_array[1][j])
                break
        fig, ax = plt.subplots(1, 3, figsize=(16, 4))
        ax[0].bar(hist_array[0][1][:-1], Rate_factor * cumulative_threshold_array[0], width=np.diff(hist_array[0][1]),
                  align="edge",
                  edgecolor="black")
        bin0_len = int(hist_array[0][1][1] - hist_array[0][1][0])
        ax[0].set_xlabel("ER/keV threshold per deposition in LAr")
        ax[0].set_ylabel(" Rate mHz")
        # ax[0].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        ax[0].set_yscale("log")
        ax[0].minorticks_on()
        # ax[0].grid(which="major", linestyle="-", linewidth=0.8, alpha=0.7)
        # ax[0].grid(which="minor", linestyle=":", linewidth=0.5, alpha=0.4)

        ax[1].bar(hist_array[1][1][:-1], Rate_factor * cumulative_threshold_array[1], width=np.diff(hist_array[1][1]),
                  align="edge",
                  edgecolor="black")
        ax[1].set_xlabel("ER/keV threshold per deposition in CF4")
        bin1_len = int(hist_array[1][1][1] - hist_array[1][1][0])
        ax[1].set_ylabel("Rate mHz")
        # ax[1].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        # ax[1].grid(True)
        ax[1].set_yscale("log")
        ax[1].minorticks_on()
        # ax[1].grid(which="major", linestyle="-", linewidth=0.8, alpha=0.7)
        # ax[1].grid(which="minor", linestyle=":", linewidth=0.5, alpha=0.4)

        ax[2].bar(hist_array[2][1][:-1], Rate_factor * cumulative_threshold_array[2], width=np.diff(hist_array[2][1]),
                  align="edge",
                  edgecolor="black")
        ax[2].set_xlabel("ER/keV threshold per deposition ")
        bin2_len = int(hist_array[2][1][1] - hist_array[2][1][0])
        ax[2].set_ylabel("Rate mHz")
        # ax[2].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        ax[2].set_yscale("log")
        # ax[2].grid(True)
        ax[2].minorticks_on()
        # ax[2].grid(which="major", linestyle="-", linewidth=0.8, alpha=0.7)
        # ax[2].grid(which="minor", linestyle=":", linewidth=0.5, alpha=0.4)

        plt.savefig(self.plot_path + "Cs_1E7_ER_perdepostion_cumulative_coldrate.pdf")
    def read_ER_CF_per_deposit_rate_cumulative(self):
        # rate factor in mHz
        Rate_factor = self.gamma_rate*1000 / (self.G4_events_gamma)
        ER_Ar = self.merged_df[self.merged_df["Volume"] == "LAr_phys"]["ER_near/eV"] / 1000
        ER_CF4 = self.merged_df[self.merged_df["Volume"] == "hydraulic_fluid_phys"]["ER_near/eV"] / 1000
        ER_sum = self.merged_df["ER_near/eV"] / 1000

        hist_array = [None]

        hist_array[0] = np.histogram(ER_CF4, bins=100,range=(0, 660))



        cumulative_threshold_array = [None]

        cumulative_threshold_array[0] = np.array([sum(hist_array[0][0][i:]) for i in range(len(hist_array[0][0]))])
        # find if compton edge exist in CF4 cumulative spectrum
        for j in range(len(cumulative_threshold_array[0])):
            if hist_array[0][1][j] > 500:
                print("500 keV edge", cumulative_threshold_array[0][j])
                break

        # find what bin first reach 150mHz
        for i in range(len(cumulative_threshold_array[0])):
            if cumulative_threshold_array[0][i]*Rate_factor<=200:
                print(cumulative_threshold_array[0][i],i,"is the threshold")
                break
        fig, ax = plt.subplots(1, 1, figsize=(5, 4))
        ax.bar(hist_array[0][1][i-2:-1], Rate_factor * cumulative_threshold_array[0][i-2:], width=np.diff(hist_array[0][1][i-2:]),
                  align="edge",
                  edgecolor="black")
        bin0_len = int(hist_array[0][1][1] - hist_array[0][1][0])
        ax.set_xlabel("ER/keV threshold per deposition in CF4")
        ax.set_ylabel(" Rate mHz")
        ax.set_yscale("log")
        # ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        ax.minorticks_on()
        # ax[0].grid(which="major", linestyle="-", linewidth=0.8, alpha=0.7)
        # ax[0].grid(which="minor", linestyle=":", linewidth=0.5, alpha=0.4)


        plt.savefig(self.plot_path + "Cs_1E5_CF4_ER_perdepostion_cumulative_coldrate.pdf")
    def read_ER_Ar_CF_1d_sum_rate(self):
        # calcualte sum of ER classified in Ar and CF4 per event
        # and sum rate is over both volume in Ar and CF4

        Rate_factor = self.gamma_rate / (3600*self.G4_events_gamma) # /h per geant run file
        ER_Ar_sum = self.merged_df[self.merged_df["Volume"]=="LAr_phys"].groupby("Event")["ER_near/eV"].sum()/1000
        ER_CF4_sum = self.merged_df[self.merged_df["Volume"]=="hydraulic_fluid_phys"].groupby("Event")["ER_near/eV"].sum()/1000

        ER_sum = self.merged_df.groupby("Event")["ER_near/eV"].sum() / 1000

        hist_array = [None] * 3

        hist_array[0] = np.histogram(ER_Ar_sum, bins=50)
        hist_array[1] = np.histogram(ER_CF4_sum, bins=50)
        hist_array[2] = np.histogram(ER_sum, bins=50)


        fig, ax = plt.subplots(1,3, figsize=(14, 4))
        ax[0].plot(hist_array[0][1][:-1], Rate_factor*hist_array[0][0])
        bin0_len = int(hist_array[0][1][1]-hist_array[0][1][0])
        ax[0].set_xlabel("ER/keV per event in LAr")
        ax[0].set_ylabel("Rate/(h*"+str(bin0_len)+" keV)")
        ax[0].minorticks_on()
        ax[0].grid(which="major", linestyle="-", linewidth=0.8, alpha=0.7)
        ax[0].grid(which="minor", linestyle=":", linewidth=0.5, alpha=0.4)

        ax[1].plot(hist_array[1][1][:-1], Rate_factor*hist_array[1][0])
        ax[1].set_xlabel("ER/keV per event in CF4")
        bin1_len = int(hist_array[1][1][1] - hist_array[1][1][0])
        ax[1].set_ylabel("Rate/(h*"+str(bin1_len)+" keV)")
        ax[1].grid(True)
        ax[1].minorticks_on()
        ax[1].grid(which="major", linestyle="-", linewidth=0.8, alpha=0.7)
        ax[1].grid(which="minor", linestyle=":", linewidth=0.5, alpha=0.4)

        ax[2].plot(hist_array[2][1][:-1], Rate_factor*hist_array[2][0])
        ax[2].set_xlabel("ER/keV per event ")
        bin2_len = int(hist_array[2][1][1] - hist_array[2][1][0])
        ax[2].set_ylabel("Rate/(h*"+str(bin2_len)+" keV)")
        ax[2].grid(True)
        ax[2].minorticks_on()
        ax[2].grid(which="major", linestyle="-", linewidth=0.8, alpha=0.7)
        ax[2].grid(which="minor", linestyle=":", linewidth=0.5, alpha=0.4)




        plt.savefig(self.plot_path + "Cs_1E5_ER_coldrate.pdf")

    def read_ER_Ar_CF_1d_sum_rate_cummulative(self):
        # calcualte sum of ER classified in Ar and CF4 per event
        # cumulative, event rate above NR threshold
        Rate_factor = self.gamma_rate / (3600*self.G4_events_gamma) # /h per geant run file
        ER_Ar_sum = self.merged_df[self.merged_df["Volume"]=="LAr_phys"].groupby("Event")["ER_near/eV"].sum()/1000
        ER_CF4_sum = self.merged_df[self.merged_df["Volume"]=="hydraulic_fluid_phys"].groupby("Event")["ER_near/eV"].sum()/1000

        ER_sum = self.merged_df.groupby("Event")["ER_near/eV"].sum() / 1000

        hist_array = [None] * 3

        hist_array[0] = np.histogram(ER_Ar_sum, bins=50)
        hist_array[1] = np.histogram(ER_CF4_sum, bins=50)
        hist_array[2] = np.histogram(ER_sum, bins=50)

        cumulative_threshold_array = [None]*3

        cumulative_threshold_array[0] = np.array([sum(hist_array[0][0][i:]) for i in range(len(hist_array[0][0]))])
        cumulative_threshold_array[1] = np.array([sum(hist_array[1][0][i:]) for i in range(len(hist_array[1][0]))])
        cumulative_threshold_array[2] = np.array([sum(hist_array[2][0][i:]) for i in range(len(hist_array[2][0]))])




        fig, ax = plt.subplots(1,3, figsize=(14, 4))
        ax[0].plot(hist_array[0][1][:-1], Rate_factor*cumulative_threshold_array[0])
        bin0_len = int(hist_array[0][1][1]-hist_array[0][1][0])
        ax[0].set_xlabel("ER/keV threshold in LAr")
        ax[0].set_ylabel("Rate/(h)")
        ax[0].minorticks_on()
        ax[0].grid(which="major", linestyle="-", linewidth=0.8, alpha=0.7)
        ax[0].grid(which="minor", linestyle=":", linewidth=0.5, alpha=0.4)

        ax[1].plot(hist_array[1][1][:-1], Rate_factor*cumulative_threshold_array[1])
        ax[1].set_xlabel("ER/keV threshold in CF4")
        bin1_len = int(hist_array[1][1][1] - hist_array[1][1][0])
        ax[1].set_ylabel("Rate/(h)")
        ax[1].grid(True)
        ax[1].minorticks_on()
        ax[1].grid(which="major", linestyle="-", linewidth=0.8, alpha=0.7)
        ax[1].grid(which="minor", linestyle=":", linewidth=0.5, alpha=0.4)

        ax[2].plot(hist_array[2][1][:-1], Rate_factor*cumulative_threshold_array[2])
        ax[2].set_xlabel("ER/keV threshold ")
        bin2_len = int(hist_array[2][1][1] - hist_array[2][1][0])
        ax[2].set_ylabel("Rate/(h)")
        ax[2].grid(True)
        ax[2].minorticks_on()
        ax[2].grid(which="major", linestyle="-", linewidth=0.8, alpha=0.7)
        ax[2].grid(which="minor", linestyle=":", linewidth=0.5, alpha=0.4)




        plt.savefig(self.plot_path + "Cs_1E5_ER_coldrate_cumulative.pdf")

    def read_ER_Ar_CF_2d_sum(self):
        # calcualte sum of ER classified in Ar and CF4 per event
        # 2d histogram


        ER_Ar_sum = self.merged_df[self.merged_df["Volume"]=="LAr_phys"].groupby("Event")["ER_near/eV"].sum()/1000
        ER_CF4_sum = self.merged_df[self.merged_df["Volume"]=="hydraulic_fluid_phys"].groupby("Event")["ER_near/eV"].sum()/1000

        evt = pd.concat([ER_Ar_sum, ER_CF4_sum], axis=1)
        evt.columns = ["ER_Ar_keV", "ER_CF4_keV"]  # <-- you create these names
        evt = evt.fillna(0)

        fig, ax = plt.subplots()

        sc = ax.hist2d(evt["ER_Ar_keV"], evt["ER_CF4_keV"], bins=50,
                       cmap="plasma", norm="log", alpha=0.7)

        ax.set_xlabel("ER_Ar/keV")
        ax.set_ylabel("ER_CF4/keV")
        # ax.set_xlim(0, 200)
        # ax.set_ylim(-100, 800)
        cbar = plt.colorbar(sc[3], ax=ax)
        cbar.set_label("Counts(log)")
        plt.savefig(self.plot_path + "Cs_1E5_ER_density.pdf")


    def gamma_rejection_rate_per_keV_vs_Setiz(self):
        # rate factor in mHz
        expfile_name = "Cold-Cs-11_17-18_exposures"
        bkgfile_name = "Background-11_26-30_exposures"

        Rate_factor = self.gamma_rate*1000 / (self.G4_events_gamma)
        ER_Ar = self.merged_df[self.merged_df["Volume"] == "LAr_phys"]["ER_near/eV"] / 1000 # in keV


        hist_array = [None]
        hist_array[0] = np.histogram(ER_Ar, bins=100,range=(0, 1200))

        # get probablity per scattering and the statistics
        cumulative_threshold_per_scatter_array = [None]

        cumulative_threshold_per_scatter_array[0] = np.array(
            [sum(hist_array[0][0][i:]) for i in range(len(hist_array[0][0]))])

        # histogram per scattering per keV
        cumulative_threshold_array = [None]
        energy_deposit_list = [hist_array[0][0][i]*hist_array[0][1][i] for i in range(len(hist_array[0][0]))]
        cumulative_threshold_array[0] = np.array([sum(energy_deposit_list[i:]) for i in range(len(energy_deposit_list))])

        source_exposure_df = self.read_exposure(expfile_name + ".txt")
        print(source_exposure_df.loc[:, 0])
        background_exposure_df = self.read_exposure(bkgfile_name + ".txt")

        Setiz = [1.3445287166423177, 1.4677096307281403, 1.6077252261931916, 1.7676644948295235, 2.163478457894038,
                 2.698514892785409, 3.038557782566206, 3.44261366411884]  # in keV
        source_pressure_list = source_exposure_df.loc[:, 0].to_list()
        print(source_pressure_list)
        bkg_pressure_list = background_exposure_df.loc[:, 0].to_list()
        print(bkg_pressure_list)
        exp_life_time = source_exposure_df.loc[:, 1].to_list()

        background_time = background_exposure_df.loc[:, 1].to_list()

        exp_life_time_sig = source_exposure_df.loc[:, 2].to_list()

        background_time_sig = background_exposure_df.loc[:, 2].to_list()
        bkg_pressure_recon_list = []
        exp_rate_list = []  # in mHz
        background_rate_list = []
        clean_rate_list = []
        exp_sigma_list = []
        background_sigma_list = []
        clean_sigma_list = []
        updated_Setiz_list = []
        for i in range(len(source_pressure_list)):
            src_pressure = source_pressure_list[i]
            source_bkg_pressure_match = True
            try:
                pressure_index = bkg_pressure_list.index(src_pressure)
                print(src_pressure, pressure_index)
                source_bkg_pressure_match = True
            except:
                print("source pressure is not found in background ", src_pressure)
                source_bkg_pressure_match = False
                # only get pressure entries that shows in both src and bkg
            if source_bkg_pressure_match:
                updated_Setiz_list.append(Setiz[i])
                exp_rate = 1000 / exp_life_time[i]
                background_rate = 1000 / background_time[pressure_index]
                clean_rate = exp_rate - background_rate
                exp_sigma = exp_life_time_sig[i] * 1000 / (exp_life_time[i]) ** 2
                back_sigma = background_time_sig[pressure_index] * 1000 / (background_time[pressure_index]) ** 2
                clean_sigma = np.sqrt(exp_sigma ** 2 + back_sigma ** 2)
                exp_rate_list.append(exp_rate)
                bkg_pressure_recon_list.append(src_pressure)
                background_rate_list.append(background_rate)
                clean_rate_list.append(clean_rate)
                exp_sigma_list.append(exp_sigma)
                background_sigma_list.append(back_sigma)
                clean_sigma_list.append(clean_sigma)

        print(exp_rate_list)

        rejection_PS_list = []
        rejection_PS_sigma_list = []

        rejection_PK_list = []
        rejection_PK_sigma_list = []

        # interpolation rate

        for j in range(len(updated_Setiz_list)):
            threshold = updated_Setiz_list[j]
            for i in range(len(hist_array[0][1])):
                if threshold>= hist_array[0][1][i]:
                    # rejection per scattering, PS meaning perscattering
                    counts = cumulative_threshold_per_scatter_array[0][i] + (threshold - hist_array[0][1][i]) * (
                                cumulative_threshold_per_scatter_array[0][i + 1] - cumulative_threshold_per_scatter_array[0][i]) / (
                                         hist_array[0][1][i + 1] - hist_array[0][1][i])
                    rate_PS = Rate_factor * (counts)

                    rate_PS_sigma = rate_PS / np.sqrt(counts)

                    rejection_PS = exp_rate_list[j] / rate_PS
                    rejection_PS_list.append(rejection_PS)
                    rejection_PS_sigma = np.sqrt(
                        (exp_sigma_list[i] / rate_PS) ** 2 + (exp_rate_list[i] * rate_PS_sigma / rate_PS ** 2) ** 2)
                    rejection_PS_sigma_list.append(rejection_PS_sigma)



                    # rejection per keV, PK meaning Per keV Per scattering
                    counts_times_keV = cumulative_threshold_array[0][i]+(threshold-hist_array[0][1][i])*(cumulative_threshold_array[0][i+1]-cumulative_threshold_array[0][i])/(hist_array[0][1][i+1]-hist_array[0][1][i])

                    rate_PK = Rate_factor*(counts_times_keV)
                    rate_PK_sigma = rate_PK / np.sqrt(counts)
                    rejection_PK = exp_rate_list[j] / rate_PK
                    rejection_PK_list.append(rejection_PK)
                    rejection_sigma = np.sqrt(
                        (exp_sigma_list[i] / rate_PK) ** 2 + (exp_rate_list[i] * rate_PK_sigma / rate_PK ** 2) ** 2)
                    rejection_PK_sigma_list.append(rejection_sigma)
                    break
        output_dict = {
            'Pressure [bara]': bkg_pressure_recon_list,
            'Updated Setiz [keV]': updated_Setiz_list,
            "Exp Rate [mHz]": exp_rate_list,
            "Bkg Rate [mHz]": background_rate_list,
            "Clean Rate [mHz]": clean_rate_list,
            "Exp Sigma [mHz]": exp_sigma_list,
            "Bkg Sigma [mHz]": background_sigma_list,
            "Clean Sigma [mHz]": clean_sigma_list,
            "Rejection Rate Scattering[mHz]": rejection_PS_list,
            "Rejection Sigma Scattering[mHz]": rejection_PS_sigma_list,
            "Rejection Rate KeV[mHz]": rejection_PK_list,
            "Rejection Sigma KeV[mHz]": rejection_PK_sigma_list}
        print(output_dict)
        df = pd.DataFrame(output_dict)
        print(df)
        save_path = os.path.join(self.plot_path, expfile_name + "_output.txt")
        df.to_csv(save_path, index=False)


        plt.savefig(self.plot_path + "Co_1E7_ER_rejection_vs_Setiz_perkeV.pdf")
    def read_exposure(self,filename):
        # Define your path (we'll use a relative path)
        file_path = os.path.join('..', 'exp_exposure', filename)

        # Read the file
        # sep='\s+' handles any number of spaces or tabs as delimiters
        df = pd.read_csv(file_path, sep='\s+', skiprows=1, header=None)

        return df

    def read_original_spectrum(self):
        list1, list2, list3 = [], [], []
        b_older = 0
        with open("./AmLiNO3_simulated_prob_MeV_raw.txt") as f:
            for line in f:
                a, b = line.split()
                list1.append(float(a))
                list2.append(float(b))
                b_older=b
                list3.append(float(b))

        return(list2, list1)


class test_csv():
    def __init__(self):
        list1 = [1,3,4.5,6.7,8.9]
        with open("/data/runzezhang/test.csv", 'w', newline='') as myfile:
            wr = csv.writer(myfile)
            wr.writerow(list1)
        with open("/data/runzezhang/test.csv", 'r') as file:
            reader = csv.reader(file)
            # Read the first row (assuming single row for simplicity)
            number_list = next(reader)
            # Convert the strings to floats
            number_list = [float(value) for value in number_list]

        print(number_list)

if __name__=="__main__":
    sn = SN()
    # test = test_csv()