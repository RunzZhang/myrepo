import pandas as pd
import matplotlib.pyplot as plt
import csv
import numpy as np
import os
import pickle
class SN():
    def __init__(self):
        # after generate new files, you need to select the capture ratio/source for different configs in read_files function.
        # then choose the correct signal/noise of with clause in read files.
        # at last change the self.name and plot_name in plot function
        # v2: change back to 2 backgrounds but with finer definitions
        # v4 kill duplicated NRERs
        self.doped_path = "/data/runzezhang/result/GR_sims/chunked_root_files_Ba_1E6_phot_track/"

        # self.base_path = "/data/runzezhang/result/GR_sims/chunked_root_files_Ba_1E5_ar_inside/"
        # self.base_path2 = "/data/runzezhang/result/GR_sims/chunked_root_files_Ba_1E5_ar_inside/"  # for gamma path


        # self.base_path = "/data/runzezhang/result/GR_sims/chunked_root_files_Ba_1E6_CF4_inside/"
        # self.base_path2 = "/data/runzezhang/result/GR_sims/chunked_root_files_Ba_1E6_CF4_inside/"

        # self.base_path = "/data/runzezhang/result/GR_sims/chunked_root_files_Ba_ar_1E6_phot/"
        # self.base_path2 = "/data/runzezhang/result/GR_sims/chunked_root_files_Ba_ar_1E6_phot/"

        self.base_path = "/data/runzezhang/result/GR_sims/chunked_root_files_low_e_1E5/"
        self.base_path2 = "/data/runzezhang/result/GR_sims/chunked_root_files_low_e_1E5/"

        # self.base_path = "/data/runzezhang/result/GR_sims/chunked_root_files_Ba_1E5_lar/"
        # self.base_path2 = "/data/runzezhang/result/GR_sims/chunked_root_files_Ba_1E5_lar/"  # for gamma path

        self.plot_path = '/data/runzezhang/result/GR_sims/plot/'
        self.false_1 = "PN_false1.csv"
        self.false_2 = "PN_false2.csv"
        self.signal = "PN_sig.csv"
        self.name1 = "Correlated ER Background"
        self.name2 = "Hard Scatter Background"
        self.name = "Backgrounds"
        self.plot_name = self.name+"PN_1E6_wt_gamma_outside.pdf"
        self.gamma_bkg_path = self.plot_path+"Co_gamma position_distribution.pdf"
        self.pho_threshold = 100

        cols = ["Event","name", "R/mm", "Z/mm", "Volume", "Process", "ER_near/eV", "Multiplicity"]

        self.df_primary_list = []
        self.df_all_list = []
        self.df_phot_list = []



        #982 statics false 1
        for i in range(1,54):
        # for i in range(1, 26):
        # for i in range(26, 54):
            try:
        # for i in range(1, 11):
                self.main_body(i)
            except Exception as e:
                print(e)
                continue
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

        self.info_primary_path = self.base_path + f"Cs_gamma_1E6_info_primary_scube_part{i}.csv"
        self.info_all_path = self.base_path + f"Cs_gamma_1E6_info_scube_all_part{i}.csv"
        self.info_phot_path = self.doped_path + f"Cs_gamma_1E6_info_scube_phot_part{i}.csv"

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

        self.gamma_rate = 2.58*1.67e4 # /s


        self.G4_events_gamma =  5E6 # only 50 chunks
        # self.G4_events_gamma = 1E8  # only 50 chunks
        self.G4_phot_gamma = 1E8
        self.ambient_bubble = 5 # /h

        temp_df_primary = pd.read_csv(self.info_primary_path)

        self.df_primary_list.append(temp_df_primary)

        temp_df_all = pd.read_csv(self.info_all_path)

        self.df_all_list.append(temp_df_all)

        temp_df_phot = pd.read_csv(self.info_phot_path)
        # print("phot temp", temp_df_phot)
        self.df_phot_list.append(temp_df_phot)

    def combine_df(self):
        self.merged_df_primary = pd.concat(self.df_primary_list, ignore_index=True)
        self.merged_df_all = pd.concat(self.df_all_list, ignore_index=True)
        # self.merged_df_phot = pd.concat(self.df_phot_list, ignore_index=True)
        self.merged_df_phot = self.df_phot_list[0]  # usually photo_list only has 1 chunked file

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
        # self.read_ER_Ar_CF_1d_sum_counts()

        # self.read_emit_spectrum()
        # self.gamma_rejection_rate_per_keV_vs_Setiz()
        # self.gamma_rejection_rate_vs_Setiz()

        # self.write_sims_results()
        self.plot_sims_results()

        #doped analyasis
        # self.read_ER_Ar_doped()
        # self.write_doped_sims_results()

        # photo process analysis
        # self.read_ER_Ar_pho_per_deposit_rate()
    def read_emit_spectrum(self):
        df_init_emit = self.merged_df[(self.merged_df["Volume"]=='calibration_Be_phys')&(self.merged_df["name"]=='gamma')&(self.merged_df["Step ID"]==1)]
        fig, ax = plt.subplots()
        ax.hist(df_init_emit["PreKinetic/MeV"]*1000, bins=1200, range=(0,600))
        ax.set_xlabel("Gamma Energy/keV")
        ax.set_ylabel("Counts")
        plt.savefig(self.plot_path + "Ba_init_spectrum.pdf")

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
        plt.savefig(self.plot_path+"Co_1E5_position.pdf")

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
        plt.savefig(self.plot_path + "Co_1E5_zslice_position_density.pdf")

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
        plt.savefig(self.plot_path+"Co_1E5_position_density.pdf")


    def read_multiplicity(self):
        multiplicity= self.merged_df.groupby("Event")["Multiplicity"].max()
        max_m = self.merged_df["Multiplicity"].max()
        print(max_m)
        bins = np.arange(1, max_m + 2)
        fig,ax = plt.subplots()
        ax.hist(multiplicity, bins= bins, align="left", rwidth=0.9, density=True)

        ax.set_xlabel("Multiplicity")
        ax.set_ylabel("Probability")
        plt.savefig(self.plot_path+"Co_1E7_multi.pdf")
    def read_ER_Ar_doped(self):
        # per energy deposition
        # the sum is per event

        ER_Ar = self.merged_df[self.merged_df["Volume"] == "LAr_phys"]["PreKinetic/MeV"] * 1000
        # the first peak
        print(self.merged_df[(self.merged_df["Volume"] == "LAr_phys")&(self.merged_df["PreKinetic/MeV"].between(0.0025,0.0040))])
        fig, ax = plt.subplots(1, 3, figsize=(14, 4))
        # ax[0].hist(ER_Ar, bins=60,range= (0,600),align="left")
        array = ax[0].hist(ER_Ar, bins=1400, range=(0, 700), align="left")
        print("counts", array[0][:20])
        print("bins", array[1][:20])
        # bins=12000, range=(0, 1200))
        ax[0].set_xlabel("photo absorption of Ar [keV] ")
        ax[0].set_ylabel("Counts")
        ax[0].set_yscale("log")
        # ax[0].set_title("Co source - 1332 keV gamma")
        ax[0].set_title("Ba source")


        plt.savefig(self.plot_path + "Ba_1E6_ar_photo.pdf")

    def read_Ar_multiplicity(self):
        multiplicity = self.merged_df[self.merged_df["Volume"]=="LAr_phys"].groupby("Event")["Multiplicity"].max()
        max_m = self.merged_df["Multiplicity"].max()
        print(max_m)
        bins = np.arange(1, max_m + 2)
        fig, ax = plt.subplots()
        ax.hist(multiplicity, bins=bins, align="left", rwidth=0.9, density=True)

        ax.set_xlabel("Argon Multiplicity")
        ax.set_ylabel("Probability")
        plt.savefig(self.plot_path + "Co_1E7_argon_multi.pdf")
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

        plt.savefig(self.plot_path + "Co_1E5_ER_all.pdf")
    def read_ER_Ar_pho_per_deposit_rate(self):
        # per energy deposition and total
        # MHz
        print("photo ",self.merged_phot_df.head(10))
        Rate_factor = self.gamma_rate * 1000 / (self.G4_events_gamma)
        ER_Ar = self.merged_phot_df[self.merged_phot_df["Volume"] == "LAr_phys"]["ER_near/eV"] / 1000


        hist_array = [None] * 3

        hist_array[0] = np.histogram(ER_Ar, bins=100, range=(0, 10))


        fig, ax = plt.subplots(1, 3, figsize=(16, 4))
        ax[0].bar(hist_array[0][1][:-1], Rate_factor * hist_array[0][0], width=np.diff(hist_array[0][1]),
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

        plt.savefig(self.plot_path + "Ba_5E6_ER_perdepostion_photo.pdf")

    def read_ER_Ar_CF_per_deposit_rate(self):
        # per energy deposition and total
        # MHz
        Rate_factor = self.gamma_rate * 1000 / (self.G4_events_gamma)
        # ER_Ar = self.merge_df_primary[self.merged_df_primary["Volume"]=="LAr_phys"]["ER_near/eV"]/1000

        ER_Ar = \
            self.merged_df_primary[
                (self.merged_df_all["Volume"] == "LAr_phys")][
                "ER_near/eV"]
        ER_Ar_photo_total = \
            self.merged_df_all[
                (self.merged_df_all["Volume"] == "LAr_phys") ][
                "ER_near/eV"]
        ER_sum = \
            self.merged_df_all[
                (self.merged_df_all["Volume"] == "LAr_phys") & (self.merged_df_all["Process"] == "compt")][
                "ER_near/eV"]

        # ER_sum = \
        #     self.merged_df_all[
        #         (self.merged_df_all["Volume"] == "LAr_phys") ][
        #         "ER_near/eV"] / 1000



        # estimate gamma rejction level
        print(len(ER_sum), len(self.merged_df_all["ER_near/eV"]))
        hist_array = [None] * 3

        hist_array[0] = np.histogram(ER_Ar, bins=1200, range=(-1, 1200))
        hist_array[1] = np.histogram(ER_Ar_photo_total, bins=1200, range=(-1, 1200))
        hist_array[2] = np.histogram(ER_sum, bins=1200, range=(-1, 1200))

        fig, ax = plt.subplots(1, 3, figsize=(16, 4))
        # ax[0].bar(hist_array[0][1][:-1], Rate_factor * hist_array[0][0], width=np.diff(hist_array[0][1]),
        #           align="edge",
        #           edgecolor="black")
        histarray0_1d = np.insert(hist_array[0][0], 0, 0)

        ax[0].plot(hist_array[0][1], histarray0_1d)
        bin0_len = int(hist_array[0][1][1] - hist_array[0][1][0])
        ax[0].set_xlabel("ER/eV per photo in LAr")
        ax[0].set_ylabel(" Count /(bin[" + str(bin0_len) + " eV])")
        # ax[0].ticklabel_format(axis="y",style="sci", scilimits=(0, 0) )
        ax[0].set_yscale("log")
        ax[0].minorticks_on()
        # ax[0].grid(which="major", linestyle="-", linewidth=0.8, alpha=0.7)
        # ax[0].grid(which="minor", linestyle=":", linewidth=0.5, alpha=0.4)

        # ax[1].bar(hist_array[1][1][:-1], Rate_factor * hist_array[1][0], width=np.diff(hist_array[1][1]),
        #           align="edge",
        #           edgecolor="black")
        histarray1_1d = np.insert(hist_array[1][0], 0, 0)
        ax[1].plot(hist_array[1][1],histarray1_1d)
        ax[1].set_xlabel("ER/eV per photo in LAr")
        bin1_len = int(hist_array[1][1][1] - hist_array[1][1][0])
        ax[1].set_ylabel("Count /([" + str(bin1_len) + " eV])")
        ax[1].set_yscale("log")
        # ax[1].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        # ax[1].grid(True)
        ax[1].minorticks_on()
        # ax[1].grid(which="major", linestyle="-", linewidth=0.8, alpha=0.7)
        # ax[1].grid(which="minor", linestyle=":", linewidth=0.5, alpha=0.4)

        # ax[2].bar(hist_array[2][1][:-1], Rate_factor * hist_array[2][0], width=np.diff(hist_array[2][1]),
        #           align="edge",
        #           edgecolor="black")
        histarray2_1d = np.insert(hist_array[2][0], 0, 0)
        ax[2].plot(hist_array[2][1], histarray2_1d)
        ax[2].set_xlabel("ER/eV per deposition ")
        bin2_len = int(hist_array[2][1][1] - hist_array[2][1][0])
        ax[2].set_ylabel("Count /([" + str(bin2_len) + " eV])")
        # ax[2].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        ax[2].set_yscale("log")
        # ax[2].grid(True)
        ax[2].minorticks_on()
        # ax[2].grid(which="major", linestyle="-", linewidth=0.8, alpha=0.7)
        # ax[2].grid(which="minor", linestyle=":", linewidth=0.5, alpha=0.4)

        plt.savefig(self.plot_path + "low_e_check.pdf")


    def read_ER_Ar_CF_per_deposit_rate_cumulative(self):
        # rate factor in mHz


        Rate_factor = self.gamma_rate*1000 / (self.G4_events_gamma)
        ER_Ar = self.merged_df[self.merged_df["Volume"] == "LAr_phys"]["ER_near/eV"] / 1000
        ER_CF4 = self.merged_df[self.merged_df["Volume"] == "hydraulic_fluid_phys"]["ER_near/eV"] / 1000
        ER_CF4 = ER_Ar # just plot the argon
        ER_sum = self.merged_df["ER_near/eV"] / 1000
        print('max(ER_Ar)', max(ER_Ar))
        print('max(ER_CF4)', max(ER_CF4))

        hist_array = [None] * 3

        hist_array[0] = np.histogram(ER_Ar, bins=100,range=(0, 1200))
        hist_array[1] = np.histogram(ER_CF4, bins=100,range=(0, 1200))
        hist_array[2] = np.histogram(ER_sum, bins=100,range=(0, 1200))
        cumulative_threshold_array = [None]*3

        cumulative_threshold_array[0] = np.array([sum(hist_array[0][0][i:]) for i in range(len(hist_array[0][0]))])
        cumulative_threshold_array[1] = np.array([sum(hist_array[1][0][i:]) for i in range(len(hist_array[1][0]))])
        cumulative_threshold_array[2] = np.array([sum(hist_array[2][0][i:]) for i in range(len(hist_array[2][0]))])
        # find first 2 bins and rate for argon
        print("argon bin", hist_array[0][1][:4])
        print("argon rate", Rate_factor * cumulative_threshold_array[0][:3])
        # find if compton edge exist in LAr cumulative spectrum
        print("argon maximum energy", max(hist_array[0][1][:-1]))
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
        ax[0].set_ylabel("ER Rate mHz")
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

        plt.savefig(self.plot_path + "Ba_1E7_ER_perdepostion_cumulative.pdf")

    def gamma_rejection_rate_vs_Setiz(self):
        # rate factor in mHz
        expfile_name = "60Co-12_15-16_exposures"
        exp_file_list = ["60Co-12_15-16_exposures_zoom"]

        bkgfile_name = "Background-11_26-30_exposures"
        Rate_factor = self.gamma_rate * 1000 / (self.G4_events_gamma)
        ER_Ar = self.merged_df[self.merged_df["Volume"] == "LAr_phys"]["ER_near/eV"] / 1000  # in keV

        for expfile_name in exp_file_list:
            hist_array = [None]
            hist_array[0] = np.histogram(ER_Ar, bins=100, range=(0, 1200))

            # get probablity per scattering and the statistics
            cumulative_threshold_per_scatter_array = [None]

            cumulative_threshold_per_scatter_array[0] = np.array(
                [sum(hist_array[0][0][i:]) for i in range(len(hist_array[0][0]))])

            # histogram per scattering per keV
            cumulative_threshold_array = [None]
            energy_deposit_list = [hist_array[0][0][i] * hist_array[0][1][i] for i in range(len(hist_array[0][0]))]
            cumulative_threshold_array[0] = np.array(
                [sum(energy_deposit_list[i:]) for i in range(len(energy_deposit_list))])

            source_exposure_df = self.read_exposure(expfile_name + ".txt")
            # print(source_exposure_df.loc[:, 0])
            print("txt read source", source_exposure_df)
            background_exposure_df = self.read_exposure(bkgfile_name + ".txt")
            print("txt read bkg", background_exposure_df)
            Seitz_pressure_list = np.arange(2.25, 6.5, 0.25)
            print(Seitz_pressure_list)
            Setiz = [np.float64(1.3445287166423177), np.float64(1.4677096307281403), np.float64(1.6077252261931916),
                     np.float64(1.7676644948295235), np.float64(1.9513368144218666), np.float64(2.163478457894038),
                     np.float64(2.41003021032974), np.float64(2.698514892785409), np.float64(3.038557782566206),
                     np.float64(3.44261366411884), np.float64(3.9269986461894346), np.float64(4.513378617028501),
                     np.float64(5.230956628311766), np.float64(6.119753114689943), np.float64(7.235636873731044),
                     np.float64(8.658238672578), np.float64(10.503756197125261)]  # in keV
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
                    pressure_index_seitz = int(np.where(Seitz_pressure_list == src_pressure)[0])
                    # pressure_index_seitz = Seitz_pressure_list.index(src_pressure)
                    print(src_pressure, pressure_index)
                    source_bkg_pressure_match = True
                except:
                    print("source pressure is not found in background ", src_pressure)
                    source_bkg_pressure_match = False
                    # only get pressure entries that shows in both src and bkg
                if source_bkg_pressure_match:
                    updated_Setiz_list.append(Setiz[pressure_index_seitz])
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
                    if threshold >= hist_array[0][1][i]:
                        # rejection per scattering, PS meaning perscattering
                        counts = cumulative_threshold_per_scatter_array[0][i] + (threshold - hist_array[0][1][i]) * (
                                cumulative_threshold_per_scatter_array[0][i + 1] -
                                cumulative_threshold_per_scatter_array[0][i]) / (
                                         hist_array[0][1][i + 1] - hist_array[0][1][i])
                        rate_PS = Rate_factor * (counts)

                        rate_PS_sigma = rate_PS / np.sqrt(counts)

                        rejection_PS = exp_rate_list[j] / rate_PS
                        rejection_PS_list.append(rejection_PS)
                        rejection_PS_sigma = np.sqrt(
                            (exp_sigma_list[i] / rate_PS) ** 2 + (exp_rate_list[i] * rate_PS_sigma / rate_PS ** 2) ** 2)
                        rejection_PS_sigma_list.append(rejection_PS_sigma)

                        # rejection per keV, PK meaning Per keV Per scattering
                        counts_times_keV = cumulative_threshold_array[0][i] + (threshold - hist_array[0][1][i]) * (
                                cumulative_threshold_array[0][i + 1] - cumulative_threshold_array[0][i]) / (
                                                   hist_array[0][1][i + 1] - hist_array[0][1][i])

                        rate_PK = Rate_factor * (counts_times_keV)
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


            save_path = os.path.join(self.plot_path, expfile_name + "_output.txt")
            df.to_csv(save_path, index=False)
            print("save to path", save_path)

    def gamma_rejection_rate_per_keV_vs_Setiz(self):
        # rate factor in mHz

        expfile_name = "60Co-12_15-16_exposures"
        exp_file_list = ["60Co-12_15-16_exposures"]
        bkgfile_name = "Background-11_26-30_exposures"

        Rate_factor = self.gamma_rate * 1000 / (self.G4_events_gamma)
        ER_Ar = self.merged_df[self.merged_df["Volume"] == "LAr_phys"]["ER_near/eV"] / 1000  # in keV

        hist_array = [None]
        # hist_array[0] = np.histogram(ER_Ar, bins=100, range=(0, 1200))
        hist_array[0] = np.histogram(ER_Ar, bins=12000, range=(0, 1200))
        # 12 keV -> 1 Setiz threshold there is no change for gamma rejection
        # we need 0.1 keV, and this gives us 4800 bins

        # transfer edge to mid point per bin

        # get probablity per scattering and the statistics
        cumulative_threshold_per_scatter_array = [None]

        cumulative_threshold_per_scatter_array[0] = np.array(
            [sum(hist_array[0][0][i:]) for i in range(len(hist_array[0][0]))])

        # histogram per scattering per keV
        cumulative_threshold_array = [None]
        energy_deposit_list = [hist_array[0][0][i] * hist_array[0][1][i] for i in range(len(hist_array[0][0]))]

        cumulative_threshold_array[0] = np.array(
            [sum(energy_deposit_list[i:]) for i in range(len(energy_deposit_list))])
        print("total count* energy Co", cumulative_threshold_per_scatter_array[0][0],cumulative_threshold_per_scatter_array[0][0]/cumulative_threshold_array[0][0])

        for expfile_name in exp_file_list:
            source_exposure_df = self.read_exposure(expfile_name + ".txt")
            # print(source_exposure_df.loc[:, 0])
            print("txt read source", source_exposure_df)
            background_exposure_df = self.read_exposure(bkgfile_name + ".txt")
            print("txt read bkg", background_exposure_df)
            Seitz_pressure_list = np.arange(2.25, 6.5, 0.25)
            print(Seitz_pressure_list)
            # keV
            Q_setiz = [1.3445287166423177, 1.4677096307281403, 1.6077252261931916, 1.7676644948295235,
                       1.9513368144218666, 2.163478457894038, 2.41003021032974, 2.698514892785409, 3.038557782566206,
                       3.44261366411884, 3.9269986461894346, 4.513378617028501, 5.230956628311766, 6.119753114689943,
                       7.235636873731044, 8.658238672578, 10.503756197125261]
            # keV
            E_ion = [0.6935170594033994, 0.7463906266235716, 0.8055756681017613, 0.8721142361762, 0.9472715210586028,
                     1.0325954589172766, 1.12999567281159, 1.2418491587814853, 1.371143722549567, 1.5216750993325374,
                     1.6983220883949837, 1.907436620175941, 2.1574068879170554, 2.459486450041476, 2.829041931492084,
                     3.2874775552634508, 3.865286624630827]
            # g / cc
            rho_l = [1.190713998903872, 1.190901456615719, 1.1910886308789765, 1.1912755227862943, 1.191462133423365,
                     1.1916484638689993, 1.1918345151951764, 1.192020288467111, 1.192205784743309, 1.1923910050756303,
                     1.1925759505093436, 1.1927606220831863, 1.1929450208294252, 1.1931291477739037, 1.1933130039361075,
                     1.193496590329213, 1.1936799079601481]
            # nm
            Rl = [5.609511255785613, 5.807717275806958, 6.020272497298806, 6.248793390358893, 6.495148643090569,
                  6.761510412830444, 7.050418580730343, 7.364861891716328, 7.70838179881121, 8.085206324999868,
                  8.500425035138141, 8.960220171784458, 9.472176544801727, 10.0457032178012, 10.69261696165673,
                  11.427965165560835, 12.271210788702833]

            compound_x = []
            for i in range(len(E_ion)):
                x = E_ion[i] / (rho_l[i] * Rl[i])
                compound_x.append(x)

            Setiz = [np.float64(1.3445287166423177), np.float64(1.4677096307281403), np.float64(1.6077252261931916),
                     np.float64(1.7676644948295235), np.float64(1.9513368144218666), np.float64(2.163478457894038),
                     np.float64(2.41003021032974), np.float64(2.698514892785409), np.float64(3.038557782566206),
                     np.float64(3.44261366411884), np.float64(3.9269986461894346), np.float64(4.513378617028501),
                     np.float64(5.230956628311766), np.float64(6.119753114689943), np.float64(7.235636873731044),
                     np.float64(8.658238672578), np.float64(10.503756197125261)]  # in keV
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
            updated_compoundx_list = []
            updated_Eion_list = []

            for i in range(len(source_pressure_list)):
                src_pressure = source_pressure_list[i]
                source_bkg_pressure_match = True
                try:
                    pressure_index = bkg_pressure_list.index(src_pressure)
                    pressure_index_seitz = int(np.where(Seitz_pressure_list == src_pressure)[0])
                    # pressure_index_seitz = Seitz_pressure_list.index(src_pressure)
                    print(src_pressure, pressure_index)
                    source_bkg_pressure_match = True
                except:
                    print("source pressure is not found in background ", src_pressure)
                    source_bkg_pressure_match = False
                    # only get pressure entries that shows in both src and bkg
                if source_bkg_pressure_match:
                    updated_Setiz_list.append(Setiz[pressure_index_seitz])
                    updated_compoundx_list.append(compound_x[pressure_index_seitz])
                    updated_Eion_list.append(E_ion[pressure_index_seitz])
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
                threshold_Eion = updated_Eion_list[j]
                for i in range(len(hist_array[0][1])):
                    if threshold >= hist_array[0][1][i]:
                        # rejection per scattering, PS meaning perscattering
                        counts = cumulative_threshold_per_scatter_array[0][i] + (threshold - hist_array[0][1][i]) * (
                                cumulative_threshold_per_scatter_array[0][i + 1] -
                                cumulative_threshold_per_scatter_array[0][i]) / (
                                         hist_array[0][1][i + 1] - hist_array[0][1][i])
                        rate_PS = Rate_factor * (counts)

                        rate_PS_sigma = rate_PS / np.sqrt(counts)

                        rejection_PS = clean_rate_list[j] / rate_PS
                        rejection_PS_list.append(rejection_PS)
                        rejection_PS_sigma = np.sqrt(
                            (clean_sigma_list[i] / rate_PS) ** 2 + (
                                        clean_rate_list[i] * rate_PS_sigma / rate_PS ** 2) ** 2)
                        rejection_PS_sigma_list.append(rejection_PS_sigma)
                    if threshold_Eion >= hist_array[0][1][i]:
                        # rejection per keV, PK meaning Per keV Per scattering

                        # counts_times_keV = cumulative_threshold_array[0][i] + (threshold_Eion - hist_array[0][1][i]) * (
                        #         cumulative_threshold_array[0][i + 1] - cumulative_threshold_array[0][i]) / (
                        #                            hist_array[0][1][i + 1] - hist_array[0][1][i]) # interpolation
                        counts_times_keV = cumulative_threshold_array[0][0]  # all energy
                        counts_Eion = cumulative_threshold_per_scatter_array[0][i] + (
                                    threshold_Eion - hist_array[0][1][i]) * (
                                              cumulative_threshold_per_scatter_array[0][i + 1] -
                                              cumulative_threshold_per_scatter_array[0][i]) / (
                                              hist_array[0][1][i + 1] - hist_array[0][1][i])

                        rate_PK = Rate_factor * (counts_times_keV)
                        rate_PK_sigma = rate_PK / np.sqrt(counts_Eion)
                        rejection_PK = clean_rate_list[j] / rate_PK
                        rejection_PK_list.append(rejection_PK)
                        rejection_sigma = np.sqrt(
                            (clean_sigma_list[i] / rate_PK) ** 2 + (
                                        clean_rate_list[i] * rate_PK_sigma / rate_PK ** 2) ** 2)
                        rejection_PK_sigma_list.append(rejection_sigma)
                        break
            output_dict = {
                'Pressure [bara]': bkg_pressure_recon_list,
                'Updated Setiz [keV]': updated_Setiz_list,
                'Eion_rl-1_rhol-1 [10GeVcm**2 g-1]': updated_compoundx_list,
                "Exp Rate [mHz]": exp_rate_list,
                "Bkg Rate [mHz]": background_rate_list,
                "Clean Rate [mHz]": clean_rate_list,
                "Exp Sigma [mHz]": exp_sigma_list,
                "Bkg Sigma [mHz]": background_sigma_list,
                "Clean Sigma [mHz]": clean_sigma_list,
                "Rejection Rate Scattering[]": rejection_PS_list,
                "Rejection Sigma Scattering[]": rejection_PS_sigma_list,
                "Rejection Rate KeV[/keV]": rejection_PK_list,
                "Rejection Sigma KeV[/keV]": rejection_PK_sigma_list}

            df = pd.DataFrame(output_dict)
            # print(df["Rejection Rate Scattering[mHz]"])
            # print(df["Rejection Rate KeV[mHz]"])
            print(df["Exp Rate [mHz]"])
            print(df["Bkg Rate [mHz]"])
            print(df["Clean Rate [mHz]"])
            save_path = os.path.join(self.plot_path, expfile_name + "_output.txt")
            df.to_csv(save_path, index=False)
        # check the shape of two array
        fig, ax = plt.subplots(1, 2, figsize=(10, 4))
        ax[0].plot(hist_array[0][1][:-1], cumulative_threshold_per_scatter_array[0])
        ax[0].set_xlabel("thershold [keV]")
        ax[0].set_ylabel("Counts")
        ax[0].set_title("Cumulative counts vs threshold")
        ax[0].set_xlim(0, 5)
        ax[0].set_ylim(3e5, 3.5e5)
        # ax[0].set_yscale("log")

        ax[1].plot(hist_array[0][1][:-1], cumulative_threshold_array[0])
        ax[1].set_xlabel("thershold [keV]")
        ax[1].set_ylabel("Counts*energy [KeV]")
        ax[1].set_title("Cumulative counts*energy vs threshold")
        ax[1].set_xlim(0, 5)
        ax[1].set_ylim(3.1e7, 3.12e7)

        for k in range(20):
            print(cumulative_threshold_array[0][k] / cumulative_threshold_per_scatter_array[0][k])
            print(hist_array[0][1][k])
        # ax[1].set_yscale("log")
        print(hist_array[0][1][1] - hist_array[0][1][0], "keV width")
        plt.savefig(self.plot_path + "Co_cumulative_counts_function.pdf")
    def plot_sims_results(self):
        # rate factor in mHz

        Rate_factor = self.gamma_rate * 1000 / (self.G4_events_gamma)
        ER_Ar_primary = self.merged_df_primary[self.merged_df_primary["Volume"] == "LAr_phys"][
                            "ER_near/eV"]  # in keV

        hist_array_primary = [None]
        # hist_array[0] = np.histogram(ER_Ar, bins=100, range=(0, 1200))
        hist_array_primary[0] = np.histogram(ER_Ar_primary, bins=10000, range=(0, 1000))
        # 12 keV -> 1 Setiz threshold there is no change for gamma rejection
        # we need 0.1 keV, and this gives us 4800 bins

        # transfer edge to mid point per bin

        # get probablity per scattering and the statistics
        cumulative_threshold_per_scatter_array_primary = [None]

        cumulative_threshold_per_scatter_array_primary[0] = np.array(
            [sum(hist_array_primary[0][0][i:]) for i in range(len(hist_array_primary[0][0]))])

        # histogram per scattering per keV
        cumulative_threshold_array_primary = [None]
        energy_deposit_list_primary = [hist_array_primary[0][0][i] * hist_array_primary[0][1][i] for i in
                                       range(len(hist_array_primary[0][0]))]

        cumulative_threshold_array_primary[0] = np.array(
            [sum(energy_deposit_list_primary[i:]) for i in range(len(energy_deposit_list_primary))])
        print("total count", cumulative_threshold_per_scatter_array_primary[0][0],
              cumulative_threshold_per_scatter_array_primary[0][0] / cumulative_threshold_array_primary[0][0])

        ER_Ar_all = self.merged_df_all[self.merged_df_all["Volume"] == "LAr_phys"][
                        "ER_near/eV"]  # in keV

        hist_array_all = [None]
        # hist_array[0] = np.histogram(ER_Ar, bins=100, range=(0, 1200))
        hist_array_all[0] = np.histogram(ER_Ar_all, bins=10000, range=(0, 1000))
        # 12 keV -> 1 Setiz threshold there is no change for gamma rejection
        # we need 0.1 keV, and this gives us 4800 bins

        # transfer edge to mid point per bin

        # get probablity per scattering and the statistics
        cumulative_threshold_per_scatter_array_all = [None]

        cumulative_threshold_per_scatter_array_all[0] = np.array(
            [sum(hist_array_all[0][0][i:]) for i in range(len(hist_array_all[0][0]))])

        # histogram per scattering per keV
        cumulative_threshold_array_all = [None]
        energy_deposit_list_all = [hist_array_all[0][0][i] * hist_array_all[0][1][i] for i in
                                   range(len(hist_array_all[0][0]))]

        cumulative_threshold_array_all[0] = np.array(
            [sum(energy_deposit_list_all[i:]) for i in range(len(energy_deposit_list_all))])
        print("total count", cumulative_threshold_per_scatter_array_all[0][0],
              cumulative_threshold_per_scatter_array_all[0][0] / cumulative_threshold_array_all[0][0])

        output_list = [Rate_factor, hist_array_all, cumulative_threshold_per_scatter_array_all[0],
                       cumulative_threshold_array_all[0],
                       hist_array_primary, cumulative_threshold_per_scatter_array_primary[0],
                       cumulative_threshold_array_primary[0]]

        fig, ax = plt.subplots(2, 4, figsize=(25, 10))
        print("counts pdf", hist_array_primary[0][0][:10])
        ax[0,0].plot(hist_array_primary[0][1][:-1],hist_array_primary[0][0])
        ax[0,0].set_xlabel("Energy [eV]")
        ax[0,0].set_ylabel("Counts")
        ax[0,0].set_xlim(0,1)

        ax[0,1].plot(hist_array_primary[0][1][:-1], cumulative_threshold_per_scatter_array_primary[0])
        ax[0,1].set_xlabel("Energy [eV]")
        ax[0,1].set_ylabel("Cumulative Counts")
        ax[0,1].set_xlim(0, 1)

        ax[0,2].plot(hist_array_primary[0][1][:-1], energy_deposit_list_primary)
        ax[0,2].set_xlabel("Energy [keV]")
        ax[0,2].set_ylabel("Energy deposit per bin [eV]")
        ax[0,2].set_xlim(0, 1)

        ax[0,3].plot(hist_array_primary[0][1][:-1], cumulative_threshold_array_primary[0])
        ax[0,3].set_xlabel("Energy [keV]")
        ax[0,3].set_ylabel("Cumulative Energy Deposit [eV]")
        ax[0,3].set_xlim(0, 1)

        ax[1,0].plot(hist_array_all[0][1][:-1], hist_array_all[0][0])
        ax[1,0].set_xlabel("Energy [eV]")
        ax[1,0].set_ylabel("Counts")
        ax[1,0].set_xlim(0, 1)

        ax[1,1].plot(hist_array_all[0][1][:-1], cumulative_threshold_per_scatter_array_all[0])
        ax[1,1].set_xlabel("Energy [eV]")
        ax[1,1].set_ylabel("Cumulative Counts")
        ax[1,1].set_xlim(0, 1)

        ax[1,2].plot(hist_array_all[0][1][:-1], energy_deposit_list_all)
        ax[1,2].set_xlabel("Energy [keV]")
        ax[1,2].set_ylabel("Energy deposit per bin [eV]")
        ax[1,2].set_xlim(0, 1)

        ax[1,3].plot(hist_array_all[0][1][:-1], cumulative_threshold_array_all[0])
        ax[1,3].set_xlabel("Energy [keV]")
        ax[1,3].set_ylabel("Cumulative Energy Deposit [eV]")
        ax[1,3].set_xlim(0, 1)

        plt.savefig(self.plot_path+"output_spectrum.pdf")






    def write_sims_results(self):
        # rate factor in mHz

        Rate_factor = self.gamma_rate * 1000 / (self.G4_events_gamma)
        ER_Ar_primary = self.merged_df_primary[self.merged_df_primary["Volume"] == "LAr_phys"][
                            "ER_near/eV"] / 1000  # in keV

        hist_array_primary = [None]
        # hist_array[0] = np.histogram(ER_Ar, bins=100, range=(0, 1200))
        hist_array_primary[0] = np.histogram(ER_Ar_primary, bins=12000, range=(0, 1200))
        # 12 keV -> 1 Setiz threshold there is no change for gamma rejection
        # we need 0.1 keV, and this gives us 4800 bins

        # transfer edge to mid point per bin

        # get probablity per scattering and the statistics
        cumulative_threshold_per_scatter_array_primary = [None]

        cumulative_threshold_per_scatter_array_primary[0] = np.array(
            [sum(hist_array_primary[0][0][i:]) for i in range(len(hist_array_primary[0][0]))])

        # histogram per scattering per keV
        cumulative_threshold_array_primary = [None]
        energy_deposit_list_primary = [hist_array_primary[0][0][i] * hist_array_primary[0][1][i] for i in
                                       range(len(hist_array_primary[0][0]))]

        cumulative_threshold_array_primary[0] = np.array(
            [sum(energy_deposit_list_primary[i:]) for i in range(len(energy_deposit_list_primary))])
        print("total count* energy Co", cumulative_threshold_per_scatter_array_primary[0][0],
              cumulative_threshold_per_scatter_array_primary[0][0] / cumulative_threshold_array_primary[0][0])

        ER_Ar_all = self.merged_df_all[self.merged_df_all["Volume"] == "LAr_phys"][
                        "ER_near/eV"] / 1000  # in keV

        hist_array_all = [None]
        # hist_array[0] = np.histogram(ER_Ar, bins=100, range=(0, 1200))
        hist_array_all[0] = np.histogram(ER_Ar_all, bins=12000, range=(0, 1200))
        # 12 keV -> 1 Setiz threshold there is no change for gamma rejection
        # we need 0.1 keV, and this gives us 4800 bins

        # transfer edge to mid point per bin

        # get probablity per scattering and the statistics
        cumulative_threshold_per_scatter_array_all = [None]

        cumulative_threshold_per_scatter_array_all[0] = np.array(
            [sum(hist_array_all[0][0][i:]) for i in range(len(hist_array_all[0][0]))])

        # histogram per scattering per keV
        cumulative_threshold_array_all = [None]
        energy_deposit_list_all = [hist_array_all[0][0][i] * hist_array_all[0][1][i] for i in
                                   range(len(hist_array_all[0][0]))]

        cumulative_threshold_array_all[0] = np.array(
            [sum(energy_deposit_list_all[i:]) for i in range(len(energy_deposit_list_all))])
        print("total count* energy Co", cumulative_threshold_per_scatter_array_all[0][0],
              cumulative_threshold_per_scatter_array_all[0][0] / cumulative_threshold_array_all[0][0])

        output_list = [Rate_factor, hist_array_all, cumulative_threshold_per_scatter_array_all[0],
                       cumulative_threshold_array_all[0],
                       hist_array_primary, cumulative_threshold_per_scatter_array_primary[0],
                       cumulative_threshold_array_primary[0]]
        # output form, rate facotr to mHz, enenrgy edges, counts above the bin edge, counts* counts above the bin edge
        with open("/data/runzezhang/result/GR_sims/Ba_output_5E6.pkl", "wb") as f:
            pickle.dump(output_list, f)


    def write_doped_sims_results(self):
        # rate factor in mHz

        Rate_factor = self.gamma_rate * 1000 / (self.G4_events_gamma)
        ER_Ar = self.merged_df_phot[self.merged_df_phot["Volume"] == "LAr_phys"]["PreKinetic/MeV"] * 1000  # in keV

        hist_array = [None]
        # hist_array[0] = np.histogram(ER_Ar, bins=100, range=(0, 1200))
        hist_array[0] = np.histogram(ER_Ar, bins=12000, range=(0, 1200))
        # 12 keV -> 1 Setiz threshold there is no change for gamma rejection
        # we need 0.1 keV, and this gives us 4800 bins

        # transfer edge to mid point per bin

        # get probablity per scattering and the statistics
        cumulative_threshold_per_scatter_array = [None]

        cumulative_threshold_per_scatter_array[0] = np.array(
            [sum(hist_array[0][0][i:]) for i in range(len(hist_array[0][0]))])

        # histogram per scattering per keV
        cumulative_threshold_array = [None]
        energy_deposit_list = [hist_array[0][0][i] * hist_array[0][1][i] for i in range(len(hist_array[0][0]))]

        cumulative_threshold_array[0] = np.array(
            [sum(energy_deposit_list[i:]) for i in range(len(energy_deposit_list))])
        print("total count* energy Co", cumulative_threshold_per_scatter_array[0][0],cumulative_threshold_per_scatter_array[0][0]/cumulative_threshold_array[0][0])

        output_list = [Rate_factor ,hist_array, cumulative_threshold_per_scatter_array[0], cumulative_threshold_array[0]]
        # output form, rate facotr to mHz, enenrgy edges, counts above the bin edge, counts* counts above the bin edge
        with open("/data/runzezhang/result/GR_sims/Ba_doped_output.pkl", "wb") as f:
            pickle.dump(output_list, f)
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


        plt.savefig(self.plot_path + "Co_1E5_CF4_ER_perdepostion_cumulative_coldrate.pdf")
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




        plt.savefig(self.plot_path + "Co_1E5_ER_coldrate.pdf")
    def read_ER_Ar_CF_1d_sum_counts(self):
        # calcualte sum counts of ER classified in Ar and CF4 per event
        # and sum rate is over both volume in Ar and CF4

        Rate_factor = self.gamma_rate / (3600*self.G4_events_gamma) # /h per geant run file
        ER_Ar_sum = self.merged_df[self.merged_df["Volume"]=="LAr_phys"].groupby("Event")["ER_near/eV"].sum()/1000
        ER_CF4_sum = self.merged_df[self.merged_df["Volume"]=="hydraulic_fluid_phys"].groupby("Event")["ER_near/eV"].sum()/1000

        ER_sum = self.merged_df.groupby("Event")["ER_near/eV"].sum() / 1000

        hist_array = [None] * 3




        fig, ax = plt.subplots(1,3, figsize=(14, 4))
        ax[0].hist(ER_Ar_sum, bins=50, range=(0,1400))
        bin_len = 1400/50
        ax[0].set_xlabel("ER/keV per event in LAr")
        ax[0].set_ylabel("Rate/(h*"+str(bin_len)+" keV)")
        ax[0].minorticks_on()
        # ax[0].grid(which="major", linestyle="-", linewidth=0.8, alpha=0.7)
        # ax[0].grid(which="minor", linestyle=":", linewidth=0.5, alpha=0.4)

        ax[1].hist(ER_CF4_sum, bins=50, range=(0,1400))
        ax[1].set_xlabel("ER/keV per event in CF4")

        ax[1].set_ylabel("Rate/(h*"+str(bin_len)+" keV)")
        # ax[1].grid(True)
        ax[1].minorticks_on()
        # ax[1].grid(which="major", linestyle="-", linewidth=0.8, alpha=0.7)
        # ax[1].grid(which="minor", linestyle=":", linewidth=0.5, alpha=0.4)

        ax[2].hist(ER_sum, bins=50, range=(0,1400))
        ax[2].set_xlabel("ER/keV per event ")
        ax[2].set_ylabel("Rate/(h*"+str(bin_len)+" keV)")
        # ax[2].grid(True)
        ax[2].minorticks_on()
        # ax[2].grid(which="major", linestyle="-", linewidth=0.8, alpha=0.7)
        # ax[2].grid(which="minor", linestyle=":", linewidth=0.5, alpha=0.4)




        plt.savefig(self.plot_path + "Co_1E5_ER_counts.pdf")

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




        plt.savefig(self.plot_path + "Co_1E5_ER_coldrate_cumulative.pdf")

    def read_ER_Ar_CF_2d_sum(self):
        # calcualte sum of ER classified in Ar and CF4 per event
        # 2d histogram


        ER_Ar_sum = self.merged_df[self.merged_df["Volume"]=="LAr_phys"].groupby("Event")["ER_near/eV"].sum()/1000
        ER_CF4_sum = self.merged_df[self.merged_df["Volume"]=="hydraulic_fluid_phys"].groupby("Event")["ER_near/eV"].sum()/1000

        evt = pd.concat([ER_Ar_sum, ER_CF4_sum], axis=1)
        evt.columns = ["ER_Ar_keV", "ER_CF4_keV"]  # <-- you create these names
        evt = evt.fillna(0)

        fig, ax = plt.subplots()
        print("max x and y", max(evt["ER_Ar_keV"]),max(evt["ER_CF4_keV"]))
        sc = ax.hist2d(evt["ER_Ar_keV"], evt["ER_CF4_keV"], bins=50,
                       cmap="plasma", norm="log", alpha=0.7)

        ax.set_xlabel("ER_Ar/keV")
        ax.set_ylabel("ER_CF4/keV")
        # ax.set_xlim(0, 200)
        # ax.set_ylim(-100, 800)
        cbar = plt.colorbar(sc[3], ax=ax)
        cbar.set_label("Counts(log)")
        plt.savefig(self.plot_path + "Co_1E5_ER_density.pdf")




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
    def read_exposure(self,filename):
        # Define your path (we'll use a relative path)
        file_path = os.path.join('..', 'exp_exposure', filename)

        # Read the file
        # sep='\s+' handles any number of spaces or tabs as delimiters
        df = pd.read_csv(file_path, sep='\s+', skiprows=1, header=None)

        return df



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