import pandas as pd
import matplotlib.pyplot as plt
import csv
import numpy as np
class SN():
    def __init__(self):
        # after generate new files, you need to select the capture ratio/source for different configs in read_files function.
        # then choose the correct signal/noise of with clause in read files.
        # at last change the self.name and plot_name in plot function
        # v2: change back to 2 backgrounds but with finer definitions
        # v4 kill duplicated NRERs
        self.base_path = "/data/runzezhang/result/TN_sims_D/chunked_root_files_Cs_1E5/"
        self.base_path2 = "/data/runzezhang/result/TN_sims_D/chunked_root_files_Cs_1E5/" # for gamma path

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
        # for i in range(1,101):
        # for i in range(1, 11):
        #     self.main_body(i)
        self.main_body(1)
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

        self.gamma_rate = 18e6 # /s


        self.G4_events_gamma =  1E5
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
        # ER distribution per row
        # self.read_ER_Ar_CF()
        # self.read_ER_Ar_CF_1d_sum()
        # self.read_ER_Ar_CF_2d_sum()
        self.read_ER_Ar_CF_1d_sum_rate()



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
        ax.set_ylabel("Counts")
        plt.savefig(self.plot_path+"Cs_1E5_multi.pdf")

    def read_ER_Ar_CF(self):
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

    def read_ER_Ar_CF_1d_sum(self):
        # calcualte sum of ER classified in Ar and CF4


        ER_Ar_sum = self.merged_df[self.merged_df["Volume"]=="LAr_phys"].groupby("Event")["ER_near/eV"].sum()/1000
        ER_CF4_sum = self.merged_df[self.merged_df["Volume"]=="hydraulic_fluid_phys"].groupby("Event")["ER_near/eV"].sum()/1000

        ER_sum = self.merged_df.groupby("Event")["ER_near/eV"].sum() / 1000

        fig, ax = plt.subplots(1,2, figsize=(9, 4))
        ax[0].hist(ER_Ar_sum, bins=50,align="left")
        ax[0].set_xlabel("ER/keV per event in LAr")
        ax[0].set_ylabel("Counts")

        ax[1].hist(ER_CF4_sum, bins=50,align="left")
        ax[1].set_xlabel("ER/keV per event in CF4")
        ax[1].set_ylabel("Counts")


        plt.savefig(self.plot_path + "Cs_1E5_ER_sum_volume.pdf")

    def read_ER_Ar_CF_1d_sum_rate(self):
        # calcualte sum of ER classified in Ar and CF4
        Rate_factor = self.gamma_rate / (3600*self.G4_events_gamma) # /h per geant run file
        ER_Ar_sum = Rate_factor*self.merged_df[self.merged_df["Volume"]=="LAr_phys"].groupby("Event")["ER_near/eV"].sum()/1000
        ER_CF4_sum = Rate_factor*self.merged_df[self.merged_df["Volume"]=="hydraulic_fluid_phys"].groupby("Event")["ER_near/eV"].sum()/1000

        ER_sum = Rate_factor*self.merged_df.groupby("Event")["ER_near/eV"].sum() / 1000

        fig, ax = plt.subplots(1,3, figsize=(14, 4))
        ax[0].hist(ER_Ar_sum, bins=50,align="left")
        ax[0].set_xlabel("ER/keV per event in LAr")
        ax[0].set_ylabel("Rate/h")

        ax[1].hist(ER_CF4_sum, bins=50,align="left")
        ax[1].set_xlabel("ER/keV per event in CF4")
        ax[1].set_ylabel("Rate/h")

        ax[2].hist(ER_sum, bins=50, align="left")
        ax[2].set_xlabel("ER/keV per event ")
        ax[2].set_ylabel("Rate/h")




        plt.savefig(self.plot_path + "Cs_1E5_ER_rate.pdf")


    def read_ER_Ar_CF_2d_sum(self):
        # calcualte sum of ER classified in Ar and CF4


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