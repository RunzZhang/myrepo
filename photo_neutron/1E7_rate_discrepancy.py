import pandas as pd
import matplotlib.pyplot as plt
import csv
import numpy as np
class SN():
    def __init__(self,gamma=False):
        # after generate new files, you need to select the capture ratio/source for different configs in read_files function.
        # then choose the correct signal/noise of with clause in read files.
        # at last change the self.name and plot_name in plot function
        # v2: change back to 2 backgrounds but with finer definitions
        # v4 kill duplicated NRERs
        self.base_path = "/data/runzezhang/result/TN_sims_D/chunked_root_files_pn_outside_1E7/"
        self.base_path2 = "/data/runzezhang/result/TN_sims_D/chunked_root_files_pn_1E6_outside_x71_nolead/" # for gamma path


        self.plot_path = '/data/runzezhang/result/TN_sims_D/plot/'
        self.false_1 = "PN_false1.csv"
        self.false_2 = "PN_false2.csv"
        self.signal = "PN_sig.csv"
        self.name1 = "Correlated ER Background"
        self.name2 = "Hard Scatter Background"
        self.name = "Backgrounds"
        self.plot_name = self.name+"PN_1E6_wt_gamma_outside.pdf"
        self.pho_threshold = 200

        self.neutron_ini_list = []
        self.neutron_ar_ke_list =[]
        self.neutron_ar_ke_alter_list = []
        self.neutron_x_list = []

        self.neutron_ini_list2 = []
        self.neutron_ar_ke_list2 = []
        self.neutron_ar_ke_alter_list2 = []
        self.neutron_x_list2 = []

        self.gamma = gamma


        #982 statics false 1
        # for i in range(1,101):
        for i in range(1, 11):
            self.main_body(i)


        self.plot_neutron_spectrum_ini_lin()
        self.plot_neutron_spectrum_Ar_lin()
        self.plot_neutron_spectrum_x_lin()


    def main_body(self,i):
        print(i)
        self.false_1 = f"PN_1E7_false1_part{i}.csv"
        self.false_2 = f"PN_1E7_false2_part{i}.csv"
        self.false_3 = f"PN_1E7_false3_part{i}.csv"
        self.false_gamma_1 = f"PN_gamma_1E7_false1_part{i}.csv"
        self.signal = f"PN_1E7_sig_part{i}.csv"



        self.ini_path = self.base_path+ f"PN_1E7_ini_part{i}.csv"
        self.ar_ke_path = self.base_path + f"PN_1E7_ke_part{i}.csv"
        self.ini_x_path = self.base_path + f"PN_1E7_inix_part{i}.csv"
        self.ar_ke_alter_path = self.base_path + f"PN_1E7_ke_part{i}.csv"



        self.ini_path2 = self.base_path2 + f"PN_1E7_ini_part{i}.csv"
        self.ar_ke_path2 = self.base_path2 + f"PN_1E7_ke_part{i}.csv"
        self.ini_x_path2 = self.base_path2 + f"PN_1E7_inix_part{i}.csv"
        self.ar_ke_alter_path2 = self.base_path2 + f"PN_1E7_ke_part{i}.csv"


        self.read_files()


        # self.read_files_s_to_N1()

# main funtion we use
    def read_files(self):
        self.original_Activity = 5 # original activity in the paper
        self.Activity = 5  # source practical activity in mivro curie for 50 bubbles/hour
        # self.Activity = 0.0416  # source activity in mivro curie
        # self.capture_ratio = 1.164E-3 # 1125eV 1.4g/cm Ar
        # self.capture_ratio = 0.121 # 400 eV 1.4g/cm3 Ar
        # self.capture_ratio = 0.267  # 350 eV
        self.capture_ratio = 0.116  # 400 eV
        # self.capture_ratio = 6.52E-3  # 700 eV
        # self.capture_ratio = 1.158E-3  # 1125 eV
        # self.rate = 435.6 #/s # CF neutron rate 9 mucurie
        # self.rate = 0.1968 #PN neutron rate /s
        # self.rate = 2.52e4  # PN neutron rate /s PNNL
        self.rate = 0.86  # PN neutron rate /s LZ 5micro Bismuth
        self.gamma_rate = 1.27e4 # PN gamma rate/s
        # self.G4_events= 1E5
        self.G4_events = 1E6
        self.G4_events_gamma =  1E6
        self.ambient_bubble = 5 # /h

        self.T = 1e-3
        self.G4_sig_time=(self.G4_events / self.rate)
        self.G4_noise_time = self.G4_events / self.rate
        self.G4_gamma_time = self.G4_events_gamma/self.gamma_rate

        # Initial amli spectrm
        with open(self.ini_path, 'r') as file:
            reader = csv.reader(file)
            # Read the first row (assuming single row for simplicity)
            number_list = next(reader)
            # Convert the strings to floats
            self.neutron_ini_raw_list = [float(value)*1e6 for value in number_list]


            # the [0] is NR number and [1:] is the photon numbers
        self.neutron_ini_list +=  self.neutron_ini_raw_list


        with open(self.ini_x_path, 'r') as file:
            reader = csv.reader(file)
            # Read the first row (assuming single row for simplicity)
            number_list = next(reader)
            # Convert the strings to floats
            self.neutron_ini_raw_xlist = [float(value)*1e6 for value in number_list]


            # the [0] is NR number and [1:] is the photon numbers
        self.neutron_x_list +=  self.neutron_ini_raw_xlist

        with open(self.ar_ke_path, 'r') as file:
            reader = csv.reader(file)
            # Read the first row (assuming single row for simplicity)
            number_list = next(reader)
            # Convert the strings to floats
            self.ar_ke_raw_list = [float(value)*1e6 for value in number_list] # in eV
            # self.noise_raw_list = [float(value)  for value in number_list]
        self.neutron_ar_ke_list += self.ar_ke_raw_list

        with open(self.ar_ke_alter_path, 'r') as file:
            reader = csv.reader(file)
            # Read the first row (assuming single row for simplicity)
            number_list = next(reader)
            # Convert the strings to floats
            self.ar_ke_alter_raw_list = [float(value)*1e6 for value in number_list] # in eV
            # self.noise_raw_list = [float(value)  for value in number_list]
        self.neutron_ar_ke_alter_list += self.ar_ke_raw_list





        with open(self.ini_path2, 'r') as file:
            reader = csv.reader(file)
            # Read the first row (assuming single row for simplicity)
            number_list = next(reader)
            # Convert the strings to floats
            self.neutron_ini_raw_list2 = [float(value)*1e6 for value in number_list]


            # the [0] is NR number and [1:] is the photon numbers
        self.neutron_ini_list2 +=  self.neutron_ini_raw_list2

        with open(self.ini_x_path2, 'r') as file:
            reader = csv.reader(file)
            # Read the first row (assuming single row for simplicity)
            number_list = next(reader)
            # Convert the strings to floats
            self.neutron_ini_raw_xlist2 = [float(value)*1e6 for value in number_list]


            # the [0] is NR number and [1:] is the photon numbers
        self.neutron_x_list2 +=  self.neutron_ini_raw_xlist2


        with open(self.ar_ke_path2, 'r') as file:
            reader = csv.reader(file)
            # Read the first row (assuming single row for simplicity)
            number_list = next(reader)
            # Convert the strings to floats
            self.ar_ke_raw_list2 = [float(value)*1e6 for value in number_list] # in eV
            # self.noise_raw_list = [float(value)  for value in number_list]
        self.neutron_ar_ke_list2 += self.ar_ke_raw_list2

        with open(self.ar_ke_alter_path2, 'r') as file:
            reader = csv.reader(file)
            # Read the first row (assuming single row for simplicity)
            number_list = next(reader)
            # Convert the strings to floats
            self.ar_ke_alter_raw_list2 = [float(value)*1e6 for value in number_list] # in eV
            # self.noise_raw_list = [float(value)  for value in number_list]
        self.neutron_ar_ke_alter_list2 += self.ar_ke_raw_list2





    def plot_neutron_spectrum_ini_lin(self):
        from matplotlib.ticker import LogLocator

        lin_bins =  np.arange(0,1e5,50)
        # counts_ini, bin_edges_ini, patches_ini = plt.hist(self.neutron_ini_list, bins=log_bins, density= True)
        # counts_ke, bin_edges_ke, patches_ke = plt.hist(self.neutron_ar_ke_list, bins=log_bins, density= True)

        counts_ini, bin_edges_ini, patches_ini = plt.hist(self.neutron_ini_list, bins=lin_bins)
        counts_ini2, bin_edges_ini2, patches_ini2 = plt.hist(self.neutron_ini_list2, bins=lin_bins)
        plt.clf()

        bins_ini = bin_edges_ini[:-1]
        bins_ini2 = bin_edges_ini2[:-1]


        Rate_ini = counts_ini*3600/self.G4_sig_time # rate in /h
        Rate_ini2 = counts_ini2 * 3600 / self.G4_sig_time  # rate in /h


        print("sum of initial rate",sum(Rate_ini)/3600,"sum of initial count", sum(counts_ini),"activity",self.rate)



        plt.plot(bins_ini, Rate_ini, drawstyle="steps-mid", label="abnormal rate")
        plt.plot(bins_ini2, Rate_ini2, drawstyle="steps-mid", label="'normal' rate")


        # plt.yscale("log")
        plt.xlabel("Energy (eV)", fontsize=16)
        plt.ylabel(r"Rate (event/hr)", fontsize=16)
        plt.title("initial energy")
        # plt.gca().xaxis.set_major_locator(LogLocator(base=10.0, numticks=15))
        plt.xlim([8e4, 1e5])
        plt.ylim([0, 35])
        plt.legend()
        plt.savefig(self.plot_path + "PN_ini_ini.pdf", bbox_inches='tight')
        plt.clf()

    def plot_neutron_spectrum_Ar_lin(self):
        from matplotlib.ticker import LogLocator

        lin_bins =  np.arange(0,1e5,50)



        counts_ke, bin_edges_ke, patches_ke = plt.hist(self.neutron_ar_ke_alter_list, bins=lin_bins)
        counts_ke2, bin_edges_ke2, patches_ke2 = plt.hist(self.neutron_ar_ke_alter_list2, bins=lin_bins)
        plt.clf()




        Rate_ke = counts_ke *3600/self.G4_sig_time  # rate in /h
        Rate_ke2 = counts_ke2 * 3600 / self.G4_sig_time  # rate in /h

        bins_ke = bin_edges_ke[:-1]
        bins_ke2 = bin_edges_ke2[:-1]





        plt.plot(bins_ke, Rate_ke, drawstyle="steps-mid", label="Abnormal rate")
        plt.plot(bins_ke2, Rate_ke2, drawstyle="steps-mid", label="'normal' rate")

        # plt.yscale("log")
        plt.xlabel("Energy (eV)", fontsize=16)
        plt.ylabel(r"Rate (event/hr)", fontsize=16)
        plt.title("first entering energy")
        # plt.gca().xaxis.set_major_locator(LogLocator(base=10.0, numticks=15))
        plt.xlim([0, 1e5])
        plt.ylim([0, 10])
        plt.legend()
        plt.savefig(self.plot_path + "PN_ini_ke_ar.pdf", bbox_inches='tight')
        plt.clf()


    def plot_neutron_spectrum_x_lin(self):
        from matplotlib.ticker import LogLocator

        lin_bins =  np.arange(-80e7,0,50)

        print(min(self.neutron_x_list),max(self.neutron_x_list))
        print(min(self.neutron_x_list2), max(self.neutron_x_list2))
        counts_ini, bin_edges_ini, patches_ini = plt.hist(self.neutron_x_list, bins=lin_bins)
        counts_ini2, bin_edges_ini2, patches_ini2 = plt.hist(self.neutron_x_list2, bins=lin_bins)

        plt.clf()



        Rate_ini = counts_ini*3600/self.G4_sig_time # rate in /h
        Rate_ini2 = counts_ini2 * 3600 / self.G4_sig_time  # rate in /h


        bins_ini = bin_edges_ini[:-1]/1e7
        bins_ini2 = bin_edges_ini2[:-1]/1e7



        plt.plot(bins_ini, Rate_ini, drawstyle="steps-mid", label="Abnorml Rate")
        plt.plot(bins_ini2, Rate_ini2, drawstyle="steps-mid", label="'Normal' Rate")


        # plt.yscale("log")
        plt.xlabel("x (cm)", fontsize=16)
        plt.ylabel(r"Rate (event/hr)", fontsize=16)
        # plt.gca().xaxis.set_major_locator(LogLocator(base=10.0, numticks=15))
        plt.xlim([-80, 0])
        # plt.ylim([1e-1, 1e4])
        plt.legend()
        plt.savefig(self.plot_path + "PN_ini_x.pdf", bbox_inches='tight')
        plt.clf()

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
    sn = SN(gamma=True)
    # test = test_csv()