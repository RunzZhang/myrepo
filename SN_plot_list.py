import pandas as pd
import matplotlib.pyplot as plt
import csv
import numpy as np
class SN():
    def __init__(self):
        # after generate new files, you need to select the capture ratio/source for different configs in read_files function.
        # then choose the correct signal/noise of with clause in read files.
        # at last change the self.name and plot_name in plot function
        self.base_path = "/data/runzezhang/result/TN_sims_D/chunked_root_files/"
        self.base_path2 = "/data/runzezhang/result/TN_sims_D/chunked_root_files/"
        self.plot_path = '/data/runzezhang/result/TN_sims_D/plot/'
        self.false_1 = "Cf_false1.csv"
        self.false_2 = "Cf_false2.csv"
        self.signal = "Cf_sig.csv"
        self.name = "Background 1 400 eV"
        self.plot_name = self.name+"1E7.png"
        self.signal_final_list = []
        self.noise_final_list =[]


        #982 statics false 1
        for i in range(1,101):
            self.main_body(i)
        self.combine_data()

    def main_body(self,i):
        print(i)
        self.false_1 = f"Cf_1E7_false1_part{i}.csv"
        self.false_2 = f"Cf_1E7_false2_part{i}.csv"
        self.signal = f"Cf_1E7_sig_part{i}.csv"

        self.false_1_path = self.base_path + self.false_1
        self.false_2_path = self.base_path + self.false_2
        self.signal_path = self.base_path + self.signal


        self.read_files()

        # self.read_files_s_to_N1()

# main funtion we use
    def read_files(self):
        self.Activity = 0.032  # source activity in mivro curie for 50 bubbles/hour
        # self.Activity = 0.0416  # source activity in mivro curie
        # self.capture_ratio = 1.164E-3 # 1125eV 1.4g/cm Ar
        # self.capture_ratio = 0.121 # 400 eV 1.4g/cm3 Ar
        # self.capture_ratio = 0.267  # 350 eV
        self.capture_ratio = 0.116  # 400 eV
        # self.capture_ratio = 6.52E-3  # 700 eV
        # self.capture_ratio = 1.158E-3  # 1125 eV
        self.rate = 435.6 #/s # CF neutron rate 9 mucurie
        # self.rate = 0.56 #AmLi neutron rate
        self.G4_events= 1E7
        self.G4_sig_time=(self.G4_events / self.rate)
        # with open(self.base_path + "Ar_photon_AmLi2.csv", 'r') as file:
        # with open(self.base_path + "Ar_photon_Cf.csv", 'r') as file:
        with open(self.signal_path, 'r') as file:
            reader = csv.reader(file)
            # Read the first row (assuming single row for simplicity)
            number_list = next(reader)
            # Convert the strings to floats
            self.sig_raw_list = [float(value) for value in number_list]
        self.signal_final_list = self.signal_final_list + self.sig_raw_list
        # self.sig_raw_df = pd.read_csv('C:\\Users\\24230\\Downloads\\Ar_photon.csv', quoting=csv.QUOTE_ALL)
        # self.sig_raw_list = self.sig_raw_df.columns.to_list()
        #
        # self.sig_raw_list = list(map(float, self.sig_raw_list))

        print("capture event number", len(self.sig_raw_list))
        # self.G4_noise_time = 1E7/self.rate
        self.G4_noise_time = self.G4_events / self.rate
        # with open("/data/runzezhang/result/TN_e_sims/scatter_spectrum_CF.csv", 'r') as file:
        # Noise 1
        # with open(self.base_path + "photon_capture_n_sing_scatterg_AmLi.csv", 'r') as file:
        # with open(self.base_path + "photon_capture_n_sing_scatterg_CF.csv", 'r') as file:
        with open(self.false_1_path, 'r') as file:
        # Noise 2
        # with open(self.base_path + "n_huge_scatterg_AmLi2.csv", 'r') as file:
        # with open(self.base_path + "n_huge_scatterg_CF2.csv", 'r') as file:
        # with open(self.false_2_path, 'r') as file:
            reader = csv.reader(file)
            # Read the first row (assuming single row for simplicity)
            number_list = next(reader)
            # Convert the strings to floats
            self.noise_raw_list = [float(value) for value in number_list]
            # self.noise_raw_list = [float(value)  for value in number_list]
        self.noise_final_list = self.noise_final_list+ self.noise_raw_list
        # self.noise_raw_df = pd.read_csv('C:\\Users\\24230\\Downloads\\scatter_spectrum.csv', quoting=csv.QUOTE_ALL)
        # self.noise_raw_list = self.noise_raw_df.columns.to_list()
        # self.noise_raw_list = list(map(float, self.noise_raw_list))
    def combine_data(self):
        plt.hist(self.noise_final_list, bins= 1000)
        plt.xlabel("photon number")
        plt.ylabel("frequency")
        plt.savefig(self.plot_path+"noise_list.png")
        print(len(self.noise_final_list))
        max_noise_photon = round(max(self.noise_final_list))
        print("max",max_noise_photon)
        # form the threshold function
        threshold_list  = []
        bin_size = round(max_noise_photon/10)# if the max noise photon is too large, then we need to modity this bc of RAM
        for i in range(0,max_noise_photon):
        # for i in range(0,round(max_noise_photon*0.1)):
            if i%bin_size==0:
                percentage = (i / max_noise_photon) * 100
                print(f"Progress: {percentage:.0f}%")
            threshold_list.append(i)
        print("ready to generate graph")
        # plot the S/N ratio picture
        self.plot_sn(threshold_list)
        self.hist_info()

    def hist_info(self):
        sig_counts, sig_bin_edges = np.histogram(self.signal_final_list, bins=100)
        sig_normalized_counts = sig_counts*self.capture_ratio/self.G4_sig_time
        sig_bin_centers = (sig_bin_edges[:-1] + sig_bin_edges[1:]) / 2
        plt.bar(sig_bin_centers, sig_normalized_counts, width=sig_bin_edges[1] - sig_bin_edges[0], color='red',label='signal')

        noise_counts, noise_bin_edges = np.histogram(self.noise_final_list, bins=100)
        noise_normalized_counts = noise_counts / self.G4_noise_time
        noise_bin_centers = (noise_bin_edges[:-1] + noise_bin_edges[1:]) / 2
        plt.bar(noise_bin_centers, noise_normalized_counts, width=noise_bin_edges[1] - noise_bin_edges[0], color='blue',label='noise')

        # Set x-label and y-label with font size
        # plt.xlabel('Value', fontsize=14)
        # plt.ylabel('Frequency (normalized)', fontsize=14)
        # plt.hist(self.signal_final_list, color="red", label='signal')
        # plt.hist(self.noise_final_list,color='blue',label='noise')

        plt.xlabel("photon detected by SiPM #", fontsize=16)
        plt.ylabel("signal/noise rate #/s", fontsize=16)
        plt.yscale('log')
        plt.legend()
        plot_name = "sn1_1E7"
        plt.savefig(self.plot_path+plot_name)
        # plt.show()

    def prepare(self, threshold): # filter the value above the threshold
        self.sig = [value for value in self.signal_final_list if value >= threshold]
        self.noise = [value for value in self.noise_final_list if value >= threshold]
        sig_len = len(self.sig)
        noise_len = len(self.noise)
        return (sig_len,noise_len)
    def plot_sn(self, threshold_list):
        signal_rate_list = []
        noise_rate_list =[]
        signal_num_list = []
        noise_num_list = []
        SN_ratio = []
        photon_n_list = []

        print("cut",max(self.noise_final_list))
        length = round(max(self.signal_final_list))
        point = []
        for i in range(length):
            photon_n_list.append(i)
            (sig_num,noise_num)= self.prepare(i)
            # change signal_number form /s to /h
            signal_rate_list.append(self.Activity*3600*sig_num*self.capture_ratio/(9*self.G4_sig_time))
            noise_rate_list.append(3600*self.Activity*noise_num/(9*self.G4_noise_time))
            signal_num_list.append(sig_num)
            noise_num_list.append(noise_num)
            if noise_num !=0:
                SN_ratio.append((sig_num*self.capture_ratio/self.G4_sig_time)/(noise_num/self.G4_noise_time))
            else:
                point.append(i)
                # print("point", point)
                # SN_ratio.append(max(SN_ratio)) # append line in the graph
                SN_ratio.append(max(SN_ratio)*1E5)  # append inf line in the graph
        for j in range(len(signal_rate_list)):
            if photon_n_list[j]>200:
                print("output",j,signal_rate_list[j],noise_rate_list[j])
                print("stat num", noise_num_list[j])
                break
        print("sig rate",max(signal_rate_list))
        if point != []:
            print("sig rate after cut", signal_rate_list[point[0]])
        print("noise stat N", len(self.noise_final_list))
        print("noise rate",max(noise_rate_list))
        print("SN",max(SN_ratio))
        print("noise uncetainty", 1.29*max(noise_rate_list)/len(self.noise_final_list))
        print("sig_stats", signal_num_list[0], "noise_stats", noise_num_list[0])

        # plt.plot(photon_n_list,signal_rate_list,color='red',label='signal')
        # plt.plot(photon_n_list,noise_rate_list,color='blue',label='noise')
        # plt.plot(photon_n_list,SN_ratio,color='green',label='ratio')
        fig, ax1 = plt.subplots()

        # Plot dataset 1 and dataset 2 on the left y-axis
        line1, = ax1.plot(photon_n_list, signal_rate_list, 'g-', label='Signal')
        line2, = ax1.plot(photon_n_list, noise_rate_list, 'b-', label='Noise')
        ax1.ticklabel_format(style='sci', scilimits=(-2, 3), axis='y')
        ax1.set_ylim([2e-5,20])

        # Set the labels and title
        ax1.set_xlabel('Photon Number Threshold (number)',fontsize = 16)
        ax1.set_ylabel('Detected Event Rate (event/hr)', color='black',fontsize = 16)
        ax1.axvline(x=200, color='black')
        ax1.set_yscale('log')

        # Create another y-axis that shares the same x-axis
        ax2 = ax1.twinx()

        # Plot dataset 3 on the right y-axis
        line3, = ax2.plot(photon_n_list, SN_ratio, 'r-', label='Signal to Noise Ratio')

        # Set the label for the second y-axis
        ax2.set_ylabel('signal to noise ratio', color='black',fontsize = 16)
        ax2.set_ylim([0.5,2.2])

        lines = [line1, line2, line3]
        # lines = [line1,  line3]
        labels = [line.get_label() for line in lines]
        fig.legend(lines, labels, loc='upper right', bbox_to_anchor=(0.9, 0.85))
        # Show the plot

        plt.title(self.name, fontsize = 16)
        plt.savefig(self.plot_path + self.plot_name)
        # plt.show()



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