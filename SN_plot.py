import pandas as pd
import matplotlib.pyplot as plt
import csv
import numpy as np
class SN():
    def __init__(self):
        self.old_read_files()
        # self.new_read_files()


    def old_read_files(self):
        self.capture_ratio = 1.164E-3
        # self.capture_ratio = 0.121
        self.rate = 435.6 #/s
        # self.rate = 0.56
        self.G4_events= 1E6
        self.G4_sig_time=(self.G4_events / self.rate)
        with open("/data/runzezhang/result/TN_e_sims/Ar_photon_CF.csv", 'r') as file:
            reader = csv.reader(file)
            # Read the first row (assuming single row for simplicity)
            number_list = next(reader)
            # Convert the strings to floats
            self.sig_raw_list = [float(value) for value in number_list]
        # self.sig_raw_df = pd.read_csv('C:\\Users\\24230\\Downloads\\Ar_photon.csv', quoting=csv.QUOTE_ALL)
        # self.sig_raw_list = self.sig_raw_df.columns.to_list()
        #
        # self.sig_raw_list = list(map(float, self.sig_raw_list))

        print( len(self.sig_raw_list))
        self.G4_noise_time = 1E7/self.rate
        with open("/data/runzezhang/result/TN_e_sims/scatter_spectrum_CF.csv", 'r') as file:
            reader = csv.reader(file)
            # Read the first row (assuming single row for simplicity)
            number_list = next(reader)
            # Convert the strings to floats
            self.noise_raw_list = [float(value) * 10 * 0.03 * 0.2 / (1000) for value in number_list]

        # self.noise_raw_df = pd.read_csv('C:\\Users\\24230\\Downloads\\scatter_spectrum.csv', quoting=csv.QUOTE_ALL)
        # self.noise_raw_list = self.noise_raw_df.columns.to_list()
        # self.noise_raw_list = list(map(float, self.noise_raw_list))
        print(len(self.noise_raw_list))
        max_noise_photon = round(max(self.noise_raw_list))
        print("max",max_noise_photon)
        # form the threshold function
        threshold_list  = []
        for i in range(0,max_noise_photon):
            threshold_list.append(i)
        self.plot_sn(threshold_list)
        # self.hist_info()



    # def new_read_files(self):
    #     self.sig_raw_df = pd.read_csv('C:\\Users\\24230\\Downloads\\Ar_photon.csv', quoting=csv.QUOTE_ALL)
    #     self.sig_raw_list = self.sig_raw_df.columns.to_list()
    #
    #     self.sig_raw_list = list(map(float, self.sig_raw_list))
    #
    #     self.noise_raw_df = pd.read_csv('C:\\Users\\24230\\Downloads\\scatter_spectrum_zip.csv', quoting=csv.QUOTE_ALL)
    #     self.noise_raw_list = self.noise_raw_df.columns.to_list()
    #
    #     self.noise_raw_list = list(map(float, self.noise_raw_list))
    #
    #     # self.plot()
    # def supress(self):
    #     self.noise_p_raw_list = self.noise_p_raw_list[:10*len(self.sig_raw_list)]
    #     with open("C:\\Users\\24230\\Downloads\\scatter_spectrum_zip.csv", 'w', newline='') as myfile:
    #         wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    #         wr.writerow(self.noise_p_raw_list)
    #     print("len",len(self.noise_p_raw_list))
    def hist_info(self):
        sig_counts, sig_bin_edges = np.histogram(self.sig_raw_list, bins=100)
        sig_normalized_counts = sig_counts*self.capture_ratio/self.G4_sig_time
        sig_bin_centers = (sig_bin_edges[:-1] + sig_bin_edges[1:]) / 2
        plt.bar(sig_bin_centers, sig_normalized_counts, width=sig_bin_edges[1] - sig_bin_edges[0], color='red',label='signal')

        noise_counts, noise_bin_edges = np.histogram(self.noise_raw_list, bins=100)
        noise_normalized_counts = noise_counts / self.G4_noise_time
        noise_bin_centers = (noise_bin_edges[:-1] + noise_bin_edges[1:]) / 2
        plt.bar(noise_bin_centers, noise_normalized_counts, width=noise_bin_edges[1] - noise_bin_edges[0], color='blue',label='noise')

        # Set x-label and y-label with font size
        # plt.xlabel('Value', fontsize=14)
        # plt.ylabel('Frequency (normalized)', fontsize=14)
        # plt.hist(self.sig_raw_list, color="red", label='signal')
        # plt.hist(self.noise_raw_list,color='blue',label='noise')

        plt.xlabel("photon detected by SiPM #", fontsize=16)
        plt.ylabel("signal/noise rate #/s", fontsize=16)
        plt.yscale('log')
        plt.legend()
        plt.show()
    def prepare(self, threshold):
        self.sig = [value for value in self.sig_raw_list if value >= threshold]
        self.noise = [value for value in self.noise_raw_list if value >= threshold]
        sig_len = len(self.sig)
        noise_len = len(self.noise)
        return (sig_len,noise_len)
    def plot_sn(self, threshold_list):
        signal_number_list = []
        noise_number_list =[]
        SN_ratio = []
        photon_n_list = []
        print("cut",max(self.noise_raw_list))
        length = round(max(self.sig_raw_list))
        for i in range(length):
            photon_n_list.append(i)
            (sig_num,noise_num)= self.prepare(i)
            signal_number_list.append(sig_num*self.capture_ratio/self.G4_sig_time)
            noise_number_list.append(noise_num/self.G4_noise_time)
            if noise_num !=0:
                SN_ratio.append((sig_num*self.capture_ratio/self.G4_sig_time)/(noise_num/self.G4_noise_time))
            else:
                SN_ratio.append(max(SN_ratio))
        print("sig rate",max(signal_number_list))
        print("noise stat N", len(self.noise_raw_list))

        # plt.plot(photon_n_list,signal_number_list,color='red',label='signal')
        # plt.plot(photon_n_list,noise_number_list,color='blue',label='noise')
        # plt.plot(photon_n_list,SN_ratio,color='green',label='ratio')
        fig, ax1 = plt.subplots()

        # Plot dataset 1 and dataset 2 on the left y-axis
        line1, = ax1.plot(photon_n_list, signal_number_list, 'g-', label='signal')
        line2, = ax1.plot(photon_n_list, noise_number_list, 'b-', label='noise')

        # Set the labels and title
        ax1.set_xlabel('photon number threshold',fontsize = 16)
        ax1.set_ylabel('detected event number', color='black',fontsize = 16)

        # Create another y-axis that shares the same x-axis
        ax2 = ax1.twinx()

        # Plot dataset 3 on the right y-axis
        line3, = ax2.plot(photon_n_list, SN_ratio, 'r-', label='signal to noise ratio')

        # Set the label for the second y-axis
        ax2.set_ylabel('signal to noise ratio', color='black',fontsize = 16)

        lines = [line1, line2, line3]
        labels = [line.get_label() for line in lines]
        fig.legend(lines, labels, loc='upper right', bbox_to_anchor=(0.9, 0.85))
        # Show the plot
        plt.show()



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