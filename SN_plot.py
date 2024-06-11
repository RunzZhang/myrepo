import pandas as pd
import matplotlib.pyplot as plt
import csv
import numpy as np
class SN():
    def __init__(self):
        self.old_read_files()
        # self.new_read_files()


    def old_read_files(self):
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

        with open("/data/runzezhang/result/TN_e_sims/scatter_spectrum_CF.csv", 'r') as file:
            reader = csv.reader(file)
            # Read the first row (assuming single row for simplicity)
            number_list = next(reader)
            # Convert the strings to floats
            self.noise_raw_list = [float(value) * 10 * 0.03 * 0.2 / 1000 for value in number_list]

        # self.noise_raw_df = pd.read_csv('C:\\Users\\24230\\Downloads\\scatter_spectrum.csv', quoting=csv.QUOTE_ALL)
        # self.noise_raw_list = self.noise_raw_df.columns.to_list()
        # self.noise_raw_list = list(map(float, self.noise_raw_list))
        print(len(self.noise_raw_list))
        max_noise_photon = round(max(self.noise_raw_list))
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
        plt.hist(self.sig_raw_list,color="red",label='sig')
        plt.hist(self.noise_raw_list,color='blue',label='noise')
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
        for i in range(len(threshold_list)):
            (sig_num,noise_num)= self.prepare(i)
            signal_number_list.append(sig_num)
            noise_number_list.append(noise_num)
            if noise_num !=0:
                SN_ratio.append(sig_num/noise_num)
            else:
                SN_ratio.append(0)
        plt.plot(threshold_list,signal_number_list,color='red',legend='signal')
        plt.plot(threshold_list,noise_number_list,color='blue',legend='noise')
        plt.plot(threshold_list,SN_ratio,color='green',legend='ratio')

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