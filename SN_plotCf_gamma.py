import pandas as pd
import matplotlib.pyplot as plt
import csv
import numpy as np
class SN():
    def __init__(self):
        # after generate new files, you need to select the capture ratio/source for different configs in read_files function.
        # then choose the correct signal/noise of with clause in read files.
        # at last change the self.name and plot_name in plot function
        self.base_path = "/data/runzezhang/result/TN_sims_D/chunked_root_files_Cf/"
        self.base_path2 = "/data/runzezhang/result/TN_sims_D/chunked_root_files_Cf/"
        self.plot_path = '/data/runzezhang/result/TN_sims_D/plot/'
        self.false_1 = "Cf_false1.csv"
        self.false_2 = "Cf_false2.csv"
        self.signal = "Cf_sig.csv"
        self.name1 = "Gamma spectrum"
        self.name2 = "Gamma despostion Rate"
        self.name = "Backgrounds"
        self.plot_name = self.name+"Cf_gamma_spectrum.pdf"
        self.signal_final_list = []
        self.noise1_final_list =[]
        self.noise2_final_list = []
        self.noise3_final_list = []
        self.noise5_final_list = []
        self.noise6_final_list = []
        self.gamma_list = []
        self.yield_rate = []



        #982 statics false 1
        for i in range(1,101):
        # for i in range(1, 2):
            self.main_body(i)
        self.plot_G()
        # self.plot_gamma()


    def main_body(self,i):
        print(i)
        self.false_1 = f"Cf_1E7_false1_part{i}.csv"
        self.false_2 = f"Cf_1E7_false2_part{i}.csv"
        self.false_3 = f"Cf_1E7_false3_part{i}.csv"
        self.false_4 = f"Cf_1E7_false4_part{i}.csv"
        self.false_5 = f"Cf_1E7_false5_part{i}.csv"
        self.false_6 = f"Cf_1E7_false6_part{i}.csv"
        self.signal = f"Cf_1E7_sig_part{i}.csv"

        self.false_1_path = self.base_path + self.false_1
        self.false_2_path = self.base_path + self.false_2
        self.false_3_path = self.base_path + self.false_3
        self.false_4_path = self.base_path + self.false_4
        self.false_5_path = self.base_path + self.false_5
        self.false_6_path = self.base_path + self.false_6
        self.signal_path = self.base_path + self.signal


        self.read_files()

        # self.read_files_s_to_N1()

# main funtion we use
    def read_files(self):



        with open(self.false_3_path, 'r') as file: # electron yielding photons per event, similar to 6, but 6 is electron energy instead of photon nums
            reader = csv.reader(file)
            # Read the first row (assuming single row for simplicity)
            number_list = next(reader)
            # Convert the strings to floats
            self.noise3_raw_list = [float(value)*1000/( 1E6 * 10 * 0.03 * 0.2) for value in number_list] # electron recoiled energy in MeV
        self.noise3_final_list  += self.noise3_raw_list
        print("F3", len(self.noise3_raw_list))



        with open(self.false_4_path, 'r') as file: # gamma energy per particle
            reader = csv.reader(file)
            # Read the first row (assuming single row for simplicity)
            number_list = next(reader)
            # Convert the strings to floats
            self.noise4_raw_list = [float(value) for value in number_list]
        self.gamma_list += self.noise4_raw_list
        print("F", len(self.noise4_raw_list))

        with open(self.false_5_path, 'r') as file: # gamma energy per event
            reader = csv.reader(file)
            # Read the first row (assuming single row for simplicity)
            number_list = next(reader)
            # Convert the strings to floats
            self.noise5_raw_list = [float(value) for value in number_list]
        self.noise5_final_list += self.noise5_raw_list
        print("F5", len(self.noise5_raw_list),self.noise5_raw_list[:10])

        with open(self.false_6_path, 'r') as file:  # electron deposit energy per event
            reader = csv.reader(file)
            # Read the first row (assuming single row for simplicity)
            number_list = next(reader)
            # Convert the strings to floats
            self.noise6_raw_list = [float(value) for value in number_list]  # electron recoiled energy in MeV
        self.noise6_final_list += self.noise6_raw_list
        print("F6", len(self.noise6_raw_list))

        print(len(self.noise3_raw_list), len(self.noise5_raw_list), len(self.noise6_raw_list))
        if len(self.noise6_raw_list)== len(self.noise5_raw_list):

            for i in range(len(self.noise5_raw_list)):
                if self.noise5_raw_list[i]==0:
                    continue
                else:
                    self.yield_rate.append(self.noise6_raw_list[i]/self.noise5_raw_list[i])
        if self.yield_rate == []: # in case error in later plot sections
            self.yield_rate.append(0)


    def plot_G(self):
        bin_num =100
        plt.hist(self.gamma_list, bins= bin_num)
        plt.xlabel("gamma energy/MeV")
        plt.ylabel("counts")
        plt.savefig(self.plot_name)


    def combine_data(self):
        # get rate vs diff threshold
        print(len(self.noise1_final_list))
        max_noise1_photon = round(max(self.noise1_final_list))
        print("max",max_noise1_photon)
        # form the threshold function
        threshold1_list  = []
        bin_size = round(max_noise1_photon/10)# if the max noise photon is too large, then we need to modity this bc of RAM
        for i in range(0,max_noise1_photon):
        # for i in range(0,round(max_noise_photon*0.1)):
            if i%bin_size==0:
                percentage = (i / max_noise1_photon) * 100
                print(f"Noise 1 Progress: {percentage:.0f}%")
            threshold1_list.append(i)

        print(len(self.noise2_final_list))
        max_noise2_photon = round(max(self.noise2_final_list))
        print("max", max_noise2_photon)
        # form the threshold function
        threshold2_list = []
        bin_size = round(
            max_noise2_photon / 10)  # if the max noise photon is too large, then we need to modity this bc of RAM
        for i in range(0, max_noise2_photon):
            # for i in range(0,round(max_noise_photon*0.1)):
            if i % bin_size == 0:
                percentage = (i / max_noise2_photon) * 100
                print(f"Noise 2 Progress: {percentage:.0f}%")
            threshold2_list.append(i)

        # plot the S/N ratio picture
        result1 = self.unit_transfer( self.noise1_final_list,threshold1_list)
        result2 = self.unit_transfer(self.noise2_final_list, threshold2_list)
        return (result1, result2)

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




    def prepare(self,noise_list, threshold): # filter the value above the threshold
        self.sig = [value for value in self.signal_final_list if value >= threshold]
        self.noise = [value for value in noise_list if value >= threshold]
        sig_len = len(self.sig)
        noise_len = len(self.noise)
        return (sig_len,noise_len)
    def unit_transfer(self, noise_list, threshold_list):
        signal_rate_list = []
        noise_rate_list =[]
        signal_num_list = []
        noise_num_list = []
        SN_ratio = []
        photon_n_list = []

        print("cut",max(noise_list))
        length = round(max(self.signal_final_list))
        point = []
        for i in range(length):
            photon_n_list.append(i)
            (sig_num,noise_num)= self.prepare(noise_list,i)
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
        print("noise stat N", len(noise_list))
        print("noise rate",max(noise_rate_list))
        print("SN",max(SN_ratio))
        print("noise uncetainty", 1.29*max(noise_rate_list)/len(noise_list))
        print("sig_stats", signal_num_list[0], "noise_stats", noise_num_list[0])

        return(signal_rate_list, photon_n_list, noise_rate_list,  SN_ratio)

        # plt.show()

    def plot_gamma(self):
        fig, (ax1, ax3) = plt.subplots(1, 2, figsize=(12, 5))  # ax1 for first plot, ax3 for second plot

        # ======== FIRST PLOT (your original one) ========
        # Plot dataset 1 and dataset 2 on the left y-axis
        line1, = ax1.hist(self.gamma_list, 'g-')
        ax1.ticklabel_format(style='sci', scilimits=(-2, 3), axis='y')
        # ax1.set_xlim([0, 600])
        # ax1.set_ylim([1e-3, 20])

        ax1.set_xlabel('Gamma sepctrum by Cf 252 inelastic scattering in Ar', fontsize=16)
        ax1.set_ylabel('Counts', color='black', fontsize=16)
        # ax1.axvline(x=200, color='black', linestyle='dotted')
        ax1.set_yscale('log')



        # Legend for first plot
        # lines_group1 = [line1, line2, line3]
        # labels_group1 = [line.get_label() for line in lines_group1]
        # ax1.legend(lines_group1, labels_group1, loc='upper right')

        ax1.set_title(self.name1, fontsize=16)

        # ======== SECOND PLOT (side-by-side) ========
        # Example plot — replace with your own data
        line4, = ax3.hist(self.yield_rate, 'g-')

        ax3.ticklabel_format(style='sci', scilimits=(-2, 3), axis='y')
        # ax3.set_xlim([0, 600])
        # ax3.set_ylim([1e-3, 20])

        ax3.set_xlabel('Yield Rate ER/Gamma each event', fontsize=16)
        ax3.set_ylabel('Counts', color='black', fontsize=16)
        ax3.axvline(x=200, color='black', linestyle='dotted')
        ax3.set_yscale('log')


        # Legend for first plot


        ax3.set_title(self.name2, fontsize=16)

        # Adjust spacing so plots don’t overlap
        plt.tight_layout()

        # Save or show
        plt.savefig(self.plot_path + self.plot_name)
    def plot_sn(self, sig1, sig2, pho1, pho2, noise1, noise2, sn1, sn2):
        fig, (ax1, ax3) = plt.subplots(1, 2, figsize=(12, 5))  # ax1 for first plot, ax3 for second plot

        # ======== FIRST PLOT (your original one) ========
        # Plot dataset 1 and dataset 2 on the left y-axis
        line1, = ax1.plot(pho1, sig1, 'g-', label='Neutron Capture Signal')
        line2, = ax1.plot(pho1, noise1, 'b-', label='Background')
        ax1.ticklabel_format(style='sci', scilimits=(-2, 3), axis='y')
        ax1.set_xlim([0, 600])
        ax1.set_ylim([1e-3, 20])

        ax1.set_xlabel('Photon Number Threshold (number)', fontsize=16)
        ax1.set_ylabel('Rate (event/hr)', color='black', fontsize=16)
        ax1.axvline(x=200, color='black', linestyle='dotted')
        ax1.set_yscale('log')

        # Create another y-axis for SNR
        ax2 = ax1.twinx()
        line3, = ax2.plot(pho1, sn1, 'r-', label='SNR')
        ax2.set_ylabel('Signal to noise ratio', color='black', fontsize=16)
        ax2.set_ylim([0.5, 2.2])

        # Legend for first plot
        lines_group1 = [line1, line2, line3]
        labels_group1 = [line.get_label() for line in lines_group1]
        ax1.legend(lines_group1, labels_group1, loc='upper right')

        ax1.set_title(self.name1, fontsize=16)

        # ======== SECOND PLOT (side-by-side) ========
        # Example plot — replace with your own data
        line4, = ax3.plot(pho2, sig2, 'g-', label='Neutron Capture Signal')
        line5, = ax3.plot(pho2, noise2, 'b-', label='Background')
        ax3.ticklabel_format(style='sci', scilimits=(-2, 3), axis='y')
        ax3.set_xlim([0, 600])
        ax3.set_ylim([1e-3, 20])

        ax3.set_xlabel('Photon Number Threshold (number)', fontsize=16)
        ax3.set_ylabel('Rate (event/hr)', color='black', fontsize=16)
        ax3.axvline(x=200, color='black', linestyle='dotted')
        ax3.set_yscale('log')

        # Create another y-axis for SNR
        ax4 = ax3.twinx()
        line6, = ax4.plot(pho2, sn2, 'r-', label='SNR')
        ax4.set_ylabel('Signal to noise ratio', color='black', fontsize=16)
        ax4.set_ylim([0.5, 2.2])

        # Legend for first plot
        lines_group2 = [line4, line5, line6]
        labels_group2 = [line.get_label() for line in lines_group2]
        ax3.legend(lines_group2, labels_group2, loc='upper right')

        ax3.set_title(self.name2, fontsize=16)

        # Adjust spacing so plots don’t overlap
        plt.tight_layout()

        # Save or show
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