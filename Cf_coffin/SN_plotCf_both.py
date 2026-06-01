import pandas as pd
import matplotlib.pyplot as plt
import csv
import numpy as np
import math, pickle
import matplotlib.ticker as ticker
import matplotlib.colors as colors
class SN():
    def __init__(self,gamma=False,full_gamma = False):
        # after generate new files, you need to select the capture ratio/source for different configs in read_files function.
        # then choose the correct signal/noise of with clause in read files.
        # at last change the self.name and plot_name in plot function
        # v2: change back to 2 backgrounds but with finer definitions
        # v4 kill duplicated NRERs

        ## change config A to B
        # replace sourcetube_A and config_A
        # self.base_path = "/data/runzezhang/result/TN_sims_D/chunked_root_files_Cf_1E7_sourcetube_B/"
        # self.base_path2 = "/data/runzezhang/result/TN_sims_D/chunked_root_files_Cf_1E7_sourcetube_B/" # for gamma path
        # pay attention to different statistics
        self.base_path = "/data/runzezhang/result/TN_sims_D/chunked_root_files_Cf_1E7_sourcetube_B/"
        self.base_path2 = "/data/runzezhang/result/TN_sims_D/chunked_root_files_Cf_1E7_sourcetube_B/"  # for different density

        # self.base_path = "/data/runzezhang/result/TN_sims_D/chunked_root_files_pn_outside_1E7/"
        # self.base_path2 = "/data/runzezhang/result/TN_sims_D/chunked_root_files_pn_1E7_outside_gamma/"  # without lead
        self.config_string = 'config_B'
        self.plot_path = '/data/runzezhang/result/TN_sims_D/plot/'
        self.TN_spectrum_path = './MC_argon_full_20250701_LSS07_2E5'
        self.false_1 = "PN_false1.csv"
        self.false_2 = "PN_false2.csv"
        self.signal = "PN_sig.csv"
        self.name1 = "Correlated ER Background"
        self.name2 = "Hard Scatter Background"
        self.name = "Backgrounds"
        self.plot_name = self.name+"PN_1E6_wt_gamma_outside.pdf"
        self.pho_threshold = 100
        self.signal_final_list = []
        self.noise1_final_list =[]
        self.noise2_final_list = []
        self.noisegamma1_final_list = []
        self.untagged_bubble_list =[]
        self.tagged_bubble_list = []
        self.neutron_ini_list = []
        self.neutron_ar_ke_list =[]
        self.neutron_ar_ke_alter_list = []
        self.gamma = gamma
        self.full_gamma = full_gamma
        self.df_phys_list = []
        self.df_geo_list = []
        self.df_energy_list = []


        #982 statics false 1
        # for i in range(1,101):
        for i in range(1, 11):
            self.main_body(i)
        # self.main_body(1)
        #check the intial neutron postions, argon volume and the intial neutron energy spectrum
        # self.check_geometry()
        #check neutron which first entering argon volum's positions and energy
        # self.neutron_spectrum_enteringLAr()
        # self.neutron_source_geometry()
        # get ssttl moderating effect and check Argon recoiled by that
        # self.coffin_phys()
        # self.source_tube_phys()
        # self.NR_spectrum_moderated_by_sstl()


        # self.NR_spectrum()
        # self.write_sims_results()
        # self.NR_multiplicity()
        # for ploting PN background tagging and SNR
        # self.untagged_bubble_rate()
        # if self.gamma:
        #     (result1, result2, resultgamma1) = self.combine_data(self.gamma)
        #     self.plot_sn_gamma_v2(result1[0], result2[0], resultgamma1[0],result1[1], result2[1], resultgamma1[1],result1[2], result2[2], resultgamma1[2], result1[3],
        #                     result2[3],resultgamma1[3])
        # else:
        #     (result1, result2) = self.combine_data(self.gamma)
        #     self.plot_sn_v2(result1[0], result2[0], result1[1], result2[1], result1[2], result2[2], result1[3],
        #                     result2[3])

        # for ploting PN spectrum and SBC detector thermalizing effect
        # self.plot_neutron_spectrum()

        # self.plot_neutron_spectrum_lin()


        #check PE density over Ar NR spectrum
        self.check_neutron_spectrum_Ar()


    def main_body(self,i):
        print(i)
        self.false_1 = f"PN_1E7_false1_part{i}.csv"
        self.false_2 = f"PN_1E7_false2_part{i}.csv"
        self.false_3 = f"PN_1E7_false3_part{i}.csv"
        self.false_gamma_1 = f"PN_gamma_1E7_false1_part{i}.csv"
        self.signal = f"PN_1E7_sig_part{i}.csv"

        self.false_1_path = self.base_path + self.false_1
        self.false_2_path = self.base_path + self.false_2
        self.false_3_path = self.base_path + self.false_3
        self.false_gamma_1_path = self.base_path2 + self.false_gamma_1
        self.signal_path = self.base_path + self.signal

        self.ini_path = self.base_path+ f"PN_1E7_ini_part{i}.csv"
        self.ar_ke_path = self.base_path + f"PN_1E7_ke_part{i}.csv"
        self.ar_ke_alter_path = self.base_path + f"PN_1E7_ke_part{i}.csv"
        self.geometry_path = self.base_path+f"PN_1E7_geo_part{i}.csv"
        self.phys_path = self.base_path + f"PN_1E7_phys_part{i}.csv"


        self.read_files()

        # self.write_sims_results()

        # self.read_files_s_to_N1()

# main funtion we use
    def read_files(self):
        self.original_Activity = 1 #
        self.Activity = 1  # unit one, the number is calculated by SBC paper
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
        self.rate = 149.952  # Cf neutron rate for SBC in 2026 Dec
        self.gamma_rate = 1.27e4 # 1.77MeV PN gamma rate/s for 5 microCurie
        self.gamma_BR = 0.0687
        # self.G4_events= 1E5
        self.G4_events = 10E5
        self.G4_events_gamma =  10E5
        self.ambient_bubble = 5 # /h

        self.T = 1e-3
        self.G4_sig_time=(self.G4_events / self.rate)
        self.G4_noise_time = self.G4_events / self.rate
        self.G4_gamma_time = self.G4_events_gamma/self.gamma_rate
        self.G4_full_gamma_time = self.G4_events_gamma*self.gamma_BR/self.gamma_rate # gamma time for whole gamma spectrum

        try:

            temp_phys_df = pd.read_csv(self.phys_path)

            self.df_phys_list.append(temp_phys_df)

            self.df_phys = pd.concat(self.df_phys_list, ignore_index=True)
        except:
            print("Fail to read source geometry")

        try:

            temp_geo_df = pd.read_csv(self.geometry_path)

            self.df_geo_list.append(temp_geo_df)

            self.df_geo = pd.concat(self.df_geo_list, ignore_index=True)
        except:
            print("Fail to read geometry")
        try:
            temp_energy_df = pd.read_csv(self.signal_path)
            # print('temp_energy_df',temp_energy_df)

            self.df_energy_list.append(temp_energy_df)

            self.df_energy = pd.concat(self.df_energy_list, ignore_index=True)
        except:
            print("failed to read energy")


        # if self.full_gamma:
        #     self.G4_gamma_time =  self.G4_full_gamma_time
        # with open(self.signal_path, 'r') as file:
        #     reader = csv.reader(file)
        #     # Read the first row (assuming single row for simplicity)
        #     number_list = next(reader)
        #     # Convert the strings to floats
        #     self.sig_raw_list = [float(value) for value in number_list]
        # self.signal_final_list = self.signal_final_list + self.sig_raw_list
        #
        # print("capture event number", len(self.sig_raw_list))
        #
        # # with open("/data/runzezhang/result/TN_e_sims/scatter_spectrum_CF.csv", 'r') as file:
        # # Noise 1,
        # with open(self.false_1_path, 'r') as file:
        #     reader = csv.reader(file)
        #     # Read the first row (assuming single row for simplicity)
        #     number_list = next(reader)
        #     # Convert the strings to floats
        #     self.noise1_raw_list = [float(value) for value in number_list[1:]]
        #     bubble_num =  number_list[0]
        #
        #
        #     # the [0] is NR number and [1:] is the photon numbers
        # self.noise1_final_list = self.noise1_final_list + self.noise1_raw_list
        # self.untagged_bubble_list.append(float(bubble_num))
        # self.tagged_bubble_list.append(len(self.noise1_raw_list))
        # # tagged number for noise1
        # # Noise 2
        # with open(self.false_2_path, 'r') as file:
        #     reader = csv.reader(file)
        #     # Read the first row (assuming single row for simplicity)
        #     number_list = next(reader)
        #     # Convert the strings to floats
        #     self.noise2_raw_list = [float(value) for value in number_list]
        #     # self.noise_raw_list = [float(value)  for value in number_list]
        # self.noise2_raw_list = list(filter(lambda x: x != 0, self.noise2_raw_list))
        # self.noise2_final_list = self.noise2_final_list + self.noise2_raw_list
        # # print(self.noise2_raw_list)
        # self.tagged_bubble_list.append(len(self.noise2_raw_list))
        # # tagged number for noise2
        #
        # with open(self.false_gamma_1_path, 'r') as file:
        #     reader = csv.reader(file)
        #     # Read the first row (assuming single row for simplicity)
        #     number_list = next(reader)
        #     # Convert the strings to floats
        #     self.noisegamma1_raw_list = [float(value) for value in number_list]
        #     # self.noise_raw_list = [float(value)  for value in number_list]
        # self.noisegamma1_final_list = self.noisegamma1_final_list + self.noisegamma1_raw_list
        #
        # # Initial amli spectrm
        # with open(self.ini_path, 'r') as file:
        #     reader = csv.reader(file)
        #     # Read the first row (assuming single row for simplicity)
        #     number_list = next(reader)
        #     # Convert the strings to floats
        #     self.neutron_ini_raw_list = [float(value)*1e6 for value in number_list]
        #
        #
        #     # the [0] is NR number and [1:] is the photon numbers
        # self.neutron_ini_list +=  self.neutron_ini_raw_list
        #
        # with open(self.ar_ke_path, 'r') as file:
        #     reader = csv.reader(file)
        #     # Read the first row (assuming single row for simplicity)
        #     number_list = next(reader)
        #     # Convert the strings to floats
        #     self.ar_ke_raw_list = [float(value)*1e6 for value in number_list] # in eV
        #     # self.noise_raw_list = [float(value)  for value in number_list]
        # self.neutron_ar_ke_list += self.ar_ke_raw_list
        #
        # with open(self.ar_ke_alter_path, 'r') as file:
        #     reader = csv.reader(file)
        #     # Read the first row (assuming single row for simplicity)
        #     number_list = next(reader)
        #     # Convert the strings to floats
        #     self.ar_ke_alter_raw_list = [float(value)*1e6 for value in number_list] # in eV
        #     # self.noise_raw_list = [float(value)  for value in number_list]
        # self.neutron_ar_ke_alter_list += self.ar_ke_raw_list


    def combine_data(self,gamma=False):
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


        # form the threshold function
        threshold2_list = []
        if len(self.noise2_final_list) !=0:
            print("noise2", len(self.noise2_final_list), self.noise2_final_list)
            max_noise2_photon = round(max(self.noise2_final_list))
            print("max", max_noise2_photon)
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
        if gamma==True: # turn on gamma read function

            print(len(self.noisegamma1_final_list))
            max_noisegamma1_photon = round(max(self.noisegamma1_final_list))
            print("max", max_noisegamma1_photon)
            # form the threshold function
            threshold_noisegamma1_list = []
            bin_size = round(
                max_noisegamma1_photon / 10)  # if the max noise photon is too large, then we need to modity this bc of RAM
            for i in range(0, max_noisegamma1_photon):
                # for i in range(0,round(max_noise_photon*0.1)):
                if i % bin_size == 0:
                    percentage = (i / max_noisegamma1_photon) * 100
                    print(f"noisegamma1  Progress: {percentage:.0f}%")
                threshold_noisegamma1_list.append(i)
            resultgamma1 = self.gamma_unit_transfer(self.noisegamma1_final_list, threshold_noisegamma1_list)
            return (result1, result2, resultgamma1)
        else:
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
        # plt.show()

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


        length = round(max(self.signal_final_list))
        point = [] # threshold cut?
        if len(noise_list)==0:
            noise_list.append(0)
            print("cut", max(noise_list))

        for i in range(length):
            photon_n_list.append(i)
            (sig_num,noise_num)= self.prepare(noise_list,i)
            # if max(noise_list)==0:
            #     print((sig_num,noise_num))
            # change signal_number form /s to /h
            signal_rate_list.append(self.Activity*3600*sig_num*self.capture_ratio/(self.original_Activity*self.G4_sig_time))
            noise_rate_list.append(3600*self.Activity*noise_num/(self.original_Activity*self.G4_noise_time))
            signal_num_list.append(sig_num)
            noise_num_list.append(noise_num)
            if sig_num ==0:
                SN_ratio.append(0)
            else:
                if noise_num !=0:
                    SN_ratio.append((sig_num*self.capture_ratio/self.G4_sig_time)/(noise_num/self.G4_noise_time))
                else:

                    point.append(i)
                # print("point", point)
                #     SN_ratio.append(max(SN_ratio)) # append line in the graph
                    SN_ratio.append(max(SN_ratio)*1E3)  # append inf line in the graph
        for j in range(len(signal_rate_list)):
            if photon_n_list[j]>self.pho_threshold:
                print("output",j,signal_rate_list[j],noise_rate_list[j])
                print("stat num", noise_num_list[j])
                break
        print("sig rate",max(signal_rate_list))
        if point != []:
            print("sig rate after cut",point ,signal_rate_list[point[0]])
        print("noise stat N", len(noise_list))
        print("noise rate",max(noise_rate_list), noise_rate_list[:10])
        self.gamma_only_num = max(noise_rate_list)

        print("SN",max(SN_ratio),SN_ratio[:20])
        print("noise uncetainty", 1.29*max(noise_rate_list)/len(noise_list))
        print("sig_stats", signal_num_list[0], "noise_stats", noise_num_list[0])

        return(signal_rate_list, photon_n_list, noise_rate_list,  SN_ratio)
    def gamma_unit_transfer(self, noise_list, threshold_list):
        signal_rate_list = []
        noise_rate_list =[]
        noise_rate_gamma_only_list =[]
        signal_num_list = []
        noise_num_list = []
        SN_ratio = []
        photon_n_list = []

        print("cut",max(noise_list))
        length = round(max(self.signal_final_list))
        point = [] # threshold cut?
        for i in range(length):
            photon_n_list.append(i)
            (sig_num,noise_num)= self.prepare(noise_list,i)
            # change signal_number form /s to /h
            signal_rate_list.append(self.Activity*3600*sig_num*self.capture_ratio/(self.original_Activity*self.G4_sig_time))
            noise_rate_list.append(self.T*self.Activity*noise_num*self.summed_bubble_only_rate/(self.original_Activity*self.G4_gamma_time))
            noise_rate_gamma_only_list.append(self.Activity*noise_num/(self.original_Activity*self.G4_gamma_time))
            # for gamma T(s)*R_gamma_photon(/s)*R_bubble(/h)
            signal_num_list.append(sig_num)
            noise_num_list.append(noise_num)
            if sig_num ==0:
                SN_ratio.append(0)
            else:
                if noise_num !=0:
                    SN_ratio.append((3600*sig_num*self.capture_ratio/self.G4_sig_time)/(noise_num*self.T*self.summed_bubble_only_rate/self.G4_gamma_time))
                else:

                    point.append(i)
                # print("point", point)
                #     SN_ratio.append(max(SN_ratio)) # append line in the graph
                    SN_ratio.append(max(SN_ratio)*1E3)  # append inf line in the graph
        for j in range(len(signal_rate_list)):
            if photon_n_list[j]>self.pho_threshold:
                print("output",j,signal_rate_list[j],noise_rate_list[j])
                print("stat num", noise_num_list[j])
                break
        print("sig rate",max(signal_rate_list))
        if point != []:
            print("sig rate after cut",point ,signal_rate_list[point[0]])
        print("noise stat N", len(noise_list))
        print("gamma accidental noise rate",max(noise_rate_list), noise_rate_list[:10])
        print("gamma rate only rate /s", max(noise_rate_gamma_only_list))

        print("SN",max(SN_ratio),SN_ratio[:20])
        print("noise uncetainty", 1.29*max(noise_rate_list)/len(noise_list))
        print("sig_stats", signal_num_list[0], "noise_stats", noise_num_list[0])

        return(signal_rate_list, photon_n_list, noise_rate_list,  SN_ratio)
    def untagged_bubble_rate(self):
        summed_bubble_num = sum(self.untagged_bubble_list)
        # considering background bubble
        summed_tagged_bubble_num = sum(self.tagged_bubble_list)

        self.summed_bubble_rate = self.Activity * 3600 * summed_bubble_num / (self.original_Activity* self.G4_sig_time)+self.ambient_bubble
        self.summed_bubble_only_rate = self.Activity * 3600 * (summed_bubble_num-summed_tagged_bubble_num) / (self.original_Activity* self.G4_sig_time)+self.ambient_bubble

        # untagged total bubble rate in /h
        # background bubble rate before tagging
        # ambient + ER correlated + huge scatter + bubble NR> 100 eV before tagging
        print("untagged total bubble rate /h (non-signal bubble rate/h)", self.summed_bubble_rate)
        # rate that scintilation is absent. Used for later coincidence gamma tagging
        # ambient + bubble only NR > 100eV before tagging
        print("bubble only event/h",self.summed_bubble_only_rate)

    def plot_sn(self, sig1, sig2,pho1, pho2, noise1, noise2,  sn1, sn2):
        fig, (ax1, ax3) = plt.subplots(1, 2, figsize=(12, 5))  # ax1 for first plot, ax3 for second plot
        x_range = [-1,430]
        left_axis_range=[3.7e-3,3e2]
        right_axis_range=[1e-1,1e3]

        # left_axis_range = [1e-1, 3e4]
        # right_axis_range = [1e-1, 1e5]

        # ======== FIRST PLOT (your original one) ========
        # Plot dataset 1 and dataset 2 on the left y-axis
        line1, = ax1.plot(pho1, sig1, 'g-', label='Neutron Capture Signal')
        line2, = ax1.plot(pho1, noise1, 'b-', label='Background')
        ax1.ticklabel_format(style='sci', scilimits=(-2, 3), axis='y')
        ax1.set_xlim(x_range)
        ax1.set_ylim(left_axis_range)
        # print("pho",pho1)
        # print("noise1", noise1)

        ax1.set_xlabel('Photon Number Threshold (number)', fontsize=16)
        ax1.set_ylabel('Rate (event/hr)', color='black', fontsize=16)
        ax1.axvline(x=self.pho_threshold, color='black', linestyle='dotted')
        ax1.set_yscale('log')
        # ax1.set_aspect('equal', adjustable="datalim")


        # Create another y-axis for SNR
        ax2 = ax1.twinx()
        line3, = ax2.plot(pho1, sn1, 'r-', label='SNR')
        # print(pho1[300:],"\n",noise1[300:],"\n",sn1[300:])
        ax2.set_ylabel('Signal to noise ratio', color='black', fontsize=16)
        # ax2.set_ylim([0, 120])
        ax2.set_ylim(right_axis_range)
        ax2.set_yscale('log')

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
        ax3.set_xlim(x_range)
        ax3.set_ylim(left_axis_range)

        ax3.set_xlabel('Photon Number Threshold (number)', fontsize=16)
        ax3.set_ylabel('Rate (event/hr)', color='black', fontsize=16)
        ax3.axvline(x=self.pho_threshold, color='black', linestyle='dotted')
        ax3.set_yscale('log')
        # ax3.set_aspect('equal')

        # Create another y-axis for SNR
        ax4 = ax3.twinx()
        line6, = ax4.plot(pho2, sn2, 'r-', label='SNR')
        ax4.set_ylabel('Signal to noise ratio', color='black', fontsize=16)
        # ax4.set_ylim([0,120])
        ax4.set_ylim(right_axis_range)
        ax4.set_yscale('log')

        # Legend for first plot
        lines_group2 = [line4, line5, line6]
        labels_group2 = [line.get_label() for line in lines_group2]
        ax3.legend(lines_group2, labels_group2, loc='upper right')

        ax3.set_title(self.name2, fontsize=16)


        # Adjust spacing so plots don’t overlap
        # fig.set_size_inches(20, 6)
        plt.tight_layout()

        # Save or show
        plt.savefig(self.plot_path + self.plot_name)
        # plt.show()

    def plot_sn_v2(self, sig1, sig2,pho1, pho2, noise1, noise2,  sn1, sn2):
        # plot sn in same graph and no SNR
        fig, ax1 = plt.subplots(1, 1, figsize=(6, 5))  # ax1 for first plot, ax3 for second plot
        x_range = [0,430]
        left_axis_range=[1e-2,2e1]

        # left_axis_range = [1e-1, 3e4]
        # right_axis_range = [1e-1, 1e5]

        # ======== FIRST PLOT (your original one) ========
        # Plot dataset 1 and dataset 2 on the left y-axis
        line1, = ax1.plot(pho1, sig1, 'r-', label='Neutron Capture Signal')
        line2, = ax1.plot(pho1, noise1, 'b-', label='Correlated ER Background')
        line5, = ax1.plot(pho2, noise2, 'g-', label='Hard Scatter Background')
        ax1.ticklabel_format(style='sci', scilimits=(-2, 3), axis='y')
        ax1.set_xlim(x_range)
        ax1.set_ylim(left_axis_range)
        # print("pho",pho1)
        # print("noise1", noise1)

        ax1.set_xlabel('Photon Number Threshold (number)', fontsize=16)
        ax1.set_ylabel('Rate (event/hr)', color='black', fontsize=16)
        ax1.axvline(x=self.pho_threshold, color='black', linestyle='dotted')
        ax1.set_yscale('log')
        # ax1.set_aspect('equal', adjustable="datalim")
        # Legend for first plot
        lines_group1 = [line1, line2, line5]
        labels_group1 = [line.get_label() for line in lines_group1]
        ax1.legend(lines_group1, labels_group1, loc='upper right')



        # Adjust spacing so plots don’t overlap
        # fig.set_size_inches(20, 6)
        plt.tight_layout()

        # Save or show
        plt.savefig(self.plot_path + self.plot_name)
        # plt.show()
    def plot_sn_gamma_v2(self, sig1, sig2,sig3, pho1, pho2, pho3,noise1, noise2, noise3, sn1, sn2,sn3):
        # plot sn in same graph and no SNR
        fig, ax1 = plt.subplots(1, 1, figsize=(6, 5))  # ax1 for first plot, ax3 for second plot
        x_range = [0,430]
        left_axis_range=[1e-3,5e0]
        print("noise2",noise2)
        # left_axis_range = [1e-1, 3e4]
        # right_axis_range = [1e-1, 1e5]

        # ======== FIRST PLOT (your original one) ========
        # Plot dataset 1 and dataset 2 on the left y-axis
        line1, = ax1.plot(pho1, sig1, 'r-', label='Neutron Capture Signal')
        line2, = ax1.plot(pho1, noise1, 'b-', label='Correlated ER Background')
        line5, = ax1.plot(pho2, noise2, 'g-', label='Hard Scatter Background')
        line6, = ax1.plot(pho3, noise3, color='brown', linestyle='-', label='Uncorrelated Gamma Background')
        ax1.ticklabel_format(style='sci', scilimits=(-2, 3), axis='y')
        ax1.set_xlim(x_range)
        ax1.set_ylim(left_axis_range)
        # print("pho",pho1)
        # print("noise1", noise1)

        ax1.set_xlabel('Photon Number Threshold (number)', fontsize=16)
        ax1.set_ylabel('Rate (event/hr)', color='black', fontsize=16)
        ax1.axvline(x=self.pho_threshold, color='black', linestyle='dotted')
        ax1.set_yscale('log')
        # ax1.set_aspect('equal', adjustable="datalim")
        # Legend for first plot
        lines_group1 = [line1, line2, line5, line6]
        labels_group1 = [line.get_label() for line in lines_group1]
        ax1.legend(lines_group1, labels_group1, loc='upper right')



        # Adjust spacing so plots don’t overlap
        # fig.set_size_inches(20, 6)
        plt.tight_layout()

        # Save or show
        print("path",self.plot_path + self.plot_name)
        plt.savefig(self.plot_path + self.plot_name)
        # plt.show()

    def plot_neutron_spectrum(self):
        from matplotlib.ticker import LogLocator
        ene_list, possibility_list = self.read_original_spectrum()
        ene_list = [value * 1e6 for value in ene_list]
        possibility_list  = [ value * self.rate for value in possibility_list]
        log_bins =  np.logspace(-3,7,100)
        # counts_ini, bin_edges_ini, patches_ini = plt.hist(self.neutron_ini_list, bins=log_bins, density= True)
        # counts_ke, bin_edges_ke, patches_ke = plt.hist(self.neutron_ar_ke_list, bins=log_bins, density= True)

        counts_ini, bin_edges_ini, patches_ini = plt.hist(self.neutron_ini_list, bins=log_bins)
        counts_ke, bin_edges_ke, patches_ke = plt.hist(self.neutron_ar_ke_list, bins=log_bins)
        counts_ke_alter, bin_edges_ke_alter, patches_ke_alter = plt.hist(self.neutron_ar_ke_alter_list, bins=log_bins)
        plt.clf()

        point = 0
        for i in range(len(counts_ini)):
            if counts_ini[i] != 0:
                point = i
                break
                # find 1st none zero counts in initial energy spectrum
        normalized_ini  = 0
        normalized_ke =0


        print("norma fact", normalized_ini, normalized_ke)
        escaping_ratio = len(self.neutron_ar_ke_list) / len(self.neutron_ini_list)

        Rate_ini = counts_ini*3600/self.G4_sig_time # rate in /h

        print("escapitng ratio",escaping_ratio,"sum of initial rate",sum(Rate_ini)/3600,"sum of initial count", sum(counts_ini),"activity",self.rate)
        Rate_ke = counts_ke *3600/self.G4_sig_time  # rate in /h
        Rate_ke_alter = counts_ke_alter * 3600 / self.G4_sig_time
        bin_factor = bin_edges_ke[1]/bin_edges_ke[0]
        bins_ini = bin_edges_ini[:-1]*bin_factor**0.5
        # bins_ini = bin_edges_ini[:-1] * bin_factor
        # print("bin edges", bin_edges_ini)
        bins_ke = bin_edges_ke[:-1]*bin_factor**0.5

        print("ini", bins_ini[:10])
        print("ke",
        bins_ke[:10])
        print("before density bins_ini, rate ini, counts ini,counts ke", Rate_ini[point:],"\n", Rate_ke[point:], "\n",
              counts_ini[point:], "\n", counts_ke[point:],"\n", bins_ini[point:],"\n", bins_ke[point:])  # print no-zero first bins

        # make the y value /h/eV
        for i in range(len(Rate_ini)):
            Rate_ini[i] = Rate_ini[i]/(bin_edges_ini[i+1]-bin_edges_ini[i])
            Rate_ke[i] = Rate_ke[i] / (bin_edges_ke[i + 1] - bin_edges_ke[i])


        plt.plot(bins_ini, Rate_ini, drawstyle="steps-mid", label="LZ PN Escaping Neutron")

        print("after density bins_ini, rate ini, counts ini,counts ke",Rate_ini[point:],"\n", Rate_ke[point:],"\n", counts_ini[point:],"\n", counts_ke[point:],"\n", bins_ini[point:],"\n", bins_ke[point:]) # print no-zero first bins
        plt.plot(bins_ke, Rate_ke, drawstyle="steps-mid", label="Neutron Entering Liquid Argon")
        # plt.plot(bins_ke, Rate_ke_alter, drawstyle="steps-mid", label="First Enter LAr Neutron lin")
        # plt.plot(ene_list, possibility_list, drawstyle="steps-mid", label="Original Spectrum dat")
        # plt.plot(bins_ini, Rate_ini, label="LZ PN Escaping Neutron")
        # plt.plot(bins_ke, Rate_ke,  label="First Enter LAr Neutron")
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("Energy (eV)", fontsize=16)
        plt.ylabel(r"Rate (event/hr/eV)", fontsize=16)
        plt.gca().xaxis.set_major_locator(LogLocator(base=10.0, numticks=15))
        plt.xlim([1e-3, 1e7])
        # plt.ylim([1e-1, 1e4])
        plt.legend()
        plt.savefig(self.plot_path + "PN_specturm_log.pdf", bbox_inches='tight')

    def plot_neutron_spectrum_lin(self):
        from matplotlib.ticker import LogLocator

        lin_bins =  np.arange(0,1e5,50)
        # counts_ini, bin_edges_ini, patches_ini = plt.hist(self.neutron_ini_list, bins=log_bins, density= True)
        # counts_ke, bin_edges_ke, patches_ke = plt.hist(self.neutron_ar_ke_list, bins=log_bins, density= True)

        counts_ini, bin_edges_ini, patches_ini = plt.hist(self.neutron_ini_list, bins=lin_bins)
        counts_ke, bin_edges_ke, patches_ke = plt.hist(self.neutron_ar_ke_list, bins=lin_bins)
        counts_ke_alter, bin_edges_ke_alter, patches_ke_alter = plt.hist(self.neutron_ar_ke_alter_list, bins=lin_bins)
        plt.clf()

        point = 0
        for i in range(len(counts_ini)):
            if counts_ini[i] != 0:
                point = i
                break
                # find 1st none zero counts in initial energy spectrum
        normalized_ini  = 0
        normalized_ke =0


        print("norma fact", normalized_ini, normalized_ke)
        escaping_ratio = len(self.neutron_ar_ke_list) / len(self.neutron_ini_list)

        Rate_ini = counts_ini*3600/self.G4_sig_time # rate in /h

        print("escapitng ratio",escaping_ratio,"sum of initial rate",sum(Rate_ini)/3600,"sum of initial count", sum(counts_ini),"activity",self.rate)
        Rate_ke = counts_ke *3600/self.G4_sig_time  # rate in /h
        Rate_ke_alter = counts_ke_alter * 3600 / self.G4_sig_time
        bin_factor = bin_edges_ke[1]-bin_edges_ke[0]
        bins_ini = bin_edges_ini[:-1]+bin_factor/2
        # bins_ini = bin_edges_ini[:-1] * bin_factor
        # print("bin edges", bin_edges_ini)
        bins_ke = bin_edges_ke[:-1]+bin_factor/2

        print("ini", bins_ini[:10])
        print("ke",
        bins_ke[:10])
        # print("before density bins_ini, rate ini, counts ini,counts ke", Rate_ini[point:],"\n", Rate_ke[point:], "\n",
        #       counts_ini[point:], "\n", counts_ke[point:],"\n", bins_ini[point:],"\n", bins_ke[point:])  # print no-zero first bins

        # make the y value /h/eV
        # for i in range(len(Rate_ini)):
        #     Rate_ini[i] = Rate_ini[i]/(bin_edges_ini[i+1]-bin_edges_ini[i])
        #     Rate_ke[i] = Rate_ke[i] / (bin_edges_ke[i + 1] - bin_edges_ke[i])


        plt.plot(bins_ini, Rate_ini, drawstyle="steps-mid", label="LZ PN Escaping Neutron")

        # print("after density bins_ini, rate ini, counts ini,counts ke",Rate_ini[point:],"\n", Rate_ke[point:],"\n", counts_ini[point:],"\n", counts_ke[point:],"\n", bins_ini[point:],"\n", bins_ke[point:]) # print no-zero first bins
        # plt.plot(bins_ke, Rate_ke, drawstyle="steps-mid", label="First Enter LAr Neutron log")
        plt.plot(bins_ke, Rate_ke_alter, drawstyle="steps-mid", label="First Enter LAr Neutron lin")
        # plt.plot(ene_list, possibility_list, drawstyle="steps-mid", label="Original Spectrum dat")
        plt.plot(bins_ini, Rate_ini, label="LZ PN Escaping Neutron")
        # plt.plot(bins_ke, Rate_ke,  label="First Enter LAr Neutron")
        # plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("Energy (eV)", fontsize=16)
        plt.ylabel(r"Rate (event/hr)", fontsize=16)
        # plt.gca().xaxis.set_major_locator(LogLocator(base=10.0, numticks=15))
        plt.xlim([0, 1e5])
        # plt.ylim([1e-1, 1e4])
        plt.legend()
        plt.savefig(self.plot_path + "PN_specturm_lin_1E6.pdf", bbox_inches='tight')
    def check_geometry(self):

        self.coffin = self.df_geo[self.df_geo["Volume"]=="cf_active_phys"]
        self.coffin["PreKinetic/keV"] = self.coffin["PreKinetic/MeV"]*1000
        self.argon = self.df_geo[self.df_geo["Volume"]=="LAr_phys"]
        self.df_geo["R/mm"]=  np.sqrt(self.df_geo["X/mm"]**2+self.df_geo["Y/mm"]**2 )
        self.coffin["R/mm"] = np.sqrt(self.coffin["X/mm"] ** 2 + self.coffin["Y/mm"] ** 2)


        ffig, ax = plt.subplots(1,2,figsize=(12,4))

        sc=ax[0].hist2d(self.df_geo["R/mm"],self.df_geo["Z/mm"],bins=50,
        cmap="plasma",norm="log",alpha=0.7)

        # ax.plot([189.95,189.95, 0.8485], [0,663.22, 714.03], color="red")
        # ax.plot([114.98,114.98, 0.75575], [0,587.01, 617.78], color="blue")
        # ax.plot([99.01,99.01, 4.34], [0,366.49, 399.82], color="red")

        ax[0].set_xlabel("R [mm]")
        ax[0].set_ylabel("Z [mm]")
        # ax[0].set_xlim(0,400)
        # ax[0].set_ylim(-100,800)
        cbar = plt.colorbar(sc[3], ax=ax[0])
        cbar.set_label("Counts(log)")

        ax[1].hist(self.coffin["PreKinetic/keV"], bins=50, alpha=0.7)

        # ax.plot([189.95,189.95, 0.8485], [0,663.22, 714.03], color="red")
        # ax.plot([114.98,114.98, 0.75575], [0,587.01, 617.78], color="blue")
        # ax.plot([99.01,99.01, 4.34], [0,366.49, 399.82], color="red")

        ax[1].set_xlabel("Energy [keV]")
        ax[1].set_ylabel("Counts")
        print(self.coffin["PreKinetic/keV"])
        plt.savefig(self.plot_path+f"Cf_1E7_position_density_{self.config_string}.pdf")
        print(self.plot_path)

    def source_tube_phys(self):

        rate_factor = 1000 * self.rate * self.Activity / (self.original_Activity * self.G4_events)


        self.active = self.df_phys[(self.df_phys["Volume"] == "cf_active_phys")&(self.df_phys["Step ID"] == 1)]

        # find center x and center Y locations
        center_x= self.active["X/mm"].mean()
        center_y = self.active["Y/mm"].mean()

        self.leaving_source = self.df_phys[(self.df_phys["name"] == "neutron") & (self.df_phys["Volume"] == "cf_source_phys")]
        self.leaving_source = self.leaving_source.loc[self.leaving_source.groupby('Event')['Step ID'].idxmax()]



        self.active["PreKinetic/keV"] = self.active["PreKinetic/MeV"] * 1000

        print("total runs", len(self.active['Event'].tolist()))

        self.leaving_source["PreKinetic/keV"] = self.leaving_source["PreKinetic/MeV"] * 1000
        self.leaving_source["PostKinetic/keV"] = self.leaving_source["PostKinetic/MeV"] * 1000
        # neutron leaving source moderated by the ssteel to below 200keV
        self.low_e_n_list = self.leaving_source[self.leaving_source["PostKinetic/keV"] < 200]["Event"].tolist()



        self.leaving_source_R = self.leaving_source
        self.leaving_source_R["R/mm"] = np.sqrt((self.leaving_source_R["X/mm"]-center_x) ** 2 + (self.leaving_source_R["Y/mm"]-center_y) ** 2)
        # slice and add legend



        ffig, ax = plt.subplots(1,2,figsize=(16,4))
        # X and Y
        sc0=ax[0].hist2d(self.leaving_source_R["R/mm"],self.leaving_source_R["Z/mm"],bins=50,
        cmap="plasma",norm="log",alpha=0.7)

        # ax.plot([189.95,189.95, 0.8485], [0,663.22, 714.03], color="red")
        # ax.plot([114.98,114.98, 0.75575], [0,587.01, 617.78], color="blue")
        # ax.plot([99.01,99.01, 4.34], [0,366.49, 399.82], color="red")

        ax[0].set_xlabel("R- $R_{center}$ [mm]")
        ax[0].set_ylabel("Z [mm]")
        # ax[0].set_xlim(-77,-75)
        # ax[0].set_ylim(-690,-680)
        cbar0 = plt.colorbar(sc0[3], ax=ax[0])
        cbar0.set_label("Counts(log)")



        ax[1].hist(self.active["PreKinetic/keV"], bins=50, alpha=0.7, label="active_region")
        ax[1].hist(self.leaving_source["PostKinetic/keV"], bins=50, alpha=0.7, label="leaving source tube")

        # ax.plot([189.95,189.95, 0.8485], [0,663.22, 714.03], color="red")
        # ax.plot([114.98,114.98, 0.75575], [0,587.01, 617.78], color="blue")
        # ax.plot([99.01,99.01, 4.34], [0,366.49, 399.82], color="red")

        ax[1].set_xlabel("Energy [keV]")
        ax[1].set_ylabel("Counts")
        ax[1].legend()

        leaving_neutron = np.histogram(self.leaving_source["PostKinetic/keV"], bins=50)
        print(leaving_neutron[0],"\n",leaving_neutron[1])
        # 200 keV
        plt.savefig(self.plot_path + f"Cf_1E7_sstl_phys_effect_{self.config_string}.pdf")
        print(self.plot_path+ f"Cf_1E7_sstl_phys_effect_{self.config_string}.pdf")


    def coffin_phys(self):

        rate_factor = 1000 * self.rate * self.Activity / (self.original_Activity * self.G4_events)

        coffin_volume_list = ["block1_phys", "block2_phys", "block3_phys",
                       "block4_phys", "block5_phys", "block6_phys", "block7_phys", "block8_phys", "block9_phys",
                       "block10_phys", "block11_phys"]
        self.active = self.df_phys[(self.df_phys["Volume"] == "cf_active_phys")&(self.df_phys["Step ID"] == 1)]

        # find center x and center Y locations
        center_x= self.active["X/mm"].mean()
        center_y = self.active["Y/mm"].mean()

        # find neutron entering argon volume event list
        lar_event_list = self.df_energy["Event"].tolist()



        # self.leaving_coffin = self.df_phys[(self.df_phys["name"] == "neutron") & (self.df_phys["Volume"].isin(coffin_volume_list))]
        # also requires neutron entering argon
        self.leaving_coffin = self.df_phys[
            (self.df_phys["name"] == "neutron") & (self.df_phys["Volume"].isin(coffin_volume_list))&(self.df_phys["Event"].isin(lar_event_list))]

        self.leaving_coffin = self.leaving_coffin.loc[self.leaving_coffin.groupby('Event')['Step ID'].idxmax()]



        self.active["PreKinetic/keV"] = self.active["PreKinetic/MeV"] * 1000

        print("total runs", len(self.active['Event'].tolist()))

        self.leaving_coffin["PreKinetic/keV"] = self.leaving_coffin["PreKinetic/MeV"] * 1000
        self.leaving_coffin["PostKinetic/keV"] = self.leaving_coffin["PostKinetic/MeV"] * 1000
        # neutron leaving source moderated by the ssteel to below 200keV
        self.low_e_n_coffin_list = self.leaving_coffin[self.leaving_coffin["PostKinetic/keV"] < 200]["Event"].tolist()



        self.leaving_coffin_R = self.leaving_coffin
        self.leaving_coffin_R["R/mm"] = np.sqrt((self.leaving_coffin_R["X/mm"]-center_x) ** 2 + (self.leaving_coffin_R["Y/mm"]-center_y) ** 2)
        # slice and add legend

        print(self.leaving_coffin_R["Y/mm"].min(), self.leaving_coffin_R["Y/mm"].max())
        self.leaving_coffin_Y_slice = self.leaving_coffin_R[(self.leaving_coffin_R["Y/mm"]>-720)&(self.leaving_coffin_R["Y/mm"]<-670)]
        self.leaving_coffin_X_slice = self.leaving_coffin_R[self.leaving_coffin_R["X/mm"].between(-96, -66)]

        if self.leaving_coffin_Y_slice.empty:
            print("Y slice empty")
            return 0
        ffig, ax = plt.subplots(1,3,figsize=(24,4))
        # X and Y
        sc0=ax[0].hist2d(self.leaving_coffin_R["X/mm"],self.leaving_coffin_R["Z/mm"],bins=50,
        cmap="plasma",norm="log",alpha=0.7)
        # sc0 = ax[0].hist2d(self.leaving_coffin_Y_slice["X/mm"], self.leaving_coffin_Y_slice["Z/mm"], bins=50,
        #                    cmap="plasma", alpha=0.7,norm=colors.LogNorm(vmin=1))

        # ax.plot([189.95,189.95, 0.8485], [0,663.22, 714.03], color="red")
        # ax.plot([114.98,114.98, 0.75575], [0,587.01, 617.78], color="blue")
        # ax.plot([99.01,99.01, 4.34], [0,366.49, 399.82], color="red")

        ax[0].set_xlabel("X [mm]")
        ax[0].set_ylabel("Z [mm]")
        # ax[0].set_xlim(-77,-75)
        # ax[0].set_ylim(-690,-680)
        cbar0 = plt.colorbar(sc0[3], ax=ax[0])
        cbar0.set_label("Counts(log)")

        sc1 = ax[1].hist2d(self.leaving_coffin_R["Y/mm"], self.leaving_coffin_R["Z/mm"], bins=50,
                           cmap="plasma", norm="log", alpha=0.7)
        # sc1 = ax[1].hist2d(self.leaving_coffin_X_slice["Y/mm"], self.leaving_coffin_X_slice["Z/mm"], bins=50,
        #                    cmap="plasma", alpha=0.7,norm=colors.LogNorm(vmin=1))

        # ax.plot([189.95,189.95, 0.8485], [0,663.22, 714.03], color="red")
        # ax.plot([114.98,114.98, 0.75575], [0,587.01, 617.78], color="blue")
        # ax.plot([99.01,99.01, 4.34], [0,366.49, 399.82], color="red")

        ax[1].set_xlabel("Y [mm]")
        ax[1].set_ylabel("Z [mm]")
        # ax[0].set_xlim(-77,-75)
        # ax[0].set_ylim(-690,-680)
        cbar1 = plt.colorbar(sc1[3], ax=ax[1])
        cbar1.set_label("Counts(log)")



        ax[2].hist(self.active["PreKinetic/keV"], bins=50, alpha=0.7, label="active_region")
        ax[2].hist(self.leaving_coffin["PostKinetic/keV"], bins=50, alpha=0.7, label="leaving coffin")

        # ax.plot([189.95,189.95, 0.8485], [0,663.22, 714.03], color="red")
        # ax.plot([114.98,114.98, 0.75575], [0,587.01, 617.78], color="blue")
        # ax.plot([99.01,99.01, 4.34], [0,366.49, 399.82], color="red")

        ax[2].set_xlabel("Energy [keV]")
        ax[2].set_ylabel("Counts")
        ax[2].legend()

        leaving_neutron = np.histogram(self.leaving_coffin["PostKinetic/keV"], bins=50)
        print(leaving_neutron[0],"\n",leaving_neutron[1])
        # 200 keV
        plt.savefig(self.plot_path + f"Cf_1E7_coffin_phys_effect_{self.config_string}.pdf")
        print(self.plot_path+ f"Cf_1E7_coffin_phys_effect_{self.config_string}.pdf")

    def neutron_spectrum_enteringLAr(self):

        rate_factor = 1000 * self.rate * self.Activity / (self.original_Activity * self.G4_events)
        self.argon = self.df_geo[self.df_geo["Volume"]=="LAr_phys"]
        # find first step in Ar and its Kinetic Energy

        self.argon["R/mm"]=  np.sqrt(self.argon["X/mm"]**2+self.argon["Y/mm"]**2 )
        self.argon["PreKinetic/keV"] = self.argon["PreKinetic/MeV"] * 1000


        #get event number for scatter and capture
        # only plot energy spectrum from scatter and capture event
        self.scatter = self.df_energy[self.df_energy["Process"].isin(['hadElastic', 'neutronInelastic'])]
        self.capture = self.df_energy[self.df_energy["Process"].isin(['nCapture'])]
        event_list = self.scatter["Event"].tolist() + self.capture["Event"].tolist()

        self.argon_intact = self.argon[self.argon['Event'].isin(event_list)]

        # slice argon
        self.argon_zslice = self.argon[(self.argon["Z/mm"]<520) & (self.argon["Z/mm"]>480)]
        #
        ar_counts,ar_edges = np.histogram(self.argon_intact["PreKinetic/keV"], bins=250,range= (0,2000))
        ar_rate  = ar_counts*rate_factor
        width = ar_edges[1]-ar_edges[0]

        ffig, ax = plt.subplots(1,3,figsize=(20,4))
        # X and Y
        sc0=ax[0].hist2d(self.argon["R/mm"],self.argon["Z/mm"],bins=50,
        cmap="plasma",norm="log",alpha=0.7)

        # ax.plot([189.95,189.95, 0.8485], [0,663.22, 714.03], color="red")
        # ax.plot([114.98,114.98, 0.75575], [0,587.01, 617.78], color="blue")
        # ax.plot([99.01,99.01, 4.34], [0,366.49, 399.82], color="red")

        ax[0].set_xlabel("R [mm]")
        ax[0].set_ylabel("Z [mm]")
        # ax[0].set_xlim(0,400)
        # ax[0].set_ylim(-100,800)
        cbar0 = plt.colorbar(sc0[3], ax=ax[0])
        cbar0.set_label("Counts(log)")

        sc1 = ax[1].hist2d(self.argon_zslice["X/mm"], self.argon_zslice["Y/mm"], bins=50,
                           cmap="plasma", norm="log", alpha=0.7)

        # ax.plot([189.95,189.95, 0.8485], [0,663.22, 714.03], color="red")
        # ax.plot([114.98,114.98, 0.75575], [0,587.01, 617.78], color="blue")
        # ax.plot([99.01,99.01, 4.34], [0,366.49, 399.82], color="red")

        ax[1].set_xlabel("X [mm]")
        ax[1].set_ylabel("Y [mm]")
        # ax[0].set_xlim(0,400)
        # ax[0].set_ylim(-100,800)
        cbar1 = plt.colorbar(sc1[3], ax=ax[1])
        cbar1.set_label("Counts(log)")



        ax[2].bar(ar_edges[:-1],ar_rate,width=width, align="edge")

        # ax.plot([189.95,189.95, 0.8485], [0,663.22, 714.03], color="red")
        # ax.plot([114.98,114.98, 0.75575], [0,587.01, 617.78], color="blue")
        # ax.plot([99.01,99.01, 4.34], [0,366.49, 399.82], color="red")

        ax[2].set_xlabel("Energy [keV]")
        ax[2].set_ylabel("Rate[mHz]")
        ax[2].set_xlim(0,2000)
        ax[2].set_yscale("log")

        plt.savefig(self.plot_path+f"Cf_1E7_position_density_first_LAr_{self.config_string}.pdf")


    def neutron_source_geometry(self):

        rate_factor = 1000 * self.rate * self.Activity / (self.original_Activity * self.G4_events)

        self.coffin = self.df_geo[self.df_geo["Volume"] == "cf_active_phys"]
        self.coffin["PreKinetic/keV"] = self.coffin["PreKinetic/MeV"] * 1000

        print(self.coffin)



        ffig, ax = plt.subplots(1,2,figsize=(16,4))
        # X and Y
        sc0=ax[0].hist2d(self.coffin["X/mm"],self.coffin["Z/mm"],bins=50,
        cmap="plasma",norm="log",alpha=0.7)

        # ax.plot([189.95,189.95, 0.8485], [0,663.22, 714.03], color="red")
        # ax.plot([114.98,114.98, 0.75575], [0,587.01, 617.78], color="blue")
        # ax.plot([99.01,99.01, 4.34], [0,366.49, 399.82], color="red")

        ax[0].set_xlabel("X [mm]")
        ax[0].set_ylabel("Z [mm]")
        # ax[0].set_xlim(0,400)
        # ax[0].set_ylim(-100,800)
        cbar0 = plt.colorbar(sc0[3], ax=ax[0])
        cbar0.set_label("Counts(log)")

        sc1 = ax[1].hist2d(self.coffin["Y/mm"], self.coffin["Z/mm"], bins=50,
                           cmap="plasma", norm="log", alpha=0.7)

        # ax.plot([189.95,189.95, 0.8485], [0,663.22, 714.03], color="red")
        # ax.plot([114.98,114.98, 0.75575], [0,587.01, 617.78], color="blue")
        # ax.plot([99.01,99.01, 4.34], [0,366.49, 399.82], color="red")

        ax[1].set_xlabel("Y [mm]")
        ax[1].set_ylabel("Z [mm]")
        # ax[0].set_xlim(0,400)
        # ax[0].set_ylim(-100,800)
        cbar1 = plt.colorbar(sc1[3], ax=ax[1])
        cbar1.set_label("Counts(log)")

        plt.savefig(self.plot_path+f"Cf_1E7_position_density_source_{self.config_string}.pdf")


    def NR_spectrum(self):
        rate_factor = 1000*self.rate*self.Activity/(self.original_Activity*self.G4_events) # /ms

        # # filter NR = 0 events
        print("initial df ", self.df_energy.head(10))
        # self.df_energy = self.df_energy[self.df_energy["Recoiled/MeV"]>0]
        self.scatter = self.df_energy[self.df_energy["Process"].isin(['hadElastic', 'neutronInelastic'])]
        self.capture = self.df_energy[self.df_energy["Process"].isin(['nCapture'])]
        max_NR_limit = max(self.scatter["Recoiled/MeV"]*1e6)
        print("maximum scatter recoil energy",max_NR_limit)
        # bin info and maybe same for both category
        # 100 ev per bin, 100 eV/ 0.100 keV 0.0001 MeV
        energy_scale = 100
        bin_num= int(max_NR_limit/100)+1
        max_bin_range= bin_num*100

        bin_range= (0,max_bin_range)
        (scatter_counts, scatter_edge) = np.histogram(self.scatter["Recoiled/MeV"]*1e6, bins=bin_num, range=bin_range)
        capture_counts = len(self.capture["Recoiled/MeV"])
        # read thermal neutron recoiled spectrum by MCMC
        self.TN_recoil_list = self.read_TN_R_spectrum()# in eV
        self.TN_recoil_keV_list = [i for i in self.TN_recoil_list]  # kev
        (capture_counts, capture_edge) = np.histogram(self.TN_recoil_keV_list, density=True,bins=bin_num, range=bin_range)
        width = capture_edge[1]-capture_edge[0]
        #0th order just a threshold

        scatter_rate_list = []
        capture_rate_list = []
        total_rate_list = []
        for threshold in scatter_edge:
            # normal
            # Efficiency_array = np.array([self.NucleationEfficiencyTrue(edge, threshold,threshold/8,threshold/8) for edge in scatter_edge] )
            # step function
            Efficiency_array = np.array(
                [self.NucleationEfficiencyTrue_Step(edge, threshold) for edge in scatter_edge])
            scatter_rate = sum(rate_factor*(Efficiency_array[1:]+Efficiency_array[:-1])*scatter_counts/2)
            # capture is different becasue density is true means is normalized also by bin width
            capture_rate = sum(rate_factor*capture_counts*(Efficiency_array[1:]+Efficiency_array[:-1])*width*capture_counts/2)
            scatter_rate_list.append(scatter_rate)
            capture_rate_list.append(capture_rate)
            total_rate_list.append(scatter_rate+capture_rate)

        print("NR edge, rate", capture_edge[:10], total_rate_list[:10])




        fig, ax = plt.subplots(1,4,figsize=(20, 4))

        ax[0].bar(scatter_edge, scatter_rate_list, width=width, align="edge")
        ax[0].set_xlabel("Energy threshold [eV]")
        ax[0].set_ylabel("Rate [mHz]")
        ax[0].set_yscale("log")
        # ax[0].set_xlim(0,2000)


        ax[1].bar(capture_edge, capture_rate_list, width=width, align="edge")
        ax[1].set_xlabel("Energy threshold [eV]")
        ax[1].set_ylabel("Rate [mHz]")
        ax[1].set_yscale("log")
        # ax[1].set_xlim(0, 2000)

        ax[2].bar(capture_edge, total_rate_list, width=width, align="edge")
        ax[2].set_xlabel("Energy threshold [eV]")
        ax[2].set_ylabel("Rate [mHz]")
        # ax[2].set_yscale("log")
        ax[2].set_xlim(0, 3500)

        # sc = ax[3].hist2d(self.scatter["PreKinetic/MeV"]*1e6, self.scatter["Recoiled/MeV"]*1e6, bins=50,
        #                   cmap="plasma", norm="log", alpha=0.7)
        #
        #
        #
        # ax[3].set_xlabel("PreKinetic [eV]")
        # ax[3].set_ylabel("Recoil [eV]")
        # # ax[0].set_xlim(0,400)
        # # ax[0].set_ylim(-100,800)
        # cbar = plt.colorbar(sc[3], ax=ax[3])
        # cbar.set_label("Counts(log)")


        plt.savefig(self.plot_path+f"Cf_1E7_energy_density_{self.config_string}.pdf")
        print(self.plot_path)


        self.NR_rate_zoomed(capture_edge, total_rate_list)


    def NR_spectrum_moderated_by_sstl(self):
        rate_factor = 1000*self.rate*self.Activity/(self.original_Activity*self.G4_events) # /ms

        # # filter NR = 0 events
        print("initial df ", self.df_energy.head(10))
        # self.df_energy = self.df_energy[self.df_energy["Recoiled/MeV"]>0]
        self.df_sstl_energy = self.df_energy[self.df_energy['Event'].isin(self.low_e_n_list)]
        self.scatter = self.df_sstl_energy[self.df_sstl_energy["Process"].isin(['hadElastic', 'neutronInelastic'])]
        self.capture = self.df_sstl_energy[self.df_sstl_energy["Process"].isin(['nCapture'])]

        max_NR_limit = max(self.scatter["Recoiled/MeV"]*1e6)
        print("maximum scatter recoil energy",max_NR_limit)
        # bin info and maybe same for both category
        # 100 ev per bin
        bin_num= int(max_NR_limit/100)+1
        max_bin_range= bin_num*100
        print("bin_num",bin_num,"max_bin_range",max_bin_range)

        bin_range= (0,max_bin_range)
        (scatter_counts, scatter_edge) = np.histogram(self.scatter["Recoiled/MeV"]*1e6, bins=bin_num, range=bin_range)
        capture_counts = len(self.capture["Recoiled/MeV"])
        # read thermal neutron recoiled spectrum by MCMC
        self.TN_recoil_list = self.read_TN_R_spectrum()# in eV
        self.TN_recoil_keV_list = [i for i in self.TN_recoil_list] #eV
        (capture_counts, capture_edge) = np.histogram(self.TN_recoil_keV_list, density=True,bins=bin_num, range=bin_range)
        width = capture_edge[1]-capture_edge[0]
        #0th order just a threshold

        scatter_rate_list = []
        capture_rate_list = []
        total_rate_list = []
        for threshold in scatter_edge:
            # normal
            # Efficiency_array = np.array([self.NucleationEfficiencyTrue(edge, threshold,threshold/8,threshold/8) for edge in scatter_edge] )
            # step function
            Efficiency_array = np.array(
                [self.NucleationEfficiencyTrue_Step(edge, threshold) for edge in scatter_edge])

            scatter_rate = sum(rate_factor*(Efficiency_array[1:]+Efficiency_array[:-1])*scatter_counts/2)
            # capture is different becasue density is true means is normalized also by bin width
            capture_rate = sum(rate_factor*capture_counts*(Efficiency_array[1:]+Efficiency_array[:-1])*width*capture_counts/2)
            scatter_rate_list.append(scatter_rate)
            capture_rate_list.append(capture_rate)
            total_rate_list.append(scatter_rate+capture_rate)

        print("NR edge, rate", capture_edge[:10], total_rate_list[:10])




        fig, ax = plt.subplots(1,4,figsize=(20, 4))

        ax[0].bar(scatter_edge, scatter_rate_list, width=width, align="edge")
        ax[0].set_xlabel("Recoiled Energy threshold [eV]")
        ax[0].set_ylabel("Rate [mHz]")
        ax[0].set_yscale("log")
        # ax[0].set_xlim(0,2000)


        ax[1].bar(capture_edge, capture_rate_list, width=width, align="edge")
        ax[1].set_xlabel("Recoiled Energy threshold [eV]")
        ax[1].set_ylabel("Rate [mHz]")
        ax[1].set_yscale("log")
        # ax[1].set_xlim(0, 2000)

        ax[2].bar(capture_edge, total_rate_list, width=width, align="edge")
        ax[2].set_xlabel("Recoiled Energy threshold [eV]")
        ax[2].set_ylabel("Rate [mHz]")
        # ax[2].set_yscale("log")
        ax[2].set_xlim(0, 3500)

        # sc = ax[3].hist2d(self.scatter["PreKinetic/MeV"]*1e6, self.scatter["Recoiled/MeV"]*1e6, bins=50,
        #                   cmap="plasma", norm="log", alpha=0.7)
        #
        #
        #
        # ax[3].set_xlabel("PreKinetic [eV]")
        # ax[3].set_ylabel("Recoil [eV]")
        # # ax[0].set_xlim(0,400)
        # # ax[0].set_ylim(-100,800)
        # cbar = plt.colorbar(sc[3], ax=ax[3])
        # cbar.set_label("Counts(log)")


        plt.savefig(self.plot_path+f"Cf_1E7_sstl_moderated_n_{self.config_string}.pdf")
        print(self.plot_path)


        self.NR_rate_zoomed(capture_edge, total_rate_list)
    def NR_multiplicity(self):
        rate_factor = 1000*self.rate*self.Activity/(self.original_Activity*self.G4_events) # /ms
        print("rate_factor",rate_factor)
        print("initial df ", self.df_energy.head(10))
        # fill capture event energy

        self.TN_recoil_list = self.read_TN_R_spectrum()
        TN_array = np.array(self.TN_recoil_list)
        capture_mask = self.df_energy["Process"].isin(['nCapture'])
        n = capture_mask.sum()

        self.df_energy.loc[capture_mask, 'Recoiled/MeV'] = TN_array[
            np.random.randint(0, len(TN_array), size=n)
        ]/1e6

        self.df_energy.loc[:,"Recoiled/eV"] = self.df_energy.loc[:,"Recoiled/MeV"]*1e6
        print('len(self.df_energy)',len(self.df_energy))
        low_NR = self.df_energy[self.df_energy["Recoiled/eV"] <0.1]
        print('low_NR.head(10)',low_NR.head(10))

        # generate filter for different threshold every 100 eV
        max_energy = max(self.df_energy["Recoiled/eV"])
        min_energy = min(self.df_energy["Recoiled/eV"])
        min_edge = round(min_energy/100,0)*100
        max_edge = (round(max_energy / 100, 0)+1) * 100
        bin_n = int((max_edge-min_edge)/100)
        print('min_edge, max_edge, bin_n',min_edge, max_edge, bin_n)

        multiplicity_list = []
        # for i in range(0,bin_n):
        for i in range(0, 30, 5):
            energy_threshold = min_edge+ i*1000 # every 1000 eV

            energy_mask = (self.df_energy["Recoiled/eV"]>= energy_threshold)
            df = self.df_energy.loc[energy_mask,:]
            # calculate the multiplicity and rate
            multiplicity = df.groupby("Event").size().tolist()
            # make histogram
            # all bins are int
            multiplicity_min = min(multiplicity)
            multiplicity_max = max(multiplicity)
            bin_num = multiplicity_max-multiplicity_min+1
            bin_range = (multiplicity_min, multiplicity_max+1)
            # print("energy threshold", energy_threshold, "min", multiplicity_min, "max", multiplicity_max, "bin_num", bin_num, "bin_range", bin_range)
            (multiplicity_counts, multiplicity_edges)= np.histogram(multiplicity, bins=bin_num, range=bin_range)
            # print("edges", multiplicity_edges)
            multiplicity_rates = multiplicity_counts*rate_factor
            multiplicity_ratios = multiplicity_counts /multiplicity_counts[0]
            multiplicity_sigma = multiplicity_counts*rate_factor/np.sqrt(multiplicity_counts)
            ratio_sigma  = (np.sqrt(multiplicity_counts)*multiplicity_counts[0]+multiplicity_counts*np.sqrt(multiplicity_counts[0]))/(multiplicity_counts[0])**2
            print('len(multiplicity_ratios),len(multiplicity_sigma),len(ratio_sigma)',len(multiplicity_ratios),len(multiplicity_sigma),len(ratio_sigma))
            multiplicity_width = multiplicity_edges[1]-multiplicity_edges[0]
            multiplicity_list.append((multiplicity_rates,multiplicity_edges,multiplicity_width,energy_threshold,multiplicity_ratios,multiplicity_sigma,ratio_sigma))

        fig, ax = plt.subplots()
        for i in range(len(multiplicity_list)):
            # ax.plot(multiplicity_list[i][1][:-1], multiplicity_list[i][0],
            #        label="threshold " + str(multiplicity_list[i][3]) + " eV")
            ax.errorbar(multiplicity_list[i][1][:-1], multiplicity_list[i][0],
                           yerr=multiplicity_list[i][5], label="threshold " + str(multiplicity_list[i][3]) + " eV")

            # ax.bar(multiplicity_list[i][1][:-1], multiplicity_list[i][0], width=multiplicity_list[i][2], align="edge", label="threshold "+str(multiplicity_list[i][3])+" eV" )
        ax.set_xlabel("Multiplicity")
        ax.set_ylabel("Rate [mHz]")
        # ax.set_yscale("log")
        ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
        ax.legend(fontsize='small')
        print("threshold range", min_edge, max_edge)
        # ax[0].set_xlim(0,2000)
        plt.savefig(self.plot_path + f"Cf_1E7_multiplicity_{self.config_string}.pdf")
        plt.clf()

        fig, ax = plt.subplots()
        for i in range(len(multiplicity_list)):
            # ax.plot(multiplicity_list[i][1][:-1], multiplicity_list[i][4],
            #        label="threshold " + str(multiplicity_list[i][3]) + " eV")
            ax.errorbar(multiplicity_list[i][1][:-1], multiplicity_list[i][4],
                        yerr=multiplicity_list[i][6], label="threshold " + str(multiplicity_list[i][3]) + " eV")

            # ax.bar(multiplicity_list[i][1][:-1], multiplicity_list[i][0], width=multiplicity_list[i][2], align="edge", label="threshold "+str(multiplicity_list[i][3])+" eV" )
        ax.set_xlabel("Multiplicity")
        ax.set_ylabel("Ratio []")

        # set_yscale("log")
        ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
        ax.legend(fontsize='small')
        plt.savefig(self.plot_path + f"Cf_1E7_ratio_{self.config_string}.pdf")

        fig, ax = plt.subplots()
        energy_th = []
        multiplicity_cut = 5
        ratio_matrix = np.empty((len(multiplicity_list), multiplicity_cut))
        ratio_err_matrix = np.empty((len(multiplicity_list), multiplicity_cut))
        rate_matrix = np.empty((len(multiplicity_list), multiplicity_cut))
        rate_err_matrix = np.empty((len(multiplicity_list), multiplicity_cut))
        # tranverse the matrix but only keep multiplicity <=4

        for i in range(len(multiplicity_list)):
            energy_th.append(multiplicity_list[i][3])
            ratio_matrix[i]= multiplicity_list[i][4][:multiplicity_cut]
            ratio_err_matrix[i] = multiplicity_list[i][6][:multiplicity_cut]
            rate_matrix[i] = multiplicity_list[i][0][:multiplicity_cut]
            rate_err_matrix[i] = multiplicity_list[i][5][:multiplicity_cut]
        ratio_t = ratio_matrix.T
        ratio_error_t = ratio_err_matrix.T
        rate_t = rate_matrix.T
        rate_error_t = rate_err_matrix.T

        for i in range(multiplicity_cut):
            # ax.plot(multiplicity_list[i][1][:-1], multiplicity_list[i][4],
            #        label="threshold " + str(multiplicity_list[i][3]) + " eV")
            ax.errorbar(energy_th, ratio_t[i],
                        yerr=ratio_error_t[i], label=f"multiplicity = {i+1}")
            print(f"multiplicity = {i+1}", energy_th,
                  ratio_t[i],ratio_error_t[i])

        ax.set_xlabel("Energy Threshold [eV]")
        ax.set_ylabel("Ratio []")

        # set_yscale("log")
        # ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
        ax.legend(fontsize='small')
        plt.savefig(self.plot_path + f"Cf_1E7_ratio_{self.config_string}_threshold.pdf")

        plt.clf()
        fig, ax = plt.subplots()
        for i in range(multiplicity_cut):
            # ax.plot(multiplicity_list[i][1][:-1], multiplicity_list[i][4],
            #        label="threshold " + str(multiplicity_list[i][3]) + " eV")
            ax.errorbar(energy_th, rate_t[i],
                        yerr=rate_error_t[i], label=f"multiplicity = {i+1}")

        ax.set_xlabel("Energy Threshold [eV]")
        ax.set_ylabel("Rate [mHz]")

        # set_yscale("log")
        # ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
        ax.legend(fontsize='small')
        plt.savefig(self.plot_path + f"Cf_1E7_rate_{self.config_string}_threshold.pdf")



        print(self.plot_path)
    def NR_rate_zoomed(self, edges, rates):
        # edges are in eV and rates in mHz
        # edges is always 100eV per bin
        fig, ax = plt.subplots()
        ax.plot(edges[:40], rates[:40])
        ax.set_xlabel("Energy threshold [eV]")
        ax.set_ylabel("Rate [mHz]")
        ax.set_ylim(0,140)
        plt.savefig(self.plot_path + "Cf_1E7_NR_zoomed.pdf")

    def write_sims_results(self):
        # rate factor in mHz

        rate_factor = 1000 * self.rate * self.Activity / (self.original_Activity * self.G4_events)  # mHz
        print("initial df ", self.df_energy.head(10))
        print("max(ER_Ar)", max(self.df_energy["Recoiled/MeV"]))
        # fill capture event energy

        self.TN_recoil_list = self.read_TN_R_spectrum()
        TN_array = np.array(self.TN_recoil_list)
        capture_mask = self.df_energy["Process"].isin(['nCapture'])
        n = capture_mask.sum()

        self.df_energy.loc[capture_mask, 'Recoiled/MeV'] = TN_array[
                                                               np.random.randint(0, len(TN_array), size=n)
                                                           ] / 1e6

        self.df_energy.loc[:, "Recoiled/eV"] = self.df_energy.loc[:, "Recoiled/MeV"] * 1e6

        NR_Ar = self.df_energy["Recoiled/eV"]   # in eV


        print("max(ER_Ar)", max(NR_Ar))
        print('self.df_energy[self.df_energy["Recoiled/eV"]>1e5]',self.df_energy[self.df_energy["Recoiled/eV"]>1e5])
        max_Ar = max(NR_Ar)
        hist_array = [None]
        hist_array[0] = np.histogram(NR_Ar, bins=int(max_Ar/100), range=(0, max_Ar))
        # hist_array[0] = np.histogram(NR_Ar, bins=350, range=(0, 3500))
        # every 100 eV per bin

        # transfer edge to mid point per bin

        # get probablity per scattering and the statistics
        cumulative_threshold_per_scatter_array = [None]

        cumulative_threshold_per_scatter_array[0] = np.array(
            [sum(hist_array[0][0][i:]) for i in range(len(hist_array[0][0]))])
        print("cumulative spectrum  in counts", cumulative_threshold_per_scatter_array[0])

        # histogram per scattering per keV
        cumulative_threshold_array = [None]
        energy_deposit_list = [hist_array[0][0][i] * hist_array[0][1][i] for i in range(len(hist_array[0][0]))]

        cumulative_threshold_array[0] = np.array(
            [sum(energy_deposit_list[i:]) for i in range(len(energy_deposit_list))])
        print("total count* energy Cf", cumulative_threshold_per_scatter_array[0][0],cumulative_threshold_per_scatter_array[0][0]/cumulative_threshold_array[0][0])

        #get first argon energy diff spectrum
        first_argon_ene_array = [None]
        #try to read
        try:#
            # also for thoes neutrons cause Ar recoil
            lar_event_list = self.df_energy["Event"].tolist()
            self.argon = self.df_geo[(self.df_geo["Volume"] == "LAr_phys")&(self.df_geo["Event"].isin(lar_event_list))]
            argon_energy = self.argon["PreKinetic/MeV"]
            argon_energy["PreKinetic/keV"] = self.argon["PreKinetic/MeV"]*1000
            # 100 per bin
            max_first_argon = max(self.argon["PreKinetic/MeV"])*1e6
            first_argon_ene_array[0] = np.histogram(argon_energy["PreKinetic/keV"],bins=int(max_first_argon/100), range=(0, max_first_argon/1000))
        except:
            print("error in neutron entering argon")
        finally:
            print("different neutron energy firstly entering argon volume")




        output_list = [rate_factor ,hist_array, cumulative_threshold_per_scatter_array[0], cumulative_threshold_array[0],first_argon_ene_array[0]]
        # output form, rate facotr to mHz, enenrgy edges, counts above the bin edge, counts* energy above the bin edge
        # with open(f"/data/runzezhang/result/TN_sims_D/Cf_output_1E7_{self.config_string}.pkl", "wb") as f:
        with open(f"/data/runzezhang/result/TN_sims_D/Cf_output_1E7_{self.config_string}.pkl", "wb") as f:
            pickle.dump(output_list, f)

    def check_neutron_spectrum_Ar(self):
        # get first argon energy diff spectrum
        first_argon_ene_array = [None]
        # try to read

        # also for those neutrons cause Ar recoil
        lar_event_list = self.df_energy["Event"].tolist()
        lar_event_list = list(dict.fromkeys(lar_event_list))

        self.argon = self.df_geo[(self.df_geo["Volume"] == "LAr_phys")&(self.df_geo["Event"].isin(lar_event_list))]
        argon_energy = self.argon
        argon_energy["PreKinetic/keV"] = argon_energy["PreKinetic/MeV"]*1000
        tagged1_event_list =argon_energy[(argon_energy["PreKinetic/keV"].values>486)&(argon_energy["PreKinetic/keV"].values<487.5)]["Event"].tolist()
        #
        # Explicitly use .loc for label/boolean alignment
        # argon_energy = argon_energy.reset_index(drop=True)
        # mask = (argon_energy["PreKinetic/keV"] > 486) & (argon_energy["PreKinetic/keV"] < 487.5)
        # tagged1_event_list = argon_energy.loc[mask, "Event"].tolist()
        print(self.df_energy[self.df_energy["Event"].isin(tagged1_event_list)])



    def NucleationEfficiencyTrue(self, r, T, sigLow, sigUp):
        if r < T:
            R = 1 / 2 * (1 + math.erf((r - T) / (sigLow * 2 ** (1 / 2))))
        else:
            R = 1 / 2 * (1 + math.erf((r - T) / (sigUp * 2 ** (1 / 2))))
        return R

    def NucleationEfficiencyTrue_Step(self, r, T):
        if r < T:
            R = 0
        else:
            R = 1
        return R
    def read_TN_R_spectrum(self):
        with open(self.TN_spectrum_path, "rb") as fp:  # Unpickling
            TN_recoil_list = pickle.load(fp)
            # print("read", self.TN_recoil_list)
        return TN_recoil_list
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
    sn = SN(gamma=True,full_gamma=False)
    # test = test_csv()