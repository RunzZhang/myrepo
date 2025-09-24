import pandas as pd
import matplotlib.pyplot as plt
import csv
import numpy as np
class SN():
    def __init__(self):
        # after generate new files, you need to select the capture ratio/source for different configs in read_files function.
        # then choose the correct signal/noise of with clause in read files.
        # at last change the self.name and plot_name in plot function
        self.base_path = "/data/runzezhang/result/TN_sims_D/chunked_root_files_AmLi_LZ_cross_check/"
        self.base_path2 = "/data/runzezhang/result/TN_sims_D/chunked_root_files_AmLi_LZ_cross_check/"
        self.plot_path = '/data/runzezhang/result/TN_sims_D/plot/'
        self.false_1 = "Cf_false1.csv"
        self.false_2 = "Cf_false2.csv"
        self.signal = "Cf_sig.csv"
        self.name = "Background 2"
        self.plot_name = self.name+"1E7.png"
        self.signal_final_list = []
        self.noise1_final_list =[]
        self.noise2_final_list = []
        self.noise3_final_list = []
        self.signal_new_final_list = []
        self.noise1_new_final_list = []
        self.noise2_new_final_list = []

        self.main_body(1)
        #982 statics false 1
        # for i in range(1,101):
        #     self.main_body(i)
        self.print_event()

    def main_body(self,i):
        print(i)
        self.false_1 = f"AmLi_1E7_false1_part{i}.csv"
        self.false_2 = f"AmLi_1E7_false2_part{i}.csv"
        self.false_3 = f"AmLi_1E7_false3_part{i}.csv"
        self.false_1_new = f"AmLi_1E7_false1_new_part{i}.csv"
        self.false_2_new = f"AmLi_1E7_false2_new_part{i}.csv"
        self.signal = f"AmLi_1E7_sig_part{i}.csv"
        self.signal_new = f"AmLi_1E7_sig_new_part{i}.csv"
        self.false_1_mid = f"AmLi_1E7_false1_mid_part{i}.csv"
        self.false_2_mid = f"AmLi_1E7_false2_mid_part{i}.csv"
        self.false_3_mid = f"AmLi_1E7_false3_mid_part{i}.csv"
        self.signal_mid = f"AmLi_1E7_sig_mid_part{i}.csv"
        self.false_1_new_mid = f"AmLi_1E7_false1_new_mid_part{i}.csv"
        self.false_2_new_mid = f"AmLi_1E7_false2_new_mid_part{i}.csv"
        self.signal_new_mid = f"AmLi_1E7_sig_new_mid_part{i}.csv"
        self.false_1_path = self.base_path + self.false_1
        self.false_2_path = self.base_path + self.false_2
        self.false_3_path = self.base_path + self.false_3
        self.false_1_new_path = self.base_path + self.false_1_new
        self.false_2_new_path = self.base_path + self.false_2_new
        self.false_1_path_mid = self.base_path + self.false_1_mid
        self.false_2_path_mid = self.base_path + self.false_2_mid
        self.false_3_path_mid = self.base_path + self.false_3_mid
        self.false_1_new_path_mid = self.base_path + self.false_1_new_mid
        self.false_2_new_path_mid = self.base_path + self.false_2_new_mid
        self.signal_path_mid = self.base_path + self.signal_mid
        self.signal_path = self.base_path + self.signal
        self.signal_new_path_mid = self.base_path + self.signal_new_mid
        self.signal_new_path = self.base_path + self.signal_new



        self.read_files()

        # self.read_files_s_to_N1()

# main funtion we use
    def read_files(self):

        with open(self.signal_path, 'r') as file:
            reader = csv.reader(file)
            # Read the first row (assuming single row for simplicity)
            number_list = next(reader)
            # Convert the strings to floats
            self.sig_raw_list = [float(value) for value in number_list]
        self.signal_final_list = self.signal_final_list + self.sig_raw_list
        self.signal_final_list = list(dict.fromkeys(self.signal_final_list))

        with open(self.false_1_path, 'r') as file:
            reader = csv.reader(file)
            # Read the first row (assuming single row for simplicity)
            number_list = next(reader)
            # Convert the strings to floats
            self.noise1_raw_list = [float(value) for value in number_list]
        self.noise1_final_list +=  self.noise1_raw_list
        self.noise1_final_list = list(dict.fromkeys(self.noise1_final_list))

        with open(self.false_2_path, 'r') as file:
            reader = csv.reader(file)
            # Read the first row (assuming single row for simplicity)
            number_list = next(reader)
            # Convert the strings to floats
            self.noise2_raw_list = [float(value) for value in number_list]
        self.noise2_final_list +=  self.noise2_raw_list
        self.noise2_final_list = list(dict.fromkeys(self.noise2_final_list))

        with open(self.false_3_path, 'r') as file:
            reader = csv.reader(file)
            # Read the first row (assuming single row for simplicity)
            number_list = next(reader)
            # Convert the strings to floats
            self.noise3_raw_list = [float(value) for value in number_list]
        self.noise3_final_list +=  self.noise3_raw_list
        self.noise3_final_list = list(dict.fromkeys(self.noise3_final_list))

        with open(self.signal_new_path, 'r') as file:
            reader = csv.reader(file)
            # Read the first row (assuming single row for simplicity)
            number_list = next(reader)
            # Convert the strings to floats
            self.sig_new_raw_list = [float(value) for value in number_list]
        self.signal_new_final_list += self.sig_new_raw_list
        self.signal_new_final_list = list(dict.fromkeys(self.signal_new_final_list))

        with open(self.false_1_new_path, 'r') as file:
            reader = csv.reader(file)
            # Read the first row (assuming single row for simplicity)
            number_list = next(reader)
            # Convert the strings to floats
            self.noise1_new_raw_list = [float(value) for value in number_list]
        self.noise1_new_final_list += self.noise1_new_raw_list
        self.noise1_new_final_list = list(dict.fromkeys(self.noise1_new_final_list))

        with open(self.false_2_new_path, 'r') as file:
            reader = csv.reader(file)
            # Read the first row (assuming single row for simplicity)
            number_list = next(reader)
            # Convert the strings to floats
            self.noise2_new_raw_list = [float(value) for value in number_list]
        self.noise2_new_final_list += self.noise2_new_raw_list
        self.noise2_new_final_list = list(dict.fromkeys(self.noise2_new_final_list))

    def print_event(self):
        print("signal list diff")
        diff_sig = set(self.signal_new_final_list) ^ set(self.signal_final_list)
        print(diff_sig)
        print("Huge scattering diff")
        diff_huge = set(self.noise2_new_raw_list) ^ set(self.noise2_final_list)
        print(diff_huge)
        print("ER scattering diff")
        self.old_scattering = self.noise1_final_list+self.noise3_final_list
        diff_scattering = set(self.old_scattering) ^ set(self.noise2_new_final_list)
        print(diff_scattering)

