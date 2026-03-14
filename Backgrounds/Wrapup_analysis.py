import pandas as pd
import matplotlib.pyplot as plt
import csv
import numpy as np
import os
import pickle
class integrated_analysis():
    def __init__(self):


        self.Co_sim_path  ='/data/runzezhang/result/TN_sims_D/Co_output.pkl'
        self.Cs_sim_path = '/data/runzezhang/result/TN_sims_D/Cs_output.pkl'

        self.main()
    def main(self):
        self.read_sims()
    def read_sims(self):
        with open(self.Co_sim_path, "rb") as f:
            self.Co_sims = pickle.load(f)
        print("self.Co_sims",self.Co_sims)
        with open(self.Cs_sim_path, "rb") as f:
            self.Cs_sims = pickle.load(f)
        print("self.Cs_sims",self.Cs_sims)



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
    IA =  integrated_analysis()
    # test = test_csv()