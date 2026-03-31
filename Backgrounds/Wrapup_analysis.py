import pandas as pd
import matplotlib.pyplot as plt
import csv
import numpy as np
import os
import pickle
from scipy.optimize import curve_fit
class integrated_analysis():
    def __init__(self):

        self.output_path = '/data/runzezhang/result/gamma_rejection/'
        self.plot_path = '/data/runzezhang/result/gamma_rejection/plot/'
        self.Co_sim_path  ='/data/runzezhang/result/TN_sims_D/Co_output_5E7.pkl'
        self.Cs_sim_path = '/data/runzezhang/result/TN_sims_D/Cs_output.pkl'
        self.Cf_simA_path = '/data/runzezhang/result/TN_sims_D/Cf_output_1E7_config_A.pkl'
        self.Cf_simB_path = '/data/runzezhang/result/TN_sims_D/Cf_output_1E7_config_B.pkl'

        self.Cs_exp_116_raw_path = ["Cold-Cs-11_17-18_exposures_mix","Cold-Cs-12_01_exposures_mix",
                                    "Cold-Cs-12_10-11_exposures_mix","Cold-Cs-1_20-21_exposures_mix"]

        self.Cs_exp_116_raw_len = len(self.Cs_exp_116_raw_path)
        self.Cs_exp_119_raw_path = ["Cold-Cs-2_2-3_exposures_zoom"]
        self.Co_exp_116_raw_path = ["60Co-12_15-16_exposures"]
        self.Co_exp_116_raw_len = len(self.Co_exp_116_raw_path)
        self.Co_exp_119_raw_path = []
        self.Co_exp_119_raw_len = len(self.Co_exp_119_raw_path)
        self.Cf_exp_116A_raw_path = ['252Cf-Coffin-A-1_7-8_exposures']
        self.Cf_exp_116A_raw_len = len(self.Cf_exp_116A_raw_path)
        self.Cf_exp_116B_raw_path = ['252Cf-Coffin-B-1_8-9_exposures']
        self.Cf_exp_116B_raw_len = len(self.Cf_exp_116B_raw_path)

        self.background_116_sorted_path = ['Background-11_7-8_exposures',"Background-11_15-17_exposures",
                                    "Background-11_26-30_exposures", "Background-12_5-8_exposures",
                                    'Background-1_12-13_exposures','Background-1_17-20_exposures']
        self.background_exp_116_raw_len = len(self.background_116_sorted_path)
        self.background_119_sorted_path = ["Background-1_30-2_2_exposures","Background-2_17-20_exposures",
                                           "Background-2_6-12_exposures","Background-2_28-3_2_exposures"]

        self.Co_exp_raw_path = self.Co_exp_116_raw_path+self.Co_exp_119_raw_path
        self.Cs_exp_raw_path = self.Cs_exp_116_raw_path+self.Cs_exp_119_raw_path
        self.backgrounds_exp_raw_path = self.background_116_sorted_path+self.background_119_sorted_path

        self.main()
    def main(self):
        # manage the path configurations
        # paths are devided by sources and temperature configurations
        # and each stage of the analysis

        self.generate_path()

        # read simulation from Geant4
        self.read_sims()

        # read the Seitz energy thresholds
        self.read_Seitz_info()

        # read experimental txt file, drop the non-sense values, and write to clean dataframe
        self.read_raw_Co_exp()
        self.read_raw_Cs_exp()
        self.read_raw_Cf_exp()

        self.read_raw_backgrounds_exp()

        # caculate the average background
        self.average_background_analysis()


        # add the Seitz energy to the exp txt files, and Seitz should have already included all temperature info, so in post-analysis
        # no demands to devide by temperature configurations
        # also caculate the clean signal and signal uncerntainty
        self.clean_signal_analysis()

        # Cf
        self.clean_NR_signal_analysis()

        # based on sims and clean signal rate, calculate gamma rejection
        # self.gamma_rejection_calculation()



        # plot
        # self.bkg_plot()
        # self.gamma_rejection_plot()
        # self.spectrums_plot()


    def generate_path(self):
        self.Co_exp_sorted_path = []
        self.Cs_exp_sorted_path = []
        self.Cf_expA_sorted_path = []
        self.Cf_expB_sorted_path = []
        self.Co_exp_rate_path = []
        self.Cs_exp_rate_path = []
        self.Cf_expA_rate_path = []
        self.Cf_expB_rate_path = []
        self.Co_exp_rejection_path = []
        self.Cs_exp_rejection_path = []
        self.Bkg_exp_sorted_path = []
        self.Bkg_average_116_path = self.output_path + "background_116_average" + ".csv"
        self.Bkg_average_119_path = self.output_path + "background_119_average" + ".csv"
        # bkg table containing seitz infos
        # separate this from above because of clean signal need to merge only on pressure column
        self.Bkg_average_116_full_info_path = self.output_path + "background_116_average_full_info" + ".csv"
        self.Bkg_average_119_full_info_path = self.output_path + "background_119_average_full_info" + ".csv"


        for exp_name in self.Co_exp_raw_path:
            self.Co_exp_sorted_path.append(self.output_path + exp_name+"_sorted.csv")
        for exp_name in self.Cs_exp_raw_path:
            self.Cs_exp_sorted_path.append(self.output_path + exp_name+"_sorted.csv")
        # Cf
        for exp_name in self.Cf_exp_116A_raw_path:
            self.Cf_expA_sorted_path.append(self.output_path + exp_name+"_sorted.csv")
        for exp_name in self.Cf_exp_116B_raw_path:
            self.Cf_expB_sorted_path.append(self.output_path + exp_name+"_sorted.csv")

        for exp_name in self.Co_exp_raw_path:
            self.Co_exp_rate_path.append(self.output_path + exp_name+"_rate.csv")
        for exp_name in self.Cs_exp_raw_path:
            self.Cs_exp_rate_path.append(self.output_path + exp_name+"_rate.csv")

        for exp_name in self.Cf_exp_116A_raw_path:
            self.Cf_expA_rate_path.append(self.output_path + exp_name+"_rate.csv")
        for exp_name in self.Cf_exp_116B_raw_path:
            self.Cf_expB_rate_path.append(self.output_path + exp_name+"_rate.csv")

        for exp_name in self.Co_exp_raw_path:
            self.Co_exp_rejection_path.append(self.output_path + exp_name+"_rejection.csv")
        for exp_name in self.Cs_exp_raw_path:
            self.Cs_exp_rejection_path.append(self.output_path + exp_name+"_rejection.csv")
        for exp_name in self.backgrounds_exp_raw_path:
            self.Bkg_exp_sorted_path.append(self.output_path + exp_name+"_sorted.csv")

    def read_sims(self):
        #output form, [rate facotr to mHz, (counts per bin,enenrgy edges), counts above the bin edge, counts* counts above the bin edge]
        with open(self.Co_sim_path, "rb") as f:
            self.Co_sims = pickle.load(f)
        print("self.Co_sims",self.Co_sims)
        with open(self.Cs_sim_path, "rb") as f:
            self.Cs_sims = pickle.load(f)
        print("self.Cs_sims",self.Cs_sims)

        with open(self.Cf_simA_path, "rb") as f:
            self.Cf_simsA = pickle.load(f)
        print("self.Cf_simsA",self.Cf_simsA)

        with open(self.Cf_simB_path, "rb") as f:
            self.Cf_simsB = pickle.load(f)
        print("self.Cf_simsB",self.Cf_simsB)

    def read_raw_Cs_exp(self):
        # read file, delete unreasonable rows and rewrite
        for i in range(len(self.Cs_exp_raw_path)):
            exposure_df = self.read_exposure(self.Cs_exp_raw_path[i] + ".txt")
            exposure_df.columns = ['Pressure [bara]', 'Lifetime [s]', 'Lifetime Error [s]','Exponential Fit 2xNLL','N.d.o.f.','Time Cut High [s]','Time Cut Low [s]']
            exposure_df = exposure_df[(exposure_df['Lifetime [s]'] <=6.92e-1) | (exposure_df['Lifetime [s]'] >= 6.94e-1)]
            exposure_df = exposure_df[
                (exposure_df['Lifetime Error [s]']/exposure_df['Lifetime [s]'] <= 0.3)]
            # add rate column
            exposure_df['Exp Rate [mHz]']= 1000/exposure_df['Lifetime [s]']
            exposure_df['Exp Rate Sigma [mHz]'] = exposure_df['Lifetime Error [s]'] * 1000 / (exposure_df['Lifetime [s]']) ** 2
            exposure_df = exposure_df.drop(columns=['Exponential Fit 2xNLL',
                                                    'N.d.o.f.', 'Time Cut High [s]', 'Time Cut Low [s]'])
            exposure_df.to_csv(self.Cs_exp_sorted_path[i],index=False)
    def read_raw_Co_exp(self):
        # read file, delete unreasonable rows and rewrite
        for i in range(len(self.Co_exp_raw_path)):
            exposure_df = self.read_exposure(self.Co_exp_raw_path[i] + ".txt")
            exposure_df.columns = ['Pressure [bara]', 'Lifetime [s]', 'Lifetime Error [s]','Exponential Fit 2xNLL','N.d.o.f.','Time Cut High [s]','Time Cut Low [s]']
            exposure_df = exposure_df[(exposure_df['Lifetime [s]'] <=6.92e-1) | (exposure_df['Lifetime [s]'] >= 6.94e-1)]
            exposure_df = exposure_df[
                (exposure_df['Lifetime Error [s]']/exposure_df['Lifetime [s]'] <= 0.3)]
            exposure_df['Exp Rate [mHz]'] = 1000 / exposure_df['Lifetime [s]']
            exposure_df['Exp Rate Sigma [mHz]'] = exposure_df['Lifetime Error [s]'] * 1000 / (
            exposure_df['Lifetime [s]']) ** 2

            exposure_df = exposure_df.drop(columns=['Exponential Fit 2xNLL',
                                                    'N.d.o.f.', 'Time Cut High [s]', 'Time Cut Low [s]'])
            exposure_df.to_csv(self.Co_exp_sorted_path[i],index=False)

    def read_raw_Cf_exp(self):
        # read file, delete unreasonable rows and rewrite
        for i in range(len(self.Cf_exp_116A_raw_path)):
            exposure_df = self.read_exposure(self.Cf_exp_116A_raw_path[i] + ".txt")
            exposure_df.columns = ['Pressure [bara]', 'Lifetime [s]', 'Lifetime Error [s]','Exponential Fit 2xNLL','N.d.o.f.','Time Cut High [s]','Time Cut Low [s]']
            exposure_df = exposure_df[(exposure_df['Lifetime [s]'] <=6.92e-1) | (exposure_df['Lifetime [s]'] >= 6.94e-1)]
            exposure_df = exposure_df[
                (exposure_df['Lifetime Error [s]']/exposure_df['Lifetime [s]'] <= 0.3)]
            exposure_df['Exp Rate [mHz]'] = 1000 / exposure_df['Lifetime [s]']
            exposure_df['Exp Rate Sigma [mHz]'] = exposure_df['Lifetime Error [s]'] * 1000 / (
            exposure_df['Lifetime [s]']) ** 2

            exposure_df = exposure_df.drop(columns=['Exponential Fit 2xNLL',
                                                    'N.d.o.f.', 'Time Cut High [s]', 'Time Cut Low [s]'])
            exposure_df.to_csv(self.Cf_expA_sorted_path[i],index=False)


        for i in range(len(self.Cf_exp_116B_raw_path)):
            exposure_df = self.read_exposure(self.Cf_exp_116B_raw_path[i] + ".txt")
            exposure_df.columns = ['Pressure [bara]', 'Lifetime [s]', 'Lifetime Error [s]','Exponential Fit 2xNLL','N.d.o.f.','Time Cut High [s]','Time Cut Low [s]']
            exposure_df = exposure_df[(exposure_df['Lifetime [s]'] <=6.92e-1) | (exposure_df['Lifetime [s]'] >= 6.94e-1)]
            exposure_df = exposure_df[
                (exposure_df['Lifetime Error [s]']/exposure_df['Lifetime [s]'] <= 0.3)]
            exposure_df['Exp Rate [mHz]'] = 1000 / exposure_df['Lifetime [s]']
            exposure_df['Exp Rate Sigma [mHz]'] = exposure_df['Lifetime Error [s]'] * 1000 / (
            exposure_df['Lifetime [s]']) ** 2

            exposure_df = exposure_df.drop(columns=['Exponential Fit 2xNLL',
                                                    'N.d.o.f.', 'Time Cut High [s]', 'Time Cut Low [s]'])
            exposure_df.to_csv(self.Cf_expB_sorted_path[i],index=False)

    def read_raw_backgrounds_exp(self):
        # read file, delete unreasonable rows and rewrite
        for i in range(len(self.backgrounds_exp_raw_path)):
            exposure_df = self.read_exposure(self.backgrounds_exp_raw_path[i] + ".txt")
            exposure_df.columns = ['Pressure [bara]', 'Lifetime [s]', 'Lifetime Error [s]', 'Exponential Fit 2xNLL',
                                   'N.d.o.f.', 'Time Cut High [s]', 'Time Cut Low [s]']
            exposure_df = exposure_df[
                (exposure_df['Lifetime [s]'] <= 6.92e-1) | (exposure_df['Lifetime [s]'] >= 6.94e-1)]
            exposure_df = exposure_df[
                (exposure_df['Lifetime Error [s]'] / exposure_df['Lifetime [s]'] <= 0.3)]
            exposure_df['Bkg Rate [mHz]'] = 1000 / exposure_df['Lifetime [s]']
            exposure_df['Bkg Rate Sigma [mHz]'] = exposure_df['Lifetime Error [s]'] * 1000 / (
            exposure_df['Lifetime [s]']) ** 2
            exposure_df = exposure_df.drop(columns=['Exponential Fit 2xNLL',
                                   'N.d.o.f.', 'Time Cut High [s]', 'Time Cut Low [s]'])
            exposure_df.to_csv(self.Bkg_exp_sorted_path[i], index=False)
    def average_background_analysis(self):
        # np array operations, keep the pressure and average value if not NA
        bkg_df_116_list = []
        for i in range(0,self.background_exp_116_raw_len+1):
            bkg_df_116 = pd.read_csv(self.Bkg_exp_sorted_path[i])
            bkg_df_116_list.append(bkg_df_116)
        combined_df_116 = pd.concat(bkg_df_116_list, ignore_index=True)
        result_df_116 = combined_df_116.groupby('Pressure [bara]').agg({
            'Lifetime [s]': 'mean',  # Simple average
            'Lifetime Error [s]': self.calculate_rss,  # Custom square root math,
            'Bkg Rate [mHz]':'mean',
            'Bkg Rate Sigma [mHz]':self.calculate_rss
        }).reset_index()
        print('result_df_116',result_df_116)
        result_df_116.to_csv(self.Bkg_average_116_path, index=False)

        # add different source uplimit

        result_df_116_full_info = pd.merge(result_df_116, self.df_energy_116_tab, on='Pressure [bara]', how="inner")
        columns_added_Cs = result_df_116_full_info.apply(self.calculate_bkg_uplimit_by_row, axis=1, args=("Cs",))
        columns_added_Co = result_df_116_full_info.apply(self.calculate_bkg_uplimit_by_row, axis=1, args=("Co",))

        result_df_116_full_info = pd.concat([result_df_116_full_info, columns_added_Cs, columns_added_Co], axis=1)
        # print('result_df_116_full_info.columns',result_df_116_full_info.columns)
        result_df_116_full_info.to_csv(self.Bkg_average_116_full_info_path, index=False)


        #119
        bkg_df_119_list = []
        for i in range(self.background_exp_116_raw_len+1,len(self.Bkg_exp_sorted_path)):
            bkg_df_119 = pd.read_csv(self.Bkg_exp_sorted_path[i])
            bkg_df_119_list.append(bkg_df_119)
        combined_df_119 = pd.concat(bkg_df_119_list, ignore_index=True)
        result_df_119 = combined_df_119.groupby('Pressure [bara]').agg({
            'Lifetime [s]': 'mean',  # Simple average
            'Lifetime Error [s]': self.calculate_rss,  # Custom square root math
            'Bkg Rate [mHz]': 'mean',
            'Bkg Rate Sigma [mHz]': self.calculate_rss
        }).reset_index()

        print('result_df_119',result_df_119)
        result_df_119.to_csv(self.Bkg_average_119_path, index=False)

        result_df_119_full_info = pd.merge(result_df_119, self.df_energy_119_tab, on='Pressure [bara]', how="inner")
        columns_added_Cs = result_df_119_full_info.apply(self.calculate_bkg_uplimit_by_row, axis=1, args=("Cs",))
        columns_added_Co = result_df_119_full_info.apply(self.calculate_bkg_uplimit_by_row, axis=1, args=("Co",))

        result_df_119_full_info = pd.concat([result_df_119_full_info, columns_added_Cs, columns_added_Co], axis=1)
        result_df_119_full_info.to_csv(self.Bkg_average_119_full_info_path, index=False)

    def clean_signal_analysis(self):
        self.df_bkg_116 = pd.read_csv(self.Bkg_average_116_path)
        # self.df_bkg_116.columns = ['Pressure [bara]','Bkg Lifetime [s]','Bkg Lifetime Error [s]','Bkg Rate [mHz]', 'Bkg Rate Sigma [mHz]']
        self.df_bkg_116 = self.df_bkg_116.rename(columns={
            'Lifetime [s]': 'Bkg Lifetime [s]',
            'Lifetime Error [s]': 'Bkg Lifetime Error [s]'
        })
        self.df_bkg_119 =  pd.read_csv(self.Bkg_average_119_path)
        # self.df_bkg_119.columns = ['Pressure [bara]', 'Bkg Lifetime [s]', 'Bkg Lifetime Error [s]','Bkg Rate [mHz]', 'Bkg Rate Sigma [mHz]']
        self.df_bkg_119 = self.df_bkg_116.rename(columns={
            'Lifetime [s]': 'Bkg Lifetime [s]',
            'Lifetime Error [s]': 'Bkg Lifetime Error [s]'
        })
        # first these are 116
        self.Cs_116_data = []
        self.Cs_119_data = []
        self.Co_116_data = []
        self.Co_119_data = []
        # print('self.Cs_exp_sorted_path',self.Cs_exp_sorted_path)
        for i in range(0,self.Cs_exp_116_raw_len):
            print("Cs 116K", self.Cs_exp_rate_path[i])
            exposure_df = pd.read_csv(self.Cs_exp_sorted_path[i])
            # merge both has the pressure value, on pressure
            merged_df = pd.merge(self.df_bkg_116,exposure_df,on='Pressure [bara]', how="inner")
            # clean rate!
            merged_df['Clean Rate [mHz]'] = merged_df['Exp Rate [mHz]'] - merged_df['Bkg Rate [mHz]']
            merged_df['Clean Rate Sigma [mHz]'] = np.sqrt(merged_df['Exp Rate Sigma [mHz]']**2 + merged_df['Bkg Rate Sigma [mHz]']**2)
            # add Seitz and Eion unit
            merged_df = pd.merge(merged_df, self.df_energy_116_tab, on='Pressure [bara]', how="inner")
            # add sims analysis to get rejection

            merged_df.to_csv(self.Cs_exp_rate_path[i], index=False)

        for i in range(self.Cs_exp_116_raw_len,len(self.Cs_exp_raw_path) ):
            print("Cs 119K", self.Cs_exp_rate_path[i])
            exposure_df = pd.read_csv(self.Cs_exp_sorted_path[i])
            # merge both has the pressure value, on pressure
            merged_df = pd.merge(self.df_bkg_119,exposure_df,on='Pressure [bara]', how="inner")
            # clean rate!
            merged_df['Clean Rate [mHz]'] = merged_df['Exp Rate [mHz]'] - merged_df['Bkg Rate [mHz]']
            merged_df['Clean Rate Sigma [mHz]'] = np.sqrt(merged_df['Exp Rate Sigma [mHz]']**2 + merged_df['Bkg Rate Sigma [mHz]']**2)
            # add Seitz and Eion unit
            merged_df = pd.merge(merged_df, self.df_energy_119_tab, on='Pressure [bara]', how="inner")
            # add sims analysis to get rejection

            merged_df.to_csv(self.Cs_exp_rate_path[i], index=False)

        for i in range(0, self.Co_exp_116_raw_len):
            print("Co 116K", self.Co_exp_rate_path[i])
            exposure_df = pd.read_csv(self.Co_exp_sorted_path[i])
            # merge both has the pressure value, on pressure
            merged_df = pd.merge(self.df_bkg_116, exposure_df, on='Pressure [bara]', how="inner")
            # clean rate!
            merged_df['Clean Rate [mHz]'] = merged_df['Exp Rate [mHz]'] - merged_df['Bkg Rate [mHz]']
            merged_df['Clean Rate Sigma [mHz]'] = np.sqrt(
                merged_df['Exp Rate Sigma [mHz]'] ** 2 + merged_df['Bkg Rate Sigma [mHz]'] ** 2)
            # add Seitz and Eion unit
            merged_df = pd.merge(merged_df, self.df_energy_116_tab, on='Pressure [bara]', how="inner")
            # add sims analysis to get rejection

            merged_df.to_csv(self.Co_exp_rate_path[i], index=False)

        for i in range(self.Co_exp_116_raw_len, len(self.Co_exp_raw_path)):
            print("Co 119K", self.Co_exp_rate_path[i])
            exposure_df = pd.read_csv(self.Co_exp_sorted_path[i])
            # merge both has the pressure value, on pressure
            merged_df = pd.merge(self.df_bkg_119, exposure_df, on='Pressure [bara]', how="inner")
            # clean rate!
            merged_df['Clean Rate [mHz]'] = merged_df['Exp Rate [mHz]'] - merged_df['Bkg Rate [mHz]']
            merged_df['Clean Rate Sigma [mHz]'] = np.sqrt(
                merged_df['Exp Rate Sigma [mHz]'] ** 2 + merged_df['Bkg Rate Sigma [mHz]'] ** 2)
            # add Seitz and Eion unit
            merged_df = pd.merge(merged_df, self.df_energy_119_tab, on='Pressure [bara]', how="inner")
            # add sims analysis to get rejection

            merged_df.to_csv(self.Co_exp_rate_path[i], index=False)

    def clean_NR_signal_analysis(self):
        self.df_bkg_116 = pd.read_csv(self.Bkg_average_116_path)
        # self.df_bkg_116.columns = ['Pressure [bara]','Bkg Lifetime [s]','Bkg Lifetime Error [s]','Bkg Rate [mHz]', 'Bkg Rate Sigma [mHz]']
        self.df_bkg_116 = self.df_bkg_116.rename(columns={
            'Lifetime [s]': 'Bkg Lifetime [s]',
            'Lifetime Error [s]': 'Bkg Lifetime Error [s]'
        })
        self.df_bkg_119 =  pd.read_csv(self.Bkg_average_119_path)
        # self.df_bkg_119.columns = ['Pressure [bara]', 'Bkg Lifetime [s]', 'Bkg Lifetime Error [s]','Bkg Rate [mHz]', 'Bkg Rate Sigma [mHz]']
        self.df_bkg_119 = self.df_bkg_116.rename(columns={
            'Lifetime [s]': 'Bkg Lifetime [s]',
            'Lifetime Error [s]': 'Bkg Lifetime Error [s]'
        })
        # first these are 116
        self.Cf_116A_data = []
        self.Cf_116B_data = []
        # print('self.Cs_exp_sorted_path',self.Cs_exp_sorted_path)
        for i in range(0,self.Cf_exp_116A_raw_len):
            print("Cf 116K A", self.Cf_expA_rate_path[i])
            exposure_df = pd.read_csv(self.Cf_expA_sorted_path[i])
            # merge both has the pressure value, on pressure
            merged_df = pd.merge(self.df_bkg_116,exposure_df,on='Pressure [bara]', how="inner")
            # clean rate!
            merged_df['Clean Rate [mHz]'] = merged_df['Exp Rate [mHz]'] - merged_df['Bkg Rate [mHz]']
            merged_df['Clean Rate Sigma [mHz]'] = np.sqrt(merged_df['Exp Rate Sigma [mHz]']**2 + merged_df['Bkg Rate Sigma [mHz]']**2)
            # add Seitz and Eion unit
            merged_df = pd.merge(merged_df, self.df_energy_116_tab, on='Pressure [bara]', how="inner")
            # add sims analysis to get rejection

            merged_df.to_csv(self.Cf_expA_rate_path[i], index=False)

        for i in range(0,self.Cf_exp_116B_raw_len):
            print("Cf 116K B", self.Cf_expB_rate_path[i])
            exposure_df = pd.read_csv(self.Cf_expB_sorted_path[i])
            # merge both has the pressure value, on pressure
            merged_df = pd.merge(self.df_bkg_116,exposure_df,on='Pressure [bara]', how="inner")
            # clean rate!
            merged_df['Clean Rate [mHz]'] = merged_df['Exp Rate [mHz]'] - merged_df['Bkg Rate [mHz]']
            merged_df['Clean Rate Sigma [mHz]'] = np.sqrt(merged_df['Exp Rate Sigma [mHz]']**2 + merged_df['Bkg Rate Sigma [mHz]']**2)
            # add Seitz and Eion unit
            merged_df = pd.merge(merged_df, self.df_energy_116_tab, on='Pressure [bara]', how="inner")
            # add sims analysis to get rejection

            merged_df.to_csv(self.Cf_expB_rate_path[i], index=False)

        fig, ax = plt.subplots(1, 2, figsize=(10, 4))
        # 1 plot to compare with original data Gray had, 2 to plot the spectrum with clean data comparasion
        for i in range(len(self.Cf_expA_rate_path)):
            exp_df = pd.read_csv(self.Cf_expA_rate_path[i])
            if i==0:
                ax[0].errorbar(exp_df['Pressure [bara]'], exp_df["Bkg Rate [mHz]"],
                               yerr=exp_df['Bkg Rate Sigma [mHz]'], label="combined bkg 116.7 K ", fmt='o', color='blue')
            ax[0].errorbar(exp_df['Pressure [bara]'],exp_df["Exp Rate [mHz]"],
                   yerr = exp_df["Exp Rate Sigma [mHz]"], label=f"Cf config A {i}",fmt = 'o')
            ax[1].errorbar(exp_df["Seitz [keV]"]*1000,exp_df["Clean Rate [mHz]"],
                   yerr = exp_df["Clean Rate Sigma [mHz]"], label=f"Cf config A {i}",fmt = 'o')

        for i in range(len(self.Cf_expB_rate_path)):
            exp_df = pd.read_csv(self.Cf_expB_rate_path[i])
            ax[0].errorbar(exp_df['Pressure [bara]'],exp_df["Exp Rate [mHz]"],
                   yerr = exp_df["Exp Rate Sigma [mHz]"], label=f"Cf config B {i}",fmt = 'o')
            ax[1].errorbar(exp_df["Seitz [keV]"]*1000,exp_df["Clean Rate [mHz]"],
                   yerr = exp_df["Clean Rate Sigma [mHz]"], label=f"Cf config B {i}",fmt = 'o')


        ax[1].plot(self.Cf_simsA[1][0][1][:-1],self.Cf_simsA[0]*self.Cf_simsA[2], label='Cf configA spectrum')
        ax[1].plot(self.Cf_simsB[1][0][1][:-1], self.Cf_simsB[0] * self.Cf_simsB[2], label='Cf configB spectrum')

        ax[0].set_xlabel('Pressure [bara]')
        ax[0].set_ylabel("Exp Rate [mHz]")
        ax[0].set_title("Check with exp rate")
        ax[0].set_xlim(1.75,5.25)
        ax[0].set_ylim(10, 100)
        ax[0].legend()

        ax[1].set_xlabel("Seitz [eV]")
        ax[1].set_ylabel("Clean Rate [mHz]")
        ax[1].set_title("Clean Rate Compare with sims")
        ax[1].set_xlim(-1, 3500)
        ax[1].set_ylim(0,600)
        ax[1].legend()

        plt.savefig(self.plot_path + "Cf_abs_rate_comparison.pdf")

    def read_Seitz_info(self):
        Seitz_pressure_list = np.arange(1.25, 6.5, 0.25)

        Seitz_116 = [0.8318354532105874, 0.8940418766838347, 0.9631954092292087, 1.0403343615139717, 1.1266932576571842,
                 1.2237481550986595, 1.3332743449138087, 1.4574206436587642, 1.5988055614249703, 1.7606432313822162,
                 1.946909785273693, 2.16256537886239, 2.4138541417933888, 2.708714138448383, 3.0573453143201323,
                 3.4730076036964754, 3.973160671933538, 4.5811194010067275, 5.328505251504318, 6.258953124220718,
                 7.43385064633319]  # in keV
        # keV
        E_ion_116 = [0.4585087383615149, 0.48723219835949405, 0.5187650745694886, 0.5534854521564034,
                 0.5918369935015713, 0.6343430463424671, 0.6816243569175929, 0.734421603748339, 0.7936241450592315,
                 0.8603071683211787, 0.9357800765380097, 1.0216501032567107, 1.1199070007748584, 1.2330370047313486,
                 1.3641782404241205, 1.5173355596283118, 1.6976822638874727, 1.9119908079547858, 2.1692593120919246,
                 2.481641649854045, 2.865860421381346]
        # g / cc
        rho_l_116 = [1.1837706571559556, 1.183965618191963, 1.1841602684182295, 1.1843546090988741, 1.1845486414895372,
                 1.1847423668374562, 1.1849357863815444, 1.18512890135247, 1.1853217129727334, 1.1855142224567425,
                 1.185706431010887, 1.1858983398336123, 1.186089950115495, 1.186281263039315, 1.1864722797801222,
                 1.186663001505315, 1.1868534293747037, 1.1870435645405828, 1.1872334081477973, 1.1874229613338119,
                 1.1876122252287777]
        # nm
        Rl_116 = [4.665037623181469, 4.8031486464831294, 4.949569693967424, 5.105073840504335, 5.270533085793664,
              5.446934770432428, 5.635401216931708, 5.8372137252300575, 6.0538416388377465, 6.286978214102999,
              6.538585000984282, 6.810947072223535, 7.106742667816778, 7.429131533695626, 7.781868241770823,
              8.169448995079783, 8.597304383423044, 9.072055318279228, 9.6018581178134, 10.19687689111723,
              10.869942000948551]

        compound_x_116 = []
        for i in range(len(E_ion_116)):
            x = E_ion_116[i] * 10 / (rho_l_116[i] * Rl_116[i])  # fit unit
            compound_x_116.append(x)

        self.dict_energy_116_tab = {"Pressure [bara]":Seitz_pressure_list,
                                    "Seitz [keV]": Seitz_116,
                                    "Eion [keV]": E_ion_116,
                                    "Eion_rl-1_rhol-1 [GeVcm**2 g-1]": compound_x_116}
        self.df_energy_116_tab = pd.DataFrame(self.dict_energy_116_tab)

        Seitz_119 = [0.4336397431016649, 0.45977540397438443, 0.4882700006652598, 0.5194078524821146,
                 0.5535164255731408, 0.5909743032392646, 0.6322208776848458, 0.677768234751849, 0.7282157681331614,
                 0.7842683220836715, 0.8467588127318099, 0.9166765658125957, 0.9952031955105735, 1.08375825146003,
                 1.1840578053726714, 1.2981902707374324, 1.4287154214193478, 1.5787948930230866, 1.7523661419564087,
                 1.9543766809995484, 2.1911035165969657]  # in keV
        # keV
        E_ion_119 = [0.25394281195661755, 0.266864943451393, 0.28081328330244465, 0.2958994998835466,
                 0.3122508371322994, 0.3300127984753024, 0.3493523708721229, 0.37046193856134035,
                 0.39356403694491504, 0.4189171955474377, 0.4468231388310798, 0.4776356895155958, 0.511771919112016,
                 0.5497261516000658, 0.592087711318551, 0.6395635940651908, 0.6930076840040227, 0.7534587276853615,
                 0.8221902660987933, 0.9007768869746815, 0.9911832613763831]
        # g / cc
        rho_l_119 = [1.1570041272719827, 1.1572302529996157, 1.1574559403955533, 1.1576811916287815, 1.1579060088506121,
                 1.1581303941948926, 1.1583543497782014, 1.1585778777000375, 1.1588009800430248, 1.1590236588730902,
                 1.1592459162396547, 1.1594677541758232, 1.1596891746985576, 1.159910179808861, 1.1601307714919509,
                 1.160350951717438, 1.1605707224394939, 1.1607900855970206, 1.1610090431138194, 1.161227596898751,
                 1.161445748845904]
        # nm
        Rl_119 = [3.704105651173001, 3.794992905503335, 3.890374002589942, 3.990590522026317, 4.096019535068052,
              4.20707837799361, 4.324230158735822, 4.447990223912374, 4.578933693884623, 4.717704425633646,
              4.865025627595007, 5.021712392413753, 5.1886868283775565, 5.366996201410651, 5.557834915522947,
              5.762571294147585, 5.98278044731765, 6.22028480464255, 6.477204741359901, 6.756021989516639,
              7.05966008057627]

        compound_x_119 = []
        for i in range(len(E_ion_119)):
            x = E_ion_119[i] * 10 / (rho_l_119[i] * Rl_119[i])  # fit unit
            compound_x_119.append(x)

        self.dict_energy_119_tab = {"Pressure [bara]": Seitz_pressure_list,
                                    "Seitz [keV]": Seitz_119,
                                    "Eion [keV]": E_ion_119,
                                    "Eion_rl-1_rhol-1 [GeVcm**2 g-1]": compound_x_119}
        self.df_energy_119_tab = pd.DataFrame(self.dict_energy_119_tab)

    def gamma_rejection_calculation(self):
        
        for i in range(len(self.Cs_exp_rate_path)):
            exp_df = pd.read_csv(self.Cs_exp_rate_path[i])
            columns_added  =exp_df.apply(self.calculate_rejection_by_row,axis=1,args=("Cs",))
            merged_df = pd.concat([exp_df, columns_added], axis=1)
            # print('Cs print(merged_df)',self.Cs_exp_rate_path[i],'\n',merged_df)
            merged_df.to_csv(self.Cs_exp_rejection_path[i], index= False)

        for i in range(len(self.Co_exp_rate_path)):
            exp_df = pd.read_csv(self.Co_exp_rate_path[i])
            columns_added  =exp_df.apply(self.calculate_rejection_by_row,axis=1,args=("Co",))
            merged_df = pd.concat([exp_df, columns_added], axis=1)
            # print('Co print(merged_df)',self.Co_exp_rate_path[i],'\n',merged_df)
            merged_df.to_csv(self.Co_exp_rejection_path[i], index= False)
            



    def read_exposure(self,filename):
        # Define your path (we'll use a relative path)
        file_path = os.path.join('..', 'exp_exposure', filename)

        # Read the file
        # sep='\s+' handles any number of spaces or tabs as delimiters
        df = pd.read_csv(file_path, sep='\s+', skiprows=1, header=None)

        return df
    def calculate_bkg_uplimit_by_row(self, row, source):

        if source == "Co":
            self.sim_list = self.Co_sims
        elif source == "Cs":
            self.sim_list = self.Cs_sims
        else:
            print("NA sources")

        self.Rate_factor = self.sim_list[0]
        # print("self.Rate_factor",self.Rate_factor)
        self.energy_edges = self.sim_list[1][0][1]
        # print("self.energy_edges", self.energy_edges)
        self.counts_cum_bin = self.sim_list[2]
        # print("self.counts_cum_bin", self.counts_cum_bin)
        self.counts_energy_cum_bin = self.sim_list[3]
        # print("self.counts_energy_cum_bin", self.counts_energy_cum_bin)

        rejection_uplimit_PS = 0
        rejection_uplimit_PK = 0
        for i in range(len(self.energy_edges)):
            if row['Seitz [keV]'] >= self.energy_edges[i]:
                # rejection per scattering, PS meaning perscattering
                counts = self.counts_cum_bin[i] + (row['Seitz [keV]'] - self.energy_edges[i]) * (
                        self.counts_cum_bin[i + 1] -
                        self.counts_cum_bin[i]) / (
                                 self.energy_edges[i + 1] - self.energy_edges[i])
                rate_PS = self.Rate_factor * (counts)

                rate_PS_sigma = rate_PS / np.sqrt(counts)


                rejection_uplimit_PS = row['Bkg Rate Sigma [mHz]']/rate_PS


                break
        for i in range(len(self.energy_edges)):
            if row['Eion [keV]'] >= self.energy_edges[i]:
                counts_times_keV = self.counts_energy_cum_bin[0]  # all energy
                counts_Eion = self.counts_cum_bin[i] + (row['Eion [keV]'] - self.energy_edges[i]) * (
                        self.counts_cum_bin[i + 1] -
                        self.counts_cum_bin[i]) / (
                                      self.energy_edges[i + 1] - self.energy_edges[i])

                rate_PK = self.Rate_factor * (counts_times_keV)
                rate_PK_sigma = rate_PK / np.sqrt(counts_Eion)

                rejection_uplimit_PK = row['Bkg Rate Sigma [mHz]'] / rate_PK
                break
        output = pd.Series({f"{source} Rejection Uplimit Scattering []":rejection_uplimit_PS,
        f"{source} Rejection Uplimit KeV [/keV]":rejection_uplimit_PK})

        return output
    def calculate_rejection_by_row(self, row, source):
        if source == "Co":
            self.sim_list = self.Co_sims
        elif source == "Cs":
            self.sim_list = self.Cs_sims
        else:
            print("NA sources")

        self.Rate_factor = self.sim_list[0]
        # print("self.Rate_factor",self.Rate_factor)
        self.energy_edges = self.sim_list[1][0][1]
        # print("self.energy_edges", self.energy_edges)
        self.counts_cum_bin = self.sim_list[2]
        # print("self.counts_cum_bin", self.counts_cum_bin)
        self.counts_energy_cum_bin = self.sim_list[3]
        # print("self.counts_energy_cum_bin", self.counts_energy_cum_bin)
        rejection_PS = 0
        rejection_PS_sigma = 0
        rejection_PK = 0
        rejection_PK_sigma = 0
        rejection_uplimit_PS = 0
        rejection_uplimit_PK = 0
        for i in range(len(self.energy_edges)):
            if row['Seitz [keV]']>= self.energy_edges[i]:
                # rejection per scattering, PS meaning perscattering
                counts = self.counts_cum_bin[i] + (row['Seitz [keV]'] - self.energy_edges[i]) * (
                        self.counts_cum_bin[i + 1] -
                        self.counts_cum_bin[i]) / (
                                 self.energy_edges[i + 1] - self.energy_edges[i])
                rate_PS = self.Rate_factor * (counts)

                rate_PS_sigma = rate_PS / np.sqrt(counts)

                rejection_PS = row['Clean Rate [mHz]'] / rate_PS
                # rejection_uplimit_PS = row['Bkg Rate Sigma [mHz]']/rate_PS

                # will be returned
                rejection_PS_sigma = np.sqrt(
                    (row['Clean Rate Sigma [mHz]'] / rate_PS) ** 2 + (row['Clean Rate [mHz]'] * rate_PS_sigma / rate_PS ** 2) ** 2)
                # will be returned
                break
        for i in range(len(self.energy_edges)):
            if row['Eion [keV]'] >= self.energy_edges[i]:
                counts_times_keV = self.counts_energy_cum_bin[0]  # all energy
                counts_Eion = self.counts_cum_bin[i] + (row['Eion [keV]'] - self.energy_edges[i]) * (
                        self.counts_cum_bin[i + 1] -
                        self.counts_cum_bin[i]) / (
                                      self.energy_edges[i + 1] - self.energy_edges[i])

                rate_PK = self.Rate_factor * (counts_times_keV)
                rate_PK_sigma = rate_PK / np.sqrt(counts_Eion)
                rejection_PK = row['Clean Rate [mHz]'] / rate_PK

                rejection_PK_sigma = np.sqrt(
                    (row['Clean Rate Sigma [mHz]'] / rate_PK) ** 2 + (row['Clean Rate [mHz]'] * rate_PK_sigma / rate_PK ** 2) ** 2)
                # rejection_uplimit_PK = row['Bkg Rate Sigma [mHz]'] / rate_PK
                break
        output = pd.Series({"Rejection Rate Scattering[]": rejection_PS,
                "Rejection Sigma Scattering[]": rejection_PS_sigma,
                "Rejection Rate KeV[/keV]": rejection_PK,
                "Rejection Sigma KeV[/keV]": rejection_PK_sigma})
        return output

    def bkg_plot(self):
        self.df_bkg_116 = pd.read_csv(self.Bkg_average_116_path)
        self.df_bkg_116 = pd.merge(self.df_bkg_116, self.df_energy_116_tab, on='Pressure [bara]', how="inner")
        print('self.df_bkg_116',self.df_bkg_116)
        self.df_bkg_119 = pd.read_csv(self.Bkg_average_119_path)
        self.df_bkg_119 = pd.merge(self.df_bkg_119, self.df_energy_119_tab, on='Pressure [bara]', how="inner")
        print('self.df_bkg_119', self.df_bkg_119)

        fig, ax = plt.subplots()
        ax.errorbar(self.df_bkg_116['Seitz [keV]'],self.df_bkg_116["Bkg Rate [mHz]"],
                   yerr = self.df_bkg_116["Bkg Rate Sigma [mHz]"],label="combined bkg 116.7 K ",fmt = 'o',color = 'r')
        ax.errorbar(self.df_bkg_119['Seitz [keV]'], self.df_bkg_119["Bkg Rate [mHz]"],
                    yerr=self.df_bkg_119["Bkg Rate Sigma [mHz]"], label="combined bkg 119.6 K ",fmt = 'o', color='b')
        ax.set_xlim(0.4,3.6)
        ax.set_ylim(5,55)
        ax.set_xlabel("Seitz [keV]")
        ax.set_ylabel("Bkg Rate [mHz]")
        ax.legend()

        plt.savefig(self.plot_path + "average_bkg_rate.pdf")

    def gamma_rejection_plot(self):
        fig, ax = plt.subplots(1, 2, figsize=(10, 4))
        fig, ax = plt.subplots(2, 1, figsize=(6, 10))
        self.fitting_list = []


        for i in range(len(self.Cs_exp_rejection_path)):
            df = pd.read_csv(self.Cs_exp_rejection_path[i])
            # print(df.columns)
            doc_label = self.Cs_exp_raw_path[i].replace('_exposures', '')
            print('doc_label',doc_label)
            # signal
            # drop 2.75,3.25, 3.75 bara pressure
            # pressure_drop_list = [2.75,3.25,3.75]
            pressure_drop_list = []
            df = df[~df['Pressure [bara]'].isin(pressure_drop_list)]
            # only positive rate
            df =  df[df['Clean Rate [mHz]']>0]

            df_fit = df[['Seitz [keV]',"Rejection Rate Scattering[]",'Eion_rl-1_rhol-1 [GeVcm**2 g-1]',"Rejection Rate KeV[/keV]"]]
            self.fitting_list.append(df_fit)



            ax[0].errorbar(df['Seitz [keV]'], df["Rejection Rate Scattering[]"],
                           yerr=df["Rejection Sigma Scattering[]"], label=doc_label, fmt='o')

            ax[1].errorbar(df['Eion_rl-1_rhol-1 [GeVcm**2 g-1]'], df["Rejection Rate KeV[/keV]"],
                           yerr=df["Rejection Sigma KeV[/keV]"], label=doc_label, fmt='o')

        for i in range(len(self.Co_exp_rejection_path)):
            df = pd.read_csv(self.Co_exp_rejection_path[i])
            # print(df.columns)
            doc_label = self.Co_exp_raw_path[i].rstrip("_exposures")
            # signal
            # drop 2.75,3.25, 3.75 bara pressure
            # pressure_drop_list = [2.75,3.25,3.75]
            pressure_drop_list = []
            df = df[~df['Pressure [bara]'].isin(pressure_drop_list)]
            df = df[df['Clean Rate [mHz]'] > 0]

            df_fit = df[['Seitz [keV]', "Rejection Rate Scattering[]", 'Eion_rl-1_rhol-1 [GeVcm**2 g-1]',
                         "Rejection Rate KeV[/keV]"]]
            self.fitting_list.append(df_fit)

            ax[0].errorbar(df['Seitz [keV]'], df["Rejection Rate Scattering[]"],
                           yerr=df["Rejection Sigma Scattering[]"], label=doc_label, fmt='o')


            ax[1].errorbar(df['Eion_rl-1_rhol-1 [GeVcm**2 g-1]'], df["Rejection Rate KeV[/keV]"],
                           yerr=df["Rejection Sigma KeV[/keV]"], label=doc_label, fmt='o')


        self.fitting_df =  pd.concat(self.fitting_list, ignore_index=True)
        [(a_fit_scatter, b_fit_scatter,x_fitted_scatter,y_fitted_scatter),(a_fit_keV, b_fit_keV,x_fitted_keV,y_fitted_keV)] = self.fitting_gamma_rejection()

        # plot the fitting function
        ax[0].plot(x_fitted_scatter,y_fitted_scatter,label = f"a,b = {a_fit_scatter:.2e} , {b_fit_scatter:.2e}", color="black")

        #gamma rejection up limit
        self.bkg_floor_plot(ax[0],"Seitz")

        ax[0].set_xlabel("Seitz [keV]")
        ax[0].set_ylabel("Gamma Rejection Per Scattering []")
        ax[0].set_title("Gamma Rejection Per Scattering ")
        # ax[0].set_ylim(1.0e-12,1.0e-2)
        # ax[0].set_xlim(0,6)
        # ax[0].set_xlim(0.8,1.5)
        # ax[0].set_ylim(1.0e-12,1.0e-2)
        ax[0].set_yscale("log")

        ax[0].legend(loc='upper right', fontsize=7)

        ax[1].plot(x_fitted_keV, y_fitted_keV, label=f"a,b = {a_fit_keV:.2e} , {b_fit_keV:.2e}", color="black")

        self.bkg_floor_plot(ax[1], "Eion")
        ax[1].set_xlabel("Eion_rl-1_rhol-1 [GeVcm**2 g-1]")
        ax[1].set_ylabel("Gamma Rejection Per keV [/keV]")
        ax[1].set_title("Gamma Rejection Per keV ")
        # ax[1].set_ylim(1.0e-14,1.0e-4)
        # ax[1].set_xlim(0.08,0.15)
        # ax[1].set_xlim(0.8,1.1)
        ax[1].set_yscale("log")
        ax[1].legend(loc='upper right', fontsize=7)

        plt.savefig(self.plot_path + "gamma_rejection.pdf")
    def bkg_floor_plot(self,ax,mode):
        self.df_bkg_116_full_info = pd.read_csv(self.Bkg_average_116_full_info_path)
        # self.df_bkg_116_full_info = pd.merge(self.df_bkg_116_full_info, self.df_energy_116_tab, on='Pressure [bara]', how="inner")
        self.df_bkg_119_full_info = pd.read_csv(self.Bkg_average_119_full_info_path)
        # self.df_bkg_119_full_info = pd.merge(self.df_bkg_119_full_info, self.df_energy_119_tab, on='Pressure [bara]', how="inner")
        # print('self.df_bkg_116_full_info.columns',self.df_bkg_116_full_info.columns)
        if mode == "Seitz":
            ax.plot(self.df_bkg_116_full_info['Seitz [keV]'], self.df_bkg_116_full_info['Cs Rejection Uplimit Scattering []'],
                       color="gray")
            # ax.plot(self.df_bkg_116_full_info['Seitz [keV]'], self.df_bkg_116_full_info['Co Rejection Uplimit Scattering []'],
            #            color="gray")

            # ax.plot(self.df_bkg_119_full_info['Seitz [keV]'], self.df_bkg_119_full_info['Cs Rejection Uplimit Scattering []'],
            #         color="gray")
            # ax.plot(self.df_bkg_119_full_info['Seitz [keV]'], self.df_bkg_119_full_info['Co Rejection Uplimit Scattering []'],
            #         color="gray")
        elif mode =="Eion":
            ax.plot(self.df_bkg_116_full_info['Eion_rl-1_rhol-1 [GeVcm**2 g-1]'], self.df_bkg_116_full_info['Cs Rejection Uplimit KeV [/keV]'],
                    color="gray")
            # ax.plot(self.df_bkg_116_full_info['Eion_rl-1_rhol-1 [GeVcm**2 g-1]'], self.df_bkg_116_full_info['Co Rejection Uplimit KeV [/keV]'],
            #         color="gray")

            # ax.plot(self.df_bkg_119_full_info['Eion_rl-1_rhol-1 [GeVcm**2 g-1]'], self.df_bkg_119_full_info['Cs Rejection Uplimit KeV [/keV]'],
            #         color="gray")
            # ax.plot(self.df_bkg_119_full_info['Eion_rl-1_rhol-1 [GeVcm**2 g-1]'], self.df_bkg_119_full_info['Co Rejection Uplimit KeV [/keV]'],
            #         color="gray")


    def spectrums_plot(self):

        self.Cs_Rate_factor = self.Cs_sims[0]
        self.Cs_energy_edges = self.Cs_sims[1][0][1]
        self.Cs_counts_cum_bin = self.Cs_sims[2]
        self.Cs_counts_energy_cum_bin = self.Cs_sims[3]

        self.Co_Rate_factor = self.Co_sims[0]
        self.Co_energy_edges = self.Co_sims[1][0][1]
        self.Co_counts_cum_bin = self.Co_sims[2]
        self.Co_counts_energy_cum_bin = self.Co_sims[3]

        fig, ax = plt.subplots(2, 2, figsize=(12, 12))
        ax[0,0].bar(self.Cs_energy_edges[:-1],self.Cs_counts_cum_bin,width=np.diff(self.Cs_energy_edges),
        align="edge",
        edgecolor="black")
        ax[0, 0].set_xlabel("ER Deposition [keV]")
        ax[0, 0].set_ylabel("Cumulative Counts []")
        ax[0, 0].set_title("Cs Cumulative Spectrum (Counts)")
        ax[0, 0].set_yscale("log")

        ax[0, 1].bar(self.Cs_energy_edges[:-1], self.Cs_counts_energy_cum_bin,width=np.diff(self.Cs_energy_edges),
        align="edge",
        edgecolor="black")
        ax[0, 1].set_xlabel("ER Deposition [keV]")
        ax[0, 1].set_ylabel("Cumulative Energy Deposited[keV]")
        ax[0, 1].set_title("Cs Cumulative Spectrum(Counts*energy)")
        ax[0, 1].set_yscale("log")

        ax[1, 0].bar(self.Co_energy_edges[:-1], self.Co_counts_cum_bin,width=np.diff(self.Co_energy_edges),
        align="edge",
        edgecolor="black")
        ax[1, 0].set_xlabel("ER Deposition [keV]")
        ax[1, 0].set_ylabel("Cumulative Counts []")
        ax[1, 0].set_title("Co Cumulative Spectrum (Counts)")
        ax[1, 0].set_yscale("log")

        ax[1, 1].bar(self.Co_energy_edges[:-1], self.Co_counts_energy_cum_bin,width=np.diff(self.Co_energy_edges),
        align="edge",
        edgecolor="black")
        ax[1, 1].set_xlabel("ER Deposition [keV]")
        ax[1, 1].set_ylabel("Cumulative Energy Deposited[keV]")
        ax[1, 1].set_title("Co Cumulative Spectrum(Counts*energy)")
        ax[1, 1].set_yscale("log")

        plt.savefig(self.plot_path + "cumulative_spectrums.pdf")

    def fitting_gamma_rejection(self):
        #
        x_per_scatter = self.fitting_df["Seitz [keV]"].values
        y_per_scatter = self.fitting_df["Rejection Rate Scattering[]"].values
        # dealing with guess
        x_min_per_scattering= min(x_per_scatter)
        x_max_per_scattering = max(x_per_scatter)
        y_min_per_scattering = min(y_per_scatter)
        y_max_per_scattering = max(y_per_scatter)
         # b is negative
        b_guess_per_scattering=-(np.log(y_max_per_scattering)-np.log(y_min_per_scattering))/(x_max_per_scattering-x_min_per_scattering)
        a_guess_scattering = y_max_per_scattering
        initial_guess_scatter = [a_guess_scattering, b_guess_per_scattering]
        popt_scatter, pcov_scatter = curve_fit(self.exp_func, x_per_scatter, y_per_scatter, p0=initial_guess_scatter)
        a_fit_scatter, b_fit_scatter= popt_scatter
        print('a_fit_scatter, b_fit_scatter',a_fit_scatter, b_fit_scatter)
        x_fitted_scatter = np.linspace(min(x_per_scatter), max(x_per_scatter), 100)
        y_fitted_scatter = self.exp_func(x_fitted_scatter, *popt_scatter)





        x_per_keV = self.fitting_df["Eion_rl-1_rhol-1 [GeVcm**2 g-1]"].values
        y_per_keV = self.fitting_df["Rejection Rate KeV[/keV]"].values
        # dealing with guess
        x_min_per_keV = min(x_per_keV)
        x_max_per_keV = max(x_per_keV)
        y_min_per_keV = min(y_per_keV)
        y_max_per_keV = max(y_per_keV)
        # b is negative
        b_guess_per_keV = -(np.log(y_max_per_keV) - np.log(y_min_per_keV)) / (
                    x_max_per_keV - x_min_per_keV)
        a_guess_scattering = y_max_per_keV
        initial_guess_keV = [a_guess_scattering, b_guess_per_keV]
        popt_keV, pcov_keV = curve_fit(self.exp_func, x_per_keV, y_per_keV, p0=initial_guess_keV)
        a_fit_keV, b_fit_keV = popt_keV
        print('a_fit_keV, b_fit_keV', a_fit_keV, b_fit_keV)
        x_fitted_keV = np.linspace(min(x_per_keV), max(x_per_keV), 100)
        y_fitted_keV = self.exp_func(x_fitted_keV, *popt_keV)

        return [(a_fit_scatter, b_fit_scatter,x_fitted_scatter,y_fitted_scatter),(a_fit_keV, b_fit_keV,x_fitted_keV,y_fitted_keV)]



    def calculate_rss(self, series):
        """Calculates sqrt(a^2 + b^2 + ...)"""
        return np.sqrt(np.sum(series ** 2))

    def exp_func(self, x, a, b):
        return a * np.exp(-b * x)





if __name__=="__main__":
    IA =  integrated_analysis()
    # test = test_csv()