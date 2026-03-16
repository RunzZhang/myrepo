import pandas as pd
import matplotlib.pyplot as plt
import csv
import numpy as np
import os
import pickle
class integrated_analysis():
    def __init__(self):

        self.output_path = '/data/runzezhang/result/gamma_rejection/'
        self.Co_sim_path  ='/data/runzezhang/result/TN_sims_D/Co_output.pkl'
        self.Cs_sim_path = '/data/runzezhang/result/TN_sims_D/Cs_output.pkl'

        self.Cs_exp_116_raw_path = ["Cold-Cs-11_17-18_exposures_mix","Cold-Cs-12_01_exposures_mix",
                                    "Cold-Cs-12_10-11_exposures_mix","Cold-Cs-1_20-21_exposures_mix"]
        self.Cs_exp_116_raw_len = len(self.Cs_exp_116_raw_path)
        self.Cs_exp_119_raw_path = ["Cold-Cs-2_2-3_exposures_zoom"]
        self.Co_exp_116_raw_path = ["60Co-12_15-16_exposures"]
        self.Co_exp_116_raw_len = len(self.Co_exp_116_raw_path)
        self.Co_exp_119_raw_path = []
        self.background_116_sorted_path = ["Background-11_26-30_exposures","QUIET-Background-1_13_exposures_zoom"]
        self.background_exp_116_raw_len = len(self.background_116_sorted_path)
        self.background_119_sorted_path = ["Background-2_6-12_exposures","Background-1_30-2_2_exposures"]

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

        # read experimental txt file, drop the non-sense values, and write to clean dataframe
        self.read_raw_Co_exp()
        self.read_raw_Cs_exp()
        self.read_raw_backgrounds_exp()

        # caculate the average background
        self.average_background_analysis()

        # read the Setiz energy thresholds
        self.read_Setiz_info()
        # add the Setiz energy to the exp txt files, and Setiz should have already included all temperature info, so in post-analysis
        # no demands to devide by temperature configurations
        # also caculate the clean signal and signal uncerntainty
        self.clean_signal_analysis()


        # based on sims and clean signal rate, calculate gamma rejection
        self.gamma_rejection_calculation()


    def generate_path(self):
        self.Co_exp_sorted_path = []
        self.Cs_exp_sorted_path = []
        self.Co_exp_rate_path = []
        self.Cs_exp_rate_path = []
        self.Co_exp_rejection_path = []
        self.Cs_exp_rejection_path = []
        self.Bkg_exp_sorted_path = []
        self.Bkg_average_116_path = self.output_path + "background_116_average" + ".csv"
        self.Bkg_average_119_path = self.output_path + "background_119_average" + ".csv"

        for exp_name in self.Co_exp_raw_path:
            self.Co_exp_sorted_path.append(self.output_path + exp_name+"_sorted.csv")
        for exp_name in self.Cs_exp_raw_path:
            self.Cs_exp_sorted_path.append(self.output_path + exp_name+"_sorted.csv")
        for exp_name in self.Co_exp_raw_path:
            self.Co_exp_rate_path.append(self.output_path + exp_name+"_rate.csv")
        for exp_name in self.Cs_exp_raw_path:
            self.Cs_exp_rate_path.append(self.output_path + exp_name+"_rate.csv")
        for exp_name in self.Co_exp_raw_path:
            self.Co_exp_rejection_path.append(self.output_path + exp_name+"_rejection.csv")
        for exp_name in self.Cs_exp_raw_path:
            self.Cs_exp_rejection_path.append(self.output_path + exp_name+"_rejection.csv")
        for exp_name in self.backgrounds_exp_raw_path:
            self.Bkg_exp_sorted_path.append(self.output_path + exp_name+"_sorted.csv")

    def read_sims(self):
        #output form, [rate facotr to mHz, enenrgy edges, counts above the bin edge, counts* counts above the bin edge]
        with open(self.Co_sim_path, "rb") as f:
            self.Co_sims = pickle.load(f)
        print("self.Co_sims",self.Co_sims)
        with open(self.Cs_sim_path, "rb") as f:
            self.Cs_sims = pickle.load(f)
        print("self.Cs_sims",self.Cs_sims)
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

    def clean_signal_analysis(self):
        self.df_bkg_116 = pd.read_csv(self.Bkg_average_116_path)
        self.df_bkg_116.columns = ['Pressure [bara]','Bkg Lifetime [s]','Bkg Lifetime Error [s]','Bkg Rate [mHz]', 'Bkg Rate Sigma [mHz]']
        self.df_bkg_119 =  pd.read_csv(self.Bkg_average_119_path)
        self.df_bkg_119.columns = ['Pressure [bara]', 'Bkg Lifetime [s]', 'Bkg Lifetime Error [s]','Bkg Rate [mHz]', 'Bkg Rate Sigma [mHz]']
        # first these are 116
        self.Cs_116_data = []
        self.Cs_119_data = []
        self.Co_116_data = []
        self.Co_119_data = []
        for i in range(0,self.Cs_exp_116_raw_len):
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

        for i in range(self.Cs_exp_116_raw_len+1,len(self.Cs_exp_raw_path) ):
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

        for i in range(self.Co_exp_116_raw_len + 1, len(self.Co_exp_raw_path)):
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


    def read_Setiz_info(self):
        Seitz_pressure_list = np.arange(1.25, 6.5, 0.25)

        Setiz_116 = [0.8318354532105874, 0.8940418766838347, 0.9631954092292087, 1.0403343615139717, 1.1266932576571842,
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
                                    "Setiz [keV]": Setiz_116,
                                    "Eion [keV]": E_ion_116,
                                    "Eion_rl-1_rhol-1 [GeVcm**2 g-1]": compound_x_116}
        self.df_energy_116_tab = pd.DataFrame(self.dict_energy_116_tab)

        Setiz_119 = [0.4336397431016649, 0.45977540397438443, 0.4882700006652598, 0.5194078524821146,
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
                                    "Setiz [keV]": Setiz_119,
                                    "Eion [keV]": E_ion_119,
                                    "Eion_rl-1_rhol-1 [GeVcm**2 g-1]": compound_x_119}
        self.df_energy_119_tab = pd.DataFrame(self.dict_energy_119_tab)

    def gamma_rejection_calculation(self):
        
        for i in range(len(self.Cs_exp_rate_path)):
            exp_df = pd.read_csv(self.Cs_exp_rate_path[i])
            columns_added  =exp_df.apply(self.calculate_rejection_by_row,axis=1,args=("Cs",))
            merged_df = pd.concat([exp_df, columns_added], axis=1)
            merged_df.to_csv(self.Cs_exp_rejection_path[i], index= False)

        for i in range(len(self.Co_exp_rate_path)):
            exp_df = pd.read_csv(self.Co_exp_rate_path[i])
            columns_added  =exp_df.apply(self.calculate_rejection_by_row,axis=1,args=("Co",))
            merged_df = pd.concat([exp_df, columns_added], axis=1)
            merged_df.to_csv(self.Co_exp_rejection_path[i], index= False)
            



    def read_exposure(self,filename):
        # Define your path (we'll use a relative path)
        file_path = os.path.join('..', 'exp_exposure', filename)

        # Read the file
        # sep='\s+' handles any number of spaces or tabs as delimiters
        df = pd.read_csv(file_path, sep='\s+', skiprows=1, header=None)

        return df
    def calculate_rejection_by_row(self, row, source):
        if source == "Co":
            self.sim_list = self.Co_sims
        elif source == "Cs":
            self.sim_list = self.Cs_sims
        else:
            print("NA sources")

        self.Rate_factor = self.sim_list[0]
        self.energy_edges = self.sim_list[1]
        self.counts_cum_bin = self.sim_list[2]
        self.counts_energy_cum_bin = self.sim_list[3]
        rejection_PS = 0
        rejection_PS_sigma = 0
        rejection_PK = 0
        rejection_PK_sigma = 0
        for i in range(len(self.energy_edges)):
            if row['Setiz [keV]']>= self.energy_edges[i]:
                # rejection per scattering, PS meaning perscattering
                counts = self.counts_cum_bin[i] + (row['Setiz [keV]'] - self.energy_edges[i]) * (
                        self.counts_cum_bin[i + 1] -
                        self.counts_cum_bin[i]) / (
                                 self.energy_edges[i + 1] - self.energy_edges[i])
                rate_PS = self.Rate_factor * (counts)

                rate_PS_sigma = rate_PS / np.sqrt(counts)

                rejection_PS = row['Clean Rate [mHz]'] / rate_PS

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
                break
        output = pd.Series({"Rejection Rate Scattering[]": rejection_PS,
                "Rejection Sigma Scattering[]": rejection_PS_sigma,
                "Rejection Rate KeV[/keV]": rejection_PK,
                "Rejection Sigma KeV[/keV]": rejection_PK_sigma})
        return output


    def calculate_rss(self, series):
        """Calculates sqrt(a^2 + b^2 + ...)"""
        return np.sqrt(np.sum(series ** 2))





if __name__=="__main__":
    IA =  integrated_analysis()
    # test = test_csv()