import pandas as pd
import matplotlib.pyplot as plt
import csv
import numpy as np
import os
import pickle
from scipy.optimize import curve_fit
import math
from scipy.stats import norm
from scipy.interpolate import interp1d
class integrated_analysis():
    def __init__(self, volume = ""):

        self.output_path = '/data/runzezhang/result/gamma_rejection/'
        self.plot_path = '/data/runzezhang/result/gamma_rejection/plot/'
        # self.Co_sim_path  ='/lzdata/runzezhang/result/GR_sims/Co_output_5E6.pkl'
        # self.Cs_sim_path = '/lzdata/runzezhang/result/GR_sims/Cs_output_5E6.pkl'
        self.base_path = "/lzdata/runzezhang/result/GR_sims/"

        self.volume_option =  volume
        if self.volume_option== ""or self.volume_option== "all":
            self.gamma_source_group = {"Cs":{"sim":{"pure_address":None,"pure_data":None, "doped_address":None,"doped_data":None},
                                             "exp":{"116K":{"raw_path":["Cold-Cs-11_17-18_exposures_mix","Cold-Cs-12_01_exposures_mix","Cold-Cs-12_10-11_exposures_mix","Cold-Cs-1_20-21_exposures_mix"],"sorted_path":[],"rate_path":[],"rejection_path":[],"plot_list":[],"plot":pd.DataFrame()},
                                                    "119K":{"raw_path":["Cold-Cs-2_2-3_exposures_zoom"],"sorted_path":[],"rate_path":[],"rejection_path":[],"plot_list":[],"plot":pd.DataFrame()}}},
                                       "Co": {"sim": {"pure_address": None, "pure_data": None, "doped_address": None,
                                                      "doped_data": None},
                                              "exp": {"116K": {"raw_path": ["60Co-Source-12_15-16_exposures_mix"],"sorted_path":[], "rate_path": [], "rejection_path": [],"plot_list":[],"plot":pd.DataFrame()},
                                                      "119K": {"raw_path": ["60Co-Source-02_06_exposures_mix"], "sorted_path":[],"rate_path": [], "rejection_path": [],"plot_list":[],"plot":pd.DataFrame()}}},
                                       "Ba": {"sim": {"pure_address": None, "pure_data": None, "doped_address": None,
                                                      "doped_data": None},
                                              "exp": {"116K": {"raw_path": ["Ba-11_19-24_exposures_mix"], "sorted_path":[],"rate_path": [], "rejection_path": [],"plot_list":[],"plot":pd.DataFrame()},
                                                      "119K": {"raw_path": [], "rate_path": [], "sorted_path":[],"rejection_path": [],"plot_list":[],"plot":pd.DataFrame()}}},
                                       "Th": {"sim": {"pure_address": None, "pure_data": None, "doped_address": None,
                                                      "doped_data": None},
                                              "exp": {"116K": {"raw_path": ["228Th-Source-11_20-21_exposures_mix"],"sorted_path":[], "rate_path": [], "rejection_path": [],"plot_list":[],"plot":pd.DataFrame()},
                                                      "119K": {"raw_path": [], "rate_path": [], "sorted_path":[],"rejection_path": [],"plot_list":[],"plot":pd.DataFrame()}}},
                                       "Hot_Cs": {"sim": {"pure_address": None, "pure_data": None, "doped_address": None,
                                                      "doped_data": None},
                                              "exp": {"116K": {"raw_path": ["Hot-Cs-11_11-12_exposures_mix"], "sorted_path":[],"rate_path": [], "rejection_path": [],"plot_list":[],"plot":pd.DataFrame()},
                                                      "119K": {"raw_path": [], "rate_path": [], "sorted_path":[],"rejection_path": [],"plot_list":[],"plot":pd.DataFrame()}}},
                                        }
        elif self.volume_option =="dome" or self.volume_option =="bulk":
            self.gamma_source_group = {
                "Cs": {"sim": {"pure_address": None, "pure_data": None, "doped_address": None, "doped_data": None},
                       "exp": {"116K": {"raw_path": ["Cold-Cs-11_17-18_fiducial_mix", "Cold-Cs-12_01_fiducial_mix",
                                                     "Cold-Cs-12_10-11_fiducial_mix", "Cold-Cs-1_20-21_fiducial_mix"],
                                        "sorted_path": [], "rate_path": [], "rejection_path": [], "plot_list": [],
                                        "plot": pd.DataFrame()},
                               "119K": {"raw_path": ["Cold-Cs-2_2-3_fiducial_mix"], "sorted_path": [],
                                        "rate_path": [], "rejection_path": [], "plot_list": [],
                                        "plot": pd.DataFrame()}}},
                "Co": {"sim": {"pure_address": None, "pure_data": None, "doped_address": None,
                               "doped_data": None},
                       "exp": {"116K": {"raw_path": ["60Co-Source-12_15-16_fiducial_mix"], "sorted_path": [],
                                        "rate_path": [], "rejection_path": [], "plot_list": [], "plot": pd.DataFrame()},
                               "119K": {"raw_path": ["60Co-Source-02_06_fiducial_mix"], "sorted_path": [],
                                        "rate_path": [], "rejection_path": [], "plot_list": [],
                                        "plot": pd.DataFrame()}}},
                "Ba": {"sim": {"pure_address": None, "pure_data": None, "doped_address": None,
                               "doped_data": None},
                       "exp": {"116K": {"raw_path": ["Ba-11_19-24_fiducial_mix"], "sorted_path": [], "rate_path": [],
                                        "rejection_path": [], "plot_list": [], "plot": pd.DataFrame()},
                               "119K": {"raw_path": [], "rate_path": [], "sorted_path": [], "rejection_path": [],
                                        "plot_list": [], "plot": pd.DataFrame()}}},
                "Th": {"sim": {"pure_address": None, "pure_data": None, "doped_address": None,
                               "doped_data": None},
                       "exp": {"116K": {"raw_path": ["228Th-Source-11_20-21_fiducial_mix"], "sorted_path": [],
                                        "rate_path": [], "rejection_path": [], "plot_list": [], "plot": pd.DataFrame()},
                               "119K": {"raw_path": [], "rate_path": [], "sorted_path": [], "rejection_path": [],
                                        "plot_list": [], "plot": pd.DataFrame()}}},
                "Hot_Cs": {"sim": {"pure_address": None, "pure_data": None, "doped_address": None,
                                   "doped_data": None},
                           "exp": {"116K": {"raw_path": ["Hot-Cs-11_11-12_fiducial_mix"], "sorted_path": [],
                                            "rate_path": [], "rejection_path": [], "plot_list": [],
                                            "plot": pd.DataFrame()},
                                   "119K": {"raw_path": [], "rate_path": [], "sorted_path": [], "rejection_path": [],
                                            "plot_list": [], "plot": pd.DataFrame()}}},
                }
        else:
            print("not availabble volume selection")


        #

        self.xe_shell_threshold = 0
        if self.volume_option == "" or self.volume_option == "all":
            self.background_group = {"116K": {"raw_path": ['Background-11_7-8_exposures',"Background-11_15-17_exposures",
                                        "Background-11_26-30_exposures", "Background-12_5-8_exposures",
                                        'Background-1_12-13_exposures','Background-1_17-20_exposures'],"sorted_path":[],"average_path":"", "full_info_path":""},
                                   "119K": { "raw_path": ["Background-1_30-2_2_exposures","Background-2_17-20_exposures",
                                               "Background-2_6-12_exposures","Background-2_28-3_2_exposures"],"sorted_path":[],"average_path":"", "full_info_path":""},
                                     "total": {"raw_path": [],
                                               "sorted_path": [],"average_path":"", "full_info_path":""}}
        else:
            self.background_group = {"116K": {"raw_path": ['Background-11_7-8_fiducial_mix',"Background-11_15-17_fiducial_mix",
                                        "Background-11_26-30_fiducial_mix", "Background-12_5-8_fiducial_mix",
                                        'Background-1_12-13_fiducial_mix','Background-1_17-20_fiducial_mix'],"sorted_path":[],"average_path":"", "full_info_path":""},
                                   "119K": { "raw_path": ["Background-1_30-2_2_fiducial_mix","Background-2_17-20-QUIET_fiducial_mix",
                                               "Background-2_6-12_fiducial_mix"],"sorted_path":[],"average_path":"", "full_info_path":""},
                                     "total": {"raw_path": [],
                                               "sorted_path": [],"average_path":"", "full_info_path":""}}
        # if self.volume_option=="":
        #     self.background_116_sorted_path = ['Background-11_7-8_exposures',"Background-11_15-17_exposures",
        #                                 "Background-11_26-30_exposures", "Background-12_5-8_exposures",
        #                                 'Background-1_12-13_exposures','Background-1_17-20_exposures']
        #     self.background_exp_116_raw_len = len(self.background_116_sorted_path)
        #     self.background_119_sorted_path = ["Background-1_30-2_2_exposures","Background-2_17-20_exposures",
        #                                        "Background-2_6-12_exposures","Background-2_28-3_2_exposures"]
        #
        #     self.backgrounds_exp_raw_path = self.background_116_sorted_path+self.background_119_sorted_path
        # else:
        #     self.background_116_sorted_path = ['Background-11_7-8_exposures', "Background-11_15-17_exposures",
        #                                        "Background-11_26-30_exposures", "Background-12_5-8_exposures",
        #                                        'Background-1_12-13_exposures', 'Background-1_17-20_exposures']
        #     self.background_exp_116_raw_len = len(self.background_116_sorted_path)
        #     self.background_119_sorted_path = ["Background-1_30-2_2_exposures", "Background-2_17-20_exposures",
        #                                        "Background-2_6-12_exposures", "Background-2_28-3_2_exposures"]
        #
        #     self.backgrounds_exp_raw_path = self.background_116_sorted_path + self.background_119_sorted_path


        self.main_v2()
    def main_v2(self):
        self.predata_process()

        self.read_Seitz_info()
        self.read_Seitz_info_C3F8()
        self.pre_background_analysis()
        # self.read_raw_backgrounds_exp()
        # self.average_background_analysis()
        #
        #
        self.bkg_subtracted_analysis()

        self.gamma_rejection_plot_v3()







    def predata_process(self):

        # # bkg data
        # self.Bkg_exp_sorted_path = []
        # self.Bkg_average_116_path = self.output_path + "background_116_average" + ".csv"
        # self.Bkg_average_119_path = self.output_path + "background_119_average" + ".csv"
        # # bkg table containing seitz infos
        # # separate this from above because of clean signal need to merge only on pressure column
        # self.Bkg_average_116_full_info_path = self.output_path + "background_116_average_full_info" + ".csv"
        # self.Bkg_average_119_full_info_path = self.output_path + "background_119_average_full_info" + ".csv"
        # for exp_name in self.backgrounds_exp_raw_path:
        #     self.Bkg_exp_sorted_path.append(self.output_path + exp_name + "_sorted.csv")



        # source data
        for source, source_config in self.gamma_source_group.items():
            # print

            source_config["sim"]["pure_address"] = self.base_path+source+self.volume_option+"_output_5E6_ERv2.pkl"
            source_config["sim"]["doped_address"] = self.base_path + source +self.volume_option+ "_doped_output_full_track.pkl"

            for temperature, temp_config in source_config["exp"].items():
                if temp_config["raw_path"] !=[]:
                    for raw_path in temp_config["raw_path"]:
                        temp_config["sorted_path"].append(self.output_path+raw_path+self.volume_option+"_sortedv2.csv")
                        temp_config["rate_path"].append(self.output_path + raw_path +self.volume_option+ "_ratev2.csv")
                        temp_config["rejection_path"].append(self.output_path + raw_path +self.volume_option+ "_rejectionv2.csv")
            # read sims
            with open(source_config["sim"]["pure_address"], "rb") as f:
                source_config["sim"]["pure_data"] = pickle.load(f)
            # print(source, "sims", source_config["sim"]["pure_data"])

            with open(source_config["sim"]["doped_address"], "rb") as f:
                source_config["sim"]["doped_data"] = pickle.load(f)
            # print(source, "sims", source_config["sim"]["doped_data"])

            # read_source exp data
            for temperature, temp_config in source_config["exp"].items():
                if temp_config["raw_path"] !=[]:
                    for raw_path_index in range(len(temp_config["raw_path"])):
                        exposure_df = self.read_exposure(temp_config["raw_path"][raw_path_index] + ".txt", volume=self.volume_option)

                        # exposure_df = exposure_df.iloc[:, :7]
                        # exposure_df.columns = ['Pressure [bara]', 'Lifetime [s]', 'Lifetime Error [s]',
                        #                        'Exponential Fit 2xNLL',
                        #                        'N.d.o.f.', 'Time Cut High [s]', 'Time Cut Low [s]']
                        exposure_df = exposure_df[
                            (exposure_df['Lifetime [s]'] <= 6.92e-1) | (exposure_df['Lifetime [s]'] >= 6.94e-1)]
                        exposure_df = exposure_df[
                            (exposure_df['Lifetime Error [s]'] / exposure_df['Lifetime [s]'] <= 0.3)]
                        # add rate column
                        exposure_df['Exp Rate [mHz]'] = 1000 / exposure_df['Lifetime [s]']
                        exposure_df['Exp Rate Sigma [mHz]'] = exposure_df['Lifetime Error [s]'] * 1000 / (
                            exposure_df['Lifetime [s]']) ** 2
                        # exposure_df = exposure_df.drop(columns=['Exponential Fit 2xNLL',
                        #                                         'N.d.o.f.', 'Time Cut High [s]', 'Time Cut Low [s]'])
                        exposure_df.to_csv(temp_config["sorted_path"][raw_path_index], index=False)


    def pre_background_analysis(self):

        self.background_group["total"]["raw_path"] = self.background_group["116K"]["raw_path"]+self.background_group["119K"]["raw_path"]
        for temp, temp_config in self.background_group.items():
            if temp != "total":
                temp_config["average_path"] = self.output_path + f"background_{temp}_average.csv"
                temp_config["full_info_path"]=self.output_path + f"background_{temp}_average_full_info" + ".csv"
                for exp in temp_config["raw_path"]:
                    temp_config["sorted_path"].append(self.output_path + exp + "_sorted.csv")
            else:
                for exp in temp_config["raw_path"]:
                    temp_config["sorted_path"].append(self.output_path + exp + "_sorted.csv")

        # calcuate the rate
        for i in range(len(self.background_group["total"]["raw_path"])):
            exposure_df = self.read_exposure(self.background_group["total"]["raw_path"][i] + ".txt", volume=self.volume_option)
            # exposure_df.columns = ['Pressure [bara]', 'Lifetime [s]', 'Lifetime Error [s]', 'Exponential Fit 2xNLL',
            #                        'N.d.o.f.', 'Time Cut High [s]', 'Time Cut Low [s]']
            exposure_df = exposure_df[
                (exposure_df['Lifetime [s]'] <= 6.92e-1) | (exposure_df['Lifetime [s]'] >= 6.94e-1)]
            exposure_df = exposure_df[
                (exposure_df['Lifetime Error [s]'] / exposure_df['Lifetime [s]'] <= 0.3)]
            exposure_df['Bkg Rate [mHz]'] = 1000 / exposure_df['Lifetime [s]']
            exposure_df['Bkg Rate Sigma [mHz]'] = exposure_df['Lifetime Error [s]'] * 1000 / (
                exposure_df['Lifetime [s]']) ** 2
            # exposure_df = exposure_df.drop(columns=['Exponential Fit 2xNLL',
            #                        'N.d.o.f.', 'Time Cut High [s]', 'Time Cut Low [s]'])
            exposure_df.to_csv(self.background_group["total"]["sorted_path"][i], index=False)


        # average

        for temp, temp_config in self.background_group.items():
            if temp != "total":
                temp_config["average_path"] = self.output_path + f"background_{temp}_average.csv"
                temp_config["full_info_path"]=self.output_path + f"background_{temp}_average_full_info" + ".csv"
                bkg_df_list = []
                for i in range(len(temp_config["sorted_path"])):

                    bkg_df = pd.read_csv(temp_config["sorted_path"][i])
                    bkg_df_list.append(bkg_df)
                combined_df = pd.concat(bkg_df_list, ignore_index=True)
                result_df = combined_df.groupby('Pressure [bara]').agg({
                    'Lifetime [s]': 'mean',  # Simple average
                    'Lifetime Error [s]': self.calculate_rss,  # Custom square root math,
                    'Bkg Rate [mHz]': 'mean',
                    'Bkg Rate Sigma [mHz]': self.calculate_rss
                }).reset_index()
                print('result_df_116', result_df)
                result_df.to_csv(temp_config["average_path"], index=False)

                # add different source uplimit

                result_df_full_info = pd.merge(result_df, self.df_energy_116_tab, on='Pressure [bara]',
                                                   how="inner")
                columns_added_Cs = result_df_full_info.apply(self.calculate_bkg_uplimit_by_row, axis=1,
                                                                 args=("Cs",))
                columns_added_Co = result_df_full_info.apply(self.calculate_bkg_uplimit_by_row, axis=1,
                                                                 args=("Co",))

                result_df_full_info = pd.concat([result_df_full_info, columns_added_Cs, columns_added_Co],
                                                    axis=1)
                # print('result_df_116_full_info.columns',result_df_116_full_info.columns)
                result_df_full_info.to_csv(temp_config["full_info_path"], index=False)


    def read_raw_backgrounds_exp(self):
        # read file, delete unreasonable rows and rewrite
        for i in range(len(self.backgrounds_exp_raw_path)):
            exposure_df = self.read_exposure(self.backgrounds_exp_raw_path[i] + ".txt", volume=self.volume_option)
            # exposure_df.columns = ['Pressure [bara]', 'Lifetime [s]', 'Lifetime Error [s]', 'Exponential Fit 2xNLL',
            #                        'N.d.o.f.', 'Time Cut High [s]', 'Time Cut Low [s]']
            exposure_df = exposure_df[
                (exposure_df['Lifetime [s]'] <= 6.92e-1) | (exposure_df['Lifetime [s]'] >= 6.94e-1)]
            exposure_df = exposure_df[
                (exposure_df['Lifetime Error [s]'] / exposure_df['Lifetime [s]'] <= 0.3)]
            exposure_df['Bkg Rate [mHz]'] = 1000 / exposure_df['Lifetime [s]']
            exposure_df['Bkg Rate Sigma [mHz]'] = exposure_df['Lifetime Error [s]'] * 1000 / (
            exposure_df['Lifetime [s]']) ** 2
            # exposure_df = exposure_df.drop(columns=['Exponential Fit 2xNLL',
            #                        'N.d.o.f.', 'Time Cut High [s]', 'Time Cut Low [s]'])
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


    def bkg_subtracted_analysis(self):
        self.df_bkg_116 = pd.read_csv(self.background_group["116K"]["average_path"])
        # self.df_bkg_116.columns = ['Pressure [bara]','Bkg Lifetime [s]','Bkg Lifetime Error [s]','Bkg Rate [mHz]', 'Bkg Rate Sigma [mHz]']
        self.df_bkg_116 = self.df_bkg_116.rename(columns={
            'Lifetime [s]': 'Bkg Lifetime [s]',
            'Lifetime Error [s]': 'Bkg Lifetime Error [s]'
        })
        self.df_bkg_119 =  pd.read_csv(self.background_group["119K"]["average_path"])
        # self.df_bkg_119.columns = ['Pressure [bara]', 'Bkg Lifetime [s]', 'Bkg Lifetime Error [s]','Bkg Rate [mHz]', 'Bkg Rate Sigma [mHz]']
        self.df_bkg_119 = self.df_bkg_116.rename(columns={
            'Lifetime [s]': 'Bkg Lifetime [s]',
            'Lifetime Error [s]': 'Bkg Lifetime Error [s]'
        })

        for source, source_config in self.gamma_source_group.items():


            # read_source exp data
            for temperature, temp_config in source_config["exp"].items():
                if temp_config["sorted_path"] !=[]:
                    for sorted_path_index in range(len(temp_config["sorted_path"])):
                        # combine the Seitz to the exp data
                        # print("Cs ",temperature, temp_config["sorted_path"])
                        exposure_df = pd.read_csv(temp_config["sorted_path"][sorted_path_index])

                        # merge both has the pressure value, on pressure
                        if temperature =="116K":
                            merged_df = pd.merge(self.df_bkg_116, exposure_df, on='Pressure [bara]', how="inner")
                        elif temperature=="119K":
                            merged_df = pd.merge(self.df_bkg_119, exposure_df, on='Pressure [bara]', how="inner")
                        else:
                            print("wrong temp")
                        # clean rate!
                        merged_df['Clean Rate [mHz]'] = merged_df['Exp Rate [mHz]'] - merged_df['Bkg Rate [mHz]']
                        merged_df['Clean Rate Sigma [mHz]'] = np.sqrt(
                            merged_df['Exp Rate Sigma [mHz]'] ** 2 + merged_df['Bkg Rate Sigma [mHz]'] ** 2)
                        # add Seitz and Eion unit
                        if temperature =="116K":
                            merged_df = pd.merge(merged_df, self.df_energy_116_tab, on='Pressure [bara]', how="inner")
                        elif temperature=="119K":
                            merged_df = pd.merge(merged_df, self.df_energy_119_tab, on='Pressure [bara]', how="inner")
                        else:
                            print("wrong temp")

                        # add sims analysis to get rejection


                        merged_df.to_csv(temp_config["rate_path"][sorted_path_index], index=False)

                        # cacluate rejection

                        exp_df = pd.read_csv(temp_config["rate_path"][sorted_path_index])

                        columns_added = exp_df.apply(self.calculate_rejection_by_row_v2, axis=1, args=(source,))
                        if source=="Cs" and sorted_path_index==3:
                            print("Cs fault", columns_added)



                        merged_df_rejection = pd.concat([exp_df, columns_added], axis=1)
                        # print('Cs print(merged_df)',self.Cs_exp_rate_path[i],'\n',merged_df)
                        merged_df_rejection.to_csv(temp_config["rejection_path"][sorted_path_index], index=False)






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


        Q_compound_x_116 = []
        for i in range(len(Seitz_116)):
            x = Seitz_116[i] * 10 / (rho_l_116[i] * Rl_116[i])  # fit unit
            Q_compound_x_116.append(x)

        self.dict_energy_116_tab = {"Pressure [bara]":Seitz_pressure_list,
                                    "Seitz [keV]": Seitz_116,
                                    "Eion [keV]": E_ion_116,
                                    "Eion_rl-1_rhol-1 [GeVcm**2 g-1]": compound_x_116, "Q_rl-1_rhol-1 [GeVcm**2 g-1]": Q_compound_x_116}
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

        Q_compound_x_119 = []
        for i in range(len(Seitz_116)):
            x = Seitz_119[i] * 10 / (rho_l_119[i] * Rl_119[i])  # fit unit
            Q_compound_x_119.append(x)

        self.dict_energy_119_tab = {"Pressure [bara]": Seitz_pressure_list,
                                    "Seitz [keV]": Seitz_119,
                                    "Eion [keV]": E_ion_119,
                                    "Eion_rl-1_rhol-1 [GeVcm**2 g-1]": compound_x_119,
                                    "Q_rl-1_rhol-1 [GeVcm**2 g-1]": Q_compound_x_119}
        self.df_energy_119_tab = pd.DataFrame(self.dict_energy_119_tab)

    def read_Seitz_info_v2(self):
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


        Q_compound_x_116 = []
        for i in range(len(Seitz_116)):
            x = Seitz_116[i] * 10 / (rho_l_116[i] * Rl_116[i])  # fit unit
            Q_compound_x_116.append(x)

        self.dict_energy_116_tab = {"Pressure [bara]":Seitz_pressure_list,
                                    "Seitz [keV]": Seitz_116,
                                    "Eion [keV]": E_ion_116,
                                    "Eion_rl-1_rhol-1 [GeVcm**2 g-1]": compound_x_116, "Q_rl-1_rhol-1 [GeVcm**2 g-1]": Q_compound_x_116}
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

        Q_compound_x_119 = []
        for i in range(len(Seitz_116)):
            x = Seitz_119[i] * 10 / (rho_l_119[i] * Rl_119[i])  # fit unit
            Q_compound_x_119.append(x)

        self.dict_energy_119_tab = {"Pressure [bara]": Seitz_pressure_list,
                                    "Seitz [keV]": Seitz_119,
                                    "Eion [keV]": E_ion_119,
                                    "Eion_rl-1_rhol-1 [GeVcm**2 g-1]": compound_x_119,
                                    "Q_rl-1_rhol-1 [GeVcm**2 g-1]": Q_compound_x_119}
        self.df_energy_119_tab = pd.DataFrame(self.dict_energy_119_tab)


    def read_Seitz_info_C3F8(self):
        Seitz_pressure_list = np.arange(0, 51, 5)

        Seitz_10 = [1.85,2.16,2.55,3.05,3.71,4.6,5.82,7.56,10.15,14.2,20.94]  # in keV
        # keV
        E_ion_10 = [1.03,1.18,1.36,1.58,1.86,2.23,2.72,3.38,4.31,5.69,7.85]


        compound_x_10 = [1.14,1.22,1.31,1.42,1.55,1.7,1.89,2.12,2.42,2.81,3.34]



        Q_compound_x_10 = [2.04, 2.23, 2.46,2.74,3.08,3.51,4.04,4.75,5.69,7.0, 8.92]


        self.dict_energy_10_tab = {"Pressure [bara]":Seitz_pressure_list,
                                    "Seitz [keV]": Seitz_10,
                                    "Eion [keV]": E_ion_10,
                                    "Eion_rl-1_rhol-1 [GeVcm**2 g-1]": compound_x_10, "Q_rl-1_rhol-1 [GeVcm**2 g-1]": Q_compound_x_10}
        self.df_energy_10_tab = pd.DataFrame(self.dict_energy_10_tab)

        Seitz_24 = [0.4, 0.45, 0.49, 0.55, 0.62, 0.69, 0.79, 0.9, 1.04, 1.22, 1.44]  # in keV
        # keV
        E_ion_24 = [0.25,0.27,0.3,0.33,0.36,0.4,0.45,0.5,0.57,0.65,0.75]


        compound_x_24 = [0.5,0.52,0.55,0.57,0.6,0.64,0.67,0.71,0.76,0.82,0.88]


        Q_compound_x_24 = [0.81, 0.86, 0.91, 0.97, 1.03, 1.1, 1.19, 1.28, 1.4, 1.53, 1.69]


        self.dict_energy_24_tab = {"Pressure [bara]": Seitz_pressure_list,
                                    "Seitz [keV]": Seitz_24,
                                    "Eion [keV]": E_ion_24,
                                    "Eion_rl-1_rhol-1 [GeVcm**2 g-1]": compound_x_24,
                                    "Q_rl-1_rhol-1 [GeVcm**2 g-1]": Q_compound_x_24}
        self.df_energy_24_tab = pd.DataFrame(self.dict_energy_24_tab)




    def read_exposure(self,filename, volume=""):
        # Define your path (we'll use a relative path)
        file_path = os.path.join('..', 'exp_exposure', filename)

        # Read the file
        # sep='\s+' handles any number of spaces or tabs as delimiters
        df = pd.read_csv(file_path, sep='\s+', skiprows=1, header=None)
        num_cols = len(df.columns)
        # print(file_path, df.head(1))
        if  num_cols ==24:
            df.columns = ['Pressure [bara]',	'Lifetime [s]',	'Lifetime Error [s]',
                          'N quality events',	'N Bulk',	'Bulk Fraction'	,'Bulk Fraction Error',
                          'Bulk Lifetime [s]',	'Bulk Lifetime Error [s]',	'N Wall',	'Wall Fraction',
                          'Wall Fraction Error',	'Wall Lifetime [s]',	'Wall Lifetime Error [s]',
                          'N Dome',	'Dome Fraction',	'Dome Fraction Error',	'Dome Lifetime [s]',
                          'Dome Lifetime Error [s]',	'N Bottom',	'Bottom Fraction',	'Bottom Fraction Error',
                          'Bottom Lifetime [s]'	,'Bottom Lifetime Error [s]']
            if volume == "":
                df = df[['Pressure [bara]', 'Lifetime [s]', 'Lifetime Error [s]']]
            elif volume =="bulk":
                df = df[['Pressure [bara]', 'Bulk Lifetime [s]', 'Bulk Lifetime Error [s]']]
                df.columns = ['Pressure [bara]', 'Lifetime [s]', 'Lifetime Error [s]']
            elif volume =="dome":
                df = df[['Pressure [bara]', 'Dome Lifetime [s]', 'Dome Lifetime Error [s]']]
                df.columns = ['Pressure [bara]', 'Lifetime [s]', 'Lifetime Error [s]']
            else:
                df = df[['Pressure [bara]', 'Lifetime [s]', 'Lifetime Error [s]']]

        elif num_cols ==7:
            df.columns = ['Pressure [bara]', 'Lifetime [s]', 'Lifetime Error [s]',
                               'Exponential Fit 2xNLL',
                               'N.d.o.f.', 'Time Cut High [s]', 'Time Cut Low [s]']
            df = df[['Pressure [bara]', 'Lifetime [s]', 'Lifetime Error [s]']]
        else:
            df = df.iloc[:, :7]
            df.columns = ['Pressure [bara]', 'Lifetime [s]', 'Lifetime Error [s]',
                          'Exponential Fit 2xNLL',
                          'N.d.o.f.', 'Time Cut High [s]', 'Time Cut Low [s]']
            df = df[['Pressure [bara]', 'Lifetime [s]', 'Lifetime Error [s]']]
            print("Column doesn't match", num_cols)


        return df
    def calculate_bkg_uplimit_by_row(self, row, source):
        try:
            self.sim_list = self.gamma_source_group[source]["sim"]["pure_data"]
            self.sim_doped_list = self.gamma_source_group[source]["sim"]["doped_data"]
        except:
            print("NA sources")

        # if source == "Co":
        #     self.sim_list = self.Co_sims
        # elif source == "Cs":
        #     self.sim_list = self.Cs_sims
        # elif source == "Cs_doped":
        #     self.sim_list = self.Cs_sims_doped
        # elif source == "Co_doped":
        #     self.sim_list = self.Co_sims_doped
        # else:
        #     print("NA sources")

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

    #source_config["sim"]["pure_address"]
    def calculate_rejection_by_row_v2(self, row, source):
        try:

            self.sim_list = self.gamma_source_group[source]["sim"]["pure_data"]
            self.sim_doped_list = self.gamma_source_group[source]["sim"]["doped_data"]
        except:
            print("NA sources")
        print(row.empty)
        if not row.empty:
            self.Rate_factor = self.sim_list[0]

            # print("self.Rate_factor",self.Rate_factor)
            self.energy_edges = self.sim_list[1][0][1]
            # print("self.energy_edges", self.energy_edges)
            self.counts_cum_bin = self.sim_list[2]
            # print("self.counts_cum_bin", self.counts_cum_bin)
            self.counts_energy_cum_bin = self.sim_list[6]
            # print("self.counts_energy_cum_bin", self.counts_energy_cum_bin)
            rejection_PS = 0
            rejection_PS_sigma = 0
            rejection_PK = 0
            rejection_PK_sigma = 0
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

                    rejection_PS = row['Clean Rate [mHz]'] / rate_PS
                    # rejection_uplimit_PS = row['Bkg Rate Sigma [mHz]']/rate_PS

                    # will be returned
                    rejection_PS_sigma = np.sqrt(
                        (row['Clean Rate Sigma [mHz]'] / rate_PS) ** 2 + (
                                    row['Clean Rate [mHz]'] * rate_PS_sigma / rate_PS ** 2) ** 2)
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
                        (row['Clean Rate Sigma [mHz]'] / rate_PK) ** 2 + (
                                    row['Clean Rate [mHz]'] * rate_PK_sigma / rate_PK ** 2) ** 2)
                    # rejection_uplimit_PK = row['Bkg Rate Sigma [mHz]'] / rate_PK
                    break
            # for xenon k shell absorption
            # only counts interactions that > 34.56 keV
            self.Rate_factor_doped = self.sim_doped_list[0]
            # print("self.Rate_factor",self.Rate_factor)
            self.energy_edges_doped = self.sim_doped_list[1][0][1]
            # print("self.energy_edges", self.energy_edges)
            self.counts_cum_bin_doped = self.sim_doped_list[2]
            # print("self.counts_cum_bin", self.counts_cum_bin)
            self.counts_energy_cum_bin_doped = self.sim_doped_list[3]
            rejection_PX = 0
            rejection_PX_sigma = 0
            energy_K = 34.56

            for i in range(len(self.energy_edges_doped)):
                # if row['Seitz [keV]']>= self.energy_edges[i]:
                #     # rejection per xenon photo absorption in k shell, PX meaning per xenon
                #     counts = self.counts_cum_bin[i] + (row['Seitz [keV]'] - self.energy_edges[i]) * (
                #             self.counts_cum_bin[i + 1] -
                #             self.counts_cum_bin[i]) / (
                #                      self.energy_edges[i + 1] - self.energy_edges[i])
                if self.energy_edges_doped[i] >= self.xe_shell_threshold:
                    counts = self.counts_cum_bin_doped[i]

                    rate_PX = self.Rate_factor_doped * (counts)

                    rate_PX_sigma = rate_PX / np.sqrt(counts)

                    rejection_PX = row['Clean Rate [mHz]'] / rate_PX
                    # rejection_uplimit_PX = row['Bkg Rate Sigma [mHz]']/rate_PX

                    # will be returned
                    rejection_PX_sigma = np.sqrt(
                        (row['Clean Rate Sigma [mHz]'] / rate_PX) ** 2 + (
                                    row['Clean Rate [mHz]'] * rate_PX_sigma / rate_PX ** 2) ** 2)
                    # will be returned
                    break
            output = pd.Series({"Rejection Rate Scattering[]": rejection_PS,
                                "Rejection Sigma Scattering[]": rejection_PS_sigma,
                                "Rejection Rate KeV[/keV]": rejection_PK,
                                "Rejection Sigma KeV[/keV]": rejection_PK_sigma,
                                "Rejection Rate Xenon Abs[]": rejection_PX,
                                "Rejection Sigma Xenon Abs[]": rejection_PX_sigma,
                                })
        else:
            output = pd.Series({"Rejection Rate Scattering[]": 0,
                                "Rejection Sigma Scattering[]": 0,
                                "Rejection Rate KeV[/keV]": 0,
                                "Rejection Sigma KeV[/keV]": 0,
                                "Rejection Rate Xenon Abs[]": 0,
                                "Rejection Sigma Xenon Abs[]": 0,
                                })
        print("source", source, "Rate[mHz]:", self.Rate_factor * self.counts_cum_bin[0], "Energy deposit [keV]",
              self.Rate_factor * self.counts_energy_cum_bin[0],
              "Xe abs Rate[mHz]: ", self.Rate_factor_doped * self.counts_cum_bin_doped[0])
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


    def gamma_rejection_plot_v3(self):

        self.fitting_list = []
        self.Cs_fitting_list = []
        self.Co_fitting_list = []
        self.Ba_fitting_list = []

        for source, source_config in self.gamma_source_group.items():
            # read_source exp data
            for temperature, temp_config in source_config["exp"].items():
                if temp_config["rejection_path"] !=[]:
                    for rejection_path_index in range(len(temp_config["rejection_path"])):
                        # combine the Seitz to the exp data
                        df = pd.read_csv(temp_config["rejection_path"][rejection_path_index])
                        print("df.columns",source,rejection_path_index,df.columns)


                        pressure_drop_list = []
                        df = df[~df['Pressure [bara]'].isin(pressure_drop_list)]
                        # only positive rate
                        df = df[df['Clean Rate [mHz]'] > 0]
                        temp_config["plot_list"].append(df)

                        df_fit = df[['Seitz [keV]', "Rejection Rate Scattering[]", 'Eion_rl-1_rhol-1 [GeVcm**2 g-1]',
                                     "Rejection Rate KeV[/keV]", 'Q_rl-1_rhol-1 [GeVcm**2 g-1]',
                                     "Rejection Rate Xenon Abs[]",
                                     'Clean Rate [mHz]', "Eion [keV]"]]

                        self.fitting_list.append(df_fit)
                        if source =="Cs":
                            self.Cs_fitting_list.append(df_fit)
                    print(temp_config["plot_list"])
                    temp_config["plot"] = pd.concat(temp_config["plot_list"], ignore_index=True)
                    temp_config["plot"] = self.concat_PT_condition(temp_config["plot"])

        # print Q vs per keV and Eion per interaction




        self.Cs_116_label = ["Cs 11/17/2025 116K", "Cs 12/01/2025 116K", "Cs 12/10/2025 116K", "Cs 01/20/2026 116K"]
        self.Cs_119_label = ["Cs 02/02/2026 119K"]
        self.Co_116_label = ["Co 12/15/2026 116K"]
        self.Co_119_label = ["Co 02/06/2026 119K"]
        self.Ba_116_label = ["Ba 11/19/2025 116K"]


        fig, ax = plt.subplots(3, 4, figsize=(40, 24))
        # fig, ax = plt.subplots(2, 1, figsize=(6, 10))

        y_config = [{"y": "Rejection Rate Scattering[]", "y_err": "Rejection Sigma Scattering[]",
                     "ylabel": "Nucleation probability (per interaction) "},
                    {"y": "Rejection Rate KeV[/keV]", "y_err": "Rejection Sigma KeV[/keV]",
                     "ylabel": "Probability per energy deposited \n (events/keV) "},
                    {"y": "Rejection Rate Xenon Abs[]", "y_err": "Rejection Sigma Xenon Abs[]",
                     "ylabel": "Nucleation probability \n (per xenon photoabsorption in K shell) "},
                    {"y": "Clean Rate [mHz]", "y_err": 'Clean Rate Sigma [mHz]',
                     "ylabel": "Background Substacted Rate [mHz]"}]
        x_config = [{"x": "Seitz [keV]", "xlabel": r"Seitz threshold [keV]"},
                    {"x": 'Eion_rl-1_rhol-1 [GeVcm**2 g-1]',
                     "xlabel": r"$E_{ion} r_l^{-1} \rho_l^{-1}$ [GeV cm$^2$ g$^{-1}$]"},
                    {"x": "Q_rl-1_rhol-1 [GeVcm**2 g-1]", "xlabel": r"$Q_{Seitz} r_l^{-1} \rho_l^{-1}$ [GeV cm$^2$ g$^{-1}$]"}]

        for i in range(4):
            for j in range(3):
                # Extract the configuration for this specific slot
                y_cfg = y_config[i]
                x_cfg = x_config[j]
                ax_ij = ax[j, i]

                for source, source_config in self.gamma_source_group.items():
                    # read_source exp data
                    for temperature, temp_config in source_config["exp"].items():
                        if not temp_config["plot"].empty:
                            ax_ij.errorbar(temp_config["plot"][x_cfg["x"]], temp_config["plot"][y_cfg["y"]],
                                           yerr=temp_config["plot"][y_cfg["y_err"]], label=str(source)+" "+str(temperature), fmt='o',
                                           markersize=8)


                ax_ij.set_xlabel(x_cfg["xlabel"],fontsize=16)
                ax_ij.set_ylabel(y_cfg["ylabel"],fontsize=16)
                ax_ij.set_yscale("log")
                ax_ij.legend(loc='lower left', fontsize=13)


        self.fitting_df = pd.concat(self.fitting_list, ignore_index=True)
        fitting_matrix = self.fitting_gamma_rejection_v3(self.fitting_df,x_config,y_config)


        # plot the fitting function
        for i in range(4):
            for j in range(3):
                # Extract the configuration for this specific slot
                y_cfg = y_config[i]
                x_cfg = x_config[j]
                ax_ij = ax[j, i]
                a_val = fitting_matrix[i][j][0]
                b_val = fitting_matrix[i][j][1]

                # ax_ij.plot(fitting_matrix[i][j][2], fitting_matrix[i][j][3],
                #            color="black")
                label_text = f"A = {a_val:.2e},\nB = {b_val:.2e}"
                # ax_ij.plot(fitting_matrix[i][j][2], fitting_matrix[i][j][3],label = label_text,
                #       color="black")
                ax_ij.plot(fitting_matrix[i][j][2], fitting_matrix[i][j][3],
                           color="black", label=label_text)
                ax_ij.legend(loc='lower left', fontsize=13)

        # plt.show()
        plt.savefig(self.plot_path + f"gamma_rejection_{self.volume_option}v3.pdf")

        plt.clf()
        # self.Qseitz_compound_xe_plot()
        # self.Ratio_plot()
        # self.time_plot()
    def Qseitz_compound_xe_plot(self):
        fig, ax = plt.subplots(1, 2, figsize=(16, 6))
        [result_Q_scatter, result_Q_keV, result_Q_xe, result_Eion_scatter, result_Eion_keV, result_Eion_xe,
         result_Q2_xe, result_Q_rate] = self.fitting_gamma_rejection_v2(self.fitting_df)

        # # for Rate Cs and Co, fit individually, only get last component
        self.fitting_df_Cs = pd.concat(self.Cs_fitting_list, ignore_index=True)
        self.fitting_df_Co = pd.concat(self.Co_fitting_list, ignore_index=True)
        result_Cs_fitting = self.fitting_gamma_rejection_v2(self.fitting_df_Cs)[-1]
        result_Co_fitting = self.fitting_gamma_rejection_v2(self.fitting_df_Co)[-1]

        ax[0].errorbar(self.df_Cs_116_plot['Q_rl-1_rhol-1 [GeVcm**2 g-1]'],
                       self.df_Cs_116_plot["Rejection Rate Xenon Abs[]"],
                       yerr=self.df_Cs_116_plot["Rejection Sigma Xenon Abs[]"], label="Cs 116K", fmt='o',markersize=8)

        # ax[1].errorbar(self.df_Co_116_plot['Q_rl-1_rhol-1 [GeVcm**2 g-1]'], self.df_Co_116_plot["Rejection Rate KeV[/keV]"],
        #                   yerr=self.df_Co_116_plot["Rejection Sigma KeV[/keV]"], label="Co 116K", fmt='o')

        ax[0].errorbar(self.df_Cs_119_plot['Q_rl-1_rhol-1 [GeVcm**2 g-1]'],
                       self.df_Cs_119_plot["Rejection Rate Xenon Abs[]"],
                       yerr=self.df_Cs_119_plot["Rejection Sigma Xenon Abs[]"], label="Cs 119K", fmt='o',markersize=8)

        # ax[1].errorbar(self.df_Co_116_plot['Q_rl-1_rhol-1 [GeVcm**2 g-1]'], self.df_Co_116_plot["Rejection Rate KeV[/keV]"],
        #                yerr=self.df_Co_116_plot["Rejection Sigma KeV[/keV]"], label="Co 116K", fmt='o')

        ax[0].errorbar(self.df_Co_116_plot['Q_rl-1_rhol-1 [GeVcm**2 g-1]'],
                       self.df_Co_116_plot["Rejection Rate Xenon Abs[]"],
                       yerr=self.df_Co_116_plot["Rejection Sigma Xenon Abs[]"], label="Co 116K", fmt='o',markersize=8)

        ax[0].errorbar(self.df_Co_119_plot['Q_rl-1_rhol-1 [GeVcm**2 g-1]'],
                       self.df_Co_119_plot["Rejection Rate Xenon Abs[]"],
                       yerr=self.df_Co_119_plot["Rejection Sigma Xenon Abs[]"], label="Co 119K", fmt='o',markersize=8)

        # ax[0].errorbar(self.df_Ba_116_plot['Q_rl-1_rhol-1 [GeVcm**2 g-1]'],
        #                self.df_Ba_116_plot["Rejection Rate Xenon Abs[]"],
        #                yerr=self.df_Ba_116_plot["Rejection Sigma Xenon Abs[]"], label="Ba 116K", fmt='o')

        ax[0].plot(self.df_Ba_116_plot['Q_rl-1_rhol-1 [GeVcm**2 g-1]'], self.df_Ba_116_plot["Rejection Rate Xenon Abs[]"],
                   label="Ba 116K 95% CL \nUpper Limit", marker='v', linestyle='None')
        a_val = result_Q2_xe[0]
        b_val = result_Q2_xe[1]
        label_text = f"A = {a_val:.2e},\nB = {b_val:.2e}"
        ax[0].plot(result_Q2_xe[2], result_Q2_xe[3], label=label_text,
                   color="black")
        # ax[0].plot(result_Q2_xe[2], result_Q2_xe[3],
        #            color="black")

        ax[0].set_xlabel(r"$Q_{Seitz} r_l^{-1} \rho_l^{-1}$ [GeV cm$^2$ g$^{-1}$]",fontsize=16 )
        ax[0].set_ylabel("Nucleation probability \n (per xenon K shell photoabsorption) ",fontsize=16)
        ax[0].set_yscale("log")
        ax[0].legend(loc='lower left', fontsize=13)




        plt.savefig(self.plot_path + "Qseitz_compound_xe_updated.pdf")

        plt.clf()
    def Ratio_plot(self):

        fig, ax = plt.subplots(1, 2, figsize=(16, 6))

        self.Cs_list = pd.concat([self.df_Cs_116_plot, self.df_Cs_119_plot], ignore_index=True)
        self.Co_list = pd.concat([self.df_Co_116_plot, self.df_Co_119_plot], ignore_index=True)
        # self.Ba_list = pd.concat([self.df_Ba_116_plot, self.df_Ba_119_plot], ignore_index=True)
        self.Ba_list = self.df_Ba_116_plot

        self.Cs_Rate_factor = self.Cs_sims[0]
        self.Cs_energy_edges = self.Cs_sims[1][0][1]
        self.Cs_counts_cum_bin = self.Cs_sims[2]
        self.Cs_counts_energy_cum_bin = self.Cs_sims[3]

        self.Cs_doped_Rate_factor = self.Cs_sims_doped[0]
        self.Cs_doped_energy_edges = self.Cs_sims_doped[1][0][1]
        self.Cs_doped_counts_cum_bin = self.Cs_sims_doped[2]
        self.Cs_doped_counts_energy_cum_bin = self.Cs_sims_doped[3]

        self.Co_Rate_factor = self.Co_sims[0]
        self.Co_energy_edges = self.Co_sims[1][0][1]
        self.Co_counts_cum_bin = self.Co_sims[2]
        self.Co_counts_energy_cum_bin = self.Co_sims[3]

        self.Co_doped_Rate_factor = self.Co_sims_doped[0]
        self.Co_doped_energy_edges = self.Co_sims_doped[1][0][1]
        self.Co_doped_counts_cum_bin = self.Co_sims_doped[2]
        self.Co_doped_counts_energy_cum_bin = self.Co_sims_doped[3]

        self.Ba_Rate_factor = self.Ba_sims[0]
        self.Ba_energy_edges = self.Ba_sims[1][0][1]
        self.Ba_counts_cum_bin = self.Ba_sims[2]
        self.Ba_counts_energy_cum_bin = self.Ba_sims[3]

        self.Ba_doped_Rate_factor = self.Ba_sims_doped[0]
        self.Ba_doped_energy_edges = self.Ba_sims_doped[1][0][1]
        self.Ba_doped_counts_cum_bin = self.Ba_sims_doped[2]
        self.Ba_doped_counts_energy_cum_bin = self.Ba_sims_doped[3]
        print("Ba, ", self.Ba_counts_cum_bin[0], "counts," , self.Ba_doped_counts_energy_cum_bin[0], "keV" )
        print("Co, ", self.Co_counts_cum_bin[0], "counts,", self.Co_doped_counts_energy_cum_bin[0], "keV")


        merged_df_Co_Cs = pd.merge(self.Cs_list, self.Co_list, on="Seitz [keV]", suffixes=('1', '2'))
        print(merged_df_Co_Cs)
        merged_df_Co_Cs['Clean Rate ratio'] = merged_df_Co_Cs['Clean Rate [mHz]1'] / merged_df_Co_Cs['Clean Rate [mHz]2']
        merged_df_Co_Cs['Clean Rate Sigma ratio'] = merged_df_Co_Cs['Clean Rate ratio'] * np.sqrt(
            (merged_df_Co_Cs['Clean Rate Sigma [mHz]1'] / merged_df_Co_Cs['Clean Rate [mHz]1']) ** 2 + (
                        merged_df_Co_Cs['Clean Rate Sigma [mHz]2'] / merged_df_Co_Cs['Clean Rate [mHz]2']) ** 2
        )


        merged_df_Ba_Cs = pd.merge(self.Cs_list, self.Ba_list, on="Seitz [keV]", suffixes=('1', '2'))
        merged_df_Ba_Cs['Clean Rate ratio'] = merged_df_Ba_Cs['Clean Rate [mHz]1'] / merged_df_Ba_Cs['Clean Rate [mHz]2']
        merged_df_Ba_Cs['Clean Rate Sigma ratio'] = merged_df_Ba_Cs['Clean Rate ratio'] * np.sqrt(
            (merged_df_Ba_Cs['Clean Rate Sigma [mHz]1'] / merged_df_Ba_Cs['Clean Rate [mHz]1']) ** 2 + (
                        merged_df_Ba_Cs['Clean Rate Sigma [mHz]2'] / merged_df_Ba_Cs['Clean Rate [mHz]2']) ** 2
        )

        ax[0].errorbar(merged_df_Co_Cs["Seitz [keV]"],
                       merged_df_Co_Cs['Clean Rate ratio'],
                       yerr=merged_df_Co_Cs['Clean Rate Sigma ratio'], label="Cs/Co", fmt='o', color = 'r')
        # ax[1].errorbar(merged_df_Ba_Cs["Seitz [keV]"],
        #                merged_df_Ba_Cs['Clean Rate ratio'],
        #                yerr=merged_df_Ba_Cs['Clean Rate Sigma ratio'], label="Cs/Ba", fmt='o', color='r')
        ax[1].plot(merged_df_Ba_Cs["Seitz [keV]"],
                       merged_df_Ba_Cs['Clean Rate ratio'], marker='v',linestyle='None', label="Cs/Ba",  color='r')

        print("Cs/Ba clean rate ratio",merged_df_Ba_Cs['Clean Rate ratio'])

        ax[0].axhline(y=(self.Cs_Rate_factor*self.Cs_counts_cum_bin[0]/(self.Co_Rate_factor*self.Co_counts_cum_bin[0])),label = 'sim per scatter',color = 'b')
        ax[0].axhline(
            y= (self.Cs_Rate_factor * self.Cs_counts_energy_cum_bin[0])/(self.Co_Rate_factor * self.Co_counts_energy_cum_bin[0]),
            label='sim per keV',color = 'g')

        print("factor of scattering", (self.Cs_Rate_factor*self.Cs_counts_cum_bin[0]/(self.Co_Rate_factor*self.Co_counts_cum_bin[0])))
        print("factor of phot", (self.Cs_Rate_factor * self.Cs_counts_energy_cum_bin[0])/(self.Co_Rate_factor * self.Co_counts_energy_cum_bin[0]))
        ax[0].axhline(
            y=(
                        self.Cs_doped_Rate_factor * self.Cs_doped_counts_cum_bin[0])/(self.Co_doped_Rate_factor * self.Co_doped_counts_cum_bin[0]),
            label='sim per xe photo',color = 'brown')

        ax[1].axhline(
            y=(self.Cs_Rate_factor * self.Cs_counts_cum_bin[0])/(self.Ba_Rate_factor * self.Ba_counts_cum_bin[0]),
            label='sim per scatter',color = 'b')
        ax[1].axhline(
            y= (self.Cs_Rate_factor * self.Cs_counts_energy_cum_bin[0]/(self.Ba_Rate_factor * self.Ba_counts_energy_cum_bin[0])),
            label='sim per keV',color = 'g')

        print("self.Cs_Rate_factor * self.Cs_counts_energy_cum_bin[0]",self.Cs_Rate_factor * self.Cs_counts_energy_cum_bin[0])
        print("self.Ba_Rate_factor * self.Ba_counts_energy_cum_bin[0]",self.Ba_Rate_factor * self.Ba_counts_energy_cum_bin[0])

        print(self.Cs_Rate_factor * self.Cs_counts_energy_cum_bin[0]/(self.Ba_Rate_factor * self.Ba_counts_energy_cum_bin[0]))
        ax[1].axhline(
            y=  (
                    self.Cs_doped_Rate_factor * self.Cs_doped_counts_cum_bin[0])/(self.Ba_doped_Rate_factor * self.Ba_doped_counts_cum_bin[0]),
            label='sim per xe photo', color= 'brown')

        ax[0].set_xlabel(r"$Q_{Seitz} [keV]$", fontsize=14)
        ax[0].set_ylabel("Cs/Co Ratio [] ", fontsize=14)
        # ax[0].set_yscale("log")
        ax[0].legend(loc='upper right', fontsize=14)

        ax[1].set_xlabel(r"$Q_{Seitz} [keV]$", fontsize=14)
        ax[1].set_ylabel("Cs/Ba Ratio [] ", fontsize=14)
        ax[1].set_yscale("log")
        ax[1].legend(loc='upper right', fontsize=14)


        plt.savefig(self.plot_path + "Ratio_compare_updated.pdf")
        
    def time_plot(self):
        import matplotlib.dates as mdates
        print(self.df_Cs_116_time_plot)

        run_dates = [
            "11/17/2025",  # Row 0: 3.00 bara
            "11/17/2025",  # Row 1: 3.50 bara
            "11/17/2025",  # Row 2: 4.00 bara
            "12/02/2025",  # Row 3: 2.25 bara
            "12/01/2025",  # Row 4: 2.50 bara (Same run batch)
            "12/01/2025",  # Row 5: 3.00 bara
            "12/01/2025",  # Row 6: 3.50 bara
            "12/01/2025",  # Row 7: 4.50 bara
            "12/11/2025",  # Row 8: 2.25 bara
            "12/10/2025",  # Row 9: 2.50 bara
            "12/10/2025",  # Row 10: 3.00 bara
            "12/10/2025",  # Row 11: 3.50 bara
            "01/20/2026"  # Row 12: 3.00 bara
        ]

        # Add the column to your existing dataframe and cast it to datetime objects
        self.df_Cs_116_time_plot['Run_Date'] = pd.to_datetime(run_dates)


        grouped = self.df_Cs_116_time_plot.groupby("Seitz [keV]")


        fig, ax = plt.subplots(figsize=(8, 6))

        # Define your target columns for the Y-axis and its error bars
        y_column = 'Clean Rate [mHz]'
        y_error = 'Clean Rate Sigma [mHz]'

        for energy, group in grouped:
            # Condition: Only plot lines that have 2 or more points
            if len(group) >= 2:
                # Sort by index to guarantee lines connect strictly in chronological order
                group_sorted = group.sort_index()

                # X-axis: Use the index (0, 1, 2...) since it matches the time sequence
                # (Note: If you created a 'Run_Date' column earlier, change this to group_sorted['Run_Date'])
                x_values = group_sorted.index

                # Plot line with markers and error bars
                ax.errorbar(group_sorted['Run_Date'], group_sorted[y_column], yerr=group_sorted[y_error],
                            fmt='o', linewidth=2, elinewidth=1.5,label=r"$Q_{Seitz}$ "+f'{energy:.2f} keV')

            # --- 3. Format the X-axis time presentation ---
            # Formats the dates on screen as MM/DD/YYYY
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%Y'))

        # Ensures matplotlib spaces out the dates nicely
        ax.xaxis.set_major_locator(mdates.AutoDateLocator())

        # Clean up layout and rotate date strings so they don't overlap
        fig.autofmt_xdate()

        # --- 4. Labels and aesthetics ---
        ax.set_xlabel("Run Date", fontsize=12)
        ax.set_ylabel("Background Substracted Rate [mHz]", fontsize=12)
        ax.legend(loc='upper right', fontsize=14)
        plt.tight_layout()
        plt.savefig(self.plot_path + "Time_stability_updated.pdf")

    def gamma_rejection_plot_PSN_v2(self):
        # print Q vs per keV and Eion per interaction

        self.fitting_list = []
        self.Cs_fitting_list = []
        self.Co_fitting_list = []
        self.Ba_fitting_list = []
        self.df_Cs_116_plot_list = []
        self.df_Cs_119_plot_list = []
        self.df_Co_116_plot_list = []
        self.df_Co_119_plot_list = []
        self.df_Ba_116_plot_list = []
        self.df_Ba_119_plot_list = []
        self.Cs_116_label = ["Cs 11/17/2025 116K", "Cs 12/01/2025 116K", "Cs 12/10/2025 116K", "Cs 01/20/2026 116K"]
        self.Cs_119_label = ["Cs 02/02/2026 119K"]
        self.Co_116_label = ["Co 12/15/2026 116K"]
        self.Co_119_label = ["Co 02/06/2026 119K"]
        self.Ba_116_label = ["Ba 11/19/2025 116K"]

        for i in range(len(self.Cs_exp_rejection_path)):
            df = pd.read_csv(self.Cs_exp_rejection_path[i])
            # print(df.columns)
            # doc_label = self.Cs_exp_raw_path[i].replace('_exposures', '')
            # doc_label = self.Cs_label[i]
            doc_label = "Cs"
            print('doc_label', doc_label)
            # signal
            # drop 2.75,3.25, 3.75 bara pressure
            # pressure_drop_list = [2.75,3.25,3.75]
            pressure_drop_list = []
            df = df[~df['Pressure [bara]'].isin(pressure_drop_list)]
            # only positive rate
            df = df[df['Clean Rate [mHz]'] > 0]
            if i <= 3:
                self.df_Cs_116_plot_list.append(df)
            else:
                self.df_Cs_119_plot_list.append(df)

            df_fit = df[['Seitz [keV]', "Rejection Rate Scattering[]", 'Eion_rl-1_rhol-1 [GeVcm**2 g-1]',
                         "Rejection Rate KeV[/keV]", 'Q_rl-1_rhol-1 [GeVcm**2 g-1]', "Rejection Rate Xenon Abs[]",
                         'Clean Rate [mHz]', "Eion [keV]"]]

            self.fitting_list.append(df_fit)
            self.Cs_fitting_list.append(df_fit)
        self.df_Cs_116_plot = pd.concat(self.df_Cs_116_plot_list, ignore_index=True)
        self.df_Cs_119_plot = pd.concat(self.df_Cs_119_plot_list, ignore_index=True)
        self.df_Cs_116_plot = self.concat_PT_condition(self.df_Cs_116_plot)




        # ax[2].errorbar(df['Q_rl-1_rhol-1 [GeVcm**2 g-1]'], df["Rejection Rate Xenon Abs[]"],
        #                yerr=df["Rejection Sigma Xenon Abs[]"], label=doc_label, fmt='o')

        for i in range(len(self.Co_exp_rejection_path)):
            df = pd.read_csv(self.Co_exp_rejection_path[i])
            # print(df.columns)
            # doc_label = self.Co_exp_raw_path[i].rstrip("_exposures")
            # doc_label = self.Co_label[i]
            doc_label = "Co"
            # signal
            # drop 2.75,3.25, 3.75 bara pressure
            # pressure_drop_list = [2.75,3.25,3.75]
            pressure_drop_list = []
            df = df[~df['Pressure [bara]'].isin(pressure_drop_list)]
            df = df[df['Clean Rate [mHz]'] > 0]

            df_fit = df[['Seitz [keV]', "Rejection Rate Scattering[]", 'Eion_rl-1_rhol-1 [GeVcm**2 g-1]',
                         "Rejection Rate KeV[/keV]", 'Q_rl-1_rhol-1 [GeVcm**2 g-1]', "Rejection Rate Xenon Abs[]",
                         'Clean Rate [mHz]', "Eion [keV]"]]

            if i <= 0:
                self.df_Co_116_plot_list.append(df)
            else:
                self.df_Co_119_plot_list.append(df)
            self.fitting_list.append(df_fit)
            self.Co_fitting_list.append(df_fit)

        self.df_Co_116_plot = pd.concat(self.df_Co_116_plot_list, ignore_index=True)
        self.df_Co_119_plot = pd.concat(self.df_Co_119_plot_list, ignore_index=True)



        for i in range(len(self.Ba_exp_rejection_path)):
            df = pd.read_csv(self.Ba_exp_rejection_path[i])
            # print(df.columns)
            # doc_label = self.Co_exp_raw_path[i].rstrip("_exposures")
            # doc_label = self.Co_label[i]
            doc_label = "Co"
            # signal
            # drop 2.75,3.25, 3.75 bara pressure
            # pressure_drop_list = [2.75,3.25,3.75]
            pressure_drop_list = []
            df = df[~df['Pressure [bara]'].isin(pressure_drop_list)]
            # df = df[df['Clean Rate [mHz]'] > 0]
            print("ba" ,df[['Clean Rate [mHz]','Clean Rate Sigma [mHz]', 'Exp Rate [mHz]','Exp Rate Sigma [mHz]','Bkg Rate [mHz]','Bkg Rate Sigma [mHz]']])

            df_fit = df[['Seitz [keV]', "Rejection Rate Scattering[]", 'Eion_rl-1_rhol-1 [GeVcm**2 g-1]',
                         "Rejection Rate KeV[/keV]", 'Q_rl-1_rhol-1 [GeVcm**2 g-1]', "Rejection Rate Xenon Abs[]",
                         'Clean Rate [mHz]',"Eion [keV]"]]

            if i <= 0:
                self.df_Ba_116_plot_list.append(df)
            else:
                self.df_Ba_119_plot_list.append(df)
            # self.fitting_list.append(df_fit)
            # self.Ba_fitting_list.append(df_fit)
        # print("Ba", self.df_Ba_116_plot_list)
        self.df_Ba_116_plot = pd.concat(self.df_Ba_116_plot_list, ignore_index=True)
        # self.df_Ba_119_plot = pd.concat(self.df_Ba_119_plot_list, ignore_index=True)


        ratio_116 = self.df_Cs_116_plot["Rejection Rate KeV[/keV]"][0]/self.df_Cs_116_plot["Rejection Rate Xenon Abs[]"][0]
        ratio_119 = self.df_Cs_119_plot["Rejection Rate KeV[/keV]"][0] / \
                    self.df_Cs_119_plot["Rejection Rate Xenon Abs[]"][0]
        print("ratio_116",ratio_116,'ratio_119',ratio_119)
        SCALE_FACTOR = (ratio_119)**(-1)  # Xenon Abs = Rate [/keV] * SCALE_FACTOR

        # Fix: changed subplots(1, 0) to subplots()
        fig, ax = plt.subplots(2,2, figsize=(15, 13))

        # Plot Cs 116K ONCE on the left axis
        ax[0,0].errorbar(
            self.df_Cs_116_plot["Seitz [keV]"],
            self.df_Cs_116_plot["Rejection Rate KeV[/keV]"],
            yerr=self.df_Cs_116_plot["Rejection Sigma KeV[/keV]"],
            label="SBC (Ar+Xe) 116K",
            fmt='o',
            markersize=8,
            color="tab:brown"  # Give datasets distinct colors
        )

        # Plot Cs 119K ONCE on the left axis
        ax[0,0].errorbar(
            self.df_Cs_119_plot["Seitz [keV]"],
            self.df_Cs_119_plot["Rejection Rate KeV[/keV]"],
            yerr=self.df_Cs_119_plot["Rejection Sigma KeV[/keV]"],
            label="SBC (Ar+Xe) 119K",
            fmt='s',
            markersize=8,
            color="tab:green"
        )

        # Set main (left) y-axis and x-axis labels
        ax[0,0].set_xlabel(r"Seitz threshold [keV]", fontsize=16)
        ax[0,0].set_xlim(0.65, 2.8)
        ax[0,0].set_ylim(1e-12, 1e-4)
        ax[0,0].set_ylabel("Probability per energy deposited (events/keV) ", fontsize=16)
        ax[0,0].set_yscale("log")
        ax[0,0].yaxis.label.set_color("red")
        ax[0,0].tick_params(axis='y', colors="red", which='both')  # 'both' colors major & minor ticks
        ax[0,0].spines['left'].set_color("red")


        # Add secondary (right) y-axis with proportional mapping
        def forward(y):
            return y * SCALE_FACTOR

        def inverse(y):
            return y / SCALE_FACTOR

        secax0 = ax[0,0].secondary_yaxis('right', functions=(forward, inverse))
        secax0.set_ylabel("Nucleation probability\n(per xenon photoabsorption in K shell) ", fontsize=16)
        secax0.yaxis.label.set_color("blue")
        secax0.tick_params(axis='y', colors="blue", which='both')
        secax0.spines['right'].set_color("blue")

        # plot the fitting lines
        self.Cs_df = pd.concat(self.Cs_fitting_list, ignore_index=True)
        [result_Q_scatter, result_Q_keV, result_Q_xe, result_Eion_scatter, result_Eion_keV, result_Eion_xe,
         result_Q2_xe, result_Q_rate] = self.fitting_gamma_rejection_v2(self.Cs_df)

        # SBC
        SBC_Q_list = result_Q_keV[2]
        SBC_keV_list = result_Q_keV[3]
        SBC_Q2_list = result_Q2_xe[2]
        SBC_Eion_list = result_Eion_keV[2]
        SBC_Q_kev_fitting = (result_Q_keV[0],result_Q_keV[1])
        SBC_Q_xe_fitting = (result_Q_xe[0], result_Q_xe[1])
        print("fitting  SBC_Q_kev_fitting", SBC_Q_kev_fitting)
        print("fitting SBC_Q_xe_fitting ", SBC_Q_xe_fitting )





        # Compute ratio list

        thermal_116_table = self.dict_energy_116_tab
        thermal_119_table = self.dict_energy_119_tab

        # result = self.interpolate_all_keys("Eion_rl-1_rhol-1 [GeVcm**2 g-1]", 0.8, thermal_116_table)
        # print(result)

        SBC_rrho_list = [SBC_Q_list[i] / SBC_Q2_list[i] for i in range(len(SBC_Q_list))]
        SBC_Eion_list = [SBC_Q_list[i] / SBC_Eion_list[i] for i in range(len(SBC_Q_list))]

        # Drexel Q2 to Xenon calculations
        SBC_fitting_len = len(SBC_Q_list)
        # Use np.linspace so length matches SBC_fitting_len exactly
        Drex_Q2_list = np.linspace(1.5, 4, SBC_fitting_len)
        Drex_phot_list = 58 * np.exp(-Drex_Q2_list / 0.2877)
        Drex_Q_116_list = self.interpolate_all_keys_vectorized("Q_rl-1_rhol-1 [GeVcm**2 g-1]", Drex_Q2_list , thermal_116_table)["Seitz [keV]"]

        # Drex_Q_119_list = [i+0.05 for i in Drex_Q_116_list]
        Drex_Q_119_list = \
        self.interpolate_all_keys_vectorized("Q_rl-1_rhol-1 [GeVcm**2 g-1]", Drex_Q2_list, thermal_119_table)[
            "Seitz [keV]"]

        # PICO Eion to keV calculations
        SBC_fitting_len = len(SBC_Q_list)
        # Use np.linspace so length matches SBC_fitting_len exactly
        PICO_Eion_list = np.linspace(0.85, 1.5, SBC_fitting_len)
        PICO_keV_list = 17e3 * np.exp(-PICO_Eion_list / 37e-3)
        PICO_Q_116_list = self.interpolate_all_keys_vectorized("Eion_rl-1_rhol-1 [GeVcm**2 g-1]", PICO_Eion_list , thermal_116_table)["Seitz [keV]"]

        # PICO_Q_119_list = [i+0.05 for i in PICO_Q_116_list]
        PICO_Q_119_list = \
        self.interpolate_all_keys_vectorized("Eion_rl-1_rhol-1 [GeVcm**2 g-1]", PICO_Eion_list, thermal_119_table)[
            "Seitz [keV]"]



        # ax[0,0].plot(
        #     Drex_Q_116_list,
        #     Drex_phot_list / SCALE_FACTOR,
        #     label="Drexel (C$_3$F$_8$+Xe) 116K",
        #
        #     color="blue"
        # )
        # ax[0,0].plot(
        #     Drex_Q_119_list,
        #     Drex_phot_list / SCALE_FACTOR,
        #     label="Drexel (C$_3$F$_8$+Xe) 119K", linestyle= '--',
        #
        #     color="blue"
        # )

        ax[0, 0].fill_betweenx(
            Drex_phot_list / SCALE_FACTOR,
            Drex_Q_116_list,
            Drex_Q_119_list,
            color="blue",
            alpha=0.3,
            label="Drexel (C$_3$F$_8$+Xe)"
        )

        # ax[0,0].plot(PICO_Q_116_list, PICO_keV_list, label="PICO C$_3$F$_8$ 116K", color="red")
        # ax[0,0].plot(PICO_Q_119_list, PICO_keV_list, label="PICO C$_3$F$_8$ 119K",  linestyle= '--', color="red")
        ax[0, 0].fill_betweenx(
            PICO_keV_list,
            PICO_Q_116_list,
            PICO_Q_119_list,
            color="red",
            alpha=0.3,
            label="PICO (C$_3$F$_8$)")



        ax[0,0].plot(SBC_Q_list, SBC_keV_list, label="SBC Best Fit", color="black")

        # fitting parameter

        # box_content0 = (f"$\\mathcal{{P}}_{{phot}} = A_{{phot}} e^{{-B_{{phot}} Q_{{Seitz}}}}$\n"
        #                 f"$A_{{phot}}$ = 0.014 K-phot$^{{-1}}$\n"    f"$B_{{phot}}$ = 4.289 keV$^{{-1}}$")
        #
        # ax[0,0].text(0.65, 0.78, box_content0,
        #         transform=ax[0,0].transAxes,
        #         fontsize=16,
        #         color='black',  # White text color
        #         verticalalignment='top',
        #         horizontalalignment='left',
        #         linespacing=1.4,  # Extra padding between lines
        #         bbox=dict(
        #             facecolor='none',  # Black background
        #             edgecolor='none',  # No border outline
        #             alpha=0.9  # Slight transparency so gridlines don't completely disappear
        #         ))
        # box_content1 = (f"$\\mathcal{{P}} = A e^{{-B E_{{ion}} / r_\\ell \\rho_\\ell}}$\n"
        #                f"A = 0.13 MeV$^{{-1}}$\n"    f"B = 8.75 keV$^{{-1}}$cm$^{{-2}}$g")
        #
        #
        #
        # box_content1 = (f"$\\mathcal{{P}}_{{edep}} = A_{{edep}} e^{{-B_{{edep}} Q_{{Seitz}}}}$\n"
        #                 f"$A_{{edep}}$ = 2.468 GeV$^{{-1}}$\n"    f"$B_{{edep}}$ = 4.289 keV$^{{-1}}$")
        #
        #
        # ax[0,0].text(0.65, 0.98, box_content1,
        #            transform=ax[0,0].transAxes,
        #            fontsize=16,
        #            color='black',  # White text color
        #            verticalalignment='top',
        #            horizontalalignment='left',
        #            linespacing=1.4,  # Extra padding between lines
        #            bbox=dict(
        #                facecolor='none',  # Black background
        #                edgecolor='none',  # No border outline
        #                alpha=0.9  # Slight transparency so gridlines don't completely disappear
        #            ))
        #
        # box_content2 = (f"Drexel (C$_3$F$_8$+Xe)")
        #
        # ax[0,0].text(0.60, 0.85, box_content2,
        #         transform=ax[0,0].transAxes,
        #         fontsize=16,
        #         color='blue',  # White text color
        #         verticalalignment='top',
        #         horizontalalignment='left',
        #         linespacing=1.4,  # Extra padding between lines
        #         bbox=dict(
        #             facecolor='none',  # Black background
        #             edgecolor='none',  # No border outline
        #             alpha=0.9  # Slight transparency so gridlines don't completely disappear
        #         ))
        #
        # box_content3 = (f"PICO C$_3$F$_8$")
        #
        # # ax_ij.legend(loc='lower left', fontsize=16, title="Source and \nTemperature", title_fontsize=16,frameon=False)
        # ax[0,0].text(0.30, 0.65, box_content3,
        #         transform=ax[0,0].transAxes,
        #         fontsize=16,
        #         color='red',  # White text color
        #         verticalalignment='top',
        #         horizontalalignment='left',
        #         linespacing=1.4,  # Extra padding between lines
        #         bbox=dict(
        #             facecolor='none',  # Black background
        #             edgecolor='none',  # No border outline
        #             alpha=0.9  # Slight transparency so gridlines don't completely disappear
        #         ))




        # PICO Eion_keV

        ax[0,0].legend(loc='lower left', fontsize=16, title=" ", title_fontsize=16,frameon=False)

       # 2nd graph that use c3F8 mapping:
        ax[0, 1].errorbar(
            self.df_Cs_116_plot["Seitz [keV]"],
            self.df_Cs_116_plot["Rejection Rate KeV[/keV]"],
            yerr=self.df_Cs_116_plot["Rejection Sigma KeV[/keV]"],
            label="SBC (Ar+Xe) 116K",
            fmt='o',
            markersize=8,
            color="tab:brown"  # Give datasets distinct colors
        )

        # Plot Cs 119K ONCE on the left axis
        ax[0, 1].errorbar(
            self.df_Cs_119_plot["Seitz [keV]"],
            self.df_Cs_119_plot["Rejection Rate KeV[/keV]"],
            yerr=self.df_Cs_119_plot["Rejection Sigma KeV[/keV]"],
            label="SBC (Ar+Xe) 119K",
            fmt='s',
            markersize=8,
            color="tab:green"
        )

        # Set main (left) y-axis and x-axis labels
        ax[0, 1].set_xlabel(r"Seitz threshold [keV]", fontsize=16)
        ax[0, 1].set_xlim(0.65, 2.8)
        ax[0, 1].set_ylim(1e-12, 1e-4)
        ax[0, 1].set_ylabel("Probability per energy deposited (events/keV) ", fontsize=16)
        ax[0, 1].set_yscale("log")
        ax[0, 1].yaxis.label.set_color("red")
        ax[0, 1].tick_params(axis='y', colors="red", which='both')  # 'both' colors major & minor ticks
        ax[0, 1].spines['left'].set_color("red")

        # Add secondary (right) y-axis with proportional mapping
        def forward(y):
            return y * SCALE_FACTOR

        def inverse(y):
            return y / SCALE_FACTOR

        secax1 = ax[0, 1].secondary_yaxis('right', functions=(forward, inverse))
        secax1.set_ylabel("Nucleation probability\n(per xenon photoabsorption in K shell) ", fontsize=16)
        secax1.yaxis.label.set_color("blue")
        secax1.tick_params(axis='y', colors="blue", which='both')
        secax1.spines['right'].set_color("blue")

        # plot the fitting lines
        self.Cs_df = pd.concat(self.Cs_fitting_list, ignore_index=True)
        [result_Q_scatter, result_Q_keV, result_Q_xe, result_Eion_scatter, result_Eion_keV, result_Eion_xe,
         result_Q2_xe, result_Q_rate] = self.fitting_gamma_rejection_v2(self.Cs_df)

        # SBC
        SBC_Q_list = result_Q_keV[2]
        SBC_keV_list = result_Q_keV[3]
        SBC_Q2_list = result_Q2_xe[2]
        SBC_Eion_list = result_Eion_keV[2]
        SBC_Q_kev_fitting = (result_Q_keV[0], result_Q_keV[1])
        SBC_Q_xe_fitting = (result_Q_xe[0], result_Q_xe[1])
        print("fitting  SBC_Q_kev_fitting", SBC_Q_kev_fitting)
        print("fitting SBC_Q_xe_fitting ", SBC_Q_xe_fitting)

        # Compute ratio list

        thermal_10_table = self.dict_energy_10_tab
        thermal_24_table = self.dict_energy_24_tab

        # result = self.interpolate_all_keys("Eion_rl-1_rhol-1 [GeVcm**2 g-1]", 0.8, thermal_10_table)
        # print(result)

        SBC_rrho_list = [SBC_Q_list[i] / SBC_Q2_list[i] for i in range(len(SBC_Q_list))]
        SBC_Eion_list = [SBC_Q_list[i] / SBC_Eion_list[i] for i in range(len(SBC_Q_list))]

        # Drexel Q2 to Xenon calculations
        SBC_fitting_len = len(SBC_Q_list)
        # Use np.linspace so length matches SBC_fitting_len exactly
        Drex_Q2_list = np.linspace(1.5, 4, SBC_fitting_len)
        Drex_phot_list = 58 * np.exp(-Drex_Q2_list / 0.2877)
        Drex_Q_10_list = \
        self.interpolate_all_keys_vectorized("Q_rl-1_rhol-1 [GeVcm**2 g-1]", Drex_Q2_list, thermal_10_table)[
            "Seitz [keV]"]

        # Drex_Q_24_list = [i+0.05 for i in Drex_Q_10_list]
        Drex_Q_24_list = \
            self.interpolate_all_keys_vectorized("Q_rl-1_rhol-1 [GeVcm**2 g-1]", Drex_Q2_list, thermal_24_table)[
                "Seitz [keV]"]

        # PICO Eion to keV calculations
        SBC_fitting_len = len(SBC_Q_list)
        # Use np.linspace so length matches SBC_fitting_len exactly
        PICO_Eion_list = np.linspace(0.85, 1.5, SBC_fitting_len)
        PICO_keV_list = 17e3 * np.exp(-PICO_Eion_list / 37e-3)
        PICO_Q_10_list = \
        self.interpolate_all_keys_vectorized("Eion_rl-1_rhol-1 [GeVcm**2 g-1]", PICO_Eion_list, thermal_10_table)[
            "Seitz [keV]"]

        # PICO_Q_24_list = [i+0.05 for i in PICO_Q_10_list]
        PICO_Q_24_list = \
            self.interpolate_all_keys_vectorized("Eion_rl-1_rhol-1 [GeVcm**2 g-1]", PICO_Eion_list, thermal_24_table)[
                "Seitz [keV]"]



        ax[0, 1].fill_betweenx(
            Drex_phot_list / SCALE_FACTOR,
            Drex_Q_10_list,
            Drex_Q_24_list,
            color="blue",
            alpha=0.3,
            label="Drexel (C$_3$F$_8$+Xe)"
        )

        ax[0, 1].fill_betweenx(
            PICO_keV_list,
            PICO_Q_10_list,
            PICO_Q_24_list,
            color="red",
            alpha=0.3,
            label="PICO (C$_3$F$_8$)")

        ax[0, 1].plot(SBC_Q_list, SBC_keV_list, label="SBC Best Fit", color="black")

        ax[0, 1].legend(loc='lower left', fontsize=16, title=" ", title_fontsize=16, frameon=False)




                     # compare to PICO
        ax[1, 0].errorbar(
            self.df_Cs_116_plot["Eion_rl-1_rhol-1 [GeVcm**2 g-1]"],
            self.df_Cs_116_plot["Rejection Rate KeV[/keV]"],
            yerr=self.df_Cs_116_plot["Rejection Sigma KeV[/keV]"],
            label="SBC (Ar+Xe) 116K",
            fmt='o',
            markersize=8,
            color="tab:brown"  # Give datasets distinct colors
        )

        # Plot Cs 119K ONCE on the left axis
        ax[1, 0].errorbar(
            self.df_Cs_119_plot["Eion_rl-1_rhol-1 [GeVcm**2 g-1]"],
            self.df_Cs_119_plot["Rejection Rate KeV[/keV]"],
            yerr=self.df_Cs_119_plot["Rejection Sigma KeV[/keV]"],
            label="SBC (Ar+Xe) 119K",
            fmt='s',
            markersize=8,
            color="tab:green"
        )
        ax[1, 0].plot(result_Eion_keV[2], result_Eion_keV[3], label="SBC Best Fit", color="black")


        # pICO result
        A = 17e3  # 0.017
        B = 37e-3  # 0.037

        # Generate x array (e.g., 100 evenly spaced points from 0.8 to 1.5)
        x_pico = np.linspace(0.8, 1.5, 100)

        # Calculate y array
        y_pico = A * np.exp(- x_pico / B)
        print("y_pico", y_pico[:10])
        ax[1, 0].plot(x_pico, y_pico,
                   color="black", linestyle='--', label="PICO Best Fit")

        # Set main (left) y-axis and x-axis labels
        ax[1, 0].set_xlabel(r"$E_{ion} r_l^{-1} \rho_l^{-1}$ [GeV cm$^2$ g$^{-1}$]", fontsize=16)
        ax[1, 0].set_xlim(0.7, 1.5)
        ax[1, 0].set_ylim(1e-12, 1e-5)
        ax[1, 0].set_ylabel("Probability per energy deposited [events/keV] ", fontsize=16)
        ax[1, 0].set_yscale("log")
        ax[1, 0].tick_params(axis='y',  which='both')  # 'both' colors major & minor ticks

        ax[1, 0].legend(loc='lower left', fontsize=16, title=" ", title_fontsize=16, frameon=False)

        # Drexel result
        ax[1, 1].errorbar(
            self.df_Cs_116_plot["Q_rl-1_rhol-1 [GeVcm**2 g-1]"],
            self.df_Cs_116_plot["Rejection Rate Xenon Abs[]"],
            yerr=self.df_Cs_116_plot["Rejection Sigma Xenon Abs[]"],
            label="SBC (Ar+Xe) 116K",
            fmt='o',
            markersize=8,
            color="tab:brown"  # Give datasets distinct colors
        )

        # Plot Cs 119K ONCE on the left axis
        ax[1, 1].errorbar(
            self.df_Cs_119_plot["Q_rl-1_rhol-1 [GeVcm**2 g-1]"],
            self.df_Cs_119_plot["Rejection Rate Xenon Abs[]"],
            yerr=self.df_Cs_119_plot["Rejection Sigma Xenon Abs[]"],
            label="SBC (Ar+Xe) 119K",
            fmt='s',
            markersize=8,
            color="tab:green"
        )
        ax[1, 1].plot(result_Q2_xe[2], result_Q2_xe[3], label="SBC Best Fit", color="black")

        A = 58  # 0.017
        B = 0.287  # 0.037

        # Generate x array (e.g., 100 evenly spaced points from 0.8 to 1.5)
        x_drexel = np.linspace(1.5, 4, 100)

        # Calculate y array
        y_drexel = A * np.exp(- x_drexel / B)
        ax[1,1].plot(x_drexel, y_drexel,
                   color="black", linestyle='--', label="Drexel Best Fit")

        ax[1, 1].set_xlabel(r"$Q_{Seitz} r_l^{-1} \rho_l^{-1}$ [GeV cm$^2$ g$^{-1}$]", fontsize=16)
        ax[1, 1].set_xlim(1.2, 3.2)
        ax[1, 1].set_ylim(1e-7, 1e0)
        ax[1, 1].set_ylabel("Nucleation probability\n[per xenon K shell photoabsorption]", fontsize=16)
        ax[1, 1].set_yscale("log")
        ax[1, 1].tick_params(axis='y', which='both')  # 'both' colors major & minor ticks
        ax[1, 1].legend(loc='lower left', fontsize=16, title=" ", title_fontsize=16, frameon=False)





        #compare to Drexel
        plt.tight_layout()
        # plt.show()
        plt.savefig(self.plot_path + "gamma_rejection_PSN_v2.pdf")

    def interpolate_all_keys_vectorized(self, target_key, target_values, data_dict):
        """
        Given a target_key and an array/list of target_values, returns a dictionary
        where each key contains the array of interpolated and extrapolated values.
        """
        if target_key not in data_dict:
            raise KeyError(f"Key '{target_key}' not found in data dictionary.")

        # Ensure inputs are NumPy arrays
        x = np.asarray(data_dict[target_key], dtype=float)
        target_values = np.asarray(target_values, dtype=float)

        # Sort x once
        sort_idx = np.argsort(x)
        x_sorted = x[sort_idx]

        result = {}
        for key, values in data_dict.items():
            if key == target_key:
                result[key] = target_values
            else:
                y_sorted = np.asarray(values, dtype=float)[sort_idx]

                # Create a linear interpolator with extrapolation enabled
                f = interp1d(x_sorted, y_sorted, kind='linear', fill_value='extrapolate')
                result[key] = f(target_values)

        return result

    def interpolate_all_keys(self, target_key, target_value, data_dict):
        """
        Given a target_key and its value target_value, interpolates and returns
        the corresponding interpolated values for ALL keys in data_dict.
        """
        if target_key not in data_dict:
            raise KeyError(f"Key '{target_key}' not found in data dictionary.")

        # x is the array corresponding to the input key
        x = np.array(data_dict[target_key])

        # np.interp requires x-coordinates to be strictly increasing
        sort_idx = np.argsort(x)
        x_sorted = x[sort_idx]

        result = {}
        for key, values in data_dict.items():
            if key == target_key:
                result[key] = float(target_value)
            else:
                y_sorted = np.array(values)[sort_idx]
                # Interpolate y at target_value based on x
                interpolated_val = np.interp(target_value, x_sorted, y_sorted)
                result[key] = float(interpolated_val)

        return result
    def gamma_rejection_plot_PSN(self):
        # print Q vs per keV and Eion per interaction

        self.fitting_list = []
        self.Cs_fitting_list = []
        self.Co_fitting_list = []
        self.Ba_fitting_list = []
        self.Hot_Cs_fitting_list = []
        self.Th_fitting_list = []
        self.df_Cs_116_plot_list = []
        self.df_Cs_119_plot_list = []
        self.df_Co_116_plot_list = []
        self.df_Co_119_plot_list = []
        self.df_Ba_116_plot_list = []
        self.df_Ba_119_plot_list = []
        self.df_Hot_Cs_116_plot_list = []
        self.df_Hot_Cs_119_plot_list = []
        self.df_Th_116_plot_list = []
        self.df_Th_119_plot_list = []

        self.Cs_116_label = ["Cs 11/17/2025 116K", "Cs 12/01/2025 116K", "Cs 12/10/2025 116K", "Cs 01/20/2026 116K"]
        self.Cs_119_label = ["Cs 02/02/2026 119K"]
        self.Co_116_label = ["Co 12/15/2026 116K"]
        self.Co_119_label = ["Co 02/06/2026 119K"]
        self.Ba_116_label = ["Ba 11/19/2025 116K"]
        self.Hot_Cs_116_label = ["Hot Cs 11/11/2025 116K"]
        self.Th_116_label = ["Th 11/20/2025 116K"]



        for i in range(len(self.Cs_exp_rejection_path)):
            df = pd.read_csv(self.Cs_exp_rejection_path[i])
            # print(df.columns)
            # doc_label = self.Cs_exp_raw_path[i].replace('_exposures', '')
            # doc_label = self.Cs_label[i]
            doc_label = "Cs"
            print('doc_label', doc_label)
            # signal
            # drop 2.75,3.25, 3.75 bara pressure
            # pressure_drop_list = [2.75,3.25,3.75]
            pressure_drop_list = []
            df = df[~df['Pressure [bara]'].isin(pressure_drop_list)]
            # only positive rate
            df = df[df['Clean Rate [mHz]'] > 0]
            if i <= 3:
                self.df_Cs_116_plot_list.append(df)
            else:
                self.df_Cs_119_plot_list.append(df)

            df_fit = df[['Seitz [keV]', "Rejection Rate Scattering[]", 'Eion_rl-1_rhol-1 [GeVcm**2 g-1]',
                         "Rejection Rate KeV[/keV]", 'Q_rl-1_rhol-1 [GeVcm**2 g-1]', "Rejection Rate Xenon Abs[]",
                         'Clean Rate [mHz]', "Eion [keV]"]]

            self.fitting_list.append(df_fit)
            self.Cs_fitting_list.append(df_fit)
        self.df_Cs_116_plot = pd.concat(self.df_Cs_116_plot_list, ignore_index=True)
        self.df_Cs_119_plot = pd.concat(self.df_Cs_119_plot_list, ignore_index=True)

        self.df_Cs_116_time_plot=self.df_Cs_116_plot
        # make Cs 116 show just as one series
        self.df_Cs_116_plot = self.concat_PT_condition(self.df_Cs_116_plot)




        # ax[2].errorbar(df['Q_rl-1_rhol-1 [GeVcm**2 g-1]'], df["Rejection Rate Xenon Abs[]"],
        #                yerr=df["Rejection Sigma Xenon Abs[]"], label=doc_label, fmt='o')

        for i in range(len(self.Co_exp_rejection_path)):
            df = pd.read_csv(self.Co_exp_rejection_path[i])
            # print(df.columns)
            # doc_label = self.Co_exp_raw_path[i].rstrip("_exposures")
            # doc_label = self.Co_label[i]
            doc_label = "Co"
            # signal
            # drop 2.75,3.25, 3.75 bara pressure
            # pressure_drop_list = [2.75,3.25,3.75]
            pressure_drop_list = []
            df = df[~df['Pressure [bara]'].isin(pressure_drop_list)]
            df = df[df['Clean Rate [mHz]'] > 0]

            df_fit = df[['Seitz [keV]', "Rejection Rate Scattering[]", 'Eion_rl-1_rhol-1 [GeVcm**2 g-1]',
                         "Rejection Rate KeV[/keV]", 'Q_rl-1_rhol-1 [GeVcm**2 g-1]', "Rejection Rate Xenon Abs[]",
                         'Clean Rate [mHz]', "Eion [keV]"]]

            if i <= 0:
                self.df_Co_116_plot_list.append(df)
            else:
                self.df_Co_119_plot_list.append(df)
            self.fitting_list.append(df_fit)
            self.Co_fitting_list.append(df_fit)

        self.df_Co_116_plot = pd.concat(self.df_Co_116_plot_list, ignore_index=True)
        self.df_Co_119_plot = pd.concat(self.df_Co_119_plot_list, ignore_index=True)



        for i in range(len(self.Ba_exp_rejection_path)):
            df = pd.read_csv(self.Ba_exp_rejection_path[i])
            # print(df.columns)
            # doc_label = self.Co_exp_raw_path[i].rstrip("_exposures")
            # doc_label = self.Co_label[i]
            doc_label = "Co"
            # signal
            # drop 2.75,3.25, 3.75 bara pressure
            # pressure_drop_list = [2.75,3.25,3.75]
            pressure_drop_list = []
            df = df[~df['Pressure [bara]'].isin(pressure_drop_list)]
            # df = df[df['Clean Rate [mHz]'] > 0]
            print("ba" ,df[['Clean Rate [mHz]','Clean Rate Sigma [mHz]', 'Exp Rate [mHz]','Exp Rate Sigma [mHz]','Bkg Rate [mHz]','Bkg Rate Sigma [mHz]']])

            df_fit = df[['Seitz [keV]', "Rejection Rate Scattering[]", 'Eion_rl-1_rhol-1 [GeVcm**2 g-1]',
                         "Rejection Rate KeV[/keV]", 'Q_rl-1_rhol-1 [GeVcm**2 g-1]', "Rejection Rate Xenon Abs[]",
                         'Clean Rate [mHz]',"Eion [keV]"]]

            if i <= 0:
                self.df_Ba_116_plot_list.append(df)
            else:
                self.df_Ba_119_plot_list.append(df)
            # self.fitting_list.append(df_fit)
            # self.Ba_fitting_list.append(df_fit)
        # print("Ba", self.df_Ba_116_plot_list)
        self.df_Ba_116_plot = pd.concat(self.df_Ba_116_plot_list, ignore_index=True)


        for i in range(len(self.Hot_Cs_exp_rejection_path)):
            df = pd.read_csv(self.Hot_Cs_exp_rejection_path[i])
            # print(df.columns)
            # doc_label = self.Co_exp_raw_path[i].rstrip("_exposures")
            # doc_label = self.Co_label[i]
            doc_label = "Co"
            # signal
            # drop 2.75,3.25, 3.75 bara pressure
            # pressure_drop_list = [2.75,3.25,3.75]
            pressure_drop_list = []
            df = df[~df['Pressure [bara]'].isin(pressure_drop_list)]
            # df = df[df['Clean Rate [mHz]'] > 0]
            print("ba" ,df[['Clean Rate [mHz]','Clean Rate Sigma [mHz]', 'Exp Rate [mHz]','Exp Rate Sigma [mHz]','Bkg Rate [mHz]','Bkg Rate Sigma [mHz]']])

            df_fit = df[['Seitz [keV]', "Rejection Rate Scattering[]", 'Eion_rl-1_rhol-1 [GeVcm**2 g-1]',
                         "Rejection Rate KeV[/keV]", 'Q_rl-1_rhol-1 [GeVcm**2 g-1]', "Rejection Rate Xenon Abs[]",
                         'Clean Rate [mHz]',"Eion [keV]"]]

            if i <= 0:
                self.df_Hot_Cs_116_plot_list.append(df)
            else:
                self.df_Hot_Cs_119_plot_list.append(df)
            # self.fitting_list.append(df_fit)
            # self.Ba_fitting_list.append(df_fit)
        # print("Ba", self.df_Ba_116_plot_list)
        self.df_Hot_Cs_116_plot = pd.concat(self.df_Hot_Cs_116_plot_list, ignore_index=True)

        # self.df_Ba_119_plot = pd.concat(self.df_Ba_119_plot_list, ignore_index=True)


        for i in range(len(self.Th_exp_rejection_path)):
            df = pd.read_csv(self.Th_exp_rejection_path[i])
            # print(df.columns)
            # doc_label = self.Co_exp_raw_path[i].rstrip("_exposures")
            # doc_label = self.Co_label[i]
            doc_label = "Co"
            # signal
            # drop 2.75,3.25, 3.75 bara pressure
            # pressure_drop_list = [2.75,3.25,3.75]
            pressure_drop_list = []
            df = df[~df['Pressure [bara]'].isin(pressure_drop_list)]
            # df = df[df['Clean Rate [mHz]'] > 0]
            print("ba" ,df[['Clean Rate [mHz]','Clean Rate Sigma [mHz]', 'Exp Rate [mHz]','Exp Rate Sigma [mHz]','Bkg Rate [mHz]','Bkg Rate Sigma [mHz]']])

            df_fit = df[['Seitz [keV]', "Rejection Rate Scattering[]", 'Eion_rl-1_rhol-1 [GeVcm**2 g-1]',
                         "Rejection Rate KeV[/keV]", 'Q_rl-1_rhol-1 [GeVcm**2 g-1]', "Rejection Rate Xenon Abs[]",
                         'Clean Rate [mHz]',"Eion [keV]"]]

            if i <= 0:
                self.df_Th_116_plot_list.append(df)
            else:
                self.df_Th_119_plot_list.append(df)
            # self.fitting_list.append(df_fit)
            # self.Th_fitting_list.append(df_fit)
        # print("Th", self.df_Th_116_plot_list)
        self.df_Th_116_plot = pd.concat(self.df_Th_116_plot_list, ignore_index=True)

        fig, ax = plt.subplots(3, 4, figsize=(40, 24))
        # fig, ax = plt.subplots(2, 1, figsize=(6, 10))

        y_config = [{"y": "Rejection Rate Scattering[]", "y_err": "Rejection Sigma Scattering[]",
                     "ylabel": "Nucleation probability (per interaction) "},
                    {"y": "Rejection Rate KeV[/keV]", "y_err": "Rejection Sigma KeV[/keV]",
                     "ylabel": "Probability per energy deposited (events/keV) "},
                    {"y": "Rejection Rate Xenon Abs[]", "y_err": "Rejection Sigma Xenon Abs[]",
                     "ylabel": "Nucleation probability (per xenon photoabsorption in K shell) "},
                    {"y": "Clean Rate [mHz]", "y_err": 'Clean Rate Sigma [mHz]',
                     "ylabel": "Background Substacted Rate [mHz]"}]
        x_config = [{"x": "Seitz [keV]", "xlabel": r"Seitz threshold [keV]"},
                    {"x": 'Eion_rl-1_rhol-1 [GeVcm**2 g-1]',
                     "xlabel": r"$E_{ion} r_l^{-1} \rho_l^{-1}$ [GeV cm$^2$ g$^{-1}$]"},
                    {"x": "Eion [keV]", "xlabel": r"$E_{ion}$"}]

        for i in range(4):
            for j in range(3):
                # Extract the configuration for this specific slot
                y_cfg = y_config[i]
                x_cfg = x_config[j]
                ax_ij = ax[j, i]

                ax_ij.errorbar(self.df_Cs_116_plot[x_cfg["x"]], self.df_Cs_116_plot[y_cfg["y"]],
                           yerr=self.df_Cs_116_plot[y_cfg["y_err"]], label="Cs 116K", fmt='o',markersize=8,color = "r")
                ax_ij.errorbar(self.df_Cs_119_plot[x_cfg["x"]], self.df_Cs_119_plot[y_cfg["y"]],
                           yerr=self.df_Cs_119_plot[y_cfg["y_err"]], label="Cs 119K", fmt='^',markersize=8,color = "orange")
                ax_ij.errorbar(self.df_Co_116_plot[x_cfg["x"]], self.df_Co_116_plot[y_cfg["y"]],
                           yerr=self.df_Co_116_plot[y_cfg["y_err"]], label="Co 116K", fmt='D',markersize=8,color = "green")
                ax_ij.errorbar(self.df_Co_119_plot[x_cfg["x"]], self.df_Co_119_plot[y_cfg["y"]],
                           yerr=self.df_Co_119_plot[y_cfg["y_err"]], label="Co 119K", fmt='s',markersize=8,color = "blue")
                # ax_ij.errorbar(self.df_Ba_116_plot[x_cfg["x"]], self.df_Ba_116_plot[y_cfg["y"]],
                #            yerr=self.df_Ba_116_plot[y_cfg["y_err"]], label="Ba 116K", fmt='o')
                ax_ij.errorbar(self.df_Hot_Cs_116_plot[x_cfg["x"]], self.df_Hot_Cs_116_plot[y_cfg["y"]],
                           yerr=self.df_Hot_Cs_116_plot[y_cfg["y_err"]], label="Hot Cs 116K", fmt='o',markersize=8,color = "purple")
                ax_ij.errorbar(self.df_Th_116_plot[x_cfg["x"]], self.df_Th_116_plot[y_cfg["y"]],
                           yerr=self.df_Th_116_plot[y_cfg["y_err"]], label="Th 116K", fmt='D',markersize=8,color = "brown")
                # ax_ij.plot(self.df_Ba_116_plot[x_cfg["x"]], self.df_Ba_116_plot[y_cfg["y"]],
                #                label="Ba 116K 95% CL \nUpper Limit", marker='v',linestyle='None')

                ax_ij.set_xlabel(x_cfg["xlabel"],fontsize=16)
                ax_ij.set_ylabel(y_cfg["ylabel"], fontsize=16)
                ax_ij.set_yscale("log")
                ax_ij.legend(loc='lower left', fontsize=16, title="Source and \nTemperature", title_fontsize=16,frameon=False)






        self.fitting_df = pd.concat(self.fitting_list, ignore_index=True)
        fitting_matrix = self.fitting_gamma_rejection_v3(self.fitting_df,x_config,y_config)


        # plot the fitting function
        for i in range(4):
            for j in range(3):
                # Extract the configuration for this specific slot
                y_cfg = y_config[i]
                x_cfg = x_config[j]
                ax_ij = ax[j, i]
                a_val = fitting_matrix[i][j][0]
                b_val = fitting_matrix[i][j][1]

                # ax_ij.plot(fitting_matrix[i][j][2], fitting_matrix[i][j][3],
                #            color="black")

                ax_ij.plot(fitting_matrix[i][j][2], fitting_matrix[i][j][3],
                      color="black", label="SBC Best Fit")
                if i==1 & j==1:




                    a_str = self.fmt_sci_tex(a_val)
                    b_str = self.fmt_sci_tex(b_val)
                    print("value",a_str,b_str)
                    # box_content = (  f"$y = A e^{{-Bx}}$\n"
                    #     f"A = ${a_str}$ keV$^{{-1}}$\n"    f"B = ${b_str}$ GeV$^{{-1}}\\cdot$cm$^{{-2}}\\cdot$g")


                    box_content = (f"$\\mathcal{{P}} = A e^{{-B E_{{ion}} / r_\\ell \\rho_\\ell}}$\n"
                                   f"A = 0.13 MeV$^{{-1}}$\n"    f"B = 8.75 GeV$^{{-1}}$cm$^{{-2}}$g" )
                    # ax_ij.legend(loc='lower left', fontsize=16, title="Source and \nTemperature", title_fontsize=16,frameon=False)
                    ax_ij.text(0.60, 0.95, box_content,
                               transform=ax_ij.transAxes,
                               fontsize=16,
                               color='black',  # White text color
                               verticalalignment='top',
                               horizontalalignment='left',
                               linespacing=1.4,  # Extra padding between lines
                               bbox=dict(
                                   facecolor='white',  # Black background
                                   edgecolor='white',  # No border outline
                                   alpha=0.9  # Slight transparency so gridlines don't completely disappear
                               ))


                # bbox = dict(boxstyle='round', facecolor='whitesmoke', alpha=0.85, edgecolor='lightgray')
                # ax_ij.legend(loc='lower left', fontsize=16, title="Source and \nTemperature", title_fontsize=16,frameon=False)

        # plt.show()
        plt.savefig(self.plot_path + "gamma_rejection_PSN.pdf")
        #
        plt.clf()
        self.Qseitz_compound_xe_plot_PSN()
    def Qseitz_compound_xe_plot_PSN(self):
        fig, ax = plt.subplots(1, 2, figsize=(16, 6))
        [result_Q_scatter, result_Q_keV, result_Q_xe, result_Eion_scatter, result_Eion_keV, result_Eion_xe,
         result_Q2_xe, result_Q_rate] = self.fitting_gamma_rejection_v2(self.fitting_df)

        # # for Rate Cs and Co, fit individually, only get last component
        self.fitting_df_Cs = pd.concat(self.Cs_fitting_list, ignore_index=True)
        self.fitting_df_Co = pd.concat(self.Co_fitting_list, ignore_index=True)
        result_Cs_fitting = self.fitting_gamma_rejection_v2(self.fitting_df_Cs)[-1]
        result_Co_fitting = self.fitting_gamma_rejection_v2(self.fitting_df_Co)[-1]

        ax[0].errorbar(self.df_Cs_116_plot['Q_rl-1_rhol-1 [GeVcm**2 g-1]'],
                       self.df_Cs_116_plot["Rejection Rate Xenon Abs[]"],
                       yerr=self.df_Cs_116_plot["Rejection Sigma Xenon Abs[]"], label="Cs 116K", fmt='o',markersize=8,color = "r")

        # ax[1].errorbar(self.df_Co_116_plot['Q_rl-1_rhol-1 [GeVcm**2 g-1]'], self.df_Co_116_plot["Rejection Rate KeV[/keV]"],
        #                   yerr=self.df_Co_116_plot["Rejection Sigma KeV[/keV]"], label="Co 116K", fmt='o')

        ax[0].errorbar(self.df_Cs_119_plot['Q_rl-1_rhol-1 [GeVcm**2 g-1]'],
                       self.df_Cs_119_plot["Rejection Rate Xenon Abs[]"],
                       yerr=self.df_Cs_119_plot["Rejection Sigma Xenon Abs[]"], label="Cs 119K", fmt='^',markersize=8,color = "orange")

        # ax[1].errorbar(self.df_Co_116_plot['Q_rl-1_rhol-1 [GeVcm**2 g-1]'], self.df_Co_116_plot["Rejection Rate KeV[/keV]"],
        #                yerr=self.df_Co_116_plot["Rejection Sigma KeV[/keV]"], label="Co 116K", fmt='o')

        ax[0].errorbar(self.df_Co_116_plot['Q_rl-1_rhol-1 [GeVcm**2 g-1]'],
                       self.df_Co_116_plot["Rejection Rate Xenon Abs[]"],
                       yerr=self.df_Co_116_plot["Rejection Sigma Xenon Abs[]"], label="Co 116K", fmt='D',markersize=8,color = "green")

        ax[0].errorbar(self.df_Co_119_plot['Q_rl-1_rhol-1 [GeVcm**2 g-1]'],
                       self.df_Co_119_plot["Rejection Rate Xenon Abs[]"],
                       yerr=self.df_Co_119_plot["Rejection Sigma Xenon Abs[]"], label="Co 119K", fmt='s',markersize=8,color = "blue")

        ax[0].errorbar(self.df_Hot_Cs_116_plot['Q_rl-1_rhol-1 [GeVcm**2 g-1]'],
                       self.df_Hot_Cs_116_plot["Rejection Rate Xenon Abs[]"],
                       yerr=self.df_Hot_Cs_116_plot["Rejection Sigma Xenon Abs[]"], label="Hot Cs 116K", fmt='s', markersize=8,
                       color="purple")

        ax[0].errorbar(self.df_Th_116_plot['Q_rl-1_rhol-1 [GeVcm**2 g-1]'],
                       self.df_Th_116_plot["Rejection Rate Xenon Abs[]"],
                       yerr=self.df_Th_116_plot["Rejection Sigma Xenon Abs[]"], label="Th 116K", fmt='D',
                       markersize=8,
                       color="brown")

        # ax[0].errorbar(self.df_Ba_116_plot['Q_rl-1_rhol-1 [GeVcm**2 g-1]'],
        #                self.df_Ba_116_plot["Rejection Rate Xenon Abs[]"],
        #                yerr=self.df_Ba_116_plot["Rejection Sigma Xenon Abs[]"], label="Ba 116K", fmt='o')

        # ax[0].plot(self.df_Ba_116_plot['Q_rl-1_rhol-1 [GeVcm**2 g-1]'], self.df_Ba_116_plot["Rejection Rate Xenon Abs[]"],
        #            label="Ba 116K 95% CL \nUpper Limit", marker='v', linestyle='None')
        a_val = result_Q2_xe[0]
        b_val = result_Q2_xe[1]
        a_str = self.fmt_sci_tex(a_val)
        b_str = self.fmt_sci_tex(b_val)

        # 2. Build string with LaTeX formatting for numbers and upright unit powers
        box_content = (
            rf"$A = {a_str}\ \mathrm{{K\text{{-}}phot}}^{{-1}}$" "\n"
            rf"$B = {b_str}\ \mathrm{{GeV}}^{{-1}}\cdot\mathrm{{cm}}^{{-2}}\cdot\mathrm{{g}}$"
        )

        ax[0].plot(result_Q2_xe[2], result_Q2_xe[3],
                   color="black")



        label_text = f"A = {a_val:.2e},\nB = {b_val:.2e}"


        ax[0].set_xlabel(r"$Q_{Seitz} r_l^{-1} \rho_l^{-1}$ [GeV cm$^2$ g$^{-1}$]",fontsize=16)
        ax[0].set_ylabel("Nucleation probability\n(per xenon K shell photoabsorption)",fontsize=16)
        ax[0].set_yscale("log")
        ax[0].legend(loc='lower left', fontsize=16,title="Source and \nTemperature",title_fontsize=16,frameon=False)
        # box_content = (    f"$y = A e^{{-Bx}}$\n" f"A = ${a_str}$ K-phot$^{{-1}}$\n"    f"B = ${b_str}$ GeV$^{{-1}}\\cdot$cm$^{{-2}}\\cdot$g")
        print("value", a_str, b_str)
        box_content = ( f"$\\mathcal{{P}} = A e^{{-B Q_{{Seitz}} / r_\\ell \\rho_\\ell}}$\n" f"A = 0.33 K-phot$^{{-1}}$\n"    f"B = 3.76 GeV$^{{-1}}$cm$^{{-2}}$g")

        ax[0].text(0.55, 0.95,box_content,
        transform=ax[0].transAxes,
        fontsize=16,
        color='black',                  # White text color
        verticalalignment='top',
        horizontalalignment='left',
        linespacing=1.4,                 # Extra padding between lines
        bbox=dict(
            facecolor='white',          # Black background
            edgecolor='white',           # No border outline
            alpha=0.9                   # Slight transparency so gridlines don't completely disappear
        ))


        # plt.show()

        plt.savefig(self.plot_path + "Qseitz_compound_xe_PSN.pdf")

        plt.clf()

    def fmt_sci_tex(self, val):
        """Converts a float to LaTeX scientific notation (e.g., 3.25 \times 10^{-1})."""
        if val == 0:
            return r"0"
        s = f"{val:.2e}"
        base, exp = s.split('e')
        return f"{base} \\times 10^{{{int(exp)}}}"
    def plot_spectrum(self):

        self.sim_list = [self.Cs_sims, self.Co_sims, self.Ba_sims]
        self.sim_doped_list = [self.Cs_sims_doped, self.Co_sims_doped, self.Ba_sims_doped]
        self.sim_tag = [r"$^{137}$Cs", r"$^{60}$Co", r"$^{133}$Ba"]
        fig, ax = plt.subplots(1, 3, figsize=(16, 5))

        for i in range(len(self.sim_list)):
            sim_result = self.sim_list[i]
            sim_doped_result = self.sim_doped_list[i]
            Rate_factor = sim_result[0]
            energy_edges_all = sim_result[1][0][1]
            energy_counts_all = sim_result[1][0][0]
            energy_edges_primary = sim_result[4][0][1]
            energy_counts_primary = sim_result[4][0][0]
            counts_energy_cum_bin_primary = sim_result[6]


            Rate_factor_doped = sim_doped_result[0]
            energy_edges_doped = sim_doped_result[1][0][1]
            energy_counts_doped = sim_doped_result[1][0][0]
            counts_cum_bin_doped = sim_doped_result[2]
            counts_energy_cum_bin_doped = sim_doped_result[3]

            edges_p5, counts_p5 = self.rebin_to_5kev(
                energy_edges_all, energy_counts_all
            )
            ax[0].step(
                edges_p5[:-1],
                Rate_factor * counts_p5,
                where="post",
                label=self.sim_tag[i],
            )

            # --- Ax[1]: Differential Energy from Cumulative Array (Re-binned) ---
            # Convert cumulative to differential spectrum
            differential_counts_primary = -np.diff(counts_energy_cum_bin_primary)
            # The diff array is shorter by 1 element, pad it or match to the edges
            edges_all_cut = energy_edges_primary[:-1]

            edges_all5, counts_all5 = self.rebin_to_5kev(
                edges_all_cut, differential_counts_primary
            )
            ax[1].step(
                edges_all5[:-1],
                Rate_factor * counts_all5,
                where="post",
                label=self.sim_tag[i],
            )

            # --- Ax[2]: Xenon-Doped Target Interactions (Re-binned) ---
            edges_d5, counts_d5 = self.rebin_to_5kev(
                energy_edges_doped, energy_counts_doped
            )
            ax[2].step(
                edges_d5[:-1],
                Rate_factor_doped * counts_d5,
                where="post",
                label=self.sim_tag[i],
            )

        # --- Polish Axis 0 ---
        # --- Polish Axis 0 ---
        ax[0].set_xlabel("Deposited Energy $E$ [keV]", fontsize=11)
        ax[0].set_ylabel("Event Rate [mHz / 5 keV]", fontsize=11)
        ax[0].set_title("Gamma Event Rate Spectrum\n(All Interaction Vertices)")
        ax[0].set_xlim(0, 1400)
        ax[0].set_yscale("log")
        ax[0].legend(frameon=True)

        # --- Polish Axis 1 ---
        ax[1].set_xlabel("Deposited Energy $E$ [keV]", fontsize=11)
        # Expressing it as energy-weighted rate makes the math clear to the reader
        ax[1].set_ylabel("Energy-Weighted Rate [keV$\cdot$mHz / 5 keV]", fontsize=11)
        ax[1].set_title("Energy Deposition Rate\n(Primary Gamma Vertices)")
        ax[1].set_xlim(0, 1400)
        ax[1].set_yscale("log")
        ax[1].legend(frameon=True)

        # --- Polish Axis 2 ---
        ax[2].set_xlabel("Deposited Energy [keV]", fontsize=11)
        ax[2].set_ylabel("Event Rate [mHz / 5 keV]", fontsize=11)
        ax[2].set_title(" Events Rate Spectrum \n(Tagged by Xenon K Shell Vacancy)")
        ax[2].set_xlim(0, 1400)
        ax[2].set_yscale("log")
        ax[2].legend(frameon=True)

        plt.savefig(self.plot_path +"spectrum_comparision_updated.pdf")

    def rebin_to_5kev(self, original_edges, original_counts, target_bin_width=5.0):
        """Aggregates arbitrary fine bins into uniform 5 keV bins."""
        max_energy = original_edges[-1]
        # Create the new uniform 5 keV bin grid
        new_edges = np.arange(0, max_energy + target_bin_width, target_bin_width)

        # Use the centers of the original bins as the data points
        bin_centers = (original_edges[:-1] + original_edges[1:]) / 2.0

        # Re-histogram the counts into the new 5 keV bins
        rebinned_counts, _ = np.histogram(
            bin_centers, bins=new_edges, weights=original_counts
        )
        return new_edges, rebinned_counts

    def concat_PT_condition(self, df):

        df_combined = df.groupby(['Pressure [bara]', 'Seitz [keV]', 'Eion_rl-1_rhol-1 [GeVcm**2 g-1]', 'Q_rl-1_rhol-1 [GeVcm**2 g-1]' , 'Eion [keV]'], as_index=False).agg({
            "Clean Rate [mHz]": 'mean',
            'Clean Rate Sigma [mHz]': lambda x: np.sqrt(np.sum(x ** 2)),
            "Rejection Rate Scattering[]": 'mean',
            "Rejection Sigma Scattering[]": lambda x: np.sqrt(np.sum(x**2)),
            "Rejection Rate KeV[/keV]": 'mean',
            "Rejection Sigma KeV[/keV]": lambda x: np.sqrt(np.sum(x ** 2)),
            "Rejection Rate Xenon Abs[]": 'mean',
            "Rejection Sigma Xenon Abs[]": lambda x: np.sqrt(np.sum(x ** 2)),
        })
        return df_combined

    def doped_gamma_rejection_plot(self):
        fig, ax = plt.subplots(2, 1, figsize=(8,14))
        # fig, ax = plt.subplots(2, 1, figsize=(6, 10))
        self.fitting_list = []

        self.Cs_label = ["Cs 11/17/2025 116K", "Cs 12/01/2025 116K", "Cs 12/10/2025 116K", "Cs 01/20/2026 116K",
                         "Cs 02/02/2026 119K"]
        self.Co_label = ["Co 12/15/2026 116K"]
        for i in range(len(self.Cs_exp_rejection_path)):
            df = pd.read_csv(self.Cs_exp_rejection_path[i])
            # print(df.columns)
            # doc_label = self.Cs_exp_raw_path[i].replace('_exposures', '')
            doc_label = self.Cs_label[i]
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

            ax[1].errorbar(df['Eion_rl-1_rhol-1 [GeVcm**2 g-1]'], df["Rejection Rate Scattering[]"],
                           yerr=df["Rejection Sigma Scattering[]"], label=doc_label, fmt='o')

        for i in range(len(self.Co_exp_rejection_path)):
            df = pd.read_csv(self.Co_exp_rejection_path[i])
            # print(df.columns)
            # doc_label = self.Co_exp_raw_path[i].rstrip("_exposures")
            doc_label = self.Co_label[i]
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


            ax[1].errorbar(df['Eion_rl-1_rhol-1 [GeVcm**2 g-1]'], df["Rejection Rate Scattering[]"],
                           yerr=df["Rejection Sigma Scattering[]"], label=doc_label, fmt='o')


        self.fitting_df =  pd.concat(self.fitting_list, ignore_index=True)
        [(a_fit_scatter, b_fit_scatter,x_fitted_scatter,y_fitted_scatter),(a_fit_keV, b_fit_keV,x_fitted_keV,y_fitted_keV)] = self.fitting_doped_gamma_rejection()

        # plot the fitting function
        # ax[0].plot(x_fitted_scatter,y_fitted_scatter,label = f"a,b = {a_fit_scatter:.2e} , {b_fit_scatter:.2e}", color="black")
        ax[0].plot(x_fitted_scatter, y_fitted_scatter,
                   color="black")

        #gamma rejection up limit
        # self.bkg_floor_plot(ax[0],"Seitz")

        ax[0].set_xlabel(r"Seitz threshold [keV]")
        ax[0].set_ylabel("Nucleation probability (per xenon photoabsorption)")
        # ax[0].set_title("Gamma Rejection Per Scattering ")
        # ax[0].set_ylim(1.0e-12,1.0e-2)
        # ax[0].set_xlim(0,6)
        # ax[0].set_xlim(0.8,1.5)
        # ax[0].set_ylim(1.0e-12,1.0e-2)
        ax[0].set_yscale("log")

        ax[0].legend(loc='lower left', fontsize=7)

        # ax[1].plot(x_fitted_keV, y_fitted_keV, label=f"a,b = {a_fit_keV:.2e} , {b_fit_keV:.2e}", color="black")
        ax[1].plot(x_fitted_keV, y_fitted_keV, color="black")

        # self.bkg_floor_plot(ax[1], "Eion")
        ax[1].set_xlabel(r"$E_{ion} r_l^{-1} \rho_l^{-1}$ [GeV cm$^2$ g$^{-1}$]")
        ax[1].set_ylabel("Nucleation probability (per xenon photoabsorption)")
        # ax[1].set_title("Gamma Rejection Per keV ")
        # ax[1].set_ylim(1.0e-14,1.0e-4)
        # ax[1].set_xlim(0.08,0.15)
        # ax[1].set_xlim(0.8,1.1)
        ax[1].set_yscale("log")
        ax[1].legend(loc='lower left', fontsize=7)

        plt.savefig(self.plot_path + "gamma_rejection_doped.pdf")
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
    def fitting_gamma_rejection_v2(self, dataframe):
        # switch Y axis. Now Q vs per kev and Eion vs per interaction
        x_Q = dataframe["Seitz [keV]"].values
        y_per_scatter = dataframe["Rejection Rate Scattering[]"].values
        # dealing with guess
        x_min_Q= min(x_Q)
        x_max_Q = max(x_Q)
        y_min_per_scattering = min(y_per_scatter)
        y_max_per_scattering = max(y_per_scatter)

        x_Eion = dataframe["Eion_rl-1_rhol-1 [GeVcm**2 g-1]"].values
        y_per_keV = dataframe["Rejection Rate KeV[/keV]"].values
        # dealing with guess
        x_min_Eion = min(x_Eion)
        x_max_Eion = max(x_Eion)
        y_min_per_keV = min(y_per_keV)
        y_max_per_keV = max(y_per_keV)

        x_Q2 = dataframe["Q_rl-1_rhol-1 [GeVcm**2 g-1]"].values
        y_per_xe = dataframe["Rejection Rate Xenon Abs[]"].values
        # dealing with guess
        x_min_Q2 = min(x_Q2)
        x_max_Q2 = max(x_Q2)
        y_min_per_xe = min(y_per_xe)
        y_max_per_xe = max(y_per_xe)

        y_rate = dataframe["Clean Rate [mHz]"].values
        y_min_rate = min(y_rate)
        y_max_rate = max(y_rate)




        result_Q_scatter  = self.fit_combination(x_Q,y_per_scatter,y_max_per_scattering,y_min_per_scattering,x_max_Q,x_min_Q)

        result_Q_keV = self.fit_combination(x_Q, y_per_keV, y_max_per_keV, y_min_per_keV, x_max_Q,
                                                x_min_Q)
        result_Q_xe = self.fit_combination(x_Q, y_per_xe, y_max_per_xe, y_min_per_xe, x_max_Q,
                                            x_min_Q)

        result_Eion_scatter = self.fit_combination(x_Eion, y_per_scatter, y_max_per_scattering, y_min_per_scattering, x_max_Eion,
                                                x_min_Eion)
        result_Eion_keV = self.fit_combination(x_Eion, y_per_keV, y_max_per_keV, y_min_per_keV, x_max_Eion,
                                            x_min_Eion)
        result_Eion_xe = self.fit_combination(x_Eion, y_per_xe, y_max_per_xe, y_min_per_xe, x_max_Eion,
                                           x_min_Eion)
        
        result_Q2_xe = self.fit_combination(x_Q2, y_per_xe, y_max_per_xe, y_min_per_xe, x_max_Q2,
                                           x_min_Q2)
        result_Q_rate  = self.fit_combination(x_Q, y_rate, y_max_rate, y_min_rate, x_max_Q,
                                           x_min_Q)
        
        
        return [result_Q_scatter,result_Q_keV,result_Q_xe,result_Eion_scatter,result_Eion_keV,result_Eion_xe,result_Q2_xe,result_Q_rate ]

        # # fit 
        # b_guess_per_scattering=-(np.log(y_max_per_scattering)-np.log(y_min_per_scattering))/(x_max_Q-x_min_Q)
        # a_guess_scattering = y_max_per_scattering
        # initial_guess_scatter = [a_guess_scattering, b_guess_per_scattering]
        # popt_scatter, pcov_scatter = curve_fit(self.exp_func, x_Q, y_per_scatter, p0=initial_guess_scatter)
        # a_fit_scatter, b_fit_scatter= popt_scatter
        # print('a_fit_scatter, b_fit_scatter',a_fit_scatter, b_fit_scatter)
        # x_fitted_scatter = np.linspace(min(x_Q), max(x_Q), 100)
        # y_fitted_scatter = self.exp_func(x_fitted_scatter, *popt_scatter)
        # 
        # 
        # 
        # 
        # 
        # 
        # # b is negative
        # b_guess_per_keV = -(np.log(y_max_per_keV) - np.log(y_min_per_keV)) / (
        #             x_max_Eion - x_min_Eion)
        # a_guess_scattering = y_max_per_keV
        # initial_guess_keV = [a_guess_scattering, b_guess_per_keV]
        # popt_keV, pcov_keV = curve_fit(self.exp_func, x_Eion, y_per_keV, p0=initial_guess_keV)
        # a_fit_keV, b_fit_keV = popt_keV
        # print('a_fit_keV, b_fit_keV', a_fit_keV, b_fit_keV)
        # x_fitted_keV = np.linspace(min(x_Eion), max(x_Eion), 100)
        # y_fitted_keV = self.exp_func(x_fitted_keV, *popt_keV)
        # 
        # 
        # # b is negative
        # b_guess_per_xe = -(np.log(y_max_per_xe) - np.log(y_min_per_xe)) / (
        #             x_max_Q2 - x_min_Q2)
        # a_guess_xe = y_max_per_xe
        # initial_guess_xe = [a_guess_xe, b_guess_per_xe]
        # popt_xe, pcov_xe = curve_fit(self.exp_func, x_Q2, y_per_xe, p0=initial_guess_xe)
        # a_fit_xe, b_fit_xe = popt_xe
        # print('a_fit_xe, b_fit_xe', a_fit_xe, b_fit_xe)
        # x_fitted_xe = np.linspace(min(x_Q2), max(x_Q2), 100)
        # y_fitted_xe = self.exp_func(x_fitted_xe, *popt_xe)

    def fitting_gamma_rejection_v3(self, dataframe , x_cfg,y_cfg):
        # switch Y axis. Now Q vs per kev and Eion vs per interaction
        result = []
        for i in range(4):
            row_result = []
            for j in range(3):
                x = dataframe[[x_cfg[j]["x"]]].values.flatten()
                y = dataframe[[y_cfg[i]["y"]]].values.flatten()
                print("shape",type(x),x.shape)
                # dealing with guess
                x_min = min(x)
                x_max = max(x)
                y_min = min(y)
                y_max = max(y)
                fit_output = self.fit_combination(x, y, y_max, y_min, x_max, x_min)
                row_result.append(fit_output)
            result.append(row_result)

        return result

    def fitting_doped_gamma_rejection(self):
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
        # for doped, y should be still per interaction, I keep the value is per interaction but the variable name
        # is per keV.
        y_per_keV = self.fitting_df["Rejection Rate Scattering[]"].values
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

    def fit_combination(self,x,y,y_max,y_min,x_max,x_min):
        b_guess_per_scattering = -(np.log(y_max) - np.log(y_min)) / (x_max - x_min)
        a_guess_scattering = (y_max+y_min)/2
        initial_guess_scatter = [a_guess_scattering, b_guess_per_scattering]
        popt_scatter, pcov_scatter = curve_fit(self.exp_func, x, y, p0=initial_guess_scatter)
        a_fit_scatter, b_fit_scatter = popt_scatter
        print('a_fit_scatter, b_fit_scatter', a_fit_scatter, b_fit_scatter)
        x_fitted_scatter = np.linspace(min(x), max(x), 100)
        y_fitted_scatter = self.exp_func(x_fitted_scatter, *popt_scatter)
        return (a_fit_scatter, b_fit_scatter, x_fitted_scatter, y_fitted_scatter)
        
    def calculate_rss(self, series):
        """Calculates sqrt(a^2 + b^2 + ...)"""
        return np.sqrt(np.sum(series ** 2))

    def exp_func(self, x, a, b):
        return a * np.exp(-b * x)

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





if __name__=="__main__":
    # IA =  integrated_analysis(volume="dome")
    IA = integrated_analysis(volume="bulk")
    # test = test_csv()