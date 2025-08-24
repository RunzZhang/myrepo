"""This file is for analyze the root files"""
"""['Transportation', 'hadElastic', 'nCapture', 'conv', 'phot', 'compt', 'Rayl', 'neutronInelastic', 'ionIoni', 
        'NoProcess', 'RadioactiveDecay', 'photonNuclear', 'Decay']
            ['neutron', 'gamma', 'Ar40', 'Ar36', 'Ar41', 'Ar37']
            ['world_phys', 'Vacuum_vessel_phys', 'Inside_vacuum_vessel_phys', 'pressure_vessel_phys', 
            'hydraulic_fluid_phys', 'HDPE_pressure_vessel_phys', 'reflector_Cu_phys', 'reflector_PTFE_phys', 
            'outer_jar_phys', 'SiPM_Holder21_phys', 'reflector_top_PTFE_phys', 'reflector_top_Cu_phys', 'LAr_phys', 
            'Camera_port_phys', 'inner_jar_phys', 'RTD_Cable_Kapton_phys', 'RTD_Cable_Cu_phys', 'Camera_System2_phys', 
            ' Iris_Holder_phys', ' Iris_phys', 'Top_Plastic_phys', 'Lens_Holder3_phys', 'Camera_System1_phys', 
            ' Iris_Holder1_phys', 'Lens_Holder4_phys', 'Side_support3_phys', 'Top_sf_phys', 'SiPM_Holder31_phys',
             'SiPM_PCB_Out42_phys', 'SiPM_Holder42_phys', 'SiPM_PCB_Inn4_phys', 'SiPM4_Inside_phys', 'SiPM_Holder51_phys', 
             'calibration_port_phys', 'calibration_Be_phys', 'calibration_air_phys', 'RTD_Cable_Cu_1_phys', 'SiPM_PCB_Out31_phys', 
             'SiPM_PCB_Inn3_phys', 'SIPM3_Si_phys', 'top_flange_phys', 'SiPM_Holder18_phys', 'SiPM1_Inside_phys', 'SIPM1_Si_phys', 
             'OJ_Spacer2_phys', 'SiPM_Holder17_phys', 'Adjustment1_phys', 'SiPM_PCB_Out38_phys', 'SiPM_Holder38_phys', 
             'SiPM3_Inside_phys', 'Hyspan_bellow_phys', 'Bottom_Flange_phys', 'Bellows_weldment_phys', ' Sensor_Holder_phys', 
             'SiPM_Holder12_phys', 'SiPM2_Inside_phys', 'Side_support1_phys', 'OJ_Spacer17_phys', 'SiPM_PCB_Out22_phys', 
             'SiPM_Holder22_phys', 'SiPM_PCB_Inn2_phys', 'SIPM2_Si_phys', 'OJ_Spacer16_phys', 'Aspheric_lens1_phys', 
             'Camera_port3_phys', 'Piezo_Cu_2_phys', 'Vacio_phys', 'Piezo_phys', 'RTD_block_phys', 'RTD_base_phys', 
             'RTD_mother2_phys', 'SiPM_Holder32_phys', 'SiPM_PCB_Out32_phys', 'PV_spool_phys', 'Base_SF_phys',
              'Piezo_Cu_1_phys', 'SiPM_Holder58_phys', 'Copper_Rod1_phys', 'PBC_phys', 'OJ_Spacer3_phys', 
              'RTD_connector_feedthru_phys', 'SiPM_Holder11_phys', 'SiPM_PCB_Inn1_phys', 'SiPM_PCB_Out11_phys', 
              'Camera_port2_phys', 'Holder1_VV_phys', 'SiPM_Holder52_phys', 'OJ_Spacer8_phys', 'Guide_rod_flange_phys', 
              'Guide_Rod2_phys', 'SiPM_PCB_Out58_phys', 'OJ_Spacer7_phys', 'SiPM_Holder41_phys', 'Sapphire1_Ssteal_phys',
               'Piezo_Cu_8_phys', 'Guide_Rod1_phys', 'SiPM_PCB_Out21_phys', 'Top_Plastic1_phys', 'Plastic_Flange_phys', 
               'Iris_brass_phys', 'Aspheric_lens_phys', 'Adjustment_phys', 'SS_Rod1_phys', 'Lens_Holder2_phys', 'Lens1_phys',
                'Sapphire1_phys', 'SiPM_Holder35_phys', 'OJ_Spacer4_phys', 'Bottom_Spacer5_phys', 'OJ_Spacer5_phys', 
                'Camera_PCB1_phys', 'Camera1_phys', 'RTD_connector_1_4_phys', 'Bottom_Spacer9_phys', 'Sapphire_Ssteal_phys',
                 'Lens_Holder_phys', 'Camera_System_phys', 'SiPM_Holder48_phys', 'Bottom_Spacer11_phys', 'Piezo_Cu_5_phys', 
                 'S_Epoxy_phys', 'RTD_mother1_phys', 'RTD_wire_clamp_phys', 'SiPM_Holder23_phys', 'SiPM_Holder33_phys', 
                 'RTD_mother3_phys', 'SiPM_PCB_Out41_phys', 'OJ_Spacer6_phys', 'OJ_Spacer1_phys', 'SiPM_Holder27_phys', 
                 'SiPM5_Inside_phys', 'SiPM_PCB_Out36_phys', 'SiPM_Holder36_phys', 'OJ_Spacer9_phys', 'SiPM_PCB_Out28_phys',
                  'SiPM_Holder28_phys', 'OJ_Spacer10_phys', 'SIPM4_Si_phys', 'SiPM_Holder37_phys', 'OJ_Spacer18_phys', 
                  'SiPM_Holder14_phys', 'SiPM_PCB_Out13_phys', 'SiPM_Holder13_phys', 'Piezo_Cu_7_phys', 'SiPM_PCB_Out18_phys',
                   'SiPM_PCB_Out51_phys', 'SiPM_PCB_Inn5_phys', 'SIPM5_Si_phys', 'SiPM_Holder16_phys', 'SiPM_PCB_Out16_phys', 
                   'SiPM_PCB_Out44_phys', 'SiPM_Holder44_phys', 'SS_Rod2_phys', 'SiPM_Holder15_phys', 'Piezo_Cu_6_phys',
                    'OJ_Spacer12_phys', 'Bottom_Spacer3_phys', 'Holder2_VV_phys', 'Sapphire_phys', 'Alignment_Ring_phys',
                     'SiPM_Holder43_phys', 'Bottom_Spacer6_phys', 'Bottom_Spacer7_phys', 'Piezo_Cu_4_phys', 'Iris_blades_phys',
                      'SiPM_PCB_Out33_phys', 'Copper_Rod2_phys', 'SiPM_PCB_Out43_phys', 'Piezo_Cu_3_phys', 'SiPM_Holder54_phys',
                       'SiPM_Holder55_phys', 'Side_support2_phys', 'Bottom_Spacer2_phys', 'SiPM_Holder53_phys', 'SiPM_Holder57_phys', 
                       'SiPM_PCB_Out57_phys', 'SiPM_PCB_Out35_phys', 'SiPM_Holder45_phys', 'SiPM_Holder56_phys', 'OJ_Spacer15_phys',
                        'OJ_Spacer14_phys', 'SiPM_PCB_Out52_phys', 'SiPM_Holder46_phys', 'SiPM_PCB_Out56_phys', 
                        'Front_Nanoguide_Holder_phys', 'SiPM_Holder47_phys', 'SiPM_Holder34_phys', 'Bottom_Spacer10_phys',
                         'Bottom_Spacer14_phys', 'Bottom_Spacer13_phys', 'Bottom_Spacer12_phys', 'SiPM_PCB_Out14_phys', 
                         'SiPM_PCB_Out12_phys', 'SiPM_PCB_Out47_phys', 'Camera_spring2_phys', 'Bottom_Spacer8_phys', 
                         'SiPM_PCB_Out53_phys', 'SiPM_PCB_Out23_phys', 'Lens_phys', 'OJ_Spacer11_phys', 'SiPM_Holder24_phys', 
                         'Bottom_Spacer17_phys', 'Bottom_Spacer18_phys', 'Holder3_VV_phys', 'SiPM_PCB_Out17_phys', 
                         'SiPM_PCB_Out24_phys', 'Bottom_Spacer1_phys', 'SiPM_PCB_Out48_phys', 'Guide_Rod3_phys',
                          'Bottom_Spacer15_phys', 'Bottom_Spacer4_phys', 'SiPM_PCB_Out37_phys', 'SiPM_PCB_Out15_phys', 
                          'SiPM_Holder26_phys', 'SiPM_PCB_Out26_phys', 'SiPM_Holder25_phys', 'SiPM_PCB_Out25_phys', 
                          'SiPM_PCB_Out45_phys', 'SiPM_PCB_Out27_phys', 'SiPM_PCB_Out34_phys', 'OJ_Spacer13_phys', 
                          'Peak_Rod2_phys', 'Rear_Nanoguide_Holder_phys', 'Screw4_phys', 'SiPM_PCB_Out46_phys', 
                          'Sensor_Covere_phys', 'Bottom_Spacer16_phys', 'SiPM_PCB_Out55_phys', 'Screw1_phys', 'Nanoguide_phys',
                           'SiPM_PCB_Out54_phys', 'Camera_spring1_phys', 'Sensor_Plate_phys', 'Peak_Rod1_phys', 'Screw3_phys', 
                           'Screw2_phys']"""
import pandas as pd
import uproot
import matplotlib.pyplot as plt
import numpy as np
import csv
# filename = "/data/runzezhang/Geant4Simulaions/g411_TN/dmx.root"

class RestructureRoot():
    def __init__(self):
        self.filepath = "/data/runzezhang/result/TN_sims_D/dmx_Cf_1E7.root"
        self.reconstruct_filepath = "/data/runzezhang/result/TN_sims_D/dmx_rcCf_1E7.csv"
        self.file = uproot.open(self.filepath)["tree"]
        print("columns: ",self.file.keys())
        #['Event', 'name', 'Parent ID', 'Track ID', 'Step ID', 'X/mm', 'Y/mm', 'Z/mm', 'Kinetic/keV', 'Recoiled/keV', 'Volume', 'Process']
        self.selected_columns = ["Event","name","Parent ID","Track ID","Step ID","X/mm","Kinetic/keV","Recoiled/keV","Volume","Process"]
        self.rows = 1000
        # self.df = self.file.arrays(self.selected_columns, library="pd").head(self.rows)
        self.df = self.file.arrays(self.selected_columns, library="pd")

        self.reconstruct()

    def reconstruct(self):
        event_number = self.df[:]["Event"].to_list()
        print(event_number[:100])
        started_point = 0
        temp_point = 1
        # find event number 1's index
        for i in range(0,len(event_number)):
            if event_number[i]==1:
                started_point = i
        print("start",started_point)
         # then if there is 0 in the event number, replace it with last none-zero event number
        for i in range(started_point, len(event_number)):
            if event_number[i]!=0:
                temp_point = event_number[i]
            else:
                event_number[i] = temp_point
        print("end",event_number[:100])

        # put the updated event_number back to data frame
        self.df.update(pd.DataFrame({'Event':event_number}))
        print(self.df[53:60])
        self.df.to_csv(self.reconstruct_filepath, sep=',', index=False, encoding='utf-8')
        # uproot.writing._dask_write.dask_write(self.df, self.reconstruct_filepath/)

class ReadRoot():
    def __init__(self):
        self.base_path = "/data/runzezhang/result/TN_box/"
        self.base_path2 = "/data/runzezhang/result/TN_box/"
        self.plot_path = '/data/runzezhang/result/TN_box/plot/'


        self.single_run()

        # self.multi_run_loop()

    def multi_run_loop(self):
        self.process = []
        self.ene = []
        self.cross_number = []
        for i in range(0,11):
            self.multi_run(i)
        print(self.process)
        self.plot_ncrystal_test()
    def multi_run(self, i):

        self.false_1 = "Cf_1E6_"+str(i)+"N_false1.csv"
        self.false_2 = "Cf_1E6_"+str(i)+"N_false2.csv"
        self.false_3 = "Cf_1E6_"+str(i)+"N_ini_false3.csv"
        self.signal = "Cf_1E6_"+str(i)+"N_sig.csv"
        self.false_1_mid = "Cf_1E6_"+str(i)+"N_false1_mid.csv"
        self.false_2_mid = "Cf_1E6_"+str(i)+"N_false2_mid.csv"
        self.false_3_mid = "Cf_1E6_"+str(i)+"N_ini_false3_mid.csv"
        self.signal_mid = "Cf_1E6_"+str(i)+"N_sig_mid.csv"
        self.false_1_path = self.base_path + self.false_1
        self.false_2_path = self.base_path + self.false_2
        self.false_3_path = self.base_path + self.false_3
        self.false_1_path_mid = self.base_path + self.false_1_mid
        self.false_2_path_mid = self.base_path + self.false_2_mid
        self.false_3_path_mid = self.base_path + self.false_3_mid
        self.signal_path_mid = self.base_path + self.signal_mid
        self.signal_path = self.base_path + self.signal
        # self.filepath = self.base_path +"dmx_lr.root"
        # self.filepath = self.base_path + "dmx_Cfneutron_Ncry_1E6.root"
        self.filepath = self.base_path +"dmx_Cfneutron_Ncry_1E6_"+str(i)+".root"

        self.file = uproot.open(self.filepath)["tree"]
        print("columns: ", self.file.keys(), len(self.file.arrays()))
        # ['Event', 'name', 'Parent ID', 'Track ID', 'Step ID', 'X/mm', 'Y/mm', 'Z/mm', 'Kinetic/keV', 'Recoiled/keV', 'Volume', 'Process']
        self.selected_columns = ["Event", "name", "Parent ID", "Track ID", "Step ID", 'X/mm', 'Y/mm', 'Z/mm', 'px/MeV',
                                 'py/MeV', 'pz/MeV', "Kinetic/keV", "Recoiled/keV", "Volume", "Process"]
        self.rows = 1000

        # self.df = self.file.arrays(self.selected_columns, library="pd").head(self.rows)
        # # process data so that it is easier to read
        # first 1000 rows
        # self.df = self.file.arrays(self.selected_columns, library="pd")
        self.df = self.file.arrays(library="pd")
        print("df", self.df.head(self.rows))
        # self.df = self.file.arrays(self.selected_columns, library="pd").head(self.rows)
        # only valid for ncrystal

        # self.df_list = []
        # for i in self.ene_dict:
        #     df_temp = self.df[""]
        #     self.modify_df_group(i)
        self.modify_df()
        # print(self.df)

        # gamma direction x -1 or +1
        # one gamma per event
        # gamma place at  [30.50, 30.52]
        # collection gamma energy
        # 1e? FISSIons
        # PDF

        # neutron spectrum on detector
        # self.neutron_momentum()
        # self.plot_neutron_momentum()

        # ncrystal test

        (n_energy, in_num, out_num)=self.ncrystal_test()
        self.ene.append(n_energy)
        self.cross_number.append(out_num/in_num)

        self.process += self.df["Process"].unique().tolist()




    def single_run(self):
        self.false_1 = "Cf_1E6_N_false1.csv"
        self.false_2 = "Cf_1E6_N_false2.csv"
        self.false_3 = "Cf_1E6_N_ini_false3.csv"
        self.signal = "Cf_1E6_N_sig.csv"
        self.false_1_mid = "Cf_1E6_N_false1_mid.csv"
        self.false_2_mid = "Cf_1E6_N_false2_mid.csv"
        self.false_3_mid = "Cf_1E6_N_ini_false3_mid.csv"
        self.signal_mid = "Cf_1E6_N_sig_mid.csv"
        self.false_1_path = self.base_path + self.false_1
        self.false_2_path = self.base_path + self.false_2
        self.false_3_path = self.base_path + self.false_3
        self.false_1_path_mid = self.base_path + self.false_1_mid
        self.false_2_path_mid = self.base_path + self.false_2_mid
        self.false_3_path_mid = self.base_path + self.false_3_mid
        self.signal_path_mid = self.base_path + self.signal_mid
        self.signal_path = self.base_path + self.signal
        # self.filepath = self.base_path +"dmx_lr.root"
        # self.filepath = self.base_path + "dmx_Cfneutron_Ncry_1E6.root"
        self.filepath = self.base_path + "dmx_Cfneutron_1E6.root"
        print(self.filepath)

        self.file = uproot.open(self.filepath)["tree"]
        print("columns: ", self.file.keys(), len(self.file.arrays()))
        # ['Event', 'name', 'Parent ID', 'Track ID', 'Step ID', 'X/mm', 'Y/mm', 'Z/mm', 'Kinetic/keV', 'Recoiled/keV', 'Volume', 'Process']
        self.selected_columns = ["Event", "name", "Parent ID", "Track ID", "Step ID", 'X/mm', 'Y/mm', 'Z/mm', 'px/MeV',
                                 'py/MeV', 'pz/MeV', "Kinetic/keV", "Recoiled/keV", "Volume", "Process"]
        self.rows = 1000

        # self.df = self.file.arrays(self.selected_columns, library="pd").head(self.rows)
        # # process data so that it is easier to read
        # first 1000 rows
        self.df = self.file.arrays(self.selected_columns, library="pd")
        print("df", self.df.head(self.rows))
        # self.df = self.file.arrays(self.selected_columns, library="pd").head(self.rows)
        # only valid for ncrystal
        # self.df_list = []
        # for i in self.ene_dict:
        #     df_temp = self.df[""]
        #     self.modify_df_group(i)
        self.modify_df()
        # print(self.df)

        # gamma direction x -1 or +1
        # one gamma per event
        # gamma place at  [30.50, 30.52]
        # collection gamma energy
        # 1e? FISSIons
        # PDF

        # neutron spectrum on detector
        self.neutron_momentum()
        self.plot_neutron_momentum()

        # ncrystal test
        # self.ncrystal_test()
        # self.plot_ncrystal_test()

        # self.gamma_event()
        # false noise 2, need to relocate directory
        # self.Huge_scatter_event()
        # signal rate, caputre in liquid argon
        # self.LAr_gamma_event()
        # single elastic scatter and capture false signal 1
        # self.single_e_n_capture_event()
        # self.FN_spectrum_v2()
        # self.plot_elastic()

        # self.find_multiplicity()
        # self.Check_inelastic()
        # test elatic and inelastic effect
        # self.bubble_rate()

    # there was some 0 in event columns, set them to corresponding value
    # for example 001002003 will be 001112223

    def reidx_event(self):
        event_number = self.df[:]["Event"].to_list()
        print(event_number[:100])
        started_point = 0
        temp_point = 1
        # find event number 1's index
        for i in range(0, len(event_number)):
            if event_number[i] == 1:
                started_point = i
        print("start", started_point)
        # then if there is 0 in the event number, replace it with last none-zero event number
        for i in range(started_point, len(event_number)):
            if event_number[i] != 0:
                temp_point = event_number[i]
            else:
                event_number[i] = temp_point
        print("end", event_number[:100], event_number[-1])

        # put the updated event_number back to data frame
        self.df.update(pd.DataFrame({'Event': event_number}))
    def string_summary(self):
        process_clean=[]
        df_process = self.df[:]["Process"].to_list()
        for element in df_process:
            if element not in process_clean:
                process_clean.append(element)
        print(process_clean)

        particle_clean = []
        df_particle = self.df[:]["name"].to_list()
        for element in df_particle:
            if element not in particle_clean:
                particle_clean.append(element)
        print(particle_clean)

        volume_clean = []
        df_volume = self.df[:]["Volume"].to_list()
        for element in df_volume:
            if element not in volume_clean:
                volume_clean.append(element)
        print(volume_clean)
    def modify_df(self):
        # change column property. Mainly this change awkuard into str
        print(self.df.dtypes)
        self.df['name'] = self.df['name'].astype(str)
        self.df['Volume'] = self.df['Volume'].astype(str)
        self.df['Process'] = self.df['Process'].astype(str)
        # this make event number correct
        self.reidx_event()

    def reidx_event_group(self,ene=0):
        event_number = self.df[:]["Event"].to_list()
        print(event_number[:100])
        started_point = 0
        temp_point = 1
        # find event number 1's index
        for i in range(0, len(event_number)):
            if event_number[i] == 1:
                started_point = i
        print("start", started_point)
        # then if there is 0 in the event number, replace it with last none-zero event number
        for i in range(started_point, len(event_number)):
            if event_number[i] != 0:
                temp_point = event_number[i]
            else:
                event_number[i] = temp_point
        print("end", event_number[:100], event_number[-1])

        # put the updated event_number back to data frame
        self.df.update(pd.DataFrame({'Event': event_number}))

    def modify_df_group(self,ene=0):
        # change column property. Mainly this change awkuard into str
        print(self.df.dtypes)
        self.df['name'] = self.df['name'].astype(str)
        self.df['Volume'] = self.df['Volume'].astype(str)
        self.df['Process'] = self.df['Process'].astype(str)
        # this make event number correct
        self.reidx_event_group(ene)

    def Capture_spectrum(self):

        self.df_Ncapture = self.df[(self.df["name"]=='neutron')&(self.df["Process"]=='nCapture')&(self.df["Volume"]!='LAr_phys')][['Event','Track ID']]
        self.df_Ncapture_check = self.df[
            (self.df["name"] == 'neutron') & (self.df["Process"] == 'nCapture') & (self.df["Volume"] != 'LAr_phys')][
            ['Event', 'Volume','Track ID']]

        print(self.df_Ncapture.head(10))
        # self.df_cap_gamma = pd.DataFrame('Event','Track ID')
        print("len",len(self.df_Ncapture.index))
        value_counts = self.df_Ncapture_check['Volume'].value_counts()
        print("occurrance",value_counts)
        # only record gamma event, whose event id same as ncap and parent id is ncap's track id.
        # change Track ID name into Parent ID so that ready for merge
        self.df_Ncapture.columns = ['Event','Parent ID']
        # select all gamma events
        self.df_cap_gamma = self.df[self.df['name'] == 'gamma' ]
        # select gamma events whose Event number is same as neutron event and parent id is neutron's track ID
        self.df_cap_gamma = pd.merge(self.df_Ncapture, self.df_cap_gamma,on=['Event','Parent ID'], how='inner')
        print(self.df_cap_gamma.head(20))
        # save these gamma event
        self.df_cap_gamma.to_csv(self.base_path +"dmx_gamma.csv", index=False)

    def LAr_Capture_spectrum(self):

        self.df_Ncapture = self.df[(self.df["name"]=='neutron')&(self.df["Process"]=='nCapture')&(self.df["Volume"]=='LAr_phys')][['Event','Track ID']]
        self.df_Ncapture_check = self.df[
            (self.df["name"] == 'neutron') & (self.df["Process"] == 'nCapture') & (self.df["Volume"] != 'LAr_phys')][
            ['Event', 'Volume','Track ID']]

        print(self.df_Ncapture.head(10),self.df_Ncapture["Event"].unique()[:50])
        # self.df_cap_gamma = pd.DataFrame('Event','Track ID')
        value_counts = self.df_Ncapture_check['Volume'].value_counts()
        print("occurrance",value_counts)
        # only record gamma event, whose event id same as ncap and parent id is ncap's track id.
        # change Track ID name into Parent ID so that ready for merge
        self.df_Ncapture.columns = ['Event','Parent ID']
        # select all gamma events
        self.df_cap_gamma = self.df[self.df['name'] == 'gamma' ]
        # select gamma events whose Event number is same as neutron event and parent id is neutron's track ID
        self.df_cap_gamma_merged = pd.merge(self.df_Ncapture, self.df_cap_gamma,on=['Event','Parent ID'], how='inner')
        print("gamma merged", len(self.df_cap_gamma_merged["Event"].unique()),
              self.df_cap_gamma_merged["Event"].unique()[:20])
        # save these gamma event
        self.df_cap_gamma_merged.to_csv(self.signal_path_mid, index=False)


    def Huge_scatter_spectrum(self):# photon generated by the scattering instead other process

        self.df_Nscatter = self.df[
            (self.df["name"] == 'neutron') & (self.df["Process"] == 'hadElastic') & (self.df["Volume"] == 'LAr_phys')][
            ['Event', 'Volume', 'Track ID', 'Parent ID']]
        self.df_Ninelastic = self.df[
            (self.df["name"] == 'neutron') & (self.df["Process"] == 'neutronInelastic') & (
                        self.df["Volume"] == 'LAr_phys')][
            ['Event', 'Volume', 'Track ID', 'Parent ID']]
        
        self.df_capture = self.df[
            (self.df["name"] == 'neutron') & (self.df["Process"] == 'nCapture') & (self.df["Volume"] == 'LAr_phys')][
            ['Event', 'Volume','Track ID']]
        self.df_capture[["Event"]].to_csv(self.base_path2 + "capture_event_list.csv", index=False)
        print("Ela",len(self.df_Nscatter["Event"].unique()))
        print("capture",self.df_capture.head(10))
        print("inelastic", len(self.df_Ninelastic["Event"].unique()),self.df_Ninelastic.head(10))

        (self.df_sing_Nscatter, self.df_multi_Nscatter) = self.find_single_n_multi(self.df_Nscatter, "Event", "Volume")

        print("sing", self.df_sing_Nscatter)
        print("multi", self.df_multi_Nscatter)
        # self.df_cap_gamma = pd.DataFrame('Event','Track ID')
        # elastic
        merged_df = pd.merge(self.df_sing_Nscatter, self.df_Ninelastic, on=['Event'], how='left', indicator=True)
        filtered_df = merged_df[merged_df['_merge'] == 'left_only'].drop(columns=['_merge'])
        # inelastic
        # merged_df = pd.merge(self.df_sing_Nscatter, self.df_Nscatter, on=['Event'], how='left', indicator=True)
        filtered_df_inela = self.df_Ninelastic
        print("merged_xor,\n", filtered_df_inela.head(10))
        #2nd filter filter out ncapture recoiled energy
        # elastic
        merged_df2 = pd.merge(filtered_df, self.df_capture, on=['Event'], how='left', indicator=True)

        filtered_df2 = merged_df2[merged_df2['_merge'] == 'left_only'].drop(columns=['_merge'])

        merged_df3 = pd.merge(filtered_df_inela, self.df_capture, on=['Event'], how='left', indicator=True)

        filtered_df3 = merged_df3[merged_df3['_merge'] == 'left_only'].drop(columns=['_merge'])

        print("merged_xor,\n", filtered_df2.head(10))

        self.LAr_recoiled = self.df[((self.df["name"] == 'Ar40') | (self.df["name"] == 'Ar36'))&(self.df["Recoiled/keV"]>0.001) ][
            ['Event']]


        # print("LAr recoiled",self.LAr_recoiled)
        # filtered df to remove nCapture event
        # elastic and inelastic
        self.LAr_n_merged = pd.merge(filtered_df2, self.LAr_recoiled, on=['Event'], how='inner')
        self.LAr_n_merged_inela = pd.merge(filtered_df3, self.LAr_recoiled, on=['Event'], how='inner')

        n_list = self.LAr_n_merged["Event"].to_list()
        self.N_check = self.df[self.df["Event"].isin(n_list) & (
                    (self.df["name"] == 'neutron') | (self.df["name"] == 'Ar40') | (self.df["name"] == 'Ar36'))]
        self.N_check.to_csv(self.false_2_path_mid, index=False)

        n_list_inela = self.LAr_n_merged_inela["Event"].to_list()
        print("Ncheck_ inel", n_list_inela)
        self.N_check_inela = self.df[self.df["Event"].isin(n_list_inela) & (
                (self.df["name"] == 'neutron') | (self.df["name"] == 'Ar40') | (self.df["name"] == 'Ar36')| (self.df["name"] == 'gamma'))]
        self.N_check_inela.to_csv(self.false_3_path_mid, index=False)

        # print(self.LAr_n_merged)
        # print("simutanous", len(self.LAr_n_merged["Event"].unique()))

        max_values = self.N_check[( (self.df["name"] == 'Ar40') | (self.df["name"] == 'Ar36'))].groupby(['Event'])["Recoiled/keV"].max().reset_index()
        print(max_values.head(20))

        # add gamma up
        self.Ar_recoiled_list = max_values["Recoiled/keV"].to_list()
        p_observed = []
        scatter_ene = [] # in eV
        for i in range(len(self.Ar_recoiled_list)):
            # 10 /keV 0.03 and 0.2 PCE and PDE
            if i > 1E-6:
                pho_num = self.Ar_recoiled_list[i] * 1E6 * 10 * 0.03 * 0.2 / (1000)
                scatter_ene.append(self.Ar_recoiled_list[i] * 1E6)
                if pho_num > 1:
                    p_observed.append(pho_num)

        num = 0
        for i in p_observed:
            if i >= 1:
                num += 1
        print("photon observed number ", num, len(p_observed))
        print("max", max(p_observed), "\n", "min", min(p_observed))
        # plt.hist(self.Ar_recoiled_list, bins=100)
        # check event 1256
        self.df_event_1542 = self.df[
            self.df["Event"] == 1542][["Event","name","Parent ID","Track ID","Step ID","X/mm","Kinetic/keV","Recoiled/keV", "Volume","Process"]]
        self.df_event_1542.to_csv("/data/runzezhang/result/TN_sims3/event1542.csv", index=False)
        with open(self.false_2_path, 'w', newline='') as myfile:
            wr = csv.writer(myfile)
            wr.writerow(p_observed)

        self.p_observed  = p_observed
        # print photon number
        # plt.hist(p_observed, bins=100)
        # plt.xlabel("Obeserved Photon per Event")

        # print scatter scatter
        plt.hist(scatter_ene, bins=100)
        plt.xscale("log")
        plt.yscale("Log")
        print("scatter number", len(scatter_ene))
        # before 6846
        # after including inelastic scattering 8681
        plt.xlabel("scatter energy per Event")
        # plt.show()
        plt.savefig(self.plot_path+"n_huge_scatter_ene_AmLi2.png")

    def inelastic_gamma(self):
        self.df_gamma_rw = pd.read_csv(self.false_3_path_mid)
        print(self.df_gamma_rw[["Kinetic/keV"]].head(20))

        self.gamma_Scint = self.df_gamma_rw[
            (self.df_gamma_rw['Volume'] == 'LAr_phys')]
        gamma_list = self.gamma_Scint["Event"].unique()
        print("gamma filter", len(gamma_list))
        self.gamma_Scint = self.keep_1st(self.gamma_Scint)
        print("scint", self.gamma_Scint)
        # print("scint2",self.gamma_Scint[self.gamma_Scint["Parent ID"]!=1])
        self.gamma_Scint_column = self.gamma_Scint[['Event', "Track ID"]]
        self.gamma_Scint_column.columns = ['Event', "Parent ID"]
        self.df_electron = self.df[(self.df['name'] == 'e-') & (self.df['Volume'] == 'LAr_phys')]
        self.df_electron = self.keep_1st(self.df_electron)
        self.df_electron_gamma = pd.merge(self.df_electron, self.gamma_Scint_column, on=['Event', 'Parent ID'],
                                          how='inner')
        print("gamma filter 2", len(self.df_electron_gamma["Event"].unique()))
        print(self.df_electron_gamma.head(10))
        # double check gamma

        summed_values = self.df_electron_gamma.groupby(['Event'])["Recoiled/keV"].sum().reset_index()
        print(summed_values.head(20))

        # add gamma up
        self.electron_recoiled_list = summed_values["Recoiled/keV"].to_list()
        p_observed = []
        for i in range(len(self.electron_recoiled_list)):
            # 40 /keV 0.03 and 0.2 PCE and PDE
            if i > 1E-6:
                p_observed.append(self.electron_recoiled_list[i] * 1E6 * 40 * 0.03 * 0.2 / (1000))

        num = 0
        for i in p_observed:
            if i >= 1:
                num += 1
        print("photon observed number ", num, len(p_observed))
        print("max", max(p_observed), "\n", "min", min(p_observed))
        # plt.hist(self.electron_recoiled_list, bins=100)
        self.p_observed += p_observed
        with open(self.false_3_path, 'w', newline='') as myfile:
            wr = csv.writer(myfile)
            wr.writerow(self.p_observed)
        plt.hist(p_observed, bins=100)
        plt.xlabel("Obeserved Photon per Event")

    def bubble_rate(self):  # photon generated by the scattering instead other process

        self.df_Nscatter = self.df[
            (self.df["name"] == 'neutron') & (self.df["Process"] == 'hadElastic') & (self.df["Volume"] == 'LAr_phys')][
            ['Event', 'Volume', 'Track ID', 'Parent ID']]
        self.df_Ninelastic = self.df[
            (self.df["name"] == 'neutron') & (self.df["Process"] == 'neutronInelastic') & (
                    self.df["Volume"] == 'LAr_phys')][
            ['Event', 'Track ID']]

        self.df_capture = self.df[
            (self.df["name"] == 'neutron') & (self.df["Process"] == 'nCapture') & (self.df["Volume"] == 'LAr_phys')][
            ['Event', 'Volume', 'Track ID']]
        self.df_capture[["Event"]].to_csv(self.base_path2 + "capture_event_list.csv", index=False)
        print("Ela", len(self.df_Nscatter["Event"].unique()))
        print("capture", self.df_capture.head(10))

        print("inelastic", len(self.df_Ninelastic["Event"].unique()), self.df_Ninelastic.head(10))

        (self.df_sing_Nscatter, self.df_multi_Nscatter) = self.find_single_n_multi(self.df_Nscatter, "Event", "Volume")

        print("sing", self.df_sing_Nscatter)
        print("multi", self.df_multi_Nscatter)
        # self.df_cap_gamma = pd.DataFrame('Event','Track ID')
        merged_df = pd.merge(self.df_sing_Nscatter, self.df_Ninelastic, on=['Event'], how='left', indicator=True)
        filtered_df = merged_df[merged_df['_merge'] == 'left_only'].drop(columns=['_merge'])
        # use all scatter events

        filtered_df = self.df_Nscatter
        # filtered_df = self.df_sing_Nscatter
        print("merged_xor,\n", filtered_df.head(10))
        # 2nd filter filter out ncapture recoiled energy
        merged_df2 = pd.merge(filtered_df, self.df_capture, on=['Event'], how='left', indicator=True)

        filtered_df2 = merged_df2[merged_df2['_merge'] == 'left_only'].drop(columns=['_merge'])
        print("merged_xor,\n", filtered_df2.head(10))

        self.LAr_recoiled = \
        self.df[((self.df["name"] == 'Ar40') | (self.df["name"] == 'Ar36')) & (self.df["Recoiled/keV"] > 0.0125)][
            ['Event']]
        # print("LAr recoiled",self.LAr_recoiled)
        # filtered df to remove nCapture event
        self.LAr_n_merged = pd.merge(filtered_df2, self.LAr_recoiled, on=['Event'], how='inner')
        # self.LAr_n_merged = pd.merge(self.df_sing_Nscatter, self.LAr_recoiled, on=['Event'], how='inner')
        n_list = self.LAr_n_merged["Event"].to_list()
        self.N_check = self.df[self.df["Event"].isin(n_list) & (
                (self.df["name"] == 'neutron') | (self.df["name"] == 'Ar40') | (self.df["name"] == 'Ar36'))]
        self.N_check.to_csv(self.false_2_path_mid, index=False)
        # print(self.LAr_n_merged)
        # print("simutanous", len(self.LAr_n_merged["Event"].unique()))

        max_values = self.N_check[((self.df["name"] == 'Ar40') | (self.df["name"] == 'Ar36'))].groupby(['Event'])[
            "Recoiled/keV"].max().reset_index()
        print(max_values.head(20))

        # add gamma up
        self.Ar_recoiled_list = max_values["Recoiled/keV"].to_list()
        p_observed = []
        scatter_ene = []  # in eV
        for i in range(len(self.Ar_recoiled_list)):
            # 10 /keV 0.03 and 0.2 PCE and PDE
            if i > 1E-6:
                pho_num = self.Ar_recoiled_list[i] * 1E6 * 10 * 0.03 * 0.2 / (1000)
                scatter_ene.append(self.Ar_recoiled_list[i] * 1E6)
                if pho_num > 1:
                    p_observed.append(pho_num)
        print("all recoil", len(self.Ar_recoiled_list), "scintilti recoil", len(p_observed))
        num = 0
        for i in p_observed:
            if i >= 1:
                num += 1
        print("photon observed number ", num, len(p_observed))
        print("max", max(p_observed), "\n", "min", min(p_observed))
        # plt.hist(self.Ar_recoiled_list, bins=100)
        # check event 1256
        self.df_event_1542 = self.df[
            self.df["Event"] == 1542][
            ["Event", "name", "Parent ID", "Track ID", "Step ID", "X/mm", "Kinetic/keV", "Recoiled/keV", "Volume",
             "Process"]]
        self.df_event_1542.to_csv("/data/runzezhang/result/TN_sims3/event1542.csv", index=False)
        with open(self.false_2_path, 'w', newline='') as myfile:
            wr = csv.writer(myfile)
            wr.writerow(p_observed)
        # print photon number
        # plt.hist(p_observed, bins=100)
        # plt.xlabel("Obeserved Photon per Event")

        # print scatter scatter
        plt.hist(scatter_ene, bins=100)
        plt.xscale("log")
        plt.yscale("Log")
        print("scatter number", len(scatter_ene))
        # before 21434

        # only elastic
        plt.xlabel("scatter energy per Event")
        # plt.show()
        plt.savefig(self.plot_path + "n_huge_scatter_ene_AmLi2.png")

    def Huge_scatter_spectrum_CF(self):  # photon generated by the scattering instead other process

        self.df_Nscatter = self.df[
            (self.df["name"] == 'neutron') & (self.df["Process"] == 'hadElastic') & (self.df["Volume"] == 'LAr_phys')][
            ['Event', 'Volume', 'Track ID', 'Parent ID']]
        self.df_Ninelastic = self.df[
            (self.df["name"] == 'neutron') & (self.df["Process"] == 'neutronInelastic') & (
                    self.df["Volume"] == 'LAr_phys')][
            ['Event', 'Track ID']]

        self.df_capture = self.df[
            (self.df["name"] == 'neutron') & (self.df["Process"] == 'nCapture') & (self.df["Volume"] == 'LAr_phys')][
            ['Event', 'Volume', 'Track ID']]
        self.df_capture[["Event"]].to_csv(self.base_path2 + "capture_event_list.csv", index=False)
        print("capture", self.df_capture.head(10))
        print("inelastic", self.df_Ninelastic.head(10))

        (self.df_sing_Nscatter, self.df_multi_Nscatter) = self.find_single_n_multi(self.df_Nscatter, "Event", "Volume")

        print("sing", self.df_sing_Nscatter)
        print("multi", self.df_multi_Nscatter)
        # self.df_cap_gamma = pd.DataFrame('Event','Track ID')
        merged_df = pd.merge(self.df_sing_Nscatter, self.df_Ninelastic, on=['Event'], how='left', indicator=True)
        filtered_df = merged_df[merged_df['_merge'] == 'left_only'].drop(columns=['_merge'])
        print("merged_xor,\n", filtered_df.head(10))
        # 2nd filter filter out ncapture recoiled energy
        merged_df2 = pd.merge(filtered_df, self.df_capture, on=['Event'], how='left', indicator=True)

        filtered_df2 = merged_df2[merged_df2['_merge'] == 'left_only'].drop(columns=['_merge'])
        print("merged_xor,\n", filtered_df2.head(10))

        self.LAr_recoiled = \
        self.df[((self.df["name"] == 'Ar40') | (self.df["name"] == 'Ar36')) & (self.df["Recoiled/keV"] > 1E-6)][
            ['Event']]
        # print("LAr recoiled",self.LAr_recoiled)
        # filtered df to remove nCapture event
        self.LAr_n_merged = pd.merge(filtered_df2, self.LAr_recoiled, on=['Event'], how='inner')
        # self.LAr_n_merged = pd.merge(self.df_sing_Nscatter, self.LAr_recoiled, on=['Event'], how='inner')
        n_list = self.LAr_n_merged["Event"].to_list()
        self.N_check = self.df[self.df["Event"].isin(n_list) & (
                (self.df["name"] == 'neutron') | (self.df["name"] == 'Ar40') | (self.df["name"] == 'Ar36'))]
        self.N_check.to_csv(self.base_path2 + "dmx_single_n_largescatter_CF_neutron_list.csv", index=False)
        # print(self.LAr_n_merged)
        # print("simutanous", len(self.LAr_n_merged["Event"].unique()))

        max_values = self.N_check[((self.df["name"] == 'Ar40') | (self.df["name"] == 'Ar36'))].groupby(['Event'])[
            "Recoiled/keV"].max().reset_index()
        print(max_values.head(20))

        # add gamma up
        self.Ar_recoiled_list = max_values["Recoiled/keV"].to_list()
        p_observed = []
        for i in range(len(self.Ar_recoiled_list)):
            # 40 /keV 0.03 and 0.2 PCE and PDE
            if i > 1E-6:
                pho_num = self.Ar_recoiled_list[i] * 1E6 * 10 * 0.03 * 0.2 / (1000)
                if pho_num > 1:
                    p_observed.append(pho_num)

        num = 0
        for i in p_observed:
            if i >= 1:
                num += 1
        print("photon observed number ", num, len(p_observed))
        print("max", max(p_observed), "\n", "min", min(p_observed))
        # plt.hist(self.Ar_recoiled_list, bins=100)
        # check event 1256
        self.df_event_390 = self.df[
            self.df["Event"] == 390][
            ["Event", "name", "Parent ID", "Track ID", "Step ID", "X/mm", "Kinetic/keV", "Recoiled/keV", "Volume",
             "Process"]]
        self.df_event_390.to_csv("/data/runzezhang/result/TN_sims3/event390.csv", index=False)
        with open("/data/runzezhang/result/TN_sims3/n_huge_scatterg_CF2.csv", 'w', newline='') as myfile:
            wr = csv.writer(myfile)
            wr.writerow(p_observed)
        plt.hist(p_observed, bins=100)
        plt.xlabel("Obeserved Photon per Event")
        # plt.show()
        plt.savefig(self.plot_path + "n_huge_scatter_CF2.png")


    def Huge_scatter_spectrum_CF_fake(self):  # photon generated by the scattering instead other process

        self.df_Nscatter = self.df[
            (self.df["name"] == 'neutron') & (self.df["Process"] == 'hadElastic') & (self.df["Volume"] == 'LAr_phys')][
            ['Event', 'Volume', 'Track ID', 'Parent ID']]
        self.df_Ninelastic = self.df[
            (self.df["name"] == 'neutron') & (self.df["Process"] == 'neutronInelastic') & (
                    self.df["Volume"] == 'LAr_phys')][
            ['Event', 'Track ID']]

        self.df_capture = self.df[
            (self.df["name"] == 'neutron') & (self.df["Process"] == 'nCapture') & (self.df["Volume"] == 'LAr_phys')][
            ['Event', 'Volume', 'Track ID']]
        self.df_capture[["Event"]].to_csv(self.base_path2 + "capture_event_list.csv", index=False)
        print("capture", self.df_capture.head(10))
        print("inelastic", self.df_Ninelastic.head(10))

        (self.df_sing_Nscatter, self.df_multi_Nscatter) = self.find_single_n_multi(self.df_Nscatter, "Event", "Volume")

        print("sing", self.df_sing_Nscatter)
        print("multi", self.df_multi_Nscatter)
        # self.df_cap_gamma = pd.DataFrame('Event','Track ID')
        merged_df = pd.merge(self.df_sing_Nscatter, self.df_Ninelastic, on=['Event'], how='left', indicator=True)
        filtered_df = merged_df[merged_df['_merge'] == 'left_only'].drop(columns=['_merge'])
        print("merged_xor,\n", filtered_df.head(10))


        self.LAr_recoiled = \
        self.df[((self.df["name"] == 'Ar40') | (self.df["name"] == 'Ar36')) & (self.df["Recoiled/keV"] > 1E-6)][
            ['Event']]
        # print("LAr recoiled",self.LAr_recoiled)
        # filtered df to remove nCapture event
        self.LAr_n_merged = pd.merge(filtered_df, self.LAr_recoiled, on=['Event'], how='inner')
        # self.LAr_n_merged = pd.merge(self.df_sing_Nscatter, self.LAr_recoiled, on=['Event'], how='inner')
        n_list = self.LAr_n_merged["Event"].to_list()
        self.N_check = self.df[self.df["Event"].isin(n_list) & (
                (self.df["name"] == 'neutron') | (self.df["name"] == 'Ar40') | (self.df["name"] == 'Ar36'))]
        self.N_check.to_csv(self.base_path2 + "dmx_single_n_largescatter_CF_neutron_fake_list.csv", index=False)
        # print(self.LAr_n_merged)
        # print("simutanous", len(self.LAr_n_merged["Event"].unique()))

        max_values = self.N_check[((self.df["name"] == 'Ar40') | (self.df["name"] == 'Ar36'))].groupby(['Event'])[
            "Recoiled/keV"].max().reset_index()
        print(max_values.head(20))

        # add gamma up
        self.Ar_recoiled_list = max_values["Recoiled/keV"].to_list()
        p_observed = []
        for i in range(len(self.Ar_recoiled_list)):
            # 40 /keV 0.03 and 0.2 PCE and PDE
            if i > 1E-6:
                pho_num = self.Ar_recoiled_list[i] * 1E6 * 10 * 0.03 * 0.2 / (1000)
                if pho_num > 1:
                    p_observed.append(pho_num)

        num = 0
        for i in p_observed:
            if i >= 1:
                num += 1
        print("photon observed number ", num, len(p_observed))
        print("max", max(p_observed), "\n", "min", min(p_observed))
        # plt.hist(self.Ar_recoiled_list, bins=100)
        # check event 1256
        self.df_event_390 = self.df[
            self.df["Event"] == 390][
            ["Event", "name", "Parent ID", "Track ID", "Step ID", "X/mm", "Kinetic/keV", "Recoiled/keV", "Volume",
             "Process"]]
        self.df_event_390.to_csv("/data/runzezhang/result/TN_sims3/event390_fake.csv", index=False)
        with open("/data/runzezhang/result/TN_sims3/n_huge_scatterg_CF2_fake.csv", 'w', newline='') as myfile:
            wr = csv.writer(myfile)
            wr.writerow(p_observed)
        plt.hist(p_observed, bins=100)
        plt.xlabel("Obeserved Photon per Event")
        # plt.show()
        plt.savefig(self.plot_path + "n_huge_scatter_CF2_fake.png")
    def Capture_n_scatter_spectrum(self): # somehow logan made the cross where the id difference is 1 like compare 4 scatter with 5 photon generation

        self.df_Ncapture = self.df[(self.df["name"]=='neutron')&(self.df["Process"]=='nCapture')&(self.df["Volume"]!='LAr_phys')][['Event','Track ID']]
        self.df_head = self.df[(self.df["Event"]==5837)|(self.df["Event"]==2906)|(self.df["Event"]==2907)][["Event","name","Parent ID","Track ID","Step ID","X/mm","Kinetic/keV","Recoiled/keV", "Volume","Process"]]
        self.df_Nscatter = self.df[
            (self.df["name"] == 'neutron') & (self.df["Process"] == 'hadElastic') & (self.df["Volume"] == 'LAr_phys')][
            ['Event', 'Volume','Track ID', 'Parent ID']]
        self.df_head.to_csv(self.base_path + "dmx_single_n_gamma_5837.csv", index=False)
        self.df_Ninelastic = self.df[
            (self.df["name"] == 'neutron') & (self.df["Process"] == 'neutronInelastic') & (self.df["Volume"] == 'LAr_phys')][
            ['Event', 'Track ID']]
        print("inelastic", self.df_Ninelastic.head(10))

        (self.df_sing_Nscatter, self.df_multi_Nscatter) = self.find_single_n_multi(self.df_Nscatter,"Event", "Volume")

        print("sing", self.df_sing_Nscatter)
        print("multi", self.df_multi_Nscatter)
        # self.df_cap_gamma = pd.DataFrame('Event','Track ID')
        print("len",len(self.df_Ncapture.index))
        merged_df = pd.merge(self.df_sing_Nscatter, self.df_Ninelastic, on=['Event'], how='left', indicator=True)
        print("merged_xor,\n", merged_df.head(10))

        # Filter the merged DataFrame to keep only rows that are in df1 but not in df2
        result_df = merged_df[merged_df['_merge'] == 'left_only'].drop(columns=['_merge','Track ID_y'])
        result_df.columns = ['Event', 'Volume','Track ID', 'Parent ID']
        print("xor", result_df.head(20))

        self.df_n = pd.merge(result_df, self.df_Ncapture, on=['Event', 'Track ID'], how='inner')
        # self.df_n = pd.merge(self.df_sing_Nscatter, self.df_Ncapture,on=['Event','Track ID'], how='inner')


        # self.LAr_recoiled = self.df[((self.df["name"]=='Ar40') | (self.df["name"]=='Ar36') )&(self.df["Recoiled/keV"]>1E-3)][['Event','Track ID',"Recoiled/keV"]]
        # self.LAr_recoiled = \
        # self.df[((self.df["name"] == 'Ar40') | (self.df["name"] == 'Ar36')) & (self.df["Recoiled/keV"] > 0.001)][
        #     ['Event']]
        self.LAr_recoiled = self.df[((self.df["name"] == 'Ar40') | (self.df["name"] == 'Ar36')) ][
            ['Event']]
        # print("LAr recoiled",self.LAr_recoiled)
        self.LAr_n_merged = pd.merge(self.df_n,self.LAr_recoiled,on=['Event'], how='inner')
        n_list = self.LAr_n_merged["Event"].to_list()
        self.N_check = self.df[self.df["Event"].isin(n_list) & ((self.df["name"]=='neutron')|(self.df["name"]=='Ar40')| (self.df["name"] == 'Ar36'))]
        self.N_check.to_csv(self.base_path +"dmx_single_n_gamma_CF_neutron_list.csv", index=False)
        print(self.LAr_n_merged)
        print("simutanous", len(self.LAr_n_merged["Event"].unique()))
        # only record gamma event, whose event id same as ncap and parent id is ncap's track id.
        # change Track ID name into Parent ID so that ready for merge
        self.df_n_slice = self.LAr_n_merged[['Event','Track ID']]
        self.df_n_slice.columns = ['Event','Parent ID']
        # select all gamma events
        self.df_cap_gamma = self.df[self.df['name'] == 'gamma' ]
        # select gamma events whose Event number is same as neutron event and parent id is neutron's track ID
        self.df_single_n_gamma = pd.merge(self.df_n_slice, self.df_cap_gamma,on=['Event','Parent ID'], how='inner')
        print(self.df_single_n_gamma.head(20))
        # save these gamma event
        self.df_single_n_gamma.to_csv(self.base_path +"dmx_single_n_gamma_CF.csv", index=False)

    def Capture_n_scatter_spectrum_loop(self):
        self.df["Kinetic diff/MeV"] = self.df["Kinetic/keV"].diff()
        self.df["Kinetic diff/MeV"] = self.df["Kinetic diff/MeV"].fillna(0)

        self.df_Ncapture = self.df[
            (self.df["name"] == 'gamma') & (self.df["Process"] == 'nCapture') & (self.df["Volume"] != 'LAr_phys')][
            ['Event', 'Track ID']]
        self.df_Nscatter = self.df[
            (self.df["name"] == 'neutron') & (self.df["Process"] == 'hadElastic') & (self.df["Volume"] == 'LAr_phys')& (self.df["Kinetic diff/MeV"] <-0.00125)][
            ['Event', 'Volume', 'Track ID', 'Parent ID']]
        self.df_Nscatter_wo = self.df[
            (self.df["name"] == 'neutron') & (self.df["Process"] == 'hadElastic') & (
                        self.df["Volume"] == 'LAr_phys') ][
            ['Event', 'Volume', 'Track ID', 'Parent ID']]
        self.df_Ncapture_step = self.df[
            (self.df["name"] == 'neutron') & (self.df["Process"] == 'nCapture') & (self.df["Volume"] != 'LAr_phys')][
            ['Event', 'Track ID','Step ID']]
        list_scatter =self.df_Nscatter["Event"].unique()
        list_capture = self.df_Ncapture_step["Event"].unique()
        print("Scaterr unique", self.df_Nscatter["Event"].unique(),len(self.df_Nscatter["Event"].unique()))
        print("capture unique", self.df_Ncapture_step["Event"].unique())
        offset_list = []
        # check if there is +-1 relationship
        # for i in range(len(list_scatter)):
        #     for j in range(len(list_capture)):
        #         if int(list_scatter[i])-int(list_capture[j])==-1:
        #             offset_list.append(int(list_scatter[i]))
        # print("offsetlist", offset_list)
        self.df_Nscatter_wo_step = self.df[
            (self.df["name"] == 'neutron') & (self.df["Process"] == 'hadElastic') & (
                    self.df["Volume"] == 'LAr_phys')][
            ['Event', 'Volume', 'Track ID', 'Parent ID','Step ID']]
        
        self.df_Ninelastic = self.df[
            (self.df["name"] == 'neutron') & (self.df["Process"] == 'neutronInelastic') & (
                        self.df["Volume"] == 'LAr_phys')][
            ['Event', 'Track ID']]

        self.compt_scatter = self.df[(self.df["name"] == 'gamma')&(self.df["Process"] == 'compt') & (
                        self.df["Volume"] == 'LAr_phys')][
            ['Event', 'Track ID']]
        self.elcap = pd.merge(self.df_Ncapture_step, self.df_Nscatter_wo_step, on=['Event'], how='inner', indicator=True)
        self.elcap_stepfilter = self.elcap[self.elcap['Step ID_x'] > self.elcap['Step ID_y']]
        print("cap wo", len(self.df_Ncapture["Event"].unique()))
        print("Ela",len(self.df_Nscatter_wo["Event"].unique()))
        print("Elcap", len(self.elcap["Event"].unique()))
        print("Elcap step", len(self.elcap_stepfilter["Event"].unique()))
        print("compt", len(self.compt_scatter["Event"].unique()))
        print("inelastic", self.df_Ninelastic.head(10))

        (self.df_sing_Nscatter, self.df_multi_Nscatter) = self.find_single_n_multi(self.df_Nscatter, "Event", "Volume")

        print("sing", self.df_sing_Nscatter, len(self.df_sing_Nscatter["Event"].unique()))
        print("multi", self.df_multi_Nscatter)
        # self.df_cap_gamma = pd.DataFrame('Event','Track ID')
        print("len", len(self.df_Ncapture.index))
        merged_df = pd.merge(self.df_sing_Nscatter, self.df_Ninelastic, on=['Event'], how='left', indicator=True)
        print("merged_xor,\n", merged_df.head(10))

        # Filter the merged DataFrame to keep only rows that are in df1 but not in df2
        result_df = merged_df[merged_df['_merge'] == 'left_only'].drop(columns=['_merge', 'Track ID_y'])
        result_df.columns = ['Event', 'Volume', 'Track ID', 'Parent ID']
        print("xor", result_df.head(20))

        self.df_n = pd.merge(result_df, self.df_Ncapture, on=['Event', 'Track ID'], how='inner')
        # self.df_n = pd.merge(self.df_sing_Nscatter, self.df_Ncapture,on=['Event','Track ID'], how='inner')

        # self.x3_n = pd.merge(self.df_Ncapture, self.compt_scatter, on=['Event'], how='inner')
        self.x3_n = pd.merge(self.df_n, self.compt_scatter, on=['Event'], how='inner')
        # self.x3_n = pd.merge(self.elcap_stepfilter, self.compt_scatter, on=['Event'], how='inner')
        self.x3_n2 = pd.merge(self.x3_n, self.df_Nscatter_wo, on=['Event'], how='inner')
        print("cross 3 check", len(self.x3_n2["Event"].unique()), self.x3_n2.head(10))

        self.last_cross = pd.merge(self.x3_n2, self.df_sing_Nscatter, on=['Event'], how='inner')
        self.last_cross_list = self.intersection(self.x3_n2["Event"].unique(),self.df_sing_Nscatter["Event"].unique())
        print("last cross check", len(self.last_cross_list), self.last_cross_list[:10])
        # # self.LAr_recoiled = self.df[((self.df["name"]=='Ar40') | (self.df["name"]=='Ar36') )&(self.df["Recoiled/keV"]>1E-3)][['Event','Track ID',"Recoiled/keV"]]
        # self.LAr_recoiled = \
        #     self.df[((self.df["name"] == 'Ar40') | (self.df["name"] == 'Ar36')) ][
        #         ['Event']]
        # # print("LAr recoiled",self.LAr_recoiled)
        # self.LAr_n_merged = pd.merge(self.df_n, self.LAr_recoiled, on=['Event'], how='inner')
        n_list = self.last_cross["Event"].to_list()
        # self.N_check = self.df[self.df["Event"].isin(self.last_cross_list) & (
        #             (self.df["Process"] == 'nCapture')|(self.df["Process"] == 'hadElastic') | (self.df["name"] == 'Ar40') | (self.df["name"] == 'Ar36')|(self.df["Process"] == 'compt'))]
        # self.N_check.to_csv(self.base_path + "dmx_single_n_gamma_CF_neutron_list_loop2.csv", index=False)
        # n_list = self.LAr_n_merged["Event"].to_list()
        self.N_check = self.df[self.df["Event"].isin(n_list) & (
                (self.df["name"] == 'neutron') | (self.df["name"] == 'Ar40') | (self.df["name"] == 'Ar36')| (self.df["name"] == 'gamma'))]
        self.N_check.to_csv(self.false_1_path_mid, index=False)
        # print(self.LAr_n_merged)
        # print("simutanous", len(self.LAr_n_merged["Event"].unique()))
        # event_list  =  self.N_check["Event"].unique()
        # print("event", len(event_list),event_list[:10])

    def intersection(self, lst1, lst2):
        lst3 = [value for value in lst1 if value in lst2]
        return lst3

    def LAr_n_single_test(self):
        self.df_20575 = self.df[(self.df["Event"]==581)|(self.df["Event"]==2906)|(self.df["Event"]==329)|(self.df["Event"]==568)]
        self.df_20575.to_csv(self.base_path + "dmx_single_n_gamma_CF_20575.csv", index=False)
    def LAr_compare(self):
        self.df_Ncapture = self.df[
            (self.df["name"] == 'neutron') & (self.df["Process"] == 'nCapture') & (self.df["Volume"] == 'LAr_phys')][
            ['Event', 'Track ID']]



        print("N capture events", len(self.df_Ncapture["Event"].unique()),self.df_Ncapture["Event"].unique()[:50])

        # only record gamma event, whose event id same as ncap and parent id is ncap's track id.
        # change Track ID name into Parent ID so that ready for merge
        self.df_Ncapture.columns = ['Event', 'Parent ID']
        # select all gamma events
        self.df_cap_gamma = self.df[self.df['name'] == 'gamma']
        print("gamma events", len(self.df_cap_gamma["Event"].unique()), self.df_cap_gamma["Event"].unique()[:20])
        # select gamma events whose Event number is same as neutron event and parent id is neutron's track ID
        self.df_cap_gamma_merged = pd.merge(self.df_Ncapture, self.df_cap_gamma, on=['Event', 'Parent ID'], how='inner')
        # print(self.df_cap_gamma.head(20))
        print("gamma merged", len(self.df_cap_gamma_merged["Event"].unique()), self.df_cap_gamma_merged["Event"].unique()[:20])
        lost_event = list(set(self.df_Ncapture["Event"].unique()) - set(self.df_cap_gamma_merged["Event"].unique()))
        print("lost", lost_event)
        print("len2", len(self.df_cap_gamma_merged.index))
    def test_merge(self):
        df_a = pd.DataFrame({ 'B':[3,5],'C':[5,7]})
        df_b = pd.DataFrame({'B': [2,3,5,3], 'C': [3,5,7,5], 'F': [5,9,10,6], 'G': [8,10,7,9]})
        print('a\n',df_a)
        print('b\n',df_b)
        merged_df = pd.merge(df_b, df_a, on=['B','C'], how='inner')

        print(merged_df)
    def gamma_event(self):
        # if already run 1st 2 steps and obtained output csv file, one can directly run 3rd function

        self.Capture_spectrum()
        # self.Gamma_spectrum()
        self.find_gamma_e()
        self.check_capture()
        # self.plot_gamma()
    def neutron_momentum(self):
        self.df.to_csv(self.false_3_path_mid, index=False)
        # z face is 1150mm
        self.df_neutron_income = self.df[
            (self.df["name"] == 'neutron') & (self.df["Parent ID"] ==0)&(self.df["Step ID"] ==1)&(self.df["Volume"] =="physSD2")][
            ['Event', 'Volume', 'Track ID', 'X/mm','Y/mm','Z/mm','px/MeV','py/MeV','pz/MeV',"Kinetic/keV",'Parent ID']]
        print("first income shoule be 1E6", len(self.df_neutron_income["Event"].to_list()),self.df_neutron_income.head(100))

        #?? escape the outerface of sapphire
        # 157.5
        self.df_neutron_outcome = self.df[
            (self.df["name"] == 'neutron')& (
                    self.df["Parent ID"] == 0)&(self.df["Volume"] =="physWorld")&(self.df["Z/mm"] <= 158)&(self.df["Z/mm"] >= 157)&(self.df["X/mm"] <= 25)&(self.df["Y/mm"] <= 25)][
            ['Event', 'Volume', 'Track ID', 'X/mm', 'Y/mm', 'Z/mm', 'px/MeV', 'py/MeV', 'pz/MeV', "Kinetic/keV",
             'Parent ID']]
        print("first outcome",len(self.df_neutron_outcome["Event"].to_list()),self.df_neutron_outcome.head(100))
        self.df_neutron_outcome = self.keep_1st(self.df_neutron_outcome)
        

        
        # self.df_neutron_outcome.to_csv(self.false_3_path_mid, index=False)

        # neutron_energy = self.df_neutron_outcome[((self.df["name"] == 'Ar40') | (self.df["name"] == 'Ar36'))].groupby(['Event'])[
        #     "Recoiled/keV"].max().reset_index()

        # neutron_energy =   self.df_neutron_income

        neutron_energy =  self.df_neutron_outcome


        # add gamma up
        self.neutron_Ek_list=neutron_energy["Kinetic/keV"].to_list()
        self.neutron_px_list= neutron_energy["px/MeV"].to_list()
        self.neutron_py_list = neutron_energy["py/MeV"].to_list()
        self.neutron_pz_list = neutron_energy["pz/MeV"].to_list()
        # for i in range(len(self.neutron_px_list)):
        #     self.neutron_Ek_list.append((self.neutron_px_list[i]**2+self.neutron_py_list[i]**2+self.neutron_pz_list[i]**2)**0.5)
        print("nenutron in 1E6 ", len(self.neutron_Ek_list))

        with open(self.false_3_path, 'w', newline='') as myfile:
            wr = csv.writer(myfile)
            wr.writerow(self.neutron_Ek_list)

    def plot_neutron_momentum(self):

        with open(self.false_3_path, 'r') as file:
            reader = csv.reader(file)
            # Read the first row (assuming single row for simplicity)
            number_list = next(reader)
            # Convert the strings to floats
            self.noise3_raw_list = [float(value)*1e6 for value in number_list]

        # sig_counts, sig_bin_edges, _ = plt.hist(self.noise3_raw_list, bins= 100)
        sig_counts, sig_bin_edges, _ = plt.hist(self.noise3_raw_list, bins=np.logspace(-5, 7, 50))
        # sig_normalized_counts = 1
        # sig_bin_centers = (sig_bin_edges[:-1] + sig_bin_edges[1:]) / 2
        # plt.bar(sig_bin_centers, sig_normalized_counts, width=sig_bin_edges[1] - sig_bin_edges[0], color='red',
        #         label='signal')

        plt.xlabel("neutron energy/MeV", fontsize=16)
        plt.ylabel("counts", fontsize=16)
        plt.yscale('log')
        plt.xscale('log')
        plt.xlim(1e-5,1e7)
        plt.legend()
        plot_name = "sn1_neutron_outcome_1E6.png"
        plt.savefig(self.plot_path + plot_name)

    def ncrystal_test(self):
        self.df.to_csv(self.false_3_path_mid, index=False)
        # z face is 1150mm
        # self.df_neutron_income = self.df[
        #     (self.df["name"] == 'neutron') & (self.df["Z/mm"] >=1240)& (self.df["Z/mm"] <=1260)& (self.df["Parent ID"] ==0)][
        #     ['Event', 'Volume', 'Track ID', 'X/mm','Y/mm','Z/mm','px/MeV','py/MeV','pz/MeV',"Kinetic/keV",'Parent ID']]
        # self.df_neutron_income = self.df[
        #     (self.df["name"] == 'neutron') & (self.df["Step ID"] == 0) & (
        #                 self.df["Parent ID"] == 0)][
        #     ['Event', 'Volume', 'Track ID', 'X/mm', 'Y/mm', 'Z/mm', 'px/MeV', 'py/MeV', 'pz/MeV', "Kinetic/keV",
        #      'Parent ID']]
        self.df_neutron_income = self.df_neutron_cross = self.df[
            (self.df["name"] == 'neutron')  & (
                    self.df["Parent ID"] == 0)&(self.df["Step ID"] == 1)&(self.df["Volume"] == "physWorld")][
            ['Event', 'Volume', 'Track ID', 'X/mm', 'Y/mm', 'Z/mm', 'px/MeV', 'py/MeV', 'pz/MeV', "Kinetic/keV",
             'Parent ID']]

        # self.df_neutron_cross = self.df[
        #     (self.df["name"] == 'neutron') & (
        #                 (self.df["Process"] == 'neutronInelastic') | (self.df["Process"] == "nCapture")) & (
        #             self.df["Parent ID"] == 0)&(self.df["Volume"] == "physSap")][
        #     ['Event', 'Volume', 'Track ID', 'X/mm', 'Y/mm', 'Z/mm', 'px/MeV', 'py/MeV', 'pz/MeV', "Kinetic/keV",
        #      'Parent ID']]

        self.df_neutron_cross = self.df[
            (self.df["name"] == 'neutron') & (
                    (self.df["Process"] == 'hadElastic') | (self.df["Process"] == "nCapture")) & (
                    self.df["Parent ID"] == 0) & (self.df["Volume"] == "physSap")][
            ['Event', 'Volume', 'Track ID', 'X/mm', 'Y/mm', 'Z/mm', 'px/MeV', 'py/MeV', 'pz/MeV', "Kinetic/keV",
             'Parent ID']]



        # test neutron just hit sapphire
        # self.df_neutron_cross = self.df[
        #     (self.df["name"] == 'neutron') &  (
        #             self.df["Parent ID"] == 0) & (self.df["Volume"] == "physSap")][
        #     ['Event', 'Volume', 'Track ID', 'X/mm', 'Y/mm', 'Z/mm', 'px/MeV', 'py/MeV', 'pz/MeV', "Kinetic/keV",
        #      'Parent ID']]
        print("first iincome", self.df_neutron_cross.head(100))
        self.df_neutron_income = self.keep_1st(self.df_neutron_income)

        neutron_energy = \
            self.df_neutron_income

        # add gamma up
        self.neutron_Ek_dic={}
        self.neutron_event_list = neutron_energy["Event"].to_list()
        self.neutron_ek_list = neutron_energy["Kinetic/keV"].to_list()
        self.neutron_px_list = neutron_energy["px/MeV"].to_list()
        self.neutron_py_list = neutron_energy["py/MeV"].to_list()
        self.neutron_pz_list = neutron_energy["pz/MeV"].to_list()
        for i in range(len(self.neutron_px_list)):
            self.neutron_Ek_dic[self.neutron_event_list[i]]= [self.neutron_ek_list[i] ,0]

        print("enutron in 1E6 ", len(self.neutron_Ek_dic) )

        self.df_neutron_cross = self.keep_1st(self.df_neutron_cross)

        neutron_final_energy = \
            self.df_neutron_cross

        # add gamma up

        self.neutron_final_event_list = neutron_final_energy["Event"].to_list()
        self.neutron_final_ek_list = neutron_final_energy["Kinetic/keV"].to_list()
        self.neutron_final_px_list = neutron_final_energy["px/MeV"].to_list()
        self.neutron_final_py_list = neutron_final_energy["py/MeV"].to_list()
        self.neutron_final_pz_list = neutron_final_energy["pz/MeV"].to_list()
        for i in range(len(self.neutron_final_px_list)):
            self.neutron_Ek_dic[self.neutron_final_event_list[i]][1] =  self.neutron_final_ek_list[i]
        print("enutron cross  in 1E6 ", len(self.neutron_final_event_list) )
        # ouotput the energy, the incident number of neutron, and captured/inelastic neutron number
        return (self.neutron_ek_list[0],len(self.neutron_Ek_dic) ,len(self.neutron_final_event_list))

    def plot_ncrystal_test(self):

        # ini_ene = [self.neutron_Ek_dic[key][0]*1e6 for key in self.neutron_Ek_dic]
        # final_ene = [self.neutron_Ek_dic[key][1]*1e6 for key in self.neutron_Ek_dic]
        # cross_list = []
        # for i in range(len(final_ene)):
        #     if final_ene[i] != 0:
        #         cross_list.append(ini_ene[i])
        #
        #
        # # print(ini_ene)
        # print("cross list",len(cross_list))
        # print("ini, max, min", max(ini_ene), min(ini_ene))
        # # sig_counts, sig_bin_edges, _ = plt.hist(self.neutron_ek_list, bins=np.logspace(-5, 7, 50))
        # sig_counts, sig_bin_edges, _ = plt.hist(ini_ene, bins=np.logspace(-5, 7, 50))
        # # sig_counts, sig_bin_edges, _ = plt.hist(cross_list, bins= np.logspace(-5, 7, 50))
        # # sig_normalized_counts = 1
        # # sig_bin_centers = (sig_bin_edges[:-1] + sig_bin_edges[1:]) / 2
        # # plt.bar(sig_bin_centers, sig_normalized_counts, width=sig_bin_edges[1] - sig_bin_edges[0], color='red',
        # #         label='signal')

        ene = [ i*1e6 for i in self.ene]
        cross = [i for i in self.cross_number]
        print("ene",ene )
        print("cross", cross)
        plt.plot(ene,cross)
        plt.xlabel("neutron energy/eV", fontsize=16)
        plt.ylabel("counts", fontsize=16)
        plt.yscale('log')
        plt.xscale('log')
        plt.xlim(1e-5, 1e7)
        plt.legend()
        plot_name = "sn1_neutron_crystal_1E6.png"
        plt.savefig(self.plot_path + plot_name)

    def Huge_scatter_event(self):
        # single scatter spectrum
        self.Huge_scatter_spectrum()
        # self.inelastic_gamma()
        # self.Huge_scatter_spectrum_CF()
        # self.Huge_scatter_spectrum_CF_fake()

    def LAr_gamma_event(self):
        # for liquid argon capture
        # if already run 1st 2 steps and obtained output csv file, one can directly run 3rd function
        self.LAr_compare()
        self.LAr_Capture_spectrum()
        self.LAr_find_gamma_e()

    def single_e_n_capture_event(self):
        # loop and without loop is just to test the algrorithms, the result should be same
        # self.LAr_n_single_test()
        # self.Capture_n_scatter_spectrum()

        self.Capture_n_scatter_spectrum_loop()
        # self.single_n_find_gamma_e_loop()
        # get the photon number per event
        # self.single_n_find_gamma_e()
    def Gamma_spectrum(self):
        self.df_gamma_rw = pd.read_csv(self.base_path +"dmx_gamma.csv")
        print(self.df_gamma_rw[["Kinetic/keV"]].head(20))
        # we need to do severalthings:
        # gamma only in LAr or CF4
        # in 1 event number, only the first series of gammas, avoiding over-countting
        self.gamma_Scint = self.df_gamma_rw[(self.df_gamma_rw['Volume']=='LAr_phys') | (self.df_gamma_rw['Volume']=='hydraulic_fluid_phys')]
        event_p =0
        track_p=[]
        parent_p = []

        for index in range(len(self.gamma_Scint.index)):
            print(index)
            if self.gamma_Scint.iloc[index]['Event']> event_p:
                event_p = self.gamma_Scint.iloc[index]['Event']
                track_p = []
                parent_p = []
            # same event
            # if parent ID appears in previous Track ID, then drop the row
            # else record the 1st track ID in same trajactory
            else:
                if self.gamma_Scint.iloc[index]['Parent ID'] in track_p:
                    self.gamma_Scint.drop(self.gamma_Scint.index[index])
                else:
                    if self.gamma_Scint.iloc[index]['Track ID'] not in track_p:
                        track_p.append(self.gamma_Scint.iloc[index]['Track ID'] )
                    elif self.gamma_Scint.iloc[index]['Track ID'] in track_p:
                        self.gamma_Scint.drop(self.gamma_Scint.index[index])

        print(self.gamma_Scint.head(20))
        self.gamma_Scint.to_csv(self.base_path +"gamma_scint2.csv", index=False)
    def find_gamma_e(self):
        self.df_gamma_rw = pd.read_csv(self.base_path + "dmx_gamma.csv")
        print(self.df_gamma_rw[["Kinetic/keV"]].head(20))
        # we need to do severalthings:
        # gamma only in LAr or CF4
        # in 1 event number, only the first series of gammas, avoiding over-countting
        self.gamma_Scint = self.df_gamma_rw[
            (self.df_gamma_rw['Volume'] == 'LAr_phys') | (self.df_gamma_rw['Volume'] == 'hydraulic_fluid_phys')]
        self.gamma_Scint = self.keep_1st(self.gamma_Scint)
        print("scint",self.gamma_Scint)
        # print("scint2",self.gamma_Scint[self.gamma_Scint["Parent ID"]!=1])
        self.gamma_Scint_column = self.gamma_Scint[['Event',"Track ID"]]
        self.gamma_Scint_column.columns = ['Event',"Parent ID"]
        self.df_electron = self.df[(self.df['name']=='e-')&(self.df['Volume']=='LAr_phys')]
        self.df_electron = self.keep_1st(self.df_electron)
        self.df_electron_gamma = pd.merge(self.df_electron,self.gamma_Scint_column,on=['Event','Parent ID'], how='inner')
        print(self.df_electron_gamma.head(10))
        # double check gamma

        summed_values = self.df_electron_gamma.groupby(['Event'])["Recoiled/keV"].sum().reset_index()
        print(summed_values.head(20))

        # add gamma up
        self.electron_recoiled_list  = summed_values["Recoiled/keV"].to_list()
        p_observed = []
        for i in range(len(self.electron_recoiled_list)):
            # 40 /keV 0.03 and 0.2 PCE and PDE
            p_observed.append(self.electron_recoiled_list[i]*1E6*40*0.03*0.2/(1000))

        num = 0
        for i in p_observed:
            if i >= 1:
                num += 1
        print("photon observed number ", num)
        print("max", max(p_observed), "\n", "min", min(p_observed))
        # plt.hist(self.electron_recoiled_list, bins=100)
        with open("/data/runzezhang/result/TN_e_sims/photon_captureout.csv", 'w', newline='') as myfile:
            wr = csv.writer(myfile)
            wr.writerow(p_observed)
        plt.hist(p_observed, bins=100)
        plt.xlabel("Obeserved Photon per Event")
        plt.show()

    def single_n_find_gamma_e(self):
        self.df_gamma_rw = pd.read_csv(self.false_1_path_mid)
        print(self.df_gamma_rw[["Kinetic/keV"]].head(20))
        # we need to do severalthings:
        # gamma only in LAr or CF4
        # in 1 event number, only the first series of gammas, avoiding over-countting
        # self.gamma_Scint = self.df_gamma_rw[
        #     (self.df_gamma_rw['Volume'] == 'LAr_phys') | (self.df_gamma_rw['Volume'] == 'hydraulic_fluid_phys')]
        self.gamma_Scint = self.df_gamma_rw[
            (self.df_gamma_rw['Volume'] == 'LAr_phys') ]
        gamma_list = self.gamma_Scint["Event"].unique()
        print("gamma filter", len(gamma_list))
        self.gamma_Scint = self.keep_1st(self.gamma_Scint)
        print("scint",self.gamma_Scint)
        # print("scint2",self.gamma_Scint[self.gamma_Scint["Parent ID"]!=1])
        self.gamma_Scint_column = self.gamma_Scint[['Event',"Track ID"]]
        self.gamma_Scint_column.columns = ['Event',"Parent ID"]
        self.df_electron = self.df[(self.df['name']=='e-')&(self.df['Volume']=='LAr_phys')]
        self.df_electron = self.keep_1st(self.df_electron)
        self.df_electron_gamma = pd.merge(self.df_electron,self.gamma_Scint_column,on=['Event','Parent ID'], how='inner')
        print("gamma filter 2",len(self.df_electron_gamma["Event"].unique()))
        print(self.df_electron_gamma.head(10))
        # double check gamma

        summed_values = self.df_electron_gamma.groupby(['Event'])["Recoiled/keV"].sum().reset_index()
        print(summed_values.head(20))

        # add gamma up
        self.electron_recoiled_list  = summed_values["Recoiled/keV"].to_list()
        p_observed = []
        for i in range(len(self.electron_recoiled_list)):
            # 40 /keV 0.03 and 0.2 PCE and PDE
            if i> 1E-6:
                p_observed.append(self.electron_recoiled_list[i]*1E6*40*0.03*0.2/(1000))

        num = 0
        for i in p_observed:
            if i >= 1:
                num += 1
        print("photon observed number ", num, len(p_observed))
        print("max", max(p_observed), "\n", "min", min(p_observed))
        # plt.hist(self.electron_recoiled_list, bins=100)
        with open(self.false_1_path, 'w', newline='') as myfile:
            wr = csv.writer(myfile)
            wr.writerow(p_observed)
        plt.hist(p_observed, bins=100)
        plt.xlabel("Obeserved Photon per Event")
        plt.show()
    def single_n_find_gamma_e_loop(self):
        self.df_gamma_rw = pd.read_csv(self.base_path2 + "dmx_single_n_gamma_AmLi_neutron_list_loop.csv")
        print(self.df_gamma_rw[["Kinetic/keV"]].head(20))
        # we need to do severalthings:
        # gamma only in LAr or CF4
        # in 1 event number, only the first series of gammas, avoiding over-countting
        self.gamma_Scint = self.df_gamma_rw[
            ((self.df_gamma_rw['Volume'] == 'LAr_phys') | (self.df_gamma_rw['Volume'] == 'hydraulic_fluid_phys'))&(self.df_gamma_rw['Process']=="compt")]
        gamma_list = self.gamma_Scint["Event"].unique()
        print("gamma filter", len(gamma_list))
        self.gamma_Scint = self.keep_1st(self.gamma_Scint)
        print("scint",self.gamma_Scint)
        # print("scint2",self.gamma_Scint[self.gamma_Scint["Parent ID"]!=1])
        self.gamma_Scint_column = self.gamma_Scint[['Event',"Track ID"]]
        self.gamma_Scint_column.columns = ['Event',"Parent ID"]
        self.df_electron = self.df[(self.df['name']=='e-')&(self.df['Volume']=='LAr_phys')]
        self.df_electron = self.keep_1st(self.df_electron)
        self.df_electron_gamma = pd.merge(self.df_electron,self.gamma_Scint_column,on=['Event','Parent ID'], how='inner')
        print("gamma filter 2",len(self.df_electron_gamma["Event"].unique()))
        print(self.df_electron_gamma.head(10))
        # double check gamma

        summed_values = self.df_electron_gamma.groupby(['Event'])["Recoiled/keV"].sum().reset_index()
        print(summed_values.head(20))

        # add gamma up
        self.electron_recoiled_list  = summed_values["Recoiled/keV"].to_list()
        p_observed = []
        for i in range(len(self.electron_recoiled_list)):
            # 40 /keV 0.03 and 0.2 PCE and PDE
            if i> 1E-6:
                p_observed.append(self.electron_recoiled_list[i]*1E6*40*0.03*0.2/(1000))

        num = 0
        for i in p_observed:
            if i >= 1:
                num += 1
        print("photon observed number ", num, len(p_observed))
        print("max", max(p_observed), "\n", "min", min(p_observed))
        # plt.hist(self.electron_recoiled_list, bins=100)
        with open(self.base_path2+"/photon_capture_n_sing_scatterg_AmLi.csv", 'w', newline='') as myfile:
            wr = csv.writer(myfile)
            wr.writerow(p_observed)
        plt.hist(p_observed, bins=100)
        plt.xlabel("Obeserved Photon per Event")
        plt.show()
    def LAr_find_gamma_e(self):
        self.df_gamma_rw = pd.read_csv(self.signal_path_mid)
        print(self.df_gamma_rw[["Kinetic/keV"]].head(20))
        # we need to do severalthings:
        # gamma only in LAr or CF4
        # in 1 event number, only the first series of gammas, avoiding over-countting

        self.gamma_Scint = self.df_gamma_rw[self.df_gamma_rw['Volume'] == 'LAr_phys']
        print("gamma scint overcount", len(self.gamma_Scint["Event"].unique()),
              self.gamma_Scint["Event"].unique()[:20])
        self.gamma_Scint = self.keep_1st(self.gamma_Scint)
        print("gamma scint", len(self.gamma_Scint["Event"].unique()),
              self.gamma_Scint["Event"].unique()[:20])
        print("scint",self.gamma_Scint)
        # print("scint2",self.gamma_Scint[self.gamma_Scint["Parent ID"]!=1])
        self.gamma_Scint_column = self.gamma_Scint[['Event',"Track ID"]]
        self.gamma_Scint_column.columns = ['Event',"Parent ID"]
        self.df_electron = self.df[(self.df['name']=='e-')&(self.df['Volume']=='LAr_phys')]
        print("electron first", len(self.df_electron["Event"].unique()),
              self.df_electron["Event"].unique()[:20])
        self.df_electron = self.keep_1st(self.df_electron)
        print("electron afterward", len(self.df_electron["Event"].unique()),
              self.df_electron["Event"].unique()[:20])
        self.df_electron_gamma_merged = pd.merge(self.df_electron,self.gamma_Scint_column,on=['Event','Parent ID'], how='inner')
        print("gamma scint merged", len(self.df_electron_gamma_merged["Event"].unique()),
              self.df_electron_gamma_merged["Event"].unique()[:20])
        lost_event = list(set(self.gamma_Scint["Event"].unique()) - set(self.df_electron_gamma_merged["Event"].unique()))
        print("lost", lost_event)
        print(self.df_electron_gamma_merged.head(10))
        # double check gamma

        # summed_values = self.df_electron_gamma_merged.groupby(['Event', 'Parent ID'])["Recoiled/keV"].sum().reset_index()
        summed_values = self.df_electron_gamma_merged.groupby(['Event'])[
            "Recoiled/keV"].sum().reset_index()
        print("summed values", len(summed_values["Event"].unique()),
              summed_values["Event"].unique()[:20])
        print(summed_values.head(20))

        # add gamma up
        self.electron_recoiled_list  = summed_values["Recoiled/keV"].to_list()
        p_observed = []
        print("recoil list", len(self.electron_recoiled_list))
        for i in range(len(self.electron_recoiled_list)):
            # 40 /keV 0.03 and 0.2 PCE and PDE
            p_observed.append(self.electron_recoiled_list[i]*1E6*40*0.03*0.2/(1000))
        print("p observed list len",len(p_observed))

        num = 0
        for i in p_observed:
            if i >= 1:
                num += 1
        print("photon observed number ", num)
        print("max", max(p_observed), "\n", "min", min(p_observed))
        # plt.hist(self.electron_recoiled_list, bins=100)
        plt.hist(p_observed, bins=100)
        print("output len",len(p_observed))
        with open(self.signal_path, 'w', newline='') as myfile:
            wr = csv.writer(myfile)
            wr.writerow(p_observed)
        plt.xlabel("Obeserved Photon per Event")
        plt.show()
    def check_capture(self):
        self.df_gamma_rw = pd.read_csv(self.base_path + "dmx_gamma.csv")

    def plot_gamma(self):
        self.gamma = pd.read_csv(self.base_path +"dmx_gamma.csv")
        # add gamma energy together for same event
        event_p =0
        energy = 0
        energy_p = []
        track_p = []
        for index in range(len(self.gamma.index)):
            print(index)
            # new event, put energy to energy_list and initialize
            if self.gamma.iloc[index]['Event']> event_p:
                event_p = self.gamma.iloc[index]['Event']
                track_p = []
                energy_p.append(energy)
                energy = 0
           # same event,
            else:
                if self.gamma.iloc[index]['Track ID'] not in track_p:
                    track_p.append(self.gamma.iloc[index]['Track ID'])
                    energy += self.gamma.iloc[index]["Kinetic/keV"]*1000000 # to ev Actullay it is Kinetic/MeV

        print("energy", energy_p[:10])

        observed_photon = []
        for i in energy_p:
            observed_photon.append(i*0.03*0.2/100)
        num = 0
        for i in observed_photon:
            if i>1:
                num +=1
        print("photon observed number ", num)
        print("max",max(observed_photon),"\n","min", min(observed_photon))
        plt.hist(observed_photon, bins=100)
        plt.xlabel("Obeserved Photon per Event")
        plt.show()
        """14664 number has photon observation >1 """

    def FN_spectrum_v2(self):

        self.df_Arrecoil = self.df[(self.df["name"]=='Ar36')|(self.df["name"]=='Ar37')|(self.df["name"]=='Ar40')|(self.df["name"]=='Ar41')][['Event','name','Parent ID','Track ID',"Recoiled/keV",'Process']]

        print(self.df_Arrecoil.head(10))
        # capture ar41 and then radiactive decay
        # 279 first elastic scatter and recoil argon and then capture by other volume
        # self.df_test_merge = self.df[(self.df["Event"]==2694)&((self.df["name"]=='Ar36')|(self.df["name"]=='Ar37')|(self.df["name"]=='Ar40')|(self.df["name"]=='Ar41')|(self.df["name"]=='neutron'))]
        # self.df_test_merge.to_csv(self.base_path +"dmx_argon_multi.csv")
        # self.df_cap_gamma = pd.DataFrame('Event','Track ID')
        print("len",len(self.df_Arrecoil.index))
        print("unique",self.df_Arrecoil['Process'].unique())
        # only record gamma event, whose event id same as ncap and parent id is ncap's track id.
        # change Track ID name into Parent ID so that ready for merge
        self.df_Ar_recoil_merge = self.df_Arrecoil[['Event','Parent ID']]
        self.df_Ar_recoil_merge.columns = ['Event','Track ID']
        # select all neutron events
        self.df_neutron = self.df[self.df['name'] == 'neutron']
        print(self.df_neutron.head(10))

        # select neutron events whose Event number is same as argon event and track id is argon's parent ID
        self.df_neutron_mom = pd.merge(self.df_neutron, self.df_Ar_recoil_merge,on=['Event','Track ID'], how='inner')
        print(self.df_neutron_mom.head(20))
        print("neutron unique", self.df_neutron_mom['Process'].unique())
        self.df_neutron_ncap = self.df_neutron_mom[(self.df_neutron_mom['Process']=='nCapture')&(self.df_neutron_mom['Volume']=='LAr_phys')]
        self.df_neutron_ncap = self.keep_1st(self.df_neutron_ncap, ["Event", "Track ID"])
        self.df_neutron_ncap.to_csv(self.base_path +"dmx_argon_ncap.csv")
        print("neutron cap", len(self.df_neutron_ncap.index))
        self.df_neutron_ela = self.df_neutron_mom[(self.df_neutron_mom['Process'] == 'hadElastic')&(self.df_neutron_mom['Volume']=='LAr_phys')]
        # multiple scattering
        # keep first ? no, depend on how much Ar nucleus are recoiled

        self.df_neutron_ela = self.keep_1st(self.df_neutron_ela,["Event", "Track ID"])

        self.df_neutron_ela_tomerge = self.df_neutron_ela[["Event", "Track ID"]]
        self.df_neutron_ela_tomerge.columns = ["Event", "Parent ID"]
        print( "\nlen single + multiple bubble", len(self.df_neutron_ela_tomerge.index))
        # merge back to find the elastic Argon recoiled energy
        self.df_argon_ela_merged = pd.merge(self.df_neutron_ela_tomerge, self.df_Arrecoil,on=['Event','Parent ID'], how='inner')
        # find the first argon recoiled

        self.df_argon_ela_merged = self.keep_1st(self.df_argon_ela_merged,["Event", "Track ID"])

        print(self.df_argon_ela_merged.head(20), "\nlen",len(self.df_argon_ela_merged.index))
        # save these gamma event
        self.df_argon_ela_merged.to_csv(self.base_path +"dmx_argon_elastic.csv", index=False)
    
    def find_multiplicity(self):
        self.df_Arrecoil = self.df[
            (self.df["name"] == 'Ar36') | (self.df["name"] == 'Ar37') | (self.df["name"] == 'Ar40') | (
                        self.df["name"] == 'Ar41')][
            ['Event', 'name', 'Parent ID', 'Track ID', "Recoiled/keV", 'Volume','Process']]

        self.df_neutron = self.df[(self.df['name'] == 'neutron')&(self.df['Volume'] == 'LAr_phys')]
        self.df_n_cap = self.df_neutron[self.df_neutron["Process"]=='nCapture']
        self.df_n_ela = self.df_neutron[self.df_neutron["Process"] == 'hadElastic']

        self.df_n_cap = self.keep_1st(self.df_n_cap,["Event", "Track ID"])
        self.ncap_columns = self.df_n_cap[["Event", "Track ID"]]
        self.ncap_columns.columns = ["Event", "Parent ID"]
        self.ar_ncap = pd.merge(self.ncap_columns, self.df_Arrecoil,on=['Event','Parent ID'], how='inner')
        # possibly contain event like elastic + capture-> only Ar41 or 37
        self.ar_ncap = self.ar_ncap[(self.ar_ncap["name"] == 'Ar37') | (self.ar_ncap["name"] == 'Ar41')]
        self.ar_ncap = self.keep_1st(self.ar_ncap,["Event", "Track ID"])
        print("cap event", len(self.ar_ncap.index), '\n', self.ar_ncap.head(10))



        self.df_n_ela = self.keep_1st(self.df_n_ela, ["Event", "Track ID"])
        self.nela_columns = self.df_n_ela[["Event", "Track ID"]]
        self.nela_columns.columns = ["Event", "Parent ID"]
        self.ar_nela = pd.merge(self.nela_columns, self.df_Arrecoil, on=['Event', 'Parent ID'], how='inner')
        self.ar_nela = self.ar_nela[(self.ar_nela["name"] == 'Ar36') | (self.ar_nela["name"] == 'Ar40')]
        self.ar_nela = self.keep_1st(self.ar_nela, ["Event", "Track ID"])

        # count multiplicity
        (self.sig_df,self.multi_df) = self.find_single_n_multi(self.ar_nela)
        print("sig", len(self.sig_df.index),'\n', self.sig_df.head(10))
        print("multi", len(self.multi_df.index),'\n', self.multi_df.head(10))
        self.multi_clean_df = self.keep_1st(self.multi_df,["Event","Parent ID"])
        print("multi_clean", len(self.multi_clean_df.index),'\n', self.multi_clean_df.head(10))

        self.ela_sig_300 = self.sig_df[self.sig_df["Recoiled/keV"] > 1E-3]
        print("300", len(self.ela_sig_300.index), '\n', self.ela_sig_300.head(10))

        # energy after applying threshold
        (self.sig_df_t, self.multi_df_t) = self.find_single_n_multi(self.ar_nela[self.ar_nela["Recoiled/keV"]>1E-3])
        print("sig_t", len(self.sig_df_t.index), '\n', self.sig_df_t.head(10))
        print("multi_t", len(self.multi_df_t.index), '\n', self.multi_df_t.head(10))
        self.multi_clean_df_t = self.keep_1st(self.multi_df_t, ["Event", "Parent ID"])
        print("multi_clean_t", len(self.multi_clean_df_t.index), '\n', self.multi_clean_df_t.head(10))






    def plot_elastic(self):
        self.df = pd.read_csv(self.base_path +"dmx_argon_elastic.csv")
        energy_list  = self.df["Recoiled/keV"].to_list()
        energy_ev = []
        energy_1kev = []
        energy_10kev = []
        for i in energy_list:
            energy_ev.append(i*1000000)# actually Recoiled/MeV
            if i*1E6> 1000:
                energy_1kev.append(i*1E6)
                if i*1E6> 10000:
                    energy_10kev.append(i*1E6)
        plt.hist(energy_ev, np.logspace(np.log10(min(energy_ev)), np.log10(max(energy_ev)), 50))
        plt.yscale("log")
        plt.xscale("log")
        print("len", len(energy_ev))
        print("1kev",len(energy_1kev))
        print("10kev", len(energy_10kev))
        plt.show()

    def find_single_n_multi(self, df,Event="Event", Parent="Parent ID"):

        # if an entry has same event and parent ID but has different Track ID
        # df = self.keep_1st(df, ['Event', 'Track ID'])

        df['combined_tuple'] = list(zip(df.iloc[:][Event], df.iloc[:][Parent]))
        multi_appearance_mask = df['combined_tuple'].duplicated(keep=False)
        sing_appearance_mask = ~df['combined_tuple'].duplicated(keep=False)
        # find the duplicated
        filtered_df_sing = df[sing_appearance_mask]
        filtered_df_multi = df[multi_appearance_mask]
        filtered_df_sing = filtered_df_sing.drop(columns=['combined_tuple'])
        filtered_df_multi = filtered_df_multi.drop(columns=['combined_tuple'])
        print("multi", filtered_df_multi.head(10))
        print("sing", filtered_df_sing.head(10))
        return (filtered_df_sing,filtered_df_multi)


    def Check_inelastic(self):
        # self.df_event168 =  self.df[self.df["Process"]=="neutronInelastic"]
        self.df_event168 =  self.df[self.df["Event"]==168]
        print(self.df_event168)
        # print(self.df_event168.head(20))
        self.df_event168.to_csv(self.base_path +"event168.csv")
    def keep_1st(self, df, columns=['Event','Track ID']):
        # Assuming df is your DataFrame and column1, column2 are the column names
        df['combined_tuple'] = list(zip(df.iloc[:][columns[0]], df.iloc[:][columns[1]]))
        first_appearance_mask = ~df['combined_tuple'].duplicated(keep='first')
        filtered_df = df[first_appearance_mask]
        filtered_df = filtered_df.drop(columns=['combined_tuple'])
        return filtered_df
    def test_event(self,event, name="dmx_argon_multi"):

        self.df_test_merge = self.df[(self.df["Event"]==event)&((self.df["name"]=='Ar36')|(self.df["name"]=='Ar37')|(self.df["name"]=='Ar40')|(self.df["name"]=='Ar41')|(self.df["name"]=='neutron'))]

        self.df_test_merge.to_csv(self.base_path +name+".csv")


if __name__ =="__main__":
    # ReR = RestructureRoot()
    RR = ReadRoot()