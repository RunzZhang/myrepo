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
# filename = "/data/runzezhang/Geant4Simulaions/g411_TN/dmx.root"

class RestructureRoot():
    def __init__(self):
        self.filepath = "/data/runzezhang/result/TN_sims/dmx.root"
        self.reconstruct_filepath = "/data/runzezhang/result/TN_sims/dmx_rc.csv"
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
        self.filepath = "/data/runzezhang/result/TN_sims/dmx.root"
        self.file = uproot.open(self.filepath)["tree"]
        print("columns: ",self.file.keys())
        #['Event', 'name', 'Parent ID', 'Track ID', 'Step ID', 'X/mm', 'Y/mm', 'Z/mm', 'Kinetic/keV', 'Recoiled/keV', 'Volume', 'Process']
        self.selected_columns = ["Event","name","Parent ID","Track ID","Step ID","X/mm","Kinetic/keV","Recoiled/keV", "Volume","Process"]
        self.rows = 1000

        # self.df = self.file.arrays(self.selected_columns, library="pd").head(self.rows)
        # # process data so that it is easier to read
        self.df = self.file.arrays(self.selected_columns, library="pd")
        self.modify_df()

        # self.gamma_event()


        # self.FN_spectrum_v2()
        self.plot_elastic()
        # self.Check_inelastic()

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
        print("end", event_number[:100])

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

    def Capture_spectrum(self):

        self.df_Ncapture = self.df[(self.df["name"]=='neutron')&(self.df["Process"]=='nCapture')&(self.df["Volume"]!='LAr_phys')][['Event','Track ID']]

        print(self.df_Ncapture.head(10))
        # self.df_Gamma = pd.DataFrame('Event','Track ID')
        print("len",len(self.df_Ncapture.index))
        # only record gamma event, whose event id same as ncap and parent id is ncap's track id.
        # change Track ID name into Parent ID so that ready for merge
        self.df_Ncapture.columns = ['Event','Parent ID']
        # select all gamma events
        self.df_Gamma = self.df[self.df['name'] == 'gamma' ]
        # select gamma events whose Event number is same as neutron event and parent id is neutron's track ID
        self.df_Gamma = pd.merge(self.df_Ncapture, self.df_Gamma,on=['Event','Parent ID'], how='inner')
        print(self.df_Gamma.head(20))
        # save these gamma event
        self.df_Gamma.to_csv("/data/runzezhang/result/TN_sims/dmx_gamma.csv", index=False)

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
        self.Gamma_spectrum()
        self.plot_gamma()
    def Gamma_spectrum(self):
        self.df_gamma_rw = pd.read_csv("/data/runzezhang/result/TN_sims/dmx_gamma.csv")
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
        self.gamma_Scint.to_csv("/data/runzezhang/result/TN_sims/gamma_scint2.csv", index=False)

    def plot_gamma(self):
        self.gamma = pd.read_csv("/data/runzezhang/result/TN_sims/gamma_scint2.csv")
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
        plt.hist(observed_photon, bins=100)
        plt.show()
        """14664 number has photon observation >1 """

    def FN_spectrum(self):

        self.df_Arrecoil = self.df[(self.df["name"]=='Ar36')|(self.df["name"]=='Ar37')|(self.df["name"]=='Ar40')|(self.df["name"]=='Ar41')][['Event','name','Parent ID','Track ID','Process']]

        print(self.df_Arrecoil.head(10))
        # self.df_Gamma = pd.DataFrame('Event','Track ID')
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
        self.df_neutron_ncap = self.df_neutron_mom[self.df_neutron_mom['Process']=='nCapture']
        self.df_neutron_ela = self.df_neutron_mom[self.df_neutron_mom['Process'] == 'hadElastic']
        # find unique event, track ID for elastic( neutron capture only once)

        event_p = 0
        track_p = []
        parent_p = []
        for index in range(len(self.df_neutron_ela.index)):
            print(index)
            if self.df_neutron_ela.iloc[index]['Event']> event_p:
                event_p = self.df_neutron_ela.iloc[index]['Event']
                track_p = []
                parent_p = []
            # same event

            # record the 1st track ID in same trajactory
            else:
                if self.df_neutron_ela.iloc[index]['Track ID'] not in track_p:
                    track_p.append(self.df_neutron_ela.iloc[index]['Track ID'])
                elif self.df_neutron_ela.iloc[index]['Track ID'] in track_p:
                    self.df_neutron_ela.drop(self.df_neutron_ela.index[index])


        self.df_neutron_ela_tomerge = self.df_neutron_ela[["Event", "Track ID"]]
        self.df_neutron_ela_tomerge.columns = ["Event", "Parent ID"]
        # merge back to find the elastic Argon recoiled energy
        self.df_argon_ela_merged = pd.merge(self.df_neutron_ela_tomerge, self.df_Arrecoil,on=['Event','Parent ID'], how='inner')
        # find the first argon recoiled
        event_p = 0
        track_p = []
        parent_p = []
        for index in range(len(self.df_argon_ela_merged.index)):
            print(index)
            if self.df_argon_ela_merged.iloc[index]['Event']> event_p:
                event_p = self.df_argon_ela_merged.iloc[index]['Event']
                track_p = []
                parent_p = []
            # same event

            # record the 1st track ID in same trajactory
            else:
                if self.df_argon_ela_merged.iloc[index]['Track ID'] not in track_p:
                    track_p.append(self.df_argon_ela_merged.iloc[index]['Track ID'])
                elif self.df_argon_ela_merged.iloc[index]['Track ID'] in track_p:
                    self.df_argon_ela_merged.drop(self.df_neutron_ela_tomerge.index[index])

        print(self.df_argon_ela_merged.head(20))
        # save these gamma event
        self.df_argon_ela_merged.to_csv("/data/runzezhang/result/TN_sims/dmx_argon_elastic.csv", index=False)



    def FN_spectrum_v2(self):

        self.df_Arrecoil = self.df[(self.df["name"]=='Ar36')|(self.df["name"]=='Ar37')|(self.df["name"]=='Ar40')|(self.df["name"]=='Ar41')][['Event','name','Parent ID','Track ID',"Recoiled/keV",'Process']]

        print(self.df_Arrecoil.head(10))
        # self.df_Gamma = pd.DataFrame('Event','Track ID')
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
        self.df_neutron_ncap = self.df_neutron_mom[self.df_neutron_mom['Process']=='nCapture']
        self.df_neutron_ela = self.df_neutron_mom[self.df_neutron_mom['Process'] == 'hadElastic']

        self.df_neutron_ela = self.keep_1st(self.df_neutron_ela,["Event", "Track ID"])

        self.df_neutron_ela_tomerge = self.df_neutron_ela[["Event", "Track ID"]]
        self.df_neutron_ela_tomerge.columns = ["Event", "Parent ID"]
        print( "\nlen", len(self.df_neutron_ela_tomerge.index))
        # merge back to find the elastic Argon recoiled energy
        self.df_argon_ela_merged = pd.merge(self.df_neutron_ela_tomerge, self.df_Arrecoil,on=['Event','Parent ID'], how='inner')
        # find the first argon recoiled

        self.df_argon_ela_merged = self.keep_1st(self.df_argon_ela_merged,["Event", "Track ID"])

        print(self.df_argon_ela_merged.head(20), "\nlen",len(self.df_argon_ela_merged.index))
        # save these gamma event
        self.df_argon_ela_merged.to_csv("/data/runzezhang/result/TN_sims/dmx_argon_elastic.csv", index=False)
    def plot_elastic(self):
        self.df = pd.read_csv("/data/runzezhang/result/TN_sims/dmx_argon_elastic.csv")
        energy_list  = self.df["Recoiled/keV"].to_list()
        energy_ev = []
        for i in energy_list:
            energy_ev.append(i*1000)# actually Recoiled/MeV
        plt.hist(energy_ev, bins=100)
        print("len", len(energy_ev))
        plt.show()

    def Check_inelastic(self):
        # self.df_event168 =  self.df[self.df["Process"]=="neutronInelastic"]
        self.df_event168 =  self.df[self.df["Event"]==168]
        print(self.df_event168)
        # print(self.df_event168.head(20))
        self.df_event168.to_csv("/data/runzezhang/result/TN_sims/event168.csv")
    def keep_1st(self, df, columns):
        # Assuming df is your DataFrame and column1, column2 are the column names
        df['combined_tuple'] = list(zip(df.iloc[:][columns[0]], df.iloc[:][columns[1]]))
        first_appearance_mask = ~df['combined_tuple'].duplicated(keep='first')
        filtered_df = df[first_appearance_mask]
        filtered_df = filtered_df.drop(columns=['combined_tuple'])
        return filtered_df


if __name__ =="__main__":
    # ReR = RestructureRoot()
    RR = ReadRoot()