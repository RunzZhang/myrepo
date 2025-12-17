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
def test_write():
    try:
        f = uproot.open("/data/runzezhang/result/TN_sims_D/chunked_root_files_test/dmx_Cf_1E7_2.root")
        tree = f["tree"]  # Or whatever your tree name is
        print("Successfully opened the file!")
        columns= ["Event","name","Parent ID","Track ID","Step ID","X/mm","PreKinetic/MeV","Recoiled/MeV","Volume","Process"]
        # columns = ["Event", "PreKinetic/MeV","Recoiled/MeV","Process"]
        # Optional: Try to read a few entries to confirm data is there
        # df_test = tree.arrays(columns, library="pd", entry_start=30633056,entry_stop=30633066)
        df_test = tree.arrays(columns, library="pd")
        # df_test = tree.arrays(["Event", "PreKinetic/MeV"], library="pd",
                              # entry_stop=10)

        print("First 10 entries:", df_test)
    except Exception as e:
        print(f"Error opening file: {e}")

def find_entries():
    """
        Opens a ROOT file and returns the number of entries in a specified TTree.
        """
    filepath= "/data/runzezhang/result/TN_sims_D/chunked_root_files_corrupted/dmx_Cf_1E7_1.root"
    tree_name="tree"
    try:
        with uproot.open(filepath) as file:
            if tree_name in file:
                tree = file[tree_name]
                return tree.num_entries
            else:
                print(f"Tree '{tree_name}' not found in {filepath}")
                return None
    except Exception as e:
        print(f"Error opening or reading {filepath}: {e}")
        return None


class RestructureRoot():
    def __init__(self):
        self.filepath = "/data/runzezhang/result/TN_sims_D/chunked_root_files/dmx_Cf_1E7.root"
        self.reconstruct_filepath = "/data/runzezhang/result/TN_sims_D/chunked_root_files/dmx_rcCf_1E7.csv"
        self.file = uproot.open(self.filepath)["tree"]
        print("columns: ",self.file.keys())
        #['Event', 'name', 'Parent ID', 'Track ID', 'Step ID', 'X/mm', 'Y/mm', 'Z/mm', 'Kinetic/MeV', 'Recoiled/MeV', 'Volume', 'Process']
        self.selected_columns = ["Event","name","Parent ID","Track ID","Step ID","X/mm","PreKinetic/MeV","Recoiled/MeV","Volume","Process"]
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
        self.base_path = "/data/runzezhang/result/TN_sims_D/chunked_root_files_Cs_1E5/"
        self.base_path2 = "/data/runzezhang/result/TN_sims_D/chunked_root_files_Cs_1E5/"
        self.plot_path = '/data/runzezhang/result/TN_sims_D/plot/'

        # self.filepath = self.base_path +"dmx_lr.root"
        self.main_body(1)
        # for i in range(1,101):
        # for i in range(1, 11):
        #     self.main_body(i)
    def main_body(self,i):
        print(i)
        self.ini_path = self.base_path+ f"Cs_gamma_1E7_ini_part{i}.csv"
        self.ar_ke_path = self.base_path+ f"Cs_gamma_1E7_ke_part{i}.csv"
        self.false_gamma_1 = f"Cs_gamma_1E7_false1_part{i}.csv"
        self.false_gamma_2 = f"Cs_gamma_1E7_false2_part{i}.csv"
        self.false_gamma_3 = f"Cs_gamma_1E7_false3_part{i}.csv"
        self.false_gamma_1_old = f"Cs_gamma_1E7_false1_old_part{i}.csv"
        self.false_gamma_2_old = f"Cs_gamma_1E7_false2_old_part{i}.csv"
        self.false_gamma_3_old = f"Cs_gamma_1E7_false3_old_part{i}.csv"
        self.false_gamma_1_new = f"Cs_gamma_1E7_false1_new_part{i}.csv"
        self.false_gamma_2_new = f"Cs_gamma_1E7_false2_new_part{i}.csv"
        self.signal = f"Cs_gamma_1E7_sig_part{i}.csv"
        self.signal_old = f"Cs_gamma_1E7_sig_old_part{i}.csv"
        self.signal_new = f"Cs_gamma_1E7_sig_new_part{i}.csv"
        self.false_gamma_1_mid = f"Cs_gamma_1E7_false1_mid_part{i}.csv"
        self.false_gamma_2_mid = f"Cs_gamma_1E7_false2_mid_part{i}.csv"
        self.false_gamma_3_mid = f"Cs_gamma_1E7_false3_mid_part{i}.csv"
        self.signal_mid = f"Cs_gamma_1E7_sig_mid_part{i}.csv"
        self.false_gamma_1_old_mid = f"Cs_gamma_1E7_false1_old_mid_part{i}.csv"
        self.false_gamma_2_old_mid = f"Cs_gamma_1E7_false2_old_mid_part{i}.csv"
        self.false_gamma_3_old_mid = f"Cs_gamma_1E7_false3_old_mid_part{i}.csv"
        self.signal_old_mid = f"Cs_gamma_1E7_sig_old_mid_part{i}.csv"
        self.false_gamma_1_new_mid = f"Cs_gamma_1E7_false1_new_mid_part{i}.csv"
        self.false_gamma_2_new_mid = f"Cs_gamma_1E7_false2_new_mid_part{i}.csv"
        self.signal_new_mid = f"Cs_gamma_1E7_sig_new_mid_part{i}.csv"
        self.info_path = self.base_path+ f"Cs_gamma_1E6_info_scube_part{i}.csv"
        self.false_gamma_1_path = self.base_path + self.false_gamma_1
        self.false_gamma_2_path = self.base_path + self.false_gamma_2
        self.false_gamma_3_path = self.base_path + self.false_gamma_3
        self.false_gamma_1_old_path = self.base_path + self.false_gamma_1_old
        self.false_gamma_2_old_path = self.base_path + self.false_gamma_2_old
        self.false_gamma_3_old_path = self.base_path + self.false_gamma_3_old
        self.false_gamma_1_new_path = self.base_path + self.false_gamma_1_new
        self.false_gamma_2_new_path = self.base_path + self.false_gamma_2_new
        self.false_gamma_1_path_mid = self.base_path + self.false_gamma_1_mid
        self.false_gamma_2_path_mid = self.base_path + self.false_gamma_2_mid
        self.false_gamma_3_path_mid = self.base_path + self.false_gamma_3_mid
        self.false_gamma_1_old_path_mid = self.base_path + self.false_gamma_1_old_mid
        self.false_gamma_2_old_path_mid = self.base_path + self.false_gamma_2_old_mid
        self.false_gamma_3_old_path_mid = self.base_path + self.false_gamma_3_old_mid
        self.false_gamma_1_new_path_mid = self.base_path + self.false_gamma_1_new_mid
        self.false_gamma_2_new_path_mid = self.base_path + self.false_gamma_2_new_mid
        self.signal_path_mid = self.base_path + self.signal_mid
        self.signal_path = self.base_path + self.signal
        self.signal_old_path_mid = self.base_path + self.signal_old_mid
        self.signal_old_path = self.base_path + self.signal_old
        self.signal_new_path_mid = self.base_path + self.signal_new_mid
        self.signal_new_path = self.base_path + self.signal_new

        self.x_range = [0, 0]
        self.y_range = [0, 0]
        self.z_range = [0, 0]

        self.filepath = self.base_path + f"dmx_Cs_1E7_{i}.root"
        self.file = uproot.open(self.filepath)["tree"]
        # print("columns: ", self.file.keys())
        # ['Event', 'name', 'Parent ID', 'Track ID', 'Step ID', 'X/mm', 'Y/mm', 'Z/mm', 'Kinetic/MeV', 'Recoiled/MeV', 'Volume', 'Process']
        self.selected_columns = ["Event", "name", "Parent ID", "Track ID", "Step ID", "X/mm",'Y/mm', 'Z/mm',"PreKinetic/MeV","PostKinetic/MeV",
                                 "Recoiled/MeV", "Volume", "Process"]

        self.bubble_threshold = 0.0001 # MeV bubble generate threshold
        self.rows = 1000

        # self.df = self.file.arrays(self.selected_columns, library="pd").head(self.rows)
        # # process data so that it is easier to read
        # first 1000 rows
        # self.df = self.file.arrays(self.selected_columns, library="pd", entry_start=0,entry_stop=10000)
        self.df = self.file.arrays(self.selected_columns, library="pd")
        self.modify_df()

        # find all ER and save ER into csv
        self.allER()
        self.ER_distribution()
        # self.gamma_ER()

        # find all NR
        # self.allNR()




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
        self.df.loc[self.df['Process']=="","Process"]="init"
        # this make event number correct
        self.reidx_event()

    def modify_gamma_process(self):
        print("modify gamma, make gamma volume to be pre-volume consistently")
        # I prefer to change the structure after filter out all events
        # all info is pre step but the process(post step)
        # so move line i process to line +1 in the same particle(track ID)

    def source_geometry(self):
        print(self.df[(self.df["name"]=="neutron")&(self.df["Step ID"]==1)&(self.df["Parent ID"]==0)]["X/mm"])
        self.x_range[0] = min(self.x_range[0],self.df[(self.df["name"]=="neutron")&(self.df["Step ID"]==1)&(self.df["Parent ID"]==0)]["X/mm"].min())
        self.x_range[1] = max(self.x_range[1], self.df[(self.df["name"]=="neutron")&(self.df["Step ID"]==1)&(self.df["Parent ID"]==0)]["X/mm"].max())
        # self.y_range[0] = min(self.y_range[0], self.df["Y/mm"].min())
        # self.y_range[1] = max(self.y_range[1], self.df["Y/mm"].max())
        # self.z_range[0] = min(self.z_range[0], self.df["Z/mm"].min())
        # self.z_range[1] = max(self.z_range[1], self.df["Z/mm"].max())
        print("x",self.x_range)
        print("y", self.y_range)
        print("z", self.z_range)

    def allNR(self):
        self.LAr_recoiled = self.df[((self.df["name"] == 'Ar40') | (self.df["name"] == 'Ar36'))&(self.df["Recoiled/MeV"]>0) ]
        self.LAr_recoiled_event_list = self.LAr_recoiled["Event"].unique().tolist()

        print(self.LAr_recoiled_event_list)

        # sample the first 10 event to see what caused the NR
        # self.LAr_NR_sample = self.df[(self.df["Event"].isin(self.LAr_recoiled_event_list))]
        # self.LAr_NR_sample.to_csv(self.base_path+"LAr_NR_sample.csv")

    def allER(self):
        # we need to do several things:
        # gamma only in LAr or CF4
        # in 1 event number, only the first series of gammas, avoiding over-countting
        # including phot - Xray excited the electrons
        # fetch the location/multiplicity/ER
        # output matrix to [event1_data, event2_data]
        # in each eventi_data, it is matrix [(x1,y1,z1),ER1], [(x2,y2,z2),ER2,,,] including the bubble multiplicity and ER and position

        print(self.df['name'].unique())
        self.gamma_Scint = self.df[
            ((self.df['Volume'] == 'LAr_phys')|(self.df['Volume'] == 'hydraulic_fluid_phys')) &( (self.df['Process'] == "compt")|(self.df['Process'] == "phot"))& (self.df['Parent ID'] == 0)]
        # record the positions and multiplicity

        gamma_list = self.gamma_Scint["Event"].unique()
        print("gamma filter", len(gamma_list))
        # self.gamma_Scint = self.keep_1st(self.gamma_Scint)
        print("scint", self.gamma_Scint)
        # find electrons are daughter of those gammas
        self.gamma_Scint_column = self.gamma_Scint[['Event', "Track ID"]]
        self.gamma_Scint_column.columns = ['Event', "Parent ID"]
        self.df_electron = self.df[(self.df['name'] == 'e-') & ((self.df['Volume'] == 'LAr_phys')|(self.df['Volume'] == 'hydraulic_fluid_phys'))]
        self.df_electron = self.keep_1st(self.df_electron)
        self.df_electron_gamma = pd.merge(self.df_electron, self.gamma_Scint_column, on=['Event', 'Parent ID'],
                                          how='inner')
        print("electron gamma",self.df_electron_gamma.head(20)) #ok

    def ER_distribution(self):



        self.tagged_gamma = self.df_electron[(self.df_electron["name"] == "e-") & (self.df_electron["Event"] != 1)]
        # double check gamma

        summed_values = self.tagged_gamma.groupby(['Event'])["Recoiled/MeV"].sum().reset_index()

        print(summed_values.head(20))

        # add gamma up
        self.electron_recoiled_list = summed_values["Recoiled/MeV"].to_list()
        # save info

        self.electron_recoiled_event_list = summed_values["Event"].to_list()
        print(self.electron_recoiled_event_list[:3])
        high_NRER = []

        self.mom_gamma = self.df[
            ((self.df['Volume'] == 'LAr_phys')|(self.df['Volume'] == 'hydraulic_fluid_phys')) &( (self.df['Process'] == "compt")|(self.df['Process'] == "phot"))& (self.df['Parent ID'] == 0)& (self.df["Event"].isin(self.electron_recoiled_event_list))]
        self.mom_gamma_group = self.mom_gamma.groupby("Event")


        self.kid_e = self.df[(self.df['name'] == 'e-') &(self.df['Step ID'] == 1)& (self.df['Parent ID'] == 1)& (self.df["Event"].isin(self.electron_recoiled_event_list))&((self.df['Volume'] == 'LAr_phys')|(self.df['Volume'] == 'hydraulic_fluid_phys'))]
        self.kid_e_group = self.kid_e.groupby("Event")

        self.mom_gamma = self.mom_gamma.copy()
        self.mom_gamma["ER_near"] = 0.0
        half = 0.05 # mm from original cube 2*2*2 mm


        for event_id, g_evt in self.mom_gamma_group:
            # electrons for same event
            try:
                e_evt = self.kid_e_group.get_group(event_id)
            except KeyError:
                continue  # no e- in this event

            if e_evt.empty:
                continue

            # positions
            g_pos = g_evt[["X/mm", "Y/mm", "Z/mm"]].to_numpy()
            e_pos = e_evt[["X/mm", "Y/mm", "Z/mm"]].to_numpy()

            # recoil energy column name: change if yours is different
            e_E = e_evt["Recoiled/MeV"].to_numpy()

            # cube cut (vectorized): inside shape = (N_gamma, N_e)
            dx = np.abs(g_pos[:, None, 0] - e_pos[None, :, 0]) <= half
            dy = np.abs(g_pos[:, None, 1] - e_pos[None, :, 1]) <= half
            dz = np.abs(g_pos[:, None, 2] - e_pos[None, :, 2]) <= half
            inside = dx & dy & dz

            # sum E per gamma point
            ER_near = inside @ e_E  # (N_gamma,)

            # write back aligned to the same rows in mom_gamma
            self.mom_gamma.loc[g_evt.index, "ER_near"] = ER_near

        first3_events = self.mom_gamma["Event"].unique()[:3]
        print(first3_events)
        print(self.mom_gamma[self.mom_gamma["Event"].isin(first3_events)])

        self.mom_gamma["Multiplicity"] = (self.mom_gamma.groupby("Event")["Step ID"].rank(method="dense", ascending=True).astype(int))
        self.mom_gamma["R/mm"] = np.sqrt(self.mom_gamma["X/mm"]**2+self.mom_gamma["Y/mm"]**2 )

        self.mom_gamma["ER_near/eV"]= self.mom_gamma["ER_near"]*1e6

        first3_events = self.mom_gamma["Event"].unique()[:3]
        print(first3_events)
        print(self.mom_gamma[self.mom_gamma["Event"].isin(first3_events)])

        self.output_df = self.mom_gamma[["Event","name", "R/mm", "Z/mm", "Volume", "Process", "ER_near/eV", "Multiplicity"]]

        self.output_df.to_csv(self.info_path, index=False)

        self.df[self.df["Event"] == 3154].to_csv(self.base_path+"LAr_ER_abnormalER.csv")
        # self.df[(self.df["Event"].isin(self.electron_recoiled_event_list))].to_csv(self.base_path+"LAr_ER_sample_preprocess_laststep_v2.csv")



    def gamma_ER(self):

        self.tagged_gamma = self.df_electron[(self.df_electron["name"] == "e-")&(self.df_electron["Event"] != 1)]
        # double check gamma

        summed_values = self.tagged_gamma.groupby(['Event'])["Recoiled/MeV"].sum().reset_index()

        print(summed_values.head(20))

        # add gamma up
        self.electron_recoiled_list = summed_values["Recoiled/MeV"].to_list()
        # save info

        self.electron_recoiled_event_list = summed_values["Event"].to_list()
        high_NRER = []

        p_observed = []  # the first digit is always the NR number
        for i in range(len(self.electron_recoiled_list)):
            # 40 /MeV 0.03 and 0.2 PCE and PDE
            if self.electron_recoiled_list[i] > 1E-6:
                p_observed.append(self.electron_recoiled_list[i] * 1E6 * 40 * 0.03 * 0.2 / (1000))
            if self.electron_recoiled_list[i] > 350 * 1000 / (1E6 * 40 * 0.03 * 0.2):
                high_NRER.append(self.electron_recoiled_event_list[i])

        print("max", max(p_observed), "\n", "min", min(p_observed))
        print("exclusive ER", high_NRER)
        # p_observe only contains ER, if one event only has NR, it still produce bubbles that we need to compress
        print("path",self.false_gamma_1_path)
        # Cs, just save the total ER in MeV
        with open(self.false_gamma_1_path, 'w', newline='') as myfile:
            wr = csv.writer(myfile)
            wr.writerow(self.electron_recoiled_list)

        # self.df[(self.df["Event"].isin(self.electron_recoiled_event_list))].to_csv(self.base_path+"LAr_ER_sample_preprocess_laststep.csv")

    def exclude_common(self,df1, df2): # exclude same ["Event"]
        common_events = set(df1["Event"]) & set(df2["Event"])

        # 2. Drop common events from both
        df1_clean = df1[~df1["Event"].isin(common_events)]
        df2_clean = df2[~df2["Event"].isin(common_events)]

        # 3. Combine
        df_combined = pd.concat([df1_clean, df2_clean])

        # 4. Sort by Event
        df_combined = df_combined.sort_values("Event").reset_index()
        return df_combined
    def intersection(self, lst1, lst2):
        lst3 = [value for value in lst1 if value in lst2]
        return lst3


    def test_merge(self):
        df_a = pd.DataFrame({ 'B':[3,5],'C':[5,7]})
        df_b = pd.DataFrame({'B': [2,3,5,3], 'C': [3,5,7,5], 'F': [5,9,10,6], 'G': [8,10,7,9]})
        print('a\n',df_a)
        print('b\n',df_b)
        merged_df = pd.merge(df_b, df_a, on=['B','C'], how='inner')

        print(merged_df)






    def plot_elastic(self):
        self.df = pd.read_csv(self.base_path +"dmx_argon_elastic.csv")
        energy_list  = self.df["Recoiled/MeV"].to_list()
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

    def distinguish_single_all(self, df, column1="Event", column2="Track ID"):

        df_LAr_bubble_NR = self.df[
            ((self.df["name"] == 'Ar40') | (self.df["name"] == 'Ar36')) & (
                    self.df["Volume"] == 'LAr_phys') & (self.df["Recoiled/MeV"] >= self.bubble_threshold)][
            ['Event', 'Volume', 'Track ID',
             'Parent ID']]  # this to judgge whether multi bubbles. if two different TrackID > bubble threhosld in one event
        df_LAr_bubble_NR = self.keep_1st(df_LAr_bubble_NR)

        # Count distinct Track IDs per Event
        track_counts = df_LAr_bubble_NR.groupby("Event")["Track ID"].nunique()

        # 1. Events with multiple Track IDs
        multi_track_events = track_counts[track_counts > 1].index
        df_multi = df_LAr_bubble_NR[df_LAr_bubble_NR["Event"].isin(multi_track_events)]

        # 2. Events with exactly 1 Track ID
        single_track_events = track_counts[track_counts == 1].index
        df_single = df_LAr_bubble_NR[df_LAr_bubble_NR["Event"].isin(single_track_events)]

        single_list = df_single["Event"].unique().tolist()
        multi_list = df_multi["Event"].unique().tolist()



        # df['combined_tuple'] = list(zip(df_LAr_bubble_NR.iloc[:][column1], df_LAr_bubble_NR.iloc[:][column2]))
        # multi_appearance_mask = df['combined_tuple'].duplicated(keep=False)
        # sing_appearance_mask = ~df['combined_tuple'].duplicated(keep=False)
        # # find the duplicated
        # filtered_df_sing = df[sing_appearance_mask]
        # filtered_df_multi = df[multi_appearance_mask]
        # filtered_df_sing = filtered_df_sing.drop(columns=['combined_tuple'])
        # filtered_df_multi = filtered_df_multi.drop(columns=['combined_tuple'])
        # print("multi", filtered_df_multi.head(10))
        # print("sing", filtered_df_sing.head(10))
        # return (filtered_df_sing, filtered_df_multi)
        return(df[df["Event"].isin(single_list)],df[df["Event"].isin(multi_list)])


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
    # test_write()

    # find corrupted file entries
    # num= find_entries()
    # print(num)