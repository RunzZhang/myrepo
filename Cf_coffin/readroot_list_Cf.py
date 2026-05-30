"""This file is for analyze the root files"""
"""['Transportation', 'hadElastic', 'nCapture', 'conv', 'phot', 'compt', 'Rayl', 'neutronInelastic', 'ionIoni', 
        'NoProcess', 'RadioactiveDecay', 'photonNuclear', 'Decay']
            ['neutron', 'gamma', 'Ar40', 'Ar36', 'Ar41', 'Ar37']
            ['world_phys', 'Vacuum_vessel_phys', 'Inside_vacuum_vessel_phys', 'pressure_vessel_phys', 
            'hydraulic_fluid_phys', 'HDPE_pressure_vessel_phys', 'reflector_Cu_phys', 'reflector_PTFE_phys', 
            'outer_jar_phys', 'SiPM_Holder21_phys', 'reflector_bare_top_PTFE_phys', 'reflector_bare_top_Cu_phys', 'LAr_phys', 
            'Camera_port_phys', 'inner_jar_phys', 'RTD_Cable_Kapton_phys', 'RTD_Cable_Cu_phys', 'Camera_System2_phys', 
            ' Iris_Holder_phys', ' Iris_phys', 'bare_top_Plastic_phys', 'Lens_Holder3_phys', 'Camera_System1_phys', 
            ' Iris_Holder1_phys', 'Lens_Holder4_phys', 'Side_support3_phys', 'bare_top_sf_phys', 'SiPM_Holder31_phys',
             'SiPM_PCB_Out42_phys', 'SiPM_Holder42_phys', 'SiPM_PCB_Inn4_phys', 'SiPM4_Inside_phys', 'SiPM_Holder51_phys', 
             'calibration_port_phys', 'calibration_Be_phys', 'calibration_air_phys', 'RTD_Cable_Cu_1_phys', 'SiPM_PCB_Out31_phys', 
             'SiPM_PCB_Inn3_phys', 'SIPM3_Si_phys', 'bare_top_flange_phys', 'SiPM_Holder18_phys', 'SiPM1_Inside_phys', 'SIPM1_Si_phys', 
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
               'Piezo_Cu_8_phys', 'Guide_Rod1_phys', 'SiPM_PCB_Out21_phys', 'bare_top_Plastic1_phys', 'Plastic_Flange_phys', 
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
        # df_test = tree.arrays(columns, library="pd", entry_start=30633056,entry_sbare_top=30633066)
        df_test = tree.arrays(columns, library="pd")
        # df_test = tree.arrays(["Event", "PreKinetic/MeV"], library="pd",
                              # entry_sbare_top=10)

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
        self.base_path = "/data/runzezhang/result/TN_sims_D/chunked_root_files_Cf_1E7_bare_top_sourcetube_B/"
        self.base_path2 = "/data/runzezhang/result/TN_sims_D/chunked_root_files_Cf_1E7_bare_top_sourcetube_B/"
        self.plot_path = '/data/runzezhang/result/TN_sims_D/plot/'

        # self.filepath = self.base_path +"dmx_lr.root"
        # self.main_body(1)
        # for i in range(1,101):
        for i in range(1, 11):
            self.main_body(i)
    def main_body(self,i):
        print(i)
        self.ini_path = self.base_path+ f"PN_1E7_ini_part{i}.csv"
        self.ar_ke_path = self.base_path+ f"PN_1E7_ke_part{i}.csv"
        self.ini_x_path = self.base_path+ f"PN_1E7_inix_part{i}.csv"
        self.false_1 = f"PN_1E7_false1_part{i}.csv"
        self.false_2 = f"PN_1E7_false2_part{i}.csv"
        self.false_3 = f"PN_1E7_false3_part{i}.csv"
        self.false_1_old = f"PN_1E7_false1_old_part{i}.csv"
        self.false_2_old = f"PN_1E7_false2_old_part{i}.csv"
        self.false_3_old = f"PN_1E7_false3_old_part{i}.csv"
        self.false_1_new = f"PN_1E7_false1_new_part{i}.csv"
        self.false_2_new = f"PN_1E7_false2_new_part{i}.csv"
        self.signal = f"PN_1E7_sig_part{i}.csv"
        self.signal_old = f"PN_1E7_sig_old_part{i}.csv"
        self.signal_new = f"PN_1E7_sig_new_part{i}.csv"
        self.false_1_mid = f"PN_1E7_false1_mid_part{i}.csv"
        self.false_2_mid = f"PN_1E7_false2_mid_part{i}.csv"
        self.false_3_mid = f"PN_1E7_false3_mid_part{i}.csv"
        self.signal_mid = f"PN_1E7_sig_mid_part{i}.csv"
        self.false_1_old_mid = f"PN_1E7_false1_old_mid_part{i}.csv"
        self.false_2_old_mid = f"PN_1E7_false2_old_mid_part{i}.csv"
        self.false_3_old_mid = f"PN_1E7_false3_old_mid_part{i}.csv"
        self.signal_old_mid = f"PN_1E7_sig_old_mid_part{i}.csv"
        self.false_1_new_mid = f"PN_1E7_false1_new_mid_part{i}.csv"
        self.false_2_new_mid = f"PN_1E7_false2_new_mid_part{i}.csv"
        self.signal_new_mid = f"PN_1E7_sig_new_mid_part{i}.csv"
        self.false_1_path = self.base_path + self.false_1
        self.false_2_path = self.base_path + self.false_2
        self.false_3_path = self.base_path + self.false_3
        self.false_1_old_path = self.base_path + self.false_1_old
        self.false_2_old_path = self.base_path + self.false_2_old
        self.false_3_old_path = self.base_path + self.false_3_old
        self.false_1_new_path = self.base_path + self.false_1_new
        self.false_2_new_path = self.base_path + self.false_2_new
        self.false_1_path_mid = self.base_path + self.false_1_mid
        self.false_2_path_mid = self.base_path + self.false_2_mid
        self.false_3_path_mid = self.base_path + self.false_3_mid
        self.false_1_old_path_mid = self.base_path + self.false_1_old_mid
        self.false_2_old_path_mid = self.base_path + self.false_2_old_mid
        self.false_3_old_path_mid = self.base_path + self.false_3_old_mid
        self.false_1_new_path_mid = self.base_path + self.false_1_new_mid
        self.false_2_new_path_mid = self.base_path + self.false_2_new_mid
        self.signal_path_mid = self.base_path + self.signal_mid
        self.signal_path = self.base_path + self.signal
        self.signal_old_path_mid = self.base_path + self.signal_old_mid
        self.signal_old_path = self.base_path + self.signal_old
        self.signal_new_path_mid = self.base_path + self.signal_new_mid
        self.signal_new_path = self.base_path + self.signal_new
        self.geometry_path = self.base_path+f"PN_1E7_geo_part{i}.csv"
        self.phys_path = self.base_path+f"PN_1E7_phys_part{i}.csv"

        self.x_range = [0, 0]
        self.y_range = [0, 0]
        self.z_range = [0, 0]

        self.filepath = self.base_path + f"dmx_PN_1E7_{i}.root"
        self.file = uproot.open(self.filepath)["tree"]
        # print("columns: ", self.file.keys())
        # ['Event', 'name', 'Parent ID', 'Track ID', 'Step ID', 'X/mm', 'Y/mm', 'Z/mm', 'Kinetic/MeV', 'Recoiled/MeV', 'Volume', 'Process']
        self.selected_columns = ["Event", "name", "Parent ID", "Track ID", "Step ID", "X/mm","Y/mm","Z/mm","PreKinetic/MeV","PostKinetic/MeV",
                                 "Recoiled/MeV", "Volume", "Process"]

        self.bubble_threshold = 0.0001 # MeV bubble generate threshold
        self.rows = 1000

        # self.df = self.file.arrays(self.selected_columns, library="pd").head(self.rows)
        # # process data so that it is easier to read
        # first 1000 rows
        # self.df = self.file.arrays(self.selected_columns, library="pd", entry_start=0,entry_sbare_top=10000)
        self.df = self.file.arrays(self.selected_columns, library="pd")
        self.modify_df()


        #find source tube

        # self.source_geometry()

        # check ssl tube effect to neutron spectrum
        self.sstl_phys()


        self.collect_NR()
        # self.check_NR()

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
    def sstl_phys(self):
        volume_list = ["cf_active_phys", "cf_source_phys", "block1_phys", "block2_phys", "block3_phys",
   "block4_phys","block5_phys","block6_phys","block7_phys","block8_phys","block9_phys","block10_phys","block11_phys"]
        # self.phys = self.df[(self.df["name"]=="neutron")&((self.df["Volume"]=="cf_active_phys")|(self.df["Volume"]=="cf_source_phys")|(self.df["Volume"]=="BPE_coffin_phys"))]
        self.phys = self.df[(self.df["name"] == "neutron") &
                    (self.df["Volume"].isin(volume_list))]
        # print(self.phys)

        self.phys.to_csv(self.phys_path, index = False)

    def source_geometry(self):
        # also include initial energy
        self.initial_position = self.df[(self.df["name"]=="neutron")&(self.df["Volume"]=="cf_active_phys")&(self.df['Step ID'] == 1)][["Event","X/mm", "Y/mm","Z/mm", "Volume","PreKinetic/MeV"]]

        print('self.initial_position',self.initial_position)
        # argon first step info
        self.first_argon = self.df[(self.df["name"]=="neutron")&(self.df["Volume"]=="LAr_phys")]
        self.first_argon = self.first_argon.loc[self.first_argon.groupby('Event')['Step ID'].idxmin()]
        self.argon_position = self.first_argon[["Event","X/mm", "Y/mm","Z/mm", "Volume","PreKinetic/MeV","PostKinetic/MeV"]]
        # all argon info
        # self.argon_position = self.df[(self.df["name"]=="neutron")&(self.df["Volume"]=="LAr_phys")][["X/mm", "Y/mm","Z/mm", "Volume","PreKinetic/MeV","PostKinetic/MeV"]]
        self.geometry =  pd.concat([self.initial_position, self.argon_position], axis=0)
        print(self.geometry)
        print(self.initial_position)
        self.geometry.to_csv(self.geometry_path, index = False)

    def collect_NR(self):
        self.NR_scatter= self.df[(self.df["name"]=="neutron")&(self.df["Volume"]=="LAr_phys")&(self.df['Process'].isin(['hadElastic', 'neutronInelastic'])) ]
        self.NR_capture = self.df[(self.df["name"]=="neutron")&(self.df["Volume"]=="LAr_phys")&(self.df['Process'].isin(['nCapture'])) ]

        self.LAr = self.df[(self.df["name"].isin(["Ar36","Ar38", "Ar40"]))&(self.df["Recoiled/MeV"]>0)]
        self.LAr = self.keep_1st(self.LAr)
        self.LAr = self.LAr.drop(columns=["Process"])

        # collect Ar recoiled by Elastic and inelastic
        self.NR_scatter =self.keep_1st(self.NR_scatter)
        self.NR_scatter_column = self.NR_scatter[['Event', "Track ID","Process"]]
        self.NR_scatter_column.columns = ['Event', "Parent ID", "Process"]
        self.LAr_scatter = pd.merge(self.LAr, self.NR_scatter_column, on=['Event', 'Parent ID'],
                                          how='inner')
        print(self.LAr_scatter)
        # duplicates = self.LAr_scatter[self.LAr_scatter.duplicated(subset='Event',keep=False)]
        # print("duplicate",duplicates)
        self.LAr_scatter = self.LAr_scatter[["Event","Parent ID", "Process","Recoiled/MeV"]]

        #capture
        self.LAr_capture = self.keep_1st(self.NR_capture)
        self.LAr_capture = self.LAr_capture[["Event","Parent ID", "Process","Recoiled/MeV"]]


        self.NR = pd.concat([self.LAr_scatter,self.LAr_capture], axis=0)
        print('self.NR',self.NR)
        # for capture, it is determined by the recoil spectrum
        self.NR.to_csv(self.signal_path, index= False)
    def check_NR(self):
        self.event206 = self.df[self.df["Event"]==206]
        self.event206.to_csv(self.base_path+"event206.csv", index= False)

        self.inelastic= self.df[(self.df['Process'].isin(['neutronInelastic'])) ]
        inelastic_list =self.inelastic["Event"].to_list()
        print(inelastic_list)
        self.inelastic_1st = self.df[self.df["Event"]==inelastic_list[0]]
        self.inelastic_1st.to_csv(self.base_path+"eventinelastic.csv", index= False)

    def keep_1st(self, df, columns=['Event','Track ID']):
        # Assuming df is your DataFrame and column1, column2 are the column names
        df['combined_tuple'] = list(zip(df.iloc[:][columns[0]], df.iloc[:][columns[1]]))
        first_appearance_mask = ~df['combined_tuple'].duplicated(keep='first')
        filtered_df = df[first_appearance_mask]
        filtered_df = filtered_df.drop(columns=['combined_tuple'])
        return filtered_df

if __name__ =="__main__":
    # ReR = RestructureRoot()
    RR = ReadRoot()
    # test_write()

    # find corrupted file entries
    # num= find_entries()
    # print(num)