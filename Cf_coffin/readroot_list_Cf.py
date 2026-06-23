"""This file is for analyze the root files"""

import pandas as pd
import uproot
import matplotlib.pyplot as plt
import numpy as np
import csv
import sys
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
        """A. Here are chunk file repository I build from chunk python code. ANd the loop below number should match the
        number of chunked files.  By default, the  path and loop number are only ones you need to change"""
        self.base_path = "/data/runzezhang/result/TN_sims_D/chunked_root_files_Cf_1E5_100ppm_sourcetube_B/"
        self.plot_path = '/data/runzezhang/result/TN_sims_D/plot/'


        # self.main_body(1)
        for i in range(1,101):
        # for i in range(1, 11):
            self.main_body(i)
    def main_body(self,i):
        """B. Output files names for next python to read. PN originally means photo-neutron but I forgot to change it when
        switch to Cf simulations.... You don't need to change these as long as it is consistent with SN_plotCf_both.py"""
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
        self.selected_columns = ["Event", "name", "Parent ID", "Track ID", "Step ID", "X/mm","Y/mm","Z/mm","PreKinetic/MeV","PostKinetic/MeV",
                                 "Recoiled/MeV", "Volume", "Process"]

        self.bubble_threshold = 0.0001 # MeV bubble generate threshold
        self.rows = 1000


        self.df = self.file.arrays(self.selected_columns, library="pd")
        self.modify_df()


        """C. there are 3 functions: 
        source_phys is get all neutron steps that passing through Coffin volumes
        first_argon is to get neutron in the active volume and also neutrons that first entering argon volume
        collet_NR is to get all neutron interaction in the argon volumes but only save argon status. """

        #find source tube

        self.first_argon()

        # check ssl tube effect to neutron spectrum
        # self.source_phys()


        # self.collect_NR()
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
    def source_phys(self):
        volume_list = ["cf_active_phys", "cf_source_phys", "block1_phys", "block2_phys", "block3_phys",
   "block4_phys","block5_phys","block6_phys","block7_phys","block8_phys","block9_phys","block10_phys","block11_phys"]
        # self.phys = self.df[(self.df["name"]=="neutron")&((self.df["Volume"]=="cf_active_phys")|(self.df["Volume"]=="cf_source_phys")|(self.df["Volume"]=="BPE_coffin_phys"))]
        self.phys = self.df[(self.df["name"] == "neutron") &
                    (self.df["Volume"].isin(volume_list))]
        self.phys.to_csv(self.phys_path, index = False)

    def first_argon(self):
        # check particles
        print("particle name",self.df["name"].unique())
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

