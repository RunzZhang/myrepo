"""This file is for analyze the root files"""

import uproot
# filename = "/data/runzezhang/Geant4Simulaions/g411_TN/dmx.root"


class ReadRoot():
    def __init__(self):
        self.filepath = "/data/runzezhang/result/TN_sims/dmx.root"
        self.file = uproot.open(self.filepath)["tree"]
        print("columns: ",self.file.keys())
        #['Event', 'name', 'Parent ID', 'Track ID', 'Step ID', 'X/mm', 'Y/mm', 'Z/mm', 'Kinetic/keV', 'Recoiled/keV', 'Volume', 'Process']
        self.selected_columns = ["Event","name","Parent ID","Track ID","Step ID","X/mm","Kinetic/keV","Recoiled/keV","Volume","Process"]
        self.rows = 1000
        self.df = self.file.arrays(self.selected_columns, library="pd").head(self.rows)
        print(self.df)
        self.process_summary()

    def process_summary(self):
        process_clean=[]
        df_process = self.df[:]["Process"].to_list()
        for element in df_process:
            if element not in process_clean:
                process_clean.append()
        print(process_clean)





if __name__ =="__main__":
    RR = ReadRoot()