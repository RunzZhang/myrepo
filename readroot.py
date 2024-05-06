"""This file is for analyze the root files"""
import pandas as pd
import uproot
# filename = "/data/runzezhang/Geant4Simulaions/g411_TN/dmx.root"

class RestructureRoot():
    def __init__(self):
        self.filepath = "/data/runzezhang/result/TN_sims/dmx.root"
        self.reconstruct_filepath = "/data/runzezhang/result/TN_sims/dmx_rc.root"
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

        uproot.writing._dask_write.dask_write(self.df, self.reconstruct_filepath)

class ReadRoot():
    def __init__(self):
        self.filepath = "/data/runzezhang/result/TN_sims/dmx.root"
        self.file = uproot.open(self.filepath)["tree"]
        print("columns: ",self.file.keys())
        #['Event', 'name', 'Parent ID', 'Track ID', 'Step ID', 'X/mm', 'Y/mm', 'Z/mm', 'Kinetic/keV', 'Recoiled/keV', 'Volume', 'Process']
        self.selected_columns = ["Event","name","Parent ID","Track ID","Step ID","X/mm","Kinetic/keV","Recoiled/keV","Volume","Process"]
        self.rows = 1000
        # self.df = self.file.arrays(self.selected_columns, library="pd").head(self.rows)
        self.df = self.file.arrays(self.selected_columns, library="pd")
        # print(self.df)
        self.process_summary()

    def process_summary(self):
        process_clean=[]
        df_process = self.df[:]["Process"].to_list()
        for element in df_process:
            if element not in process_clean:
                process_clean.append(element)
        print(process_clean)





if __name__ =="__main__":
    ReR = RestructureRoot()
    # RR = ReadRoot()