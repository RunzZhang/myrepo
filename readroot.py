"""This file is for analyze the root files"""

import uproot
# # filename = "/data/runzezhang/Geant4Simulaions/g411_TN/dmx.root"
# filename = "/data/runzezhang/result/TN_sims/dmx.root"
# file = uproot.open(filename)["tree"]
# print(file)
# print(file.keys())
#
# df = file.arrays(["Event","name","Parent ID","Track ID","Step ID","X/mm","Kinetic/keV","Recoiled/keV","Volume","Process"], library="pd")
# df = df.head(1000)
# # df = file.arrays(["Event", "x"], library="pd")
# Capture = df['Event'].tolist()
# print(df)

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





if __name__ =="__main__":
    RR = ReadRoot()