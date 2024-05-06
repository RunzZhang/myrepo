import uproot
# filename = "/data/runzezhang/Geant4Simulaions/g411_TN/dmx.root"
filename = "/data/runzezhang/result/TN_sims/dmx.root"
file = uproot.open(filename)["tree"]
print(file)
print(file.keys())
# "Event");
#   man->CreateNtupleDColumn("name");
#   man->CreateNtupleDColumn("Parent ID");
#   man->CreateNtupleDColumn("Track ID");
#   man->CreateNtupleDColumn("Step ID");
#   man->CreateNtupleDColumn("X/mm");
#   man->CreateNtupleDColumn("Y/mm");
#   man->CreateNtupleDColumn("Z/mm");
#   man->CreateNtupleSColumn("Kinetic/keV");
#   man->CreateNtupleDColumn("Recoiled/keV");
#   man->CreateNtupleDColumn("Volume");
#   man->CreateNtupleDColumn("Process")
df = file.arrays(["Event","name","Parent ID","Track ID","Step ID","X/mm","Kinetic/keV","Recoiled/keV","Volume","Process"], library="pd")
# df = file.arrays(["Event", "x"], library="pd")
Capture = df['Event'].tolist()
print(df)
# Elastic = df['x'].tolist()