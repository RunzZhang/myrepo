import uproot
filename = "/data/runzezhang/Geant4Simulaions/g411_TN/dmx.root"
file = uproot.open(filename)["tree"]
print(file.keys())
df = file.arrays(["Event", "x"], library="pd")
Capture = df['Event'].tolist()
# Elastic = df['x'].tolist()