import uproot
filename = "Sap5cmx7.5cm3cmPad10cm25cmAmLi.root"
file = uproot.open(filename)["tree"]
print(file.keys())
df = file.arrays(["Hit", "x"], library="pd")
Capture = df['Hit'].tolist()
Elastic = df['x'].tolist()