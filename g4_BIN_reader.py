import ROOT

read_file= ROOT.TFile.Open("/data/runzezhang/result/dmx.root","READ")
tree = read_file.Get("Scintillation Hits Info")
print(tree)