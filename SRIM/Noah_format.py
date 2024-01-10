import pickle
import numpy as np
class regen_txt():
    def __init__(self):
        self.address ="/data/runzezhang/result/SRIM_MC/MC_argon_full_20231206_full"
        # in eV
        self.recoil_list = []
        self.txt = ""
        self.outputMCaddress = "/data/runzezhang/result/SRIM_MC/MC_argon_full_20231206_full_Noahformat.txt"
        self.load()
        self.rewrite()
    def load(self):
        with open(self.address, "rb") as fp:  # Unpickling
            self.recoil_list = pickle.load(fp)
            print("read", self.recoil_list)
    def rewrite(self):
        for i in range(len(self.recoil_list)):
            self.txt += str(i)+" "+str(i)+ " 100 0 "+ str(self.recoil_list[i])+" 1 0\n"
        f = open(self.outputMCaddress, "a")
        f.write(self.txt)
        f.close()


if __name__ =="__main__":
    RT = regen_txt()


