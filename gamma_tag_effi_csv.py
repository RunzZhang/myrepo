#translate G4 output file into list

import matplotlib.pyplot as plt
import numpy as np
import os, pickle
import pandas as pd
from scipy.stats import poisson
def add_paired_number(low,high,lst):

    if type(low)==int and type(high)==int and low< high and type(lst)==list:
        for i in range(low, high+1):
            lst.append(i)

    # print(lst)
    return lst

class thermal_neutron_calibration():
    def __init__(self, parent=None, join='Informacion.txt'):
        # super().__init__(parent)
        print(os.getcwd())
        self.base = os.getcwd()

        # self.Infoaddress = 'SDInformacion_428.csv'
        self.Infoaddress = 'Informacion_July18.csv'
        self.fullInfoaddress = os.path.join(self.base, self.Infoaddress)
        self.event_N_list=[10770, 4138, 31596, 19599, 184883, 36379, 12001, 20605, 63306, 18871, 77445, 6640, 61150, 10754, 33638, 17239, 36696, 12199, 681, 30010, 13433, 180585, 14845, 168179, 18587, 36118, 22727, 3391, 3940, 132822, 19929, 97370, 27002, 10366, 21280, 58203, 122465, 109169, 41456, 120842, 45476, 234387, 14913, 151374, 53623, 85910, 10745, 51086, 16029, 16346, 70420, 30837, 49398, 14526, 72860, 18395, 27804, 23181, 54640, 39489, 72224, 6459, 154248, 154201, 20895, 3746, 53307, 15409, 33371, 27986, 67163, 3990, 10852, 30495, 78646, 11452, 15448, 53706, 11159, 54955, 120026, 44144, 1921, 5115, 40063, 79976, 111068, 31326, 6149, 67146, 57312, 10125, 5808, 68404, 126896, 31276, 18004, 42277, 157481, 1673, 13889, 56421, 58008, 10307, 71019, 19727, 17448, 83074, 15429, 135887, 51898, 5810, 188228, 48475, 68204, 53617, 4494, 28698, 18316, 15549, 33127, 118269, 22288, 33815, 32871, 47666, 4763, 5580, 74754, 44461, 31737, 56270, 123791, 16539, 182010, 12749, 74974, 20430, 39842, 12051, 94590, 3212, 103773, 20432, 39015, 29579, 15080, 19926, 90670, 33317, 11109, 42349, 25934, 127228, 191929, 10482, 202987, 146158, 25993, 30321, 35106, 20844, 177576, 15802, 66153, 16, 1638, 115106, 22579, 142347, 74394, 39437, 15143, 69784, 25929, 73292, 235823, 185926, 28580, 5250, 155819, 29517, 170550, 65070, 6971, 17166, 39569, 25959, 64095, 224287, 119794, 56720, 34792, 4128, 18337, 30392, 12639, 15912, 56353, 54562, 10196, 134733, 37370, 12158, 106751, 21164, 22458, 73973, 87274, 24693, 21273, 23067, 15333, 106208, 57016, 49503, 33659, 60051, 26369, 44044, 13596, 42957, 49004, 57921, 59107, 37630, 13090, 37309, 73074, 15704, 156388, 148142, 223784, 42500, 28533, 171746, 155699, 30407, 47864, 26650, 107387, 134121, 66853, 4673, 10150, 226900, 28623, 45129, 49317, 10856, 53763, 9863, 10766, 72062, 16588, 247514, 34028, 5553, 8195, 91205, 54460, 18677, 133637, 9242, 72026, 97392, 28752, 25325, 3472, 16164, 156286, 36401, 125523, 53858, 106400, 33734, 43372, 55371, 216964, 45585, 10844, 32401, 70383, 11397, 16737, 13865, 159180, 124130, 69926, 67675, 28259, 174453, 29260, 73637, 8339, 12067, 33549, 20746, 36678, 39984, 104405, 39692, 58171, 70960, 24733, 52403, 136300, 17104, 29518, 275484, 12008, 49417, 65758, 8244, 62931, 124367, 113459, 105345, 81854, 41789, 183537, 14469, 30137, 73919, 6276, 19366, 61899, 127903, 37400, 17446, 59698, 30481, 8833, 13063, 18429, 70402, 20680, 127804, 17458, 845, 67915, 198469, 44029, 22891, 114606, 52277, 53595, 117289, 32562, 14876, 41514, 41748, 37199, 0, 37997, 40493, 75593, 117840, 34649, 54341, 212343, 72282, 14719, 0, 218205, 101, 52939, 68127, 105804, 127558, 84602, 93306, 89319, 39735, 73663, 171670, 75147, 2404, 58671, 28, 38311, 16109, 50264, 37904, 13978, 39019, 92328, 81181, 17177, 21097, 102808, 10703, 44959, 100361, 50313, 58106, 33631, 14015, 42339, 76241, 93672, 38118, 13043, 70563, 17224, 16976, 57410, 97208, 25576, 69258, 83190, 20407, 34256, 125718, 16674, 12889, 52535, 19761, 29212, 2, 18659, 74598, 23641, 3933, 33651, 64391, 108283, 18652, 88932, 41439, 51962, 111429, 84109, 66842, 92332, 36979, 11351, 135845, 66692, 154341, 198293, 4873, 76161, 20552, 64810, 11800, 17029, 36736, 64210, 0, 32457, 35879, 82394, 28613, 16557, 1246, 32828, 650, 56233, 23496, 47120, 184267, 18677, 27481, 4876, 78791, 11377, 29143, 3783, 4577, 62929, 142784, 50550, 79881, 60959, 39577, 17645, 119563, 62615, 13414, 72337, 22523, 7175, 94345, 37250, 29172, 57858, 129526, 28454, 102739, 57965, 40492, 317264, 227452, 18, 65439, 74546, 24058, 7338, 1181, 30032, 50408, 24238, 24833, 136843, 83929, 6641, 5806, 106764, 33935, 10328, 158633, 21934, 13155, 211262, 90178, 16966, 85400, 11091, 61206, 164416, 81023, 169109, 120312, 30656, 22009, 58811, 51182, 127359, 125936, 36289, 25634, 45009, 26890, 43053, 126764, 98946, 3780, 199403, 169570, 61254, 35745, 87906, 9357, 10907, 65966, 58509, 42031, 12554, 5395, 30631, 42240, 8644, 19501, 19141, 24507, 111694, 83550, 66905, 43655, 60632, 95969, 19347, 21983, 26994, 44733, 29469, 28335, 58633, 10419, 163866, 34655, 25142, 60872, 86356, 17541, 15858, 61971, 9255, 83854, 113289, 14187, 25615, 13603, 54409, 72374, 67317, 9745, 38246, 14719, 58466, 18247, 15306, 55211, 49459, 4975, 69587, 147111, 74866, 25697, 12996, 99711, 12141, 74219, 13764, 15220, 14519, 55099, 119347, 193797, 180165, 7579, 6784, 205522, 67008, 187554, 22762, 154889, 21938, 26295, 27708, 75983, 12708, 33020, 64158, 36325, 15910, 24856, 17137, 51821, 64825, 76165, 16345, 167819, 19788, 19836, 43355, 76315, 107529, 98743, 3564, 14463, 44542, 19406, 26077, 46987, 14691, 8397, 171687, 71308, 8300, 39067, 50857, 41647, 7553, 64411, 25020, 7217, 78869, 14323, 165299, 31921, 57089, 38013, 196094, 262514, 21135, 16486, 37252, 52417, 5819, 13556, 53078, 26620, 56798, 54491, 26241, 6629, 27228, 195461, 19836, 54118, 10589, 94722, 65778, 14808, 74641, 29120, 30055, 101514, 20099, 185787, 76646, 38168, 34908, 27474, 26482, 29854, 149365, 21550, 35650, 40750, 29572, 6828, 73158, 5945, 55354, 37560, 85794, 47472, 34608, 62436, 56146, 264381, 8716, 57103, 22887, 5406, 105271, 71846, 84859, 70903, 18965, 48671, 47414, 6996, 53592, 104271, 15321, 127106, 10810, 157772, 8600, 25457, 5018, 31647, 5602, 11305, 85311, 36899, 69611, 10135, 95591, 103363, 113906, 61996, 14314, 25843, 41520, 82153, 78511, 30769, 76769, 108779, 3641, 36335, 26833, 143639, 101350, 169187, 462, 84136, 37419, 29795, 14395, 22882, 63238, 33142, 201318, 91781, 85818, 90956, 16720, 15946, 36314, 83452, 3612, 25550, 14605, 67126, 30222, 64617, 18946, 24656, 13314, 238921, 16699, 31529, 21021, 21813, 6580, 8457, 140243, 89287, 94182, 35343, 133676, 2518, 26714, 40905, 27681, 71046, 53859, 35255, 26162, 5596, 42633, 89, 130916, 42857, 45581, 8, 70220, 23442, 136451, 37029, 44506, 14605, 18, 93983, 71678, 61162, 85193, 36086, 16843, 19936, 29492, 110149, 18818, 31466, 73254, 96591, 162854, 20346, 70782, 25021, 178068, 16702, 60884, 5522, 105981, 6795, 170866, 67982, 53564, 76920, 16750, 18655, 9917, 18244, 20310, 172219, 44732, 60682, 1862, 63542, 77261, 59222, 13793, 5284, 26998, 46755, 117533, 17, 120563, 84980, 58898, 126098, 55181, 14656, 12, 142884, 0, 26002, 43696, 141767, 41229, 49626, 24095, 52794, 21975, 17593, 8339, 17414, 16366, 22023, 39464, 163525, 21805, 4531, 12, 24672, 10754, 98680, 24377, 38465, 2275, 5, 61404, 15372, 54024, 53320, 43857, 3191, 26006, 37963, 14298, 56141, 28357, 45337, 35350, 14609, 43578, 37049, 10198, 54767, 44703, 26987, 1243, 88673, 54255, 74028, 28642, 801, 11262, 83233, 52160, 24667, 69640, 156557, 142250, 63563, 25544, 32506, 79240, 11550, 42854, 27885, 60918, 9568, 203882, 105682, 5033, 26720, 2287, 16583, 38227, 19836, 14068, 5454, 20508, 70931, 100492, 80470, 36917, 28063, 76457, 66706, 90701, 16322, 19307, 21284, 142751, 14008, 69802, 16738, 53997, 13949, 33134, 23237, 138867, 40433, 37648, 28817, 75036, 29333, 11416, 34364, 117340, 69376, 59093, 49198, 11295, 3595, 46840, 58093, 63115, 28020, 198799, 40051, 24717, 11692, 8180, 21374, 54467, 20460, 13948, 14238, 79865, 132247]


        # self.InfoPmtaddress = "PmtInformation_428.csv"
        # self.fullInfoPmtaddress = os.path.join(self.base, self.InfoPmtaddress)
        # self.config = 'config_422.png'
        # self.fullconfigaddress = os.path.join(self.base, self.config)
        # self.outputfile = 'result_tn_422.csv'
        # self.gammafile = 'gammaresult_April1.csv'
        # ONe entry like this event number : (incident energy(float, ev), outgoing(bool), Q Value (float) in MeV)




    def read_Information_checkparent(self):
        df = pd.read_csv(self.fullInfoaddress)
        # print(df.head(5))
        index = df.index
        idx_N = len(index)

        event_df=df[df["Process"]=="nCapture"]["Event"]
        event_list= event_df.to_list()
        event_bool_longlst=[]
        true_number = 0
        false_number = 0

        # for event_n in event_list:
        for event_n in event_list[:1000]:
            # True means electron parent can be opticalphoton, False means no
            event_bool = False
            temp_electron_list=[]
            temp_photon_list=[]
            temp_electron_df = df[(df["Event"]==event_n) & (df["Particle"]=="e-")]
            temp_photon_df = df[(df["Event"] == event_n) & (df["Particle"] == "opticalphoton")]
            # print(temp_electron_df.head(5))
            # print(temp_photon_df.head(5))
            idx_e_lst = temp_electron_df.index.to_list()
            idx_e_N = len(idx_e_lst)
            idx_photon_lst = temp_photon_df.index.to_list()
            idx_photon_N = len(idx_photon_lst)
            # print(idx_e_lst[:5])
            for idx in idx_photon_lst:
                # temp_photon_list=add_paired_number(temp_photon_df.iloc[idx]["startid"],temp_photon_df.iloc[idx]["endid"],temp_photon_list)
                # print("start",df.iloc[idx]["startid"])
                # print("end",df.iloc[idx]["endid"])
                # print(add_paired_number(int(df.iloc[idx]["startid"]), int(df.iloc[idx]["endid"]), temp_photon_list))
                temp_photon_list = add_paired_number(int(df.iloc[idx]["startid"]), int(df.iloc[idx]["endid"]), temp_photon_list)
            # print(temp_photon_list)
            for idx in idx_e_lst:
                # temp_electron_list.append(temp_electron_df.iloc[idx]["ParentID"])
                temp_electron_list.append(int(df.iloc[idx]["ParentID"]))
            # print(temp_electron_list)
            for id in temp_electron_list:
                if id in temp_photon_list:
                    if len(df[(df["Event"]==event_n)& (df["TrackID"]== id)].head(1))==0:
                        # print(True)
                        # print("electron",df[(df["Event"]==event_n)& (df["TrackID"]== id)].head(1))
                        event_bool=event_bool | True
                        # print("id",id)
                    else:
                        # print(len(df[(df["Event"]==event_n)& (df["TrackID"]== id)].head(1))==0)
                        # print(False)
                        event_bool = event_bool | False

                else:
                    event_bool = event_bool | False
            print(event_n,event_bool)
            event_bool_longlst.append(event_bool)
            self.event_N_list.append(len(temp_photon_list))
        for i in event_bool_longlst:
            if i:
                true_number += 1
            else:
                false_number +=1
        print("true",true_number)
        print("false",false_number)
        print("total",len(event_bool_longlst))

    def read_photon_number(self):
        print(self.event_N_list)
        print("min",min(self.event_N_list))
        print("max",max(self.event_N_list))
        print("mean",sum(self.event_N_list)/len(self.event_N_list))
        bins =100


        fig = plt.figure()
        ax = fig.add_subplot()
        ax.hist(self.event_N_list, bins= bins)

        ax.plot()
        ax.set_xlabel("photon_number")
        ax.set_ylabel("counts/bin")
        # ax.set_xlim(left=7,right=12)
        # ax.semilogy()


        plt.show()


    def possion_plot(self):
        fig, ax = plt.subplots(1, 1)
        mu = 54825
        mean, var, skew, kurt = poisson.stats(mu, moments='mvsk')
        x = np.arange(poisson.ppf(0.01, mu),
                      poisson.ppf(0.99, mu))
        ax.plot(x, poisson.pmf(x, mu), 'bo', ms=8, label='poisson pmf')
        plt.show()

    def read_photon_energy(self):
        # df = pd.read_csv(self.fullInfoaddress, nrows = 100000)
        df = pd.read_csv(self.fullInfoaddress)
        bins=500
        print(df.head(5))
        index = df.index
        idx_N = len(index)
        energy_df= df[df["Particle"]=="opticalphoton"]["Energy_Cinetica"]
        energy_list=(energy_df*10**6).to_list()


        fig = plt.figure()
        ax = fig.add_subplot()
        ax.hist(energy_list, bins= bins)

        ax.plot()
        ax.set_xlabel("energy/eV")
        ax.set_ylabel("counts/bin")
        ax.set_xlim(left=7,right=12)
        # ax.semilogy()

        ax2=ax.twinx()
        ax2.set_ylabel("Comparative Scintillation Strength")
        ax2.plot(self.LAr_PP,self.LAr_SCINT, color='orange')
        ax2.set_ylim(bottom=0,top=1750)
        plt.show()

    def read_photon_dep_capture_energy(self):
        # df = pd.read_csv(self.fullInfoaddress, nrows = 100000)
        capture_energy=[]
        df = pd.read_csv(self.fullInfoaddress)
        bins= 80
        print(df.head(5))
        index = df.index
        idx_N = len(index)
        for line in range(idx_N):
            try:
                if (df.iloc[line]["Process"] == "nCapture") & (df.iloc[line + 1]['Process'] == "ionIoni"):
                    capture_energy.append(df.iloc[line+1]["Energy_Cinetica"] * 10 ** 6)
            except:
                continue

        energy_df= df[df["Process"]=="ionIoni"]["Energy_Cinetica"]
        energy_list=(energy_df*10**6).to_list()
        print("total collection counts:",len(energy_list))


        print(capture_energy)
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.hist(capture_energy, bins= bins)

        ax.plot()
        ax.set_xlabel("energy/eV")
        ax.set_ylabel("counts/bin")
        # ax.set_xlim(left=7,right=12)
        # ax.semilogy()

        # ax2=ax.twinx()
        # ax2.set_ylabel("Comparative Scintillation Strength")
        # ax2.plot(self.LAr_PP,self.LAr_SCINT, color='orange')
        # ax2.set_ylim(bottom=0,top=1750)
        plt.show()
    def read_photon_dep_scatter_energy(self):
        # df = pd.read_csv(self.fullInfoaddress, nrows = 100000)
        capture_energy = []
        df = pd.read_csv(self.fullInfoaddress)
        bins = 80
        print(df.head(5))
        index = df.index
        idx_N = len(index)
        for line in range(idx_N):
            try:
                if (df.iloc[line-1]["Process"] != "nCapture") & (df.iloc[line]['Process'] == "ionIoni"):
                    capture_energy.append(df.iloc[line]["Energy_Cinetica"] * 10 ** 6)
            except:
                continue

        # energy_df= df[df["Process"]=="ionIoni"]["Energy_Cinetica"]
        # energy_list=(energy_df*10**6).to_list()

        print(capture_energy)
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.hist(capture_energy, bins=bins)

        ax.plot()
        ax.set_xlabel("energy/eV")
        ax.set_ylabel("counts/bin")
        # ax.set_xlim(left=7,right=12)
        # ax.semilogy()

        # ax2=ax.twinx()
        # ax2.set_ylabel("Comparative Scintillation Strength")
        # ax2.plot(self.LAr_PP,self.LAr_SCINT, color='orange')
        # ax2.set_ylim(bottom=0,top=1750)
        plt.show()

    def read_photon_dep_energy(self):
        # df = pd.read_csv(self.fullInfoaddress, nrows = 100000)

        df = pd.read_csv(self.fullInfoaddress)
        bins = 80
        print(df.head(5))
        index = df.index
        idx_N = len(index)

        energy_df= df[df["Process"]=="ionIoni"]["Energy_Cinetica"]
        energy_list=(energy_df*10**6).to_list()

        print(energy_list)
        fig,ax = plt.subplots(2)
        count_list=[]
        n, bins_out, patches = ax[0].hist(energy_list, bins=bins)
        for i in range(len(bins_out)-1):
            middle_point = (bins_out[i]+bins_out[i+1])/2
            count_list.append(middle_point/6.9)
        n_total = 0
        pdf = []
        for i in n:
            n_total = n_total + i
        for i in n:
            pdf.append(i/n_total)

        print("counts",count_list)
        print("pdf",pdf)

        expect = 0
        for i in range(len(count_list)):
            expect = expect + count_list[i]*pdf[i]
        tag_effi = expect*0.04


        ax[0].plot()
        ax[0].set_xlabel("energy/eV")
        ax[0].set_ylabel("counts/bin")
        # ax.set_xlim(left=7,right=12)
        # ax.semilogy()

        ax[1].plot(count_list, pdf)
        ax[1].set_xlabel("N of photons")
        ax[1].set_ylabel("possbility")

        print("Expectation", expect)
        print("Effi", tag_effi)

        plt.show()


    def read_N_ini_energy(self):
        # df = pd.read_csv(self.fullInfoaddress, nrows = 100000)
        capture_energy = []
        df = pd.read_csv(self.fullInfoaddress)
        bins = 40
        print(df.head(5))
        index = df.index
        idx_N = len(index)
        for line in range(idx_N):
            try:
                if (df.iloc[line]["Event"] == df.iloc[line-1]["Event"]+1):
                    capture_energy.append(df.iloc[line]["Energy_Cinetica"] * 10 ** 6)
            except:
                continue

        # energy_df= df[df["Process"]=="ionIoni"]["Energy_Cinetica"]
        # energy_list=(energy_df*10**6).to_list()

        print(capture_energy)
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.hist(capture_energy, bins=bins)

        ax.plot()
        ax.set_xlabel("energy/eV")
        ax.set_ylabel("counts/bin")
        # ax.set_xlim(left=7,right=12)
        # ax.semilogy()

        # ax2=ax.twinx()
        # ax2.set_ylabel("Comparative Scintillation Strength")
        # ax2.plot(self.LAr_PP,self.LAr_SCINT, color='orange')
        # ax2.set_ylim(bottom=0,top=1750)
        plt.show()

    def read_tag_G_efficiency(self):
        df_SD = pd.read_csv(self.fullInfoaddress)
        df_Pmt = pd.read_csv(self.fullInfoPmtaddress)

        SD_list = df_SD["Evento"].to_list()
        Pmt_list = df_Pmt["m"].to_list()

        SD_len = len(SD_list)
        Pmt_len = len(Pmt_list)
        j0=0
        counts=0

        for i in range(SD_len):
            for j in range(j0,Pmt_len,1):
                if SD_list[i]==Pmt_list[j]:
                    j0=j
                    counts += 1
                    break

        print("Ar Capture", SD_len, "Pmt Counts", Pmt_len, "Pmt Double Check", counts, "Eff",counts/SD_len)






if __name__=="__main__":
    # get hits number and plot positions
    tnc= thermal_neutron_calibration()
    # tnc.read_Information_checkparent()
    tnc.read_photon_number()
    # tnc.possion_plot()


    # print(add_paired_number(1,2,[2,3]))