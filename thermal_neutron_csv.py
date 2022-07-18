#translate G4 output file into list

import matplotlib.pyplot as plt
import numpy as np
import os, pickle
import pandas as pd

class thermal_neutron_calibration():
    def __init__(self, parent=None, join='Informacion.txt'):
        # super().__init__(parent)
        print(os.getcwd())
        self.base = os.getcwd()

        # self.Infoaddress = 'SDInformacion_428.csv'
        self.Infoaddress = 'Informacion_0615.csv'
        self.fullInfoaddress = os.path.join(self.base, self.Infoaddress)
        self.InfoPmtaddress = "PmtInformation_428.csv"
        self.fullInfoPmtaddress = os.path.join(self.base, self.InfoPmtaddress)
        self.config = 'config_422.png'
        self.fullconfigaddress = os.path.join(self.base, self.config)
        self.outputfile = 'result_tn_422.csv'
        self.gammafile = 'gammaresult_April1.csv'
        # ONe entry like this event number : (incident energy(float, ev), outgoing(bool), Q Value (float) in MeV)
        self.INFOmatrix = {}
        self.Q_list=[]
        self.abnormal_event=[]
        self.A = 39.962 #g/mole
        self.density = 1.779*10**(-3) # g/cm3
        self.L = 2 #cm
        self.NA = 6.02*10**23 # /mole
        #sigma_N =10**24*self.A*out/(denstity*L*NA*in) in milibarn
        self.coefi = 10**24*self.A/(self.density*self.L*self.NA)
        self.y=np.logspace(-2, -1, num = 40)
        self.sigma_matrix = []
        for i in self.y:
            self.sigma_matrix.append((i,{"Total":0, "Out":0}))
        # print(self.sigma_matrix)
        self.gamma_list=[]

        # https: // www - nds.iaea.org / exfor / servlet / E4sGetIntSection?SectID = 148438 & req = 2250 & e4up = 0 & PenSectID = 6570496 & pen = 0
        # Library: JENDL - 3.3

        self.endfE= [1.265500*10**(-2),  1.581625*10**(-2),  1.697880*10**(-2),
                     1.897750*10**(-2),  2.230510*10**(-2),  2.530000*10**(-2),
                     3.062612*10**(-2),  3.295760*10**(-2),  3.595223*10**(-2),
                     4.361020*10**(-2),  4.660445*10**(-2),  5.725668*10**(-2),
                     6.491530*10**(-2),  6.790890*10**(-2),  8.622030*10**(-2),
                     8.921335*10**(-2),  1.105178*10**(-1)]
        self.endfsig = [9.174450*10**2,   8.209190*10**2,  7.849240*10**2,
                        7.477450*10**2,   6.846480*10**2,  6.431150*10**2,
                        5.878920*10**2,   5.633990*10**2,  5.429490*10**2,
                        4.897090*10**2,   4.776170*10**2,  4.335600*10**2,
                        4.014440*10**2,   3.941080*10**2,  3.483120*10**2,
                        3.439990*10**2,   3.124790*10**2]


        self.endfE1 = [.012928235,  .014032614,  .015242812,
                       .016544912,  .017971775,  .019506993,
                       .021189312,  .022999384,  .024798654,
                       .025300000,  .027481919,  .029829529,
                       .032402083,  .035169992,  .038203117,
                       .041466574,  .045042727,  .048890449,
                       .053106852,  .057643442,  .062614718,
                       .067963507,  .073824805,  .080131202,
                       .087041865,  .094477313,  .102625211]
        self.endfsig1 = [.923960966*10**3,  .886862090*10**3,  .850933662*10**3,
                         .816767141*10**3,  .783677319*10**3,  .752212948*10**3,
                         .721740691*10**3,  .692764496*10**3,  .667165590*10**3,
                         .660526127*10**3,  .633767036*10**3,  .608320515*10**3,
                         .583678875*10**3,  .560256516*10**3,  .537581534*10**3,
                         .516017376*10**3,  .495135423*10**3,  .475255633*10**3,
                         .455990922*10**3,  .437664425*10**3,  .419915151*10**3,
                         .403043834*10**3,  .386717603*10**3,  .371206768*10**3,
                         .356188252*10**3,  .341894480*10**3,  .328034297*10**3]


        self.LAr_PP = [1.54 , 7.04 , 8.15 , 8.37 , 8.55 ,
                            8.73 , 8.85 , 8.98 , 9.11 , 9.25 ,
                            9.39 , 9.46 , 9.53 , 9.61 , 9.68 ,
                            9.76 , 9.84 , 9.99 , 10.16 , 10.32 ,
                            10.50 , 10.68 , 10.87 , 11.16 , 11.92 ,
                            15.49 , 20.66 ]

        self.LAr_SCINT = [0.0, 0.0, 85.9, 207.2, 361.1, 570.5,
                               745.1, 928.5, 1111.5, 1278.4, 1413.0, 1463.3,
                               1500.1, 1523.1, 1530.4, 1523.1, 1500.1,
                               1413.0, 1278.4, 1111.5, 928.5, 745.1, 574.3,
                               361.0, 85.9, 0.0, 0.0]



    def read_Information_ncap(self):
        df = pd.read_csv(self.fullInfoaddress)
        print(df.head(5))
        index = df.index
        idx_N = len(index)
        df_PV = df[df["Volume"]=='pressure_vessel_phys'][["x","y","z"]]
        df_OJ = df[df["Volume"] == 'outer_jar_phys'][["x", "y", "z"]]
        df_LAr = df[df["Volume"] == 'LAr_phys'][["x", "y", "z"]]
        df_VV = df[df["Volume"] == 'Vacuum_vessel_phys'][["x", "y", "z"]]
        df_HDPE = df[df["Volume"] == 'HDPE_pressure_vessel_phys'][["x", "y", "z"]]
        df_HF = df[df["Volume"] == 'hydraulic_fluid_phys'][["x", "y", "z"]]

        PV_x = df_PV["x"].to_list()
        PV_y = df_PV["y"].to_list()
        PV_z = df_PV["z"].to_list()

        OJ_x = df_OJ["x"].to_list()
        OJ_y = df_OJ["y"].to_list()
        OJ_z = df_OJ["z"].to_list()

        LAr_x = df_LAr["x"].to_list()
        LAr_y = df_LAr["y"].to_list()
        LAr_z = df_LAr["z"].to_list()

        VV_x = df_VV["x"].to_list()
        VV_y = df_VV["y"].to_list()
        VV_z = df_VV["z"].to_list()

        HDPE_x = df_HDPE["x"].to_list()
        HDPE_y = df_HDPE["y"].to_list()
        HDPE_z = df_HDPE["z"].to_list()

        HF_x = df_HF["x"].to_list()
        HF_y = df_HF["y"].to_list()
        HF_z = df_HF["z"].to_list()



        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        lenPV = len(PV_x)
        lenOJ = len(OJ_x)
        lenLAr = len(LAr_x)
        lenVV = len(VV_x)
        lenHDPE = len(HDPE_x)
        lenHF = len(HF_x)
        print(lenPV,lenOJ,lenLAr,lenVV,lenHDPE, lenHF)

        minn =1000
        ax.scatter(VV_x[:minn], VV_y[:minn], VV_z[:minn], marker='.', c='yellow', label='Vaccum Vessel')
        ax.scatter(PV_x[:minn], PV_y[:minn], PV_z[:minn], marker='.', c='blue', label='Pressure Vessel')
        ax.scatter(HDPE_x[:minn], HDPE_y[:minn], HDPE_z[:minn], marker='.', c='green', label='HDPE Pressure Vessel')
        ax.scatter(HF_x[:minn], HF_y[:minn], HF_z[:minn], marker='.', c='maroon', label='Hydraulic Fluid')
        ax.scatter(OJ_x[:minn], OJ_y[:minn], OJ_z[:minn], marker='.', c='red', label='Outer Jar')
        ax.scatter(LAr_x[:minn], LAr_y[:minn], LAr_z[:minn], marker='.', c='orange', label='LAr')



        ax.scatter(-400, 0, 500, marker='o', c='black', s=100,label = 'Neutron Source')

        ax.legend(loc=(0.7,0.85))
        ax.set_xlabel('X axis/mm')
        ax.set_ylabel('Y axis/mm')
        ax.set_zlabel('Z axis/mm')
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)
        # plt.savefig(self.fullconfigaddress)
        print('PV', lenPV/idx_N)
        print('OJ', lenOJ / idx_N)
        print('LAr', lenLAr / idx_N)
        print('VV', lenVV / idx_N)
        print('HF', lenHF / idx_N)
        print('HDPE', lenHDPE / idx_N)
        sum = (lenPV+lenOJ+lenLAr+lenVV+lenHF+lenHDPE)
        print("sum", sum, 'total', 300000, sum/300000)
        # print("other", 1-sum)
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
    # tnc.read_Information_ncap()
    # tnc.read_photon_energy()
    # tnc.read_photon_dep_capture_energy()
    # tnc.read_N_ini_energy()
    # tnc.read_photon_dep_scatter_energy()
    tnc.read_photon_dep_energy()
    # tnc.read_tag_G_efficiency()
    # tnc.plot_Q(True)
    # tnc.plot_cross(read=True)
    # tnc.read_gamma_Information()
    # tnc.plot_gamma(True)