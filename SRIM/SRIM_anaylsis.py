import matplotlib.pyplot as plt
import numpy as np
import os, pickle, random
import sys
import pandas as pd
from statistics import stdev
from scipy.optimize import curve_fit

pd.options.display.max_seq_items = None
# pd.set_option('display.max_rows', None)
# pd.set_option('display.max_columns', None)

random.seed(10)
skip = 0
Data = 0
avo = 6.0221408e+23
amu = 0
KE = 0
IonData = ''
file_lines = []
e_v = 1.602*10**(-19)
K_125 = 1.5414423840161205e-05
K_500 = 1.5414372244542094e-05
PI = np.pi

def li_func(x, k, b):
    return k * x + b


def li_func2(x, k):
    return k * x



class SRIM_EXY():
    def __init__(self):
        super().__init__()
        self.file_name = 'EXYZArgon1keV.txt'
        self.file_name_edit = self.file_name[0:-4] + '_edit.txt'
        self.displacement = []  # Want to record the displacement in each step (in meters)
        self.displacement_1d = []
        self.Ek_diff = []  # Kiniteic energy differenct between each step
        self.Ek_diff_1d = []
        self.El_diff = []  # Electronic energy loss between each step
        self.El_diff_1d = []
        self.vel = []
        self.vel_1d = []
        self.vi_1d = []
        self.time = []
        self.time_1d = []
        self.ERs = []
        self.ERs_1d = []

        self.data_ini()
        self.table = SRIM_Table()
        self.secondary_variable()
        self.scatterplot(self.E_compare_diff_1d, self.Energy,  "E_diff/eV", "Ek/keV")
        print("sqrt E",self.sqrt_E[:100])
        # self.plotEvX_fit( self.sqrt_E, self.ElStop)



        # self.scatterplot( self.ERs_1d, self.displacement_1d, "ER/eV", "x/A")
        # self.scatterplot(self.ERs_1d, self.time_1d, "ER/eV", "t/s")
        # self.scatterplot( self.ElStop, self.time_1d, "Eloss rate", "t/s")
        # self.time_displacement_hist(self.time_1d, self.displacement_1d, self.ERs_1d, self.vel_1d)
        # self.scatterplot(self.Elrecoil, self.Energy, "Erecoil", "Energy")
        # self.selected_hist(0.2,0.0)





    def data_ini(self):
        with open(self.file_name, 'r') as file:  # Searches for the row in which data starts and
            done = 0  # looks for ion data in the file (we want the atomic mass)
            for num, line in enumerate(file, 1):
                if '0000001' in line and done == 0:
                    # print('found at line:', num)
                    skip = num - 1
                    done = 1
                if 'Ion Data' in line:
                    Data = num

        with open(self.file_name, 'r') as file:  # Adds a space after first column to allow
            n = 0  # pandas to use double-space as separator
            for line in file:
                if n >= skip:
                    file_lines.append(''.join([line.strip()[:7], ' ', line.strip()[7:], '\n']))
                else:
                    # file_lines.append(''.join([line.strip(), '\n']))
                    file_lines.append(''.join([line]))
                n += 1

        with open(self.file_name, 'r') as file:  # Adds the ion data string to a variable
            n = 0
            for line in file:
                if n == Data:
                    IonData = line.rstrip()
                    IonData = (' '.join(IonData.split())).split()
                n += 1
        self.amu = float(IonData[1])  # Set the atomic mass variable
        self.mass = self.amu / avo * 0.001  # calculate mass of 1 atom in kg
        # masskev = amu/avo * 0.001 * 5.6095886e32 # mass in keV
        print(self.amu)
        print("mass",self.mass)

        # with open(self.file_name_edit, 'w') as file: # Write the editied file for pandas to open
        #   file.writelines(file_lines)

        self.daf = pd.read_csv(self.file_name, sep='  | ', skiprows=skip, header=None, engine=('python'))
        # print(self.daf.head(5))
        self.daf.rename(columns={0: 'Ion Number'}, inplace=True)
        self.daf.rename(columns={1: 'Energy (keV)'}, inplace=True)
        self.daf.rename(columns={2: 'x (A)'}, inplace=True)
        self.daf.rename(columns={3: 'y (A)'}, inplace=True)
        self.daf.rename(columns={4: 'z (A)'}, inplace=True)
        self.daf.rename(columns={5: 'Electronic Stop (eV/A)'}, inplace=True)
        self.daf.rename(columns={6: 'Energy Lost to Last Recoil (eV)'}, inplace=True)

        self.Posx = self.daf['x (A)'].tolist()
        print(self.Posx)

        self.Posy = self.daf['y (A)'].tolist()
        self.Posz = self.daf['z (A)'].tolist()
        self.Energy = self.daf['Energy (keV)'].tolist()
        # self.Energy = [1.60218e-16 * x for x in self.Energy]  # Writing all energy in Joules
        self.ElStop = self.daf['Electronic Stop (eV/A)'].tolist()
        self.Elrecoil = self.daf['Energy Lost to Last Recoil (eV)'].tolist()
        self.IonNum = self.daf['Ion Number'].tolist()
        self.Events = max(self.IonNum)
        self.KEInitial = (self.daf.iloc[0][1])  # Set kinetic energy variable (in keV)
        self.vi_1d = [np.sqrt(2*e_v*1000*self.Energy[i]/self.mass) for i in range(len(self.Energy))]
        self.sqrt_E = [np.sqrt(1000*self.Energy[i]) for i in range(len(self.Energy))]
        # print(KE)
        print(len(self.Posx), len(self.Energy), len(self.ElStop), len(self.Elrecoil), len(self.IonNum))
        print(self.Events)
        min_value = 5
        for value in self.Energy:
            if value <min_value and value != 0:
                min_value = value
        print("min energy kev",min_value)
    def secondary_variable(self):



        (self.displacement_1d,self.displacement)= self.fetch_data(len(self.Posx), 0)
        (self.Ek_diff_1d,self.Ek_diff)= self.fetch_data(len(self.Posx), 1)
        (self.El_diff_1d,self.El_diff)= self.fetch_data(len(self.Posx), 2)
        (self.vel_1d, self.vel)= self.fetch_data(len(self.Posx), 3)
        (self.time_1d, self.time)= self.fetch_data(len(self.Posx), 4)
        (self.ERs_1d ,self.ERs) = self.fetch_data(len(self.Posx), 5)
        (self.E_compare_diff_1d, self.E_compare_diff) = self.fetch_data(len(self.Posx), 6)

    def fetch_data(self, length, mode):
        # print(self.Posx[:20])
        list_2d = []
        list_1d =[]

        event_pointer = 1
        temp = []
        step_pointer = 0
        Ed_limit_pointer = []
        Ek_limit_pointer = []
        Num_limit_pointer = []
        for j in range(length):  # all length should be same
            # for j in range(16):
            if self.IonNum[j] == event_pointer:
                if step_pointer == 0:  # row 1
                    temp.append(0)
                    list_1d.append(0)
                else:  # add data
                    # mode0 displacement(m), mode1 kinetic energy loss(eV), mode2 electronic energy loss(eV), mode 3 average velocity(m/s)
                    # mode 4 time(s), # mode 5 energy loss by recoiled collision # mode 6 kinetic energy loss difference between experimental and theory
                    if mode ==0:
                        value = (1e-10) * (
                            np.sqrt((self.Posx[j] - self.Posx[j - 1]) ** 2 + (self.Posy[j] - self.Posy[j - 1]) ** 2 + (
                                    self.Posz[j] - self.Posz[j - 1]) ** 2))
                    elif mode ==1 :
                        value = (self.Energy[j] - self.Energy[j - 1])*1000
                    elif mode ==2:
                        value = self.ElStop[j] * (
                            np.sqrt((self.Posx[j] - self.Posx[j - 1]) ** 2 + (self.Posy[j] - self.Posy[j - 1]) ** 2 + (
                                    self.Posz[j] - self.Posz[j - 1]) ** 2))
                    elif mode ==3:
                        value1  = np.sqrt(2*self.Energy[j-1]*1000*e_v/(self.mass))
                        value2 = np.sqrt(2*self.Energy[j]*1000*e_v/(self.mass))
                        value = (value1+value2)/2
                        # print(value)
                    elif mode ==4:
                        distance = (1e-10) * (
                            np.sqrt((self.Posx[j] - self.Posx[j - 1]) ** 2 + (self.Posy[j] - self.Posy[j - 1]) ** 2 + (
                                    self.Posz[j] - self.Posz[j - 1]) ** 2))
                        value1 = np.sqrt(2 * self.Energy[j - 1] * 1000 * e_v / (self.mass))
                        value2 = np.sqrt(2 * self.Energy[j] * 1000 * e_v / (self.mass))
                        vel = (value1 + value2) / 2
                        value = distance/vel
                    elif mode ==5:
                        value  = -(self.Energy[j] - self.Energy[j - 1]) * 1000-self.ElStop[j] * (
                            np.sqrt((self.Posx[j] - self.Posx[j - 1]) ** 2 + (self.Posy[j] - self.Posy[j - 1]) ** 2 + (
                                    self.Posz[j] - self.Posz[j - 1]) ** 2))
                    elif mode ==6:
                        value  = (self.Energy[j] - self.Energy[j - 1]) * 1000-((self.ElStop[j]+self.ElStop[j-1])/2+self.table.output((self.Energy[j]+self.Energy[j-1])*1000/2) )* (
                            np.sqrt((self.Posx[j] - self.Posx[j - 1]) ** 2 + (self.Posy[j] - self.Posy[j - 1]) ** 2 + (
                                    self.Posz[j] - self.Posz[j - 1]) ** 2))
                        if value < -990:
                            # rule out bc step distance is too small
                            Ek_limit_pointer.append(self.Energy[j])
                            Num_limit_pointer.append(self.IonNum[j])
                            # find the row number of -990
                            Ed_limit_pointer.append(j)

                    temp.append(value)
                    list_1d.append(value)
                step_pointer += 1
            elif self.IonNum[j] == event_pointer + 1:
                # if it is in new event, should jump to next ion number and initialize pointers
                list_2d.append(temp)

                # print(j,list_2d)
                step_pointer = 0
                event_pointer += 1
                temp = []
                temp.append(0)
                list_1d.append(0)
                step_pointer += 1
            else:
                print("else")
        list_2d.append(temp)

        print("mode",mode, len(list_1d), "\n ",list_1d[:30], list_2d[:3])
        if mode ==6:
            print("lower limit 990", Ed_limit_pointer, '\n', Ek_limit_pointer,"\n", Num_limit_pointer)
        return (list_1d, list_2d)


    def plot_EvsX(self, E_list, v_list):  # both 1D
        start = 0
        end = -1
        # plot part of the list
        E_list = E_list[start:end]
        v_list = v_list[start:end]

        plt.plot(v_list, E_list)
        plt.show()

    def plotEvX_fit(self, velocity_dataset, E_lost_dataset):
        popt, pcov = curve_fit(li_func2, velocity_dataset, E_lost_dataset)

        perr = np.sqrt(np.diag(pcov))
        print("popt", popt, "\n", "perr", perr)
        print(popt[0])
        fit_data = []
        for i in velocity_dataset:
            fit_data.append(li_func2(i, popt[0]))

        plt.plot(velocity_dataset, fit_data, 'r-', label="linear fit")
        plt.xlabel("velocity m/s")
        plt.ylabel("dE/dx eV/A")
        # plt.xlim(0,500)
        plt.legend()

        plt.show()

    def time_displacement_hist(self, t, x, ER, velocity):
        bin_n = 500
        figure, axis = plt.subplots(4)

        axis[0].hist(t, bins=bin_n, label="time")
        axis[1].hist(x, bins=bin_n, label="displacement")
        axis[2].hist(ER, bins=bin_n, label="ER")
        axis[3].hist(velocity, bins=bin_n, label="v")

        plt.show()

    def scatterplot(self, y, x, ylabel="y", xlabel="x"):
        plt.scatter(x, y)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.show()

    def Eloss(self, Ediff, E_i, El):
        print(len(Ediff))
        print(len(E_i))
        print(len(El))
        Ek_diff = []
        for i in range(len(El)):
            Ek_diff.append(El[i])
            # Ek_diff.append(Ediff[i]-El[i])
            # Ek_diff.append(Ediff[i] - 0)
        plt.scatter(E_i, Ek_diff)
        plt.xlabel("Ei/eV")
        plt.ylabel("E_kdiff/j")
        plt.show()

    def selected_hist(self, uplim, lowlim):
        self.selected_Erecoil = []
        self.selected_Energy = []
        for i in range(len(self.Energy)):
            if lowlim < self.Energy[i] and self.Energy[i] < uplim:
                self.selected_Erecoil.append(self.Elrecoil[i])
                self.selected_Energy.append(self.Energy[i])


        plt.hist(self.selected_Erecoil)
        plt.xlabel("energy/eV")
        plt.show()




class SRIM_Table():
    def __init__(self):
        super().__init__()
        # self.file_name = 'argon_140keV(gas).txt'
        self.file_name = 'argon_tn(gas).txt'
        # self.energy_startpoint = "10.00 keV"
        self.energy_startpoint ="9.99999 eV"
        self.file_name_edit = self.file_name[0:-4] + '_edit.txt'


        self.step = 200  # in ev

        self.Z_tp = 18
        self.C_tf = (9 * PI ** 2 / (2 ** 7)) ** (1 / 3)
        self.A_tp = 41
        self.Z = (2 * self.Z_tp ** 0.23) ** 2
        self.EB = 13.6  # in eV
        self.a0 = 5.291 * 10 ** (-11)  # in m
        # self.a0 = 0.5291  # in A
        # self.energy = range(self.min, self.max, self.step)
        self.F1 = self.C_tf * 0.5 * (0.5 / 13.6) / (self.Z_tp ** 2 * self.Z ** (0.5))

        print("F1",self.F1)

        self.au = 0.8853 * self.a0 / (2 * self.Z_tp ** 2) # in m
        self.au = self.au*100 # in cm
        print("au", self.au)
        self.gam = 4 / 2 ** 2
        self.ev = 1.60218e-19
        self.Tar_Den = 2.1361E+22 # atoms/cm3
        self.SN = [] # shouldbe ev/cm
        self.LSS_factor = 1.2656 * (10 ** (-2))



        self.data_ini()
        self.fetch_data()
        # self.compare_theory()
        # self.plot_data()

    def data_ini(self):

        with open(self.file_name, 'r') as file:  # Searches for the row in which data starts and
            done = 0  # looks for ion data in the file (we want the atomic mass)
            for num, line in enumerate(file, 1):
                if self.energy_startpoint in line and done == 0:
                    # print('found at line:', num)
                    skip = num - 1
                    done = 1
                if 'Ion Data' in line:
                    Data = num

        with open(self.file_name, 'r') as file:  # Adds a space after first column to allow
            n = 0  # pandas to use double-space as separator
            for line in file:
                if n >= skip:
                    file_lines.append(''.join([line.strip()[:10], ' ', line.strip()[10:], '\n']))
                else:
                    # file_lines.append(''.join([line.strip(), '\n']))
                    file_lines.append(''.join([line]))
                n += 1

        # with open(self.file_name, 'r') as file:  # Adds the ion data string to a variable
        #     n = 0
        #     for line in file:
        #         if n == Data:
        #             IonData = line.rstrip()
        #             IonData = (' '.join(IonData.split())).split()
        #         n += 1
        self.amu = 39.962  # Set the atomic mass variable
        self.mass = self.amu / avo * 0.001  # calculate mass of 1 atom in kg
        # masskev = amu/avo * 0.001 * 5.6095886e32 # mass in keV
        # print(self.amu)
        # print(self.mass)

        # with open(self.file_name_edit, 'w') as file: # Write the editied file for pandas to open
        #   file.writelines(file_lines)

        self.daf = pd.read_csv(self.file_name, delim_whitespace=True, skiprows=skip, header=None, engine=('python'))
        # print("ini",self.daf)

        self.end_of_file = 0
        for i in range(self.daf.index.size): # delete the aferwards info

            try:
                print(float(self.daf.iloc[i][0]))
            except:
                self.end_of_file = i
                break
        self.daf = self.daf.head(self.end_of_file)
        # print("edit",self.daf)

        # print(self.daf.head(5))
        self.daf.rename(columns={0: 'Ion Energy'}, inplace=True)
        self.daf.rename(columns={1: 'Energy Unit'}, inplace=True)
        self.daf.rename(columns={2: 'dE/dx Elec.'}, inplace=True)
        self.daf.rename(columns={3: 'dE/dx Nuclear'}, inplace=True)
        self.daf.rename(columns={4: 'Projected Range'}, inplace=True)
        self.daf.rename(columns={6: 'Longitudinal Straggling'}, inplace=True)
        self.daf.rename(columns={8: 'Lateral Straggling'}, inplace=True)

        self.Ion_ene = self.daf['Ion Energy'].tolist()
        self.Units = self.daf['Energy Unit'].tolist()
        self.E_loss = self.daf['dE/dx Elec.'].tolist()
        self.N_loss = self.daf['dE/dx Nuclear'].tolist()
        self.Projected_range = self.daf['Projected Range'].tolist()
        self.Longitudinal_stra = self.daf['Longitudinal Straggling'].tolist()
        self.Lateral_stra = self.daf['Lateral Straggling'].tolist()

        print(self.Ion_ene)


    def fetch_data(self):
        self.vel = []
        # prpotional to velocity
        for i in range(len(self.Ion_ene)):
            if self.Units[i]=="eV":
                self.Ion_ene[i]=float(self.Ion_ene[i])
            elif self.Units[i]=="keV":
                self.Ion_ene[i]=float(self.Ion_ene[i])*1000

            else:
                self.Ion_ene[i] = float(self.Ion_ene[i])
            self.vel.append(np.sqrt(self.Ion_ene[i]))
            self.E_loss[i] = float(self.E_loss[i])

            # self.N_loss[i] = float(self.N_loss[i])*7.0574E-02 # in Mev/(mg/cm2)
            self.N_loss[i] = float(self.N_loss[i] )*(1.2656*10**(-2))# in LSS
            self.Projected_range[i] =float(self.Projected_range[i])
            self.Longitudinal_stra[i] = float(self.Longitudinal_stra[i])
            self.Lateral_stra[i] = float(self.Lateral_stra[i])

    def compare_theory(self):
        self.ep_list = []
        self.ep2_list =[]
        self.sn = []

        for i in range(len(self.Ion_ene)):
            energy = self.Ion_ene[i]
            sn_value = self.sn_func(self.F1*energy)
            self.ep_list.append(self.F1*energy)
            self.ep2_list.append(self.au*0.5*self.Ion_ene[i]/(self.Z_tp**2*self.ev**2))
            self.sn.append(sn_value)
            if energy != 0:
                # value = PI * self.au ** 2 * self.gam * energy * sn_value / (energy*self.F1) # in
                # value = PI * self.au ** 2 * self.gam * energy * sn_value
                value = 8.462*10**(-15)*self.Z_tp**2*0.5*sn_value/(2*self.Z_tp**2)

            else:
                value = 0
            self.SN.append(value*self.Tar_Den*10**(-8)) # change ev/cm to eV/A
            # self.SN.append(value)  # change ev/cm to eV/A
        print("ep",self.ep_list)
        print("ep2", self.ep2_list)
        print("sn",self.sn)

    def sn_func(self, ep):
        return np.log(1+1.1383*ep)/(2*(ep+0.001321*ep**0.21226+0.19593*ep**0.5))

    def plot_data(self):
        print("theo before", self.sn)
        # self.force_fit()
        self.force_fit_LSS()
        plt.plot(self.Ion_ene, self.N_loss, label="experimental data")
        plt.plot(self.Ion_ene, self.sn, label="theoretical curve")
        plt.xlabel("Recoiled energy/eV")
        plt.ylabel("LSS")
        print("exp", self.N_loss)
        print("theo afterwards", self.sn)

        plt.legend()
        plt.show()

    def force_fit(self):
        self.alpha = 137
        self.ratio = []
        for i in range(len(self.N_loss)):
            ratio = self.N_loss[i]/self.SN[i]
            self.ratio.append(ratio)
        self.mean_ratio = sum(self.ratio)/len(self.ratio)
        for i in range(len(self.SN)):
            self.SN[i]= self.alpha*self.SN[i]
        print("mean ratio", self.mean_ratio)
        # I think the lost factor is 1/alpha, which is 137, the mean ration is 127
        print("after edit theo", self.SN)

    def force_fit_LSS(self):
        self.ratio = []
        for i in range(len(self.N_loss)):
            ratio = self.N_loss[i]/self.sn[i]
            self.ratio.append(ratio)
        self.mean_ratio = sum(self.ratio)/len(self.ratio)
        for i in range(len(self.sn)):
            self.sn[i]= self.mean_ratio*self.sn[i]
        print("mean ratio", self.mean_ratio)
        # I think the lost factor is 1/alpha, which is 137, the mean ration is 127
        print("after edit theo", self.sn)

    def output(self, E):
        # E in eV
        sn_value = self.sn_func(self.F1*E)
        # inter_v1 = 8.462*10**(-15)*self.Z_tp**2*0.5*sn_value/(2*self.Z_tp**2)
        inter_v1 = (1/self.LSS_factor) * sn_value
        dedx = inter_v1
        # 137 is atomic fine factor
        return dedx
    # in eV/A
    #







class Model_Test():
    def __init__(self):
        self.Ek_uplim = 500
        self.ER_list = []
        self.Ek_list =[]
        self.Gaussian_sim(20000)
        # self.exp_sim(1000)
        # self.possion_sim(100000)

    def Gaussian_sim(self,N):
        p_list = []
        for i in range(N):
            Ek = np.random.uniform()*self.Ek_uplim
            ER = np.random.normal(0,(self.Ek_uplim-Ek)/20)
            p = np.random.normal(0,3)

            if ER >0:
                self.ER_list.append(ER)
                self.Ek_list.append(Ek)
                # p_list.append(p)
        print(self.ER_list)
        plt.scatter(self.Ek_list,self.ER_list)
        # plt.scatter(self.Ek_list, p_list)
        plt.xlabel("Ek/eV")
        plt.ylabel("ER/eV")
        plt.show()
    def exp_sim(self,N):
        for i in range(N):
            Ek = np.random.uniform()*self.Ek_uplim
            ER = np.random.exponential((self.Ek_uplim-Ek))
            if Ek < 490:
                self.ER_list.append(ER)
                self.Ek_list.append(Ek)
        print(self.ER_list[100])
        plt.scatter(self.Ek_list,self.ER_list)
        plt.xlabel("Ek/eV")
        plt.ylabel("ER/eV")
        plt.show()
    def possion_sim(self,N):
        for i in range(N):
            Ek = np.random.uniform()*self.Ek_uplim
            ER = self.Ek_uplim*np.random.poisson(0.5)
            if Ek < 490:
                self.ER_list.append(ER)
                self.Ek_list.append(Ek)
        plt.scatter(self.Ek_list,self.ER_list)
        plt.xlabel("Ek/eV")
        plt.ylabel("ER/eV")
        plt.show()

class SRIM_scatter():
    def __init__(self):
        # self.min  = 0
        # self.max = 140*1000 # in ev
        # self.step  = 20000 # in ev
        self.min = 0
        self.max = 140 * 10**5  # in ev
        self.step = 200  # in ev

        self.Z_tp= 18
        self.C_tf = (9*PI**2/(2**7))**(1/3)
        self.A_tp =41
        self.Z = (2*self.Z_tp**0.23)**2
        self.EB = 13.6 # in eV
        self.a0 = 5.291*10**(-11) # in m
        # self.a0 = 0.5291  # in A
        self.energy = range(self.min, self.max, self.step)
        self.F1 = self.C_tf*0.5* (0.5/13.6)/(self.Z_tp**2*self.Z**(0.5))
        print(0.3/self.F1)
        self.ep_list=[ self.C_tf*0.5* (i*0.5/13.6)/(self.Z_tp**2*self.Z**(0.5)) for i in range(self.min,self.max,self.step)]
        self.sn = [np.log(1+1.1383*e)/(2*(e+0.001321*e**0.21226+0.19593*e**0.5)) for e in self.ep_list]
        # self.sn = [np.log(e) / (2 * e) for e in self.ep_list]
        self.au = 0.8853*self.a0/(2*self.Z_tp**2)
        self.gam = 4/2**2
        self.SN = []
        for i in range(len(self.energy)):
            if self.ep_list[i] != 0:
                value  = PI* self.au**2 * self.gam * self.energy[i] *self.sn[i]/self.ep_list[i]
                value = PI * self.au ** 2 * self.gam * self.energy[i] * self.sn[i]
            else:
                value = 0
            self.SN.append(value)
        # print(self.sn)
        # plt.plot(self.energy, self.SN)
        plt.plot(self.ep_list, self.sn)
        plt.show()
        print(self.ep_list)
if __name__ == "__main__":
    srim_result  = SRIM_EXY()
    # srim_data = SRIM_Table()
    # test = model_test()
    # scatter = SRIM_scatter()


