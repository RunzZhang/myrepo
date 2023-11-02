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

def li_func(x, k, b):
    return k * x + b


def li_func2(x, k):
    return k * x



class SRIM():
    def __init__(self):
        super().__init__()
        self.file_name = 'EXYZArgon0.5keV.txt'
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
        self.secondary_variable()



        # self.scatterplot( self.ERs_1d, self.displacement_1d, "ER/eV", "x/A")
        # self.scatterplot(self.ERs_1d, self.time_1d, "ER/eV", "t/s")
        # self.scatterplot( self.ElStop, self.time_1d, "Eloss rate", "t/s")
        # self.time_displacement_hist(self.time_1d, self.displacement_1d, self.ERs_1d, self.vel_1d)
        self.scatterplot(self.Elrecoil, self.Energy, "Erecoil", "Energy")
        # self.selected_hist(0.2,0.1)




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
        print(self.mass)

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
        # print(KE)
        print(len(self.Posx), len(self.Energy), len(self.ElStop), len(self.Elrecoil), len(self.IonNum))
        print(self.Events)
    def secondary_variable(self):



        (self.displacement_1d,self.displacement)= self.fetch_data(len(self.Posx), 0)
        (self.Ek_diff_1d,self.Ek_diff)= self.fetch_data(len(self.Posx), 1)
        (self.El_diff_1d,self.El_diff)= self.fetch_data(len(self.Posx), 2)
        (self.vel_1d, self.vel)= self.fetch_data(len(self.Posx), 3)
        (self.time_1d, self.time)= self.fetch_data(len(self.Posx), 4)
        (self.ERs_1d ,self.ERs) = self.fetch_data(len(self.Posx), 5)

    def fetch_data(self, length, mode):
        # print(self.Posx[:20])
        list_2d = []
        list_1d =[]

        event_pointer = 1
        temp = []
        step_pointer = 0
        for j in range(length):  # all length should be same
            # for j in range(16):
            if self.IonNum[j] == event_pointer:
                if step_pointer == 0:  # row 1
                    temp.append(0)
                    list_1d.append(0)
                else:  # add data
                    # mode0 displacement(m), mode1 kinetic energy loss(eV), mode2 electronic energy loss(eV), mode 3 average velocity(m/s)
                    # mode 4 time(s), # mode 5 energy loss by recoiled collision
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
        for i in self.velocity_dataset:
            fit_data.append(li_func2(i, popt[0]))

        plt.plot(self.velocity_dataset, fit_data, 'r-', label="linear fit")
        plt.xlabel("velocity m/s")
        plt.ylabel("dE/dx eV/A")
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


class model_test():
    def __init__(self):
        self.Ek_uplim = 500
        self.ER_list = []
        self.Ek_list =[]
        # self.Gaussian_sim(20000)
        self.exp_sim(1000)
        # self.possion_sim(100000)

    def Gaussian_sim(self,N):
        for i in range(N):
            Ek = np.random.uniform()*self.Ek_uplim
            ER = np.random.normal(0,(self.Ek_uplim-Ek)/5)
            if ER >0:
                self.ER_list.append(ER)
                self.Ek_list.append(Ek)
        plt.scatter(self.Ek_list,self.ER_list)
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


if __name__ == "__main__":
    # SM =  SRIM()
    test = model_test()


