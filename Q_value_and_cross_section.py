#translate G4 output file into list

import matplotlib.pyplot as plt
import numpy as np
import os, pickle
import pandas as pd

class thermal_neutron_calibration():
    def __init__(self, parent=None, join='Q_sigma_Info.txt'):
        # super().__init__(parent)
        print(os.getcwd())
        self.base = os.getcwd()

        self.Infoaddress = 'Q_sigma_Info_Mar30.txt'
        self.fullInfoaddress = os.path.join(self.base, self.Infoaddress)
        self.outputfile = 'result_Mar30.bin'
        self.gammafile = 'gammaresult_Mar30.bin'
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
        print(self.sigma_matrix)
        self.gamma_list=[]



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




    def read_gamma_Information(self):
        file = open(self.fullInfoaddress, 'r')


        event_pointer=-1
        track_pointer = 1
        track_history =[]


        lines = file.readlines()
        for line in lines:
            line_list = line.split()
            try:
                # print(line_list[0])
                if line_list[0]=="Event":
                    continue
                if event_pointer == int(line_list[0]):
                    pass

                elif event_pointer +1 == int(line_list[0]) :
                    event_pointer += 1
                    track_history= []
                    pass
                elif event_pointer +2 == int(line_list[0]) :
                    print(event_pointer, 'is lost. Skip and continue.')
                    event_pointer += 2
                    track_history= []
                    pass
                else:
                    print("Pointer Not Matched! In line ",int(line_list[0]))
                    break
                    # pass




                if line_list[1]=="gamma" and line_list[3]=='1':

                    if line_list[2] in track_history:
                        pass
                    else:
                        if line_list[5] == "MeV":
                            self.gamma_list.append(float(line_list[4]))

                        elif line_list[5] == "keV":
                            self.gamma_list.append(float(line_list[4])/1000)

                        elif line_list[5] == "eV":
                            self.gamma_list.append(float(line_list[4]) / 1000000)
                        else:
                            print("Unit not recognized")
                            break
                        track_history.append(line_list[2])





            except Exception as e:
                print(e)

        # print(self.INFOmatrix)
        with open(self.gammafile, 'wb') as fp:
            pickle.dump(self.gamma_list, fp)

    def read_Information(self):
        file = open(self.fullInfoaddress, 'r')

        event_pointer = -1
        track_pointer = 1
        track_history = []

        lines = file.readlines()
        for line in lines:
            line_list = line.split()
            try:
                # print(line_list[0])
                if line_list[0] == "Event":
                    continue
                if event_pointer == int(line_list[0]):
                    pass

                elif event_pointer + 1 == int(line_list[0]):
                    event_pointer += 1
                    track_history = []
                    pass
                elif event_pointer + 2 == int(line_list[0]):
                    print(event_pointer, 'is lost. Skip and continue.')
                    event_pointer += 2
                    track_history = []
                    pass
                else:
                    print("Pointer Not Matched! In line ", int(line_list[0]))
                    break
                    # pass

                if int(line_list[0]) in self.INFOmatrix.keys():

                    if line_list[1] == "gamma" and line_list[3] == '1':
                        self.INFOmatrix[event_pointer]['OG'] = False
                        if line_list[2] in track_history:
                            pass
                        else:
                            if line_list[5] == "MeV":
                                self.INFOmatrix[event_pointer]['Q'] = self.INFOmatrix[event_pointer]['Q'] + float(line_list[4])

                            elif line_list[5] == "keV":
                                self.INFOmatrix[event_pointer]['Q'] = self.INFOmatrix[event_pointer]['Q'] + float(line_list[4]) / 1000
                            elif line_list[5] == "eV":
                                self.INFOmatrix[event_pointer]['Q'] = self.INFOmatrix[event_pointer]['Q'] + float(line_list[4]) / 1000000
                            else:
                                print("Unit not recognized")
                                break
                            track_history.append(line_list[2])
                    if line_list[1] == "neutron" and line_list[-1] == 'Outgoing':
                        self.INFOmatrix[event_pointer]['OG'] = True





                else:
                    if line_list[1] == "neutron" and line_list[-1] == 'Incident':
                        if line_list[5] == "MeV":
                            self.INFOmatrix[event_pointer] = {'IE': float(line_list[4]) * 1000000, 'OG': False, 'Q': 0}
                            # self.INFOmatrix[event_pointer]['IE'] = self.INFOmatrix[event_pointer]['Q'] + float(line_list[4])

                        elif line_list[5] == "keV":
                            self.INFOmatrix[event_pointer] = {'IE': float(line_list[4]) * 1000, 'OG': False, 'Q': 0}
                            # self.INFOmatrix[event_pointer]['Q'] = self.INFOmatrix[event_pointer]['Q'] + float(line_list[4]) / 1000
                        elif line_list[5] == "eV":
                            self.INFOmatrix[event_pointer] = {'IE': float(line_list[4]), 'OG': False, 'Q': 0}
                            # self.INFOmatrix[event_pointer]['Q'] = self.INFOmatrix[event_pointer]['Q'] + float(line_list[4]) / 1000000
                        else:
                            print("Unit not recognized")
                            break
                        # if line_list[5] == "MeV":
                        #     self.INFOmatrix[event_pointer] = {'IE': float(line_list[4])*1000000, 'OG': True, 'Q': 0}
                        #     # self.INFOmatrix[event_pointer]['IE'] = self.INFOmatrix[event_pointer]['Q'] + float(line_list[4])
                        #
                        # elif line_list[5] == "keV":
                        #     self.INFOmatrix[event_pointer] = {'IE': float(line_list[4])*1000, 'OG': True, 'Q': 0}
                        #     # self.INFOmatrix[event_pointer]['Q'] = self.INFOmatrix[event_pointer]['Q'] + float(line_list[4]) / 1000
                        # elif line_list[5] == "eV":
                        #     self.INFOmatrix[event_pointer] = {'IE': float(line_list[4]), 'OG': True, 'Q': 0}
                        #     # self.INFOmatrix[event_pointer]['Q'] = self.INFOmatrix[event_pointer]['Q'] + float(line_list[4]) / 1000000
                        # else:
                        #     print("Unit not recognized")
                        #     break

                    else:
                        print(line_list)
                        print("Error! First Step is not Neutron!")
                        break
            except Exception as e:
                print(e)

        # print(self.INFOmatrix)
        with open(self.outputfile, 'wb') as fp:
            pickle.dump(self.INFOmatrix, fp)

    def read_Information_alln(self):
        file = open(self.fullInfoaddress, 'r')

        event_pointer = -1
        track_pointer = 1
        track_history = []

        lines = file.readlines()
        for line in lines:
            line_list = line.split()
            try:
                # print(line_list[0])
                if line_list[0] == "Event":
                    continue
                if event_pointer == int(line_list[0]):
                    pass

                elif event_pointer + 1 == int(line_list[0]):
                    event_pointer += 1
                    track_history = []
                    pass
                else:
                    print("Pointer Not Matched! In line ", int(line_list[0]))
                    break

                if int(line_list[0]) in self.INFOmatrix.keys():


                    if line_list[1] == "gamma" and line_list[3] == '1':
                        self.INFOmatrix[event_pointer]['OG'] = False
                        if line_list[2] in track_history:
                            pass
                        else:
                            if line_list[5] == "MeV":
                                self.INFOmatrix[event_pointer]['Q'] = self.INFOmatrix[event_pointer]['Q'] + float(line_list[4])

                            elif line_list[5] == "keV":
                                self.INFOmatrix[event_pointer]['Q'] = self.INFOmatrix[event_pointer]['Q'] + float(line_list[4]) / 1000
                            elif line_list[5] == "eV":
                                self.INFOmatrix[event_pointer]['Q'] = self.INFOmatrix[event_pointer]['Q'] + float(line_list[4]) / 1000000
                            else:
                                print("Unit not recognized")
                                break
                            track_history.append(line_list[2])

                else:
                    if line_list[1] == "neutron":
                        self.INFOmatrix[event_pointer] = {'IE': line_list[4] + " " + line_list[5], 'OG': True, 'Q': 0}
                    else:
                        print(line_list)

                        print("Error! First Step is not Neutron!")
                        break
            except Exception as e:
                print(e)

            with open(self.outputfile, 'wb') as fp:
                pickle.dump(self.INFOmatrix, fp)

    def plot_gamma(self, read=False):
        # print(self.INFOmatrix)

        if read:
            with open(self.gammafile, 'rb') as fp:
                b = pickle.load(fp)
                self.gamma_list = b

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('Q Value in G4 Calibration')
        ax.set_xlabel('Q Value/MeV')
        ax.set_ylabel('entries/bin')
        plt.hist(self.gamma_list, bins=1000, range=(0, 6.5), log=True)
        plt.show()

        # print(self.INFOmatrix)
    def plot_Q(self, read =False):
        # print(self.INFOmatrix)

        if read:
            with open(self.outputfile, 'rb') as fp:
                b = pickle.load(fp)
                self.INFOmatrix = b


        for key in self.INFOmatrix:
            self.Q_list.append(self.INFOmatrix[key]["Q"])
            if self.INFOmatrix[key]["Q"]<6 and self.INFOmatrix[key]["Q"]>0:
                self.abnormal_event.append(key)
            # if self.INFOmatrix[key]["Q"] !=0:
            #     print(True)
        print(self.abnormal_event)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('Q Value in G4 Calibration')
        ax.set_xlabel('Q Value/MeV')
        ax.set_ylabel('entries/bin')
        plt.hist(self.Q_list,bins=70, range=(5.8,6.5) ,log= True)
        plt.show()



                # print(int(line_list[0]))
                # raw_counts = raw_counts + 1
                #
                # if line_list[-1] in sipm_matrix.keys():
                #     sipm_matrix[line_list[-1]] = sipm_matrix[line_list[-1]] + 1
                # else:
                #     sipm_matrix[line_list[-1]] = 1
                #
                # if line_list[-1] in ['SIPM1_Si_phys', 'SIPM2_Si_phys', 'SIPM3_Si_phys', 'SIPM4_Si_phys', 'SIPM5_Si_phys']:
                #     count = count + 1
                #     if line_list[7] == 'um':
                #         x.append(float(line_list[6]) / 1000)
                #     elif line_list[7] == 'mm':
                #         x.append(float(line_list[6]))
                #     elif line_list[7] == 'cm':
                #         x.append(float(line_list[6]) * 10)
                #     else:
                #         print('x unit', line_list[7])
                #
                #     if line_list[9] == 'um':
                #         y.append(float(line_list[8]) / 1000)
                #     elif line_list[9] == 'mm':
                #         y.append(float(line_list[8]))
                #     elif line_list[9] == 'cm':
                #         y.append(float(line_list[8]) * 10)
                #     else:
                #         pass
                #
                #     if line_list[11] == 'um':
                #         z.append(float(line_list[10]) / 1000)
                #     elif line_list[11] == 'mm':
                #         z.append(float(line_list[10]))
                #     elif line_list[11] == 'cm':
                #         z.append(float(line_list[10]) * 10)
                #     else:
                #         pass


        return 0
    def plot_cross(self, read = False):
        print(self.outputfile)
        if read:
            with open(self.outputfile, 'rb') as fp:
                b = pickle.load(fp)
                self.INFOmatrix = b
                print(len(self.INFOmatrix))
        else:
            print("You have to first run read_information() function before plot_cross if read = False. Otherwise might cause error!")
        Incident_matrix = []
        Outgoing_matrix = []
        for key in self.INFOmatrix:
            if self.INFOmatrix[key]['OG']:
                Incident_matrix.append(self.INFOmatrix[key]['IE'])
                Outgoing_matrix.append(self.INFOmatrix[key]['IE'])
            else:
                Incident_matrix.append(self.INFOmatrix[key]['IE'])
        print(Incident_matrix, '\n', Outgoing_matrix)

        fig, axs = plt.subplots(2, 2, figsize=(9, 3), sharex=True)
        (n0, bins0, patches0) = axs[0, 0].hist(Incident_matrix, bins=np.logspace(-2, -1, num=40))
        axs[0, 0].set_title('Incident Neutron Energy Distribution')
        axs[0, 0].grid(which='both')
        axs[0, 0].set_xlabel('Incident Neutron Energy/eV')
        axs[0, 0].set_ylabel('Numbers/bin')
        axs[0, 0].semilogx()
        (n1, bins1, patches1) = axs[0, 1].hist(Outgoing_matrix, bins=np.logspace(-2, -1, num=40))
        axs[0, 1].set_title('Outgoing Neutron Energy Distribution')
        axs[0, 1].grid(which='both')
        axs[0, 1].set_xlabel('Outgoing Neutron Energy/eV')
        axs[0, 1].set_ylabel('Numbers/bin')
        axs[0, 1].semilogx()
        plt.xscale('log')
        # axs[0,0].xscale('log')
        # axs[0,1].xscale('log')
        # axs[0].hist(Incident_matrix)
        # axs[1].hist(Outgoing_matrix)

        print('n0', len(n0), n0)
        print('n1', len(n1), n1)
        print('bins', len(bins0), bins0)

        # create new bins for x value in the following steps:
        bin_x= []
        for i in range(len(bins0)-1):
            bin_x.append((bins0[i]+bins0[i+1])/2)
        # print('patches0',len(patches0),patches0)
        ratio = []
        error = []

        for n in range(len(n0)):
            ratio_el = self.coefi * (n0[n] - n1[n]) / n0[n]
            error_el =  (self.coefi - ratio_el) * (1 / (n0[n]) ** (1 / 2) + 1 / (n0[n]) ** (1 / 2))
            # error_el = self.coefi * (1 - ratio_el) * (1 / (n0[n]) ** (1 / 2) + 1 / (n0[n]) ** (1 / 2))/ratio_el
            ratio.append(ratio_el)
            error.append(error_el)
        print('ratios', len(ratio), ratio)
        # axs[1,0].plot(bins0[:39],ratio)
        axs[1, 0].errorbar(bin_x, ratio, yerr=error, ecolor='red', label = 'Simulation Data')
        # axs[1, 0].set_xticks(np.arange(0.1,0.01,step=0.01))
        # axs[1, 0].set_xticklabels([2 * 10 ** (-2), 3 * 10 ** (-2)])
        axs[1, 0].set_title('Cross Section of Neutron-Ar Capture')
        axs[1, 0].grid(which='both')
        axs[1, 0].set_xlabel('Incident Neutron Energy/eV')
        axs[1, 0].set_ylabel('Cross Section/mb')

        axs[1, 0].semilogx()

        axs[1,0].plot(self.endfE1, self.endfsig1, label = "ENDF DATA" , color = 'orange')
        axs[1,0].legend()

        plt.show()



if __name__=="__main__":
    # get hits number and plot positions

    tnc= thermal_neutron_calibration()
    # tnc.read_Information()
    # tnc.plot_Q(True)
    tnc.plot_cross(read=True)
    # tnc.read_gamma_Information()
    # tnc.plot_gamma(True)