#translate G4 output file into list

import matplotlib.pyplot as plt
import numpy as np
import os, pickle
import pandas as pd


#This plot the Informacion 2022-01119, which checked the 3 different process of neutron capture and recoiled energy
class thermal_neutron_calibration():
    def __init__(self, parent=None, join='Q_sigma_Info.txt'):
        # super().__init__(parent)
        print(os.getcwd())
        self.base = os.getcwd()


        self.fullInfoaddress = "/data/runzezhang/result/Q_sigma_Info_0418.csv"
        self.captureoutaddress="/data/runzezhang/result/Informacion_20220119_ncout.csv"
        self.data_dic="/data/runzezhang/result/Informacion_20220119_dic.pickle"
        self.Line5582_dic = "/data/runzezhang/5582_202200418_dic.pickle"
        self.Line4745_dic = "/data/runzezhang/result/4745_202200418_dic.pickle"
        self.Line3700_dic = "/data/runzezhang/result/3700_202200418_dic.pickle"
        self.Line1186_dic = "/data/runzezhang/result/1186_202200418_dic.pickle"
        self.G5n1_input = "/data/runzezhang/result/multigamma/Informacion20230706_5n1.csv"
        self.G4n2_input = "/data/runzezhang/result/multigamma/Informacion20230706_4n2.csv"
        self.G3n2n1_input = "/data/runzezhang/result/multigamma/Informacion20230706_3n2n1.csv"
        self.G5n1_output = "/data/runzezhang/result/multigamma/Informacion20230706_5n1.pickle"
        self.G4n2_output = "/data/runzezhang/result/multigamma/Informacion20230706_4n2.pickle"
        self.G3n2n1_output = "/data/runzezhang/result/multigamma/Informacion20230706_3n2n1.pickle"
        self.Multi_G_dic = "/data/runzezhang/result/multigamma/multi_G_dic.pickle"


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

        self.gamma_list=[]
        self.outputfile="/data/runzezhang/result/spectrum_418.pickle"
        self.Qoutputfile = "/data/runzezhang/result/Q_418.pickle"
        self.Routputfile = "/data/runzezhang/result/Recoil_418.pickle"


    def read_Information(self):
        # df = pd.read_csv(self.fullInfoaddress)
        df = pd.read_csv(self.fullInfoaddress)
        print(df.head(5))
        Capture_event=[]

        for idx in range(len(df.index)):
            print(idx)
        # for idx in range(1000):

            try:
                if df.iloc[idx]['Particle'] == "Ar41":
                    if df.iloc[idx]['Event'] not in Capture_event:
                        Capture_event.append(df.iloc[idx]['Event'])

        #         if event_pointer == df.iloc[idx]['Event']:
        #             pass
        # #
                # else:
                #     print(idx)
                #     event_pointer = df.iloc[idx]['Event']
                #     track_history = []
                #     pass
                #
                #
                # if event_pointer in self.INFOmatrix.keys():
                #
                #     if df.iloc[idx]['particle name'] == "gamma" and df.iloc[idx]['Parent ID'] == 1:
                #         self.INFOmatrix[event_pointer][1] = False
                #         # make sure the gamma energy is the initial one with same Track ID and Parent ID.
                #
                #         if df.iloc[idx]['Track ID'] in track_history:
                #
                #             pass
                #         else:
                #             self.INFOmatrix[event_pointer][2] = self.INFOmatrix[event_pointer][2] + df.iloc[idx]['Kinetic E']
                #             track_history.append(df.iloc[idx]['Track ID'])
                #
                #     if df.iloc[idx]['particle name'] == "neutron" and df.iloc[idx]['Physics Process'] == 'Outgoing':
                #         self.INFOmatrix[event_pointer][1] = True
                #
                # else:
                #     if df.iloc[idx]['particle name'] == "neutron" and df.iloc[idx]['Physics Process'] == "Incident":
                #
                #         self.INFOmatrix[event_pointer] = [df.iloc[idx]['Kinetic E'],False,0]
                #
                #     else:
                #         print(df.iloc[idx])
                #         print("Error! First Step is not Neutron!")
                #         break
            except Exception as e:
                print(e)

        df_capture=df[df["Event"].isin(Capture_event)]

        print(df_capture.head(5))
        df_capture.to_csv(path_or_buf=self.captureoutaddress,sep=',',index = False)
        #
        # print(self.INFOmatrix)
        # print(self.INFOmatrix)
        # wdf = pd.DataFrame.from_dict(self.INFOmatrix, orient = 'index', columns=['IE', 'OG', 'Q'])
        # wdf.to_csv(self.outputfile, sep=',')






        event_pointer = -1
        track_pointer = 1
        track_history = []

        # for idx in range(len(df.index)):
        # # for idx in range(1000):
        #
        #     try:
        #         if event_pointer == df.iloc[idx]['Event']:
        #             pass
        # #
        #         else:
        #             print(idx)
        #             event_pointer = df.iloc[idx]['Event']
        #             track_history = []
        #             pass
        #
        #
        #         if event_pointer in self.INFOmatrix.keys():
        #
        #             if df.iloc[idx]['particle name'] == "gamma" and df.iloc[idx]['Parent ID'] == 1:
        #                 self.INFOmatrix[event_pointer][1] = False
        #                 # make sure the gamma energy is the initial one with same Track ID and Parent ID.
        #
        #                 if df.iloc[idx]['Track ID'] in track_history:
        #
        #                     pass
        #                 else:
        #                     self.INFOmatrix[event_pointer][2] = self.INFOmatrix[event_pointer][2] + df.iloc[idx]['Kinetic E']
        #                     track_history.append(df.iloc[idx]['Track ID'])
        #
        #             if df.iloc[idx]['particle name'] == "neutron" and df.iloc[idx]['Physics Process'] == 'Outgoing':
        #                 self.INFOmatrix[event_pointer][1] = True
        #
        #         else:
        #             if df.iloc[idx]['particle name'] == "neutron" and df.iloc[idx]['Physics Process'] == "Incident":
        #
        #                 self.INFOmatrix[event_pointer] = [df.iloc[idx]['Kinetic E'],False,0]
        #
        #             else:
        #                 print(df.iloc[idx])
        #                 print("Error! First Step is not Neutron!")
        #                 break
        #     except Exception as e:
        #         print(e)
        # #
        # # print(self.INFOmatrix)
        # print(self.INFOmatrix)
        # wdf = pd.DataFrame.from_dict(self.INFOmatrix, orient = 'index', columns=['IE', 'OG', 'Q'])
        # wdf.to_csv(self.outputfile, sep=',')
        #

    def read_Scheduled_info(self):
        # still read
        # df = pd.read_csv(self.captureoutaddress, nrows= 20000)
        df = pd.read_csv(self.captureoutaddress)
        print(df.head(5))
        Event_list = df['Event'].unique()
        e_n = 0
        e_p = 0
        p_n = 0
        t_n = 0
        base_time= 0
        time_list = []
        time_modified = []
        # Event|step number|edp|ek|trackID|Time
        ar_info=[]
        gamma_info =[]
        step = 0
        updated_dic={}
        step_p=0
        gamma_step = 0

        for idx in range(len(df.index)):
            if e_p< df.iloc[idx]['Event']:
                # print("start")
                #new event and it should start with ar41
                e_p = df.iloc[idx]['Event']
                step_p = 0
                if df.iloc[idx]['Particle'] == 'Ar41':
                    # print("output")
                    p_n = df.iloc[idx]['ParentID']
                    #Ek| Edop
                    # updated_dic[df.iloc[idx]['Event']]={'Ar41':{step_p:[df.iloc[idx]['Energy_Cinetica'], df.iloc[idx]['Ek'], df.iloc[idx]['TrackID'],
                    #      df.iloc[idx]['Time']]}}updated_dic[df.iloc[idx]['Event']]['Ar41'] =  {step_p: [df.iloc[idx]['Energy_Cinetica'], df.iloc[idx]['Ek'], df.iloc[idx]['TrackID'],
                    #                                           df.iloc[idx]['Time']]}
                    updated_dic[df.iloc[idx]['Event']] = {'Ar41': None, "gamma": []}
                    updated_dic[df.iloc[idx]['Event']]['Ar41'] =  {step_p: [df.iloc[idx]['Energy_Cinetica'], df.iloc[idx]['Ek'], df.iloc[idx]['TrackID'],
                                          df.iloc[idx]['Time']]}

                    base_time = df.iloc[idx]['Time']
                    step_p = step_p + 1

            elif e_p == df.iloc[idx]['Event']:
                # print("same event")
                #first step record it
                if df.iloc[idx]['Particle'] == 'gamma':
                        if df.iloc[idx]['ParentID'] == p_n and gamma_step != df.iloc[idx]['TrackID']:
                            # updated_dic[df.iloc[idx]['Event']]={'gamma': {step_p: [df.iloc[idx]['Energy_Cinetica'], df.iloc[idx]['Ek'], df.iloc[idx]['TrackID'],
                            #                       df.iloc[idx]['Time']-base_time]}}
                            gamma_step = df.iloc[idx]['TrackID']
                            try:
                                updated_dic[df.iloc[idx]['Event']]['gamma'][gamma_step]=  [df.iloc[idx]['Energy_Cinetica'], df.iloc[idx]['Ek'], df.iloc[idx]['TrackID'],df.iloc[idx]['Time']-base_time]
                            except:
                                updated_dic[df.iloc[idx]['Event']] = {'Ar41': None, "gamma": None}
                                updated_dic[df.iloc[idx]['Event']]['gamma'][gamma_step] = [
                                    df.iloc[idx]['Energy_Cinetica'], df.iloc[idx]['Ek'], df.iloc[idx]['TrackID'],
                                    df.iloc[idx]['Time'] - base_time]

                            # print(updated_dic)
                            # updated_dic[df.iloc[idx]['Event']]['gamma'] = 0
                           # # print([df.iloc[idx]['Energy_Cinetica'], df.iloc[idx]['Ek'], df.iloc[idx]['TrackID'],df.iloc[idx]['Time'] - base_time])

            else:
                print(df.iloc[idx]['Event'],"?")

        print(updated_dic)
        with  open(self.data_dic, 'wb') as handle:
            pickle.dump(updated_dic, handle, protocol=pickle.HIGHEST_PROTOCOL)
            print("finished write!")


        return True

    def read_dic(self):
        with open(self.data_dic, 'rb') as handle:
            b = pickle.load(handle)
        print(b)

    def check_Q_value(self):
        with open(self.data_dic, 'rb') as handle:
            b = pickle.load(handle)
            gamma_list=[]
            gamma = 0

        for event in b:
            gamma = 0
            for key in b[event]['gamma']:
                if key > 0:
                    gamma = b[event]['gamma'][key][0] + gamma
            gamma_list.append(gamma)

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('Q Value in G4 Calibration')
        ax.set_xlabel('Q Value/MeV')
        ax.set_ylabel('entries/bin')
        plt.hist(gamma_list, bins=70, range=(0 , 6.5 ), log=True)
        plt.show()


        print(gamma_list)


    def check_time_value(self):
        with open(self.data_dic, 'rb') as handle:
            b = pickle.load(handle)
            time_list=[]
            time = 0

        for event in b:
            time = 0
            for key in b[event]['gamma']:
                if key > 0:
                    time_list.append(b[event]['gamma'][key][3])


        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('Time in G4 Calibration')
        ax.set_xlabel('Time/ns')
        ax.set_ylabel('entries/bin')
        plt.hist(time_list, bins=150, range=(1*10**(-6) , 1.75 ), log=True)
        ax.set_xscale('log')
        plt.show()


        print(time_list)



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
            df=pd.read_csv(self.outputfile)

        self.Q_list= df['Q'].to_list()
        print(self.Q_list)


        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('Q Value in G4 Calibration')
        ax.set_xlabel('Q Value/MeV')
        ax.set_ylabel('entries/bin')
        plt.hist(self.Q_list,bins=70, range=(5.8*10**6,6.5*10**6) ,log= True)
        plt.show()



        return 0
    def plot_cross(self, read = False):
        print(self.outputfile)
        if read:
            df=pd.read_csv(self.outputfile)
        else:
            print("You have to first run read_information() function before plot_cross if read = False. Otherwise might cause error!")
        Incident_matrix = df['IE'].to_list()
        Outgoing_matrix = df[df['OG']==True]['IE'].to_list()

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

    def plot_endf(self):
        fig, axs = plt.subplots(figsize=(9, 3))
        axs.plot(self.endfE1, self.endfsig1,color = 'orange')
        axs.semilogx()
        axs.set_xlim([2*10**(-2),0.1])
        axs.get_xaxis().set_visible(False)
        plt.show()

    def read_Information_spectrum(self):
        # df = pd.read_csv(self.fullInfoaddress)
        df = pd.read_csv(self.fullInfoaddress)
        print(df.head(5))
        index_process=[round(len(df.index)/10),round(len(df.index)/5),round(len(df.index)/2),round(len(df.index)*3/4)]
        event_pointer=0
        track_pointer =[0]
        parent_pointer = [0]
        event_pointer_ar = 0
        track_pointer_ar = [0]
        parent_pointer_ar = [0]
        Capture_event=[]
        Q_pointer =[0]
        Q_list =[]
        Recoil_energy=[]
        recoil_energy_instep=[0]

        for idx in range(len(df.index)):
        # # for idx in range(10):
            if idx in index_process:
                print(idx*100/len(df.index),"%")
        # for idx in range(1000):

            if df.iloc[idx]['particle name'] == "gamma":
                if df.iloc[idx]['Event'] == event_pointer:
                    #old event
                    if df.iloc[idx]['Parent ID']==1:
                        # ar41's gamma
                        if df.iloc[idx]['Track ID'] not in track_pointer:
                            track_pointer.append(df.iloc[idx]['Track ID'])
                            Capture_event.append(df.iloc[idx]['Kinetic E']/1000000)
                            Q_pointer.append(df.iloc[idx]['Kinetic E']/1000000)
                else:
                    #new event
                    event_pointer = df.iloc[idx]['Event']
                    track_pointer = [0]
                    parent_pointer = [0]
                    Q_list.append(sum(Q_pointer))
                    Q_pointer=[]

                    # same as old event
                    if df.iloc[idx]['Parent ID']==1:
                        # ar41's gamma
                        if df.iloc[idx]['Track ID'] not in track_pointer:
                            track_pointer.append(df.iloc[idx]['Track ID'])
                            Capture_event.append(df.iloc[idx]['Kinetic E']/1000000)
                            Q_pointer.append(df.iloc[idx]['Kinetic E'] / 1000000)
            elif df.iloc[idx]['particle name'] == "Ar41":
                # record the max the recoil energy? or first step?
                if df.iloc[idx]['Event'] == event_pointer_ar:
                    #old event
                    if df.iloc[idx]['Parent ID'] == 1:
                        recoil_energy_instep.append(df.iloc[idx]['Recoil E'])


                else:
                    #new event
                    event_pointer_ar = df.iloc[idx]['Event']
                    track_pointer_ar = [0]
                    parent_pointer_ar = [0]

                    Recoil_energy.append(abstract_recoil_des(recoil_energy_instep))
                    recoil_energy_instep = [0]


                    # same as old event
                    if df.iloc[idx]['Parent ID']==1:
                        recoil_energy_instep.append(df.iloc[idx]['Recoil E'])




        with open(self.outputfile, 'wb') as f:
            pickle.dump(Capture_event, f)
        with open(self.Qoutputfile, 'wb') as f:
            pickle.dump(Q_list, f)
        with open(self.Routputfile, 'wb') as f:
            pickle.dump(Recoil_energy, f)
        print("finished saving process")

    def print_spectrum(self):
        Capture_event=[]
        Q_value=[]
        Recoil_value=[]
        with open(self.outputfile, 'rb') as f:
            Capture_event = pickle.load(f)
        with open(self.Qoutputfile, 'rb') as f:
            Q_value = pickle.load(f)
        with open(self.Routputfile, 'rb') as f:
            Recoil_value = pickle.load(f)
        print(Capture_event[0:5])
        print(Q_value[0:5])
        print(Recoil_value[0:5])
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('spectrum')
        ax.set_xlabel('Energy/MeV')
        ax.set_ylabel('entries/bin')
        plt.hist(Capture_event, bins=1000, range=(0, 6.5))
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        ax2.set_title('spectrum')
        ax2.set_xlabel('Energy/MeV')
        ax2.set_ylabel('entries/bin')
        plt.hist(Q_value, bins=1000, range=(0, 6.5))
        fig3 = plt.figure()
        ax3 = fig3.add_subplot(111)
        ax3.set_title('spectrum')
        ax3.set_xlabel('Energy/MeV')
        ax3.set_ylabel('entries/bin')
        plt.hist(Recoil_value, bins=1000, range=(0, 1000))
        plt.show()

    def read_Information_subspectrum(self):
        # df = pd.read_csv(self.fullInfoaddress)
        df = pd.read_csv(self.fullInfoaddress)
        print(df.head(5))
        index_process=[round(len(df.index)/10),round(len(df.index)/5),round(len(df.index)/2),round(len(df.index)*3/4)]
        event_pointer=0
        track_pointer =[0]
        parent_pointer = [0]
        event_pointer_ar = 0
        track_pointer_ar = [0]
        parent_pointer_ar = [0]
        Capture_event=[]
        Q_pointer =[0]
        Line5582 =[]
        Line4745 = []
        Line3700=[]
        Line1186 = []
        Q_list =[]
        Recoil_energy=[]
        recoil_energy_instep=[0]

        for idx in range(len(df.index)):
        # # for idx in range(10):
            if idx in index_process:
                print(idx*100/len(df.index),"%")
        # for idx in range(1000):

            if df.iloc[idx]['particle name'] == "gamma":
                if df.iloc[idx]['Event'] == event_pointer:
                    #old event
                    if df.iloc[idx]['Parent ID']==1:
                        # ar41's gamma
                        if df.iloc[idx]['Track ID'] not in track_pointer:
                            track_pointer.append(df.iloc[idx]['Track ID'])
                            Capture_event.append(df.iloc[idx]['Kinetic E']/1000000)
                            Q_pointer.append(df.iloc[idx]['Kinetic E']/1000000)
                else:
                    #new event
                    event_pointer = df.iloc[idx]['Event']
                    track_pointer = [0]
                    parent_pointer = [0]
                    #check the line
                    for i in Q_pointer:
                        if (i<5.592000 and i>5.572000):
                            #put in line5582
                            Line5582=copy_list(Line5582,Q_pointer)
                        elif (i<4.755000 and i>4.735000):
                            #put in line5582
                            Line4745=copy_list(Line4745,Q_pointer)
                        elif (i<3.710000 and i>3.690000):
                            #put in line5582
                            Line3700=copy_list(Line3700,Q_pointer)
                        elif (i<1.196000 and i>1.176000):
                            #put in line5582
                            Line1186=copy_list(Line1186,Q_pointer)

                    Q_list.append(sum(Q_pointer))
                    Q_pointer=[]

                    # same as old event
                    if df.iloc[idx]['Parent ID']==1:
                        # ar41's gamma
                        if df.iloc[idx]['Track ID'] not in track_pointer:
                            track_pointer.append(df.iloc[idx]['Track ID'])
                            Capture_event.append(df.iloc[idx]['Kinetic E']/1000000)
                            Q_pointer.append(df.iloc[idx]['Kinetic E'] / 1000000)





        with open(self.Line5582_dic, 'wb') as f:
            pickle.dump(Line5582, f)
        with open(self.Line4745_dic, 'wb') as f:
            pickle.dump(Line4745, f)
        with open(self.Line3700_dic, 'wb') as f:
            pickle.dump(Line3700, f)
        with open(self.Line1186_dic, 'wb') as f:
            pickle.dump(Line1186, f)
        with open(self.Qoutputfile,"wb" ) as f:
            pickle.dump(Q_list, f)

        print("finished saving process")

    def read_Information_multiG(self):
        df5n1 = pd.read_csv(self.G5n1_input)
        df4n2 = pd.read_csv(self.G4n2_input)
        df3n2n1 = pd.read_csv(self.G3n2n1_input)
        print(df5n1.head(5))
        print(df4n2.head(5))
        print(df3n2n1.head(5))
        index5n1_process=[round(len(df5n1.index)/10),round(len(df5n1.index)/5),round(len(df5n1.index)/2),round(len(df5n1.index)*3/4)]
        index4n2_process = [round(len(df4n2.index) / 10), round(len(df4n2.index) / 5), round(len(df4n2.index) / 2),
                            round(len(df4n2.index) * 3 / 4)]
        index3n2n1_process = [round(len(df3n2n1.index) / 10), round(len(df3n2n1.index) / 5), round(len(df3n2n1.index) / 2),
                            round(len(df3n2n1.index) * 3 / 4)]
        event_pointer=0
        track_pointer =[0]
        Capture_event=[]
        Q5n1_list =[]
        Q4n2_list = []
        Q3n2n1_list = []
        Recoil_energy=[]
        recoil_energy_instep=[0]
        Q_pointer = []


        for idx in range(len(df5n1.index)):
            if idx in index5n1_process:
                print("5n1",idx*100/len(df5n1.index),"%")
        # for idx in range(1000):

            if df5n1.iloc[idx]['Particle'] == "e-":
                if df5n1.iloc[idx]['Event'] == event_pointer:
                    #old event
                    if df5n1.iloc[idx]['TrackID'] not in track_pointer:
                            track_pointer.append(df5n1.iloc[idx]['TrackID'])
                            Capture_event.append(df5n1.iloc[idx]['Ek']*1000000)
                            Q_pointer.append(df5n1.iloc[idx]['Ek'] * 1000000)

                else:
                    #new event
                    event_pointer = df5n1.iloc[idx]['Event']
                    track_pointer = [0]
                    parent_pointer = [0]

                    Q5n1_list.append(sum(Q_pointer))
                    Q_pointer=[]

                    # same as old event
                    if df5n1.iloc[idx]['TrackID'] not in track_pointer:
                            track_pointer.append(df5n1.iloc[idx]['TrackID'])
                            Capture_event.append(df5n1.iloc[idx]['Ek']*1000000)
                            Q_pointer.append(df5n1.iloc[idx]['Ek'] * 1000000)
        #4n1 gamma group
        for idx in range(len(df4n2.index)):
            if idx in index4n2_process:
                print("4n2", idx * 100 / len(df4n2.index), "%")
            # for idx in range(1000):

            if df4n2.iloc[idx]['Particle'] == "e-":
                if df4n2.iloc[idx]['Event'] == event_pointer:
                    # old event
                    if df4n2.iloc[idx]['TrackID'] not in track_pointer:
                        track_pointer.append(df4n2.iloc[idx]['TrackID'])
                        Capture_event.append(df4n2.iloc[idx]['Ek'] * 1000000)
                        Q_pointer.append(df4n2.iloc[idx]['Ek'] * 1000000)

                else:
                    # new event
                    event_pointer = df4n2.iloc[idx]['Event']
                    track_pointer = [0]
                    parent_pointer = [0]

                    Q4n2_list.append(sum(Q_pointer))
                    Q_pointer = []

                    # same as old event
                    if df4n2.iloc[idx]['TrackID'] not in track_pointer:
                        track_pointer.append(df4n2.iloc[idx]['TrackID'])
                        Capture_event.append(df4n2.iloc[idx]['Ek'] * 1000000)
                        Q_pointer.append(df4n2.iloc[idx]['Ek'] * 1000000)
        #3n2n1 group
        for idx in range(len(df3n2n1.index)):
            if idx in index3n2n1_process:
                print("3n2n1", idx * 100 / len(df3n2n1.index), "%")
            # for idx in range(1000):

            if df3n2n1.iloc[idx]['Particle'] == "e-":
                if df3n2n1.iloc[idx]['Event'] == event_pointer:
                    # old event
                    if df3n2n1.iloc[idx]['TrackID'] not in track_pointer:
                        track_pointer.append(df3n2n1.iloc[idx]['TrackID'])
                        Capture_event.append(df3n2n1.iloc[idx]['Ek'] * 1000000)
                        Q_pointer.append(df3n2n1.iloc[idx]['Ek'] * 1000000)

                else:
                    # new event
                    event_pointer = df3n2n1.iloc[idx]['Event']
                    track_pointer = [0]
                    parent_pointer = [0]

                    Q3n2n1_list.append(sum(Q_pointer))
                    Q_pointer = []

                    # same as old event
                    if df3n2n1.iloc[idx]['TrackID'] not in track_pointer:
                        track_pointer.append(df3n2n1.iloc[idx]['TrackID'])
                        Capture_event.append(df3n2n1.iloc[idx]['Ek'] * 1000000)
                        Q_pointer.append(df3n2n1.iloc[idx]['Ek'] * 1000000)




        with open(self.G5n1_output, 'wb') as f:
            pickle.dump(Q5n1_list, f)
        with open(self.G4n2_output, 'wb') as f:
            pickle.dump(Q4n2_list, f)
        with open(self.G3n2n1_output, 'wb') as f:
            pickle.dump(Q3n2n1_list, f)


        print("finished saving process")

    def print_MGspectrum(self):
        lenth = 25
        bin_n =300
        Line5n1 = []
        Line4n2 =[]
        Line3n2n1 =[]

        Q_list =[]
        with open(self.G5n1_output, 'rb') as f:
            Line5n1 = pickle.load(f)
        with open(self.G4n2_output, 'rb') as f:
            Line4n2 = pickle.load(f)
        with open(self.G3n2n1_output, 'rb') as f:
            Line3n2n1 = pickle.load(f)


        print(Line5n1[0:lenth])
        print(Line4n2[0:lenth])
        print(Line3n2n1[0:lenth])
        # print(Line1186[0:lenth])
        # print(Q_list[0:lenth])
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('spectrum')
        ax.set_xlabel('Energy/eV')
        ax.set_ylabel('entries/bin')
        plt.hist(Line5n1, bins=bin_n, range=(0, 40000))
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        ax2.set_title('spectrum')
        ax2.set_xlabel('Energy/eV')
        ax2.set_ylabel('entries/bin')
        plt.hist(Line4n2, bins=bin_n, range=(0, 40000))
        fig3 = plt.figure()
        ax3 = fig3.add_subplot(111)
        ax3.set_title('spectrum')
        ax3.set_xlabel('Energy/eV')
        ax3.set_ylabel('entries/bin')
        plt.hist(Line3n2n1, bins=bin_n, range=(0, 40000))
        # fig4 = plt.figure()
        # ax4 = fig4.add_subplot(111)
        # ax4.set_title('spectrum')
        # ax4.set_xlabel('Energy/MeV')
        # ax4.set_ylabel('entries/bin')
        # plt.hist(Line1186, bins=bin_n, range=(0, 6.5))
        # fig5 = plt.figure()
        # ax5 = fig4.add_subplot(111)
        # ax5.set_title('spectrum')
        # ax5.set_xlabel('Energy/MeV')
        # ax5.set_ylabel('entries/bin')
        # plt.hist(Q_list, bins=bin_n, range=(0, 6.5))

        plt.show()

    def print_MGspectrum_LY(self):
        lenth = 25
        bin_n =300
        Line5n1 = []
        Line4n2 =[]
        Line3n2n1 =[]

        Q_list =[]
        with open(self.G5n1_output, 'rb') as f:
            Line5n1 = pickle.load(f)
        with open(self.G4n2_output, 'rb') as f:
            Line4n2 = pickle.load(f)
        with open(self.G3n2n1_output, 'rb') as f:
            Line3n2n1 = pickle.load(f)

        for i in Line5n1:
            i = i*40/1000
        for i in Line4n2:
            i = i*40/1000
        for i in Line3n2n1:
            i = i*40/1000


        print(Line5n1[0:lenth])
        print(Line4n2[0:lenth])
        print(Line3n2n1[0:lenth])
        # print(Line1186[0:lenth])
        # print(Q_list[0:lenth])
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('spectrum')
        ax.set_xlabel('photonN')
        ax.set_ylabel('entries/bin')
        plt.hist(Line5n1, bins=bin_n, range=(0, 40000))
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        ax2.set_title('spectrum')
        ax2.set_xlabel('photonN')
        ax2.set_ylabel('entries/bin')
        plt.hist(Line4n2, bins=bin_n, range=(0, 40000))
        fig3 = plt.figure()
        ax3 = fig3.add_subplot(111)
        ax3.set_title('spectrum')
        ax3.set_xlabel('photonN')
        ax3.set_ylabel('entries/bin')
        plt.hist(Line3n2n1, bins=bin_n, range=(0, 40000))
        # fig4 = plt.figure()
        # ax4 = fig4.add_subplot(111)
        # ax4.set_title('spectrum')
        # ax4.set_xlabel('Energy/MeV')
        # ax4.set_ylabel('entries/bin')
        # plt.hist(Line1186, bins=bin_n, range=(0, 6.5))
        # fig5 = plt.figure()
        # ax5 = fig4.add_subplot(111)
        # ax5.set_title('spectrum')
        # ax5.set_xlabel('Energy/MeV')
        # ax5.set_ylabel('entries/bin')
        # plt.hist(Q_list, bins=bin_n, range=(0, 6.5))

        plt.show()

    def print_spectrum(self):
        Capture_event=[]
        Q_value=[]
        Recoil_value=[]
        with open(self.outputfile, 'rb') as f:
            Capture_event = pickle.load(f)
        with open(self.Qoutputfile, 'rb') as f:
            Q_value = pickle.load(f)
        with open(self.Routputfile, 'rb') as f:
            Recoil_value = pickle.load(f)


        print(Capture_event[0:5])
        print(Q_value[0:5])
        print(Recoil_value[0:5])
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('spectrum')
        ax.set_xlabel('Energy/MeV')
        ax.set_ylabel('entries/bin')
        plt.hist(Capture_event, bins=1000, range=(0, 6.5))
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        ax2.set_title('spectrum')
        ax2.set_xlabel('Energy/MeV')
        ax2.set_ylabel('entries/bin')
        plt.hist(Q_value, bins=1000, range=(0, 6.5))
        fig3 = plt.figure()
        ax3 = fig3.add_subplot(111)
        ax3.set_title('spectrum')
        ax3.set_xlabel('Energy/MeV')
        ax3.set_ylabel('entries/bin')
        plt.hist(Recoil_value, bins=1000, range=(0, 1000))

        plt.show()

    def print_subspectrum(self):
        lenth = 25
        bin_n =300
        Line5582 = []
        Line4745 =[]
        Line3700 =[]
        Line1186 = []
        Q_list =[]
        with open(self.Line5582_dic, 'rb') as f:
            Line5582 = pickle.load(f)
        with open(self.Line4745_dic, 'rb') as f:
            Line4745 = pickle.load(f)
        with open(self.Line3700_dic, 'rb') as f:
            Line3700 = pickle.load(f)
        with open(self.Line1186_dic, 'rb') as f:
            Line1186 = pickle.load(f)
        with open(self.Qoutputfile, 'rb') as f:
            Q_list= pickle.load(f)

        print(Line5582[0:lenth])
        print(Line4745[0:lenth])
        print(Line3700[0:lenth])
        print(Line1186[0:lenth])
        print(Q_list[0:lenth])
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('spectrum')
        ax.set_xlabel('Energy/MeV')
        ax.set_ylabel('entries/bin')
        plt.hist(Line5582, bins=bin_n, range=(0, 6.5))
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        ax2.set_title('spectrum')
        ax2.set_xlabel('Energy/MeV')
        ax2.set_ylabel('entries/bin')
        plt.hist(Line4745, bins=bin_n, range=(0, 6.5))
        fig3 = plt.figure()
        ax3 = fig3.add_subplot(111)
        ax3.set_title('spectrum')
        ax3.set_xlabel('Energy/MeV')
        ax3.set_ylabel('entries/bin')
        plt.hist(Line3700, bins=bin_n, range=(0, 6.5))
        fig4 = plt.figure()
        ax4 = fig4.add_subplot(111)
        ax4.set_title('spectrum')
        ax4.set_xlabel('Energy/MeV')
        ax4.set_ylabel('entries/bin')
        plt.hist(Line1186, bins=bin_n, range=(0, 6.5))
        fig5 = plt.figure()
        ax5 = fig4.add_subplot(111)
        ax5.set_title('spectrum')
        ax5.set_xlabel('Energy/MeV')
        ax5.set_ylabel('entries/bin')
        plt.hist(Q_list, bins=bin_n, range=(0, 6.5))

        plt.show()
    def check_recoil_interaction(self):
        # df = pd.read_csv(self.fullInfoaddress)
        df = pd.read_csv(self.fullInfoaddress)
        print(df.head(5))
        index_process=[round(len(df.index)/10),round(len(df.index)/5),round(len(df.index)/2),round(len(df.index)*3/4)]
        event_pointer=0
        track_pointer =[0]
        parent_pointer = [0]
        event_pointer_ar = 0
        track_pointer_ar = [0]
        parent_pointer_ar = [0]
        Capture_event=[]
        Q_pointer =[0]
        Q_list =[]
        Recoil_energy=[]
        recoil_energy_instep=[0]
        ar_action_length=[]

        for idx in range(len(df.index)):
        # # for idx in range(10):
            if idx in index_process:
                print(idx*100/len(df.index),"%")
        # for idx in range(1000):


            elif df.iloc[idx]['particle name'] == "Ar41":
                # record the max the recoil energy? or first step?
                if df.iloc[idx]['Event'] == event_pointer_ar:
                    #old event
                    if df.iloc[idx]['Parent ID'] == 1:
                        recoil_energy_instep.append(df.iloc[idx]['Recoil E'])


                else:
                    #new event
                    event_pointer_ar = df.iloc[idx]['Event']
                    track_pointer_ar = [0]
                    parent_pointer_ar = [0]

                    Recoil_energy.append(abstract_recoil_des(recoil_energy_instep))
                    ar_action_length.append(len(recoil_energy_instep))
                    recoil_energy_instep = [0]


                    # same as old event
                    if df.iloc[idx]['Parent ID']==1:
                        recoil_energy_instep.append(df.iloc[idx]['Recoil E'])

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('spectrum')
        ax.set_xlabel('Energy/MeV')
        ax.set_ylabel('entries/bin')
        plt.hist(ar_action_length, bins=10, range=(0, 10))
        plt.show()
        print(ar_action_length)

        print("finished saving process")

def abstract_recoil(list):
    diff_list=[]
    for i in range(len(list)):
        if i<=len(list)-2:
            if list[i+1]>list[i]:
                diff_list.append(list[i+1]-list[i])
            else:
                diff_list.append(0)
    recoil=sum(diff_list)
    return recoil

def abstract_recoil_des(list):
    diff_list=[]
    for i in range(len(list)):
        if i<=len(list)-2:
            if list[i+1]<list[i]:
                diff_list.append(list[i]-list[i+1])
            else:
                diff_list.append(0)
    diff_list.append(list[-1])
    recoil=sum(diff_list)
    return recoil

def copy_list(list1, list2):
    for i in list2:
        list1.append(i)
    return list1


if __name__=="__main__":
    # get hits number and plot positions
    # print(abstract_recoil_des([1,3,5,4,7,6,5]))
    tnc= thermal_neutron_calibration()
    # tnc.read_Information_subspectrum()
    # tnc.read_Information_multiG()
    # tnc.print_MGspectrum_LY()
    # tnc.print_subspectrum()
    # tnc.check_recoil_interaction()
    # print(copy_list([1],[123,321,343]))
    # tnc.read_Information_spectrum()
    # tnc.print_spectrum()
    # tnc.read_Information()
    # tnc.plot_Q(True)
    # tnc.read_Scheduled_info()
    # tnc.read_dic()
    # tnc.check_Q_value()
    # tnc.check_time_value()
    # tnc.plot_cross(read=True)
    # tnc.read_gamma_Information()
    tnc.plot_gamma(True)
    tnc.plot_endf()