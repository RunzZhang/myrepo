# to run a new event, please check the event number, output file name, energy chain and time chain
import matplotlib.pyplot as plt
import scipy,random, pickle, Eq_sol
import sympy
from sympy import cos, sin, nsolve, Symbol
import numpy as np
from scipy.optimize import fsolve
#variable declaration
#energy in MeV, MeV to SI

random.seed(10)
PI = scipy.pi
# import warnings
# warnings.filterwarnings("ignore")

class MC_sim_full_argon():
    def __init__(self):

        # non normalized
        #[level energy, lifetime ,{next level_i: [branch ratio, gamma energy, possibility boundary]}]
        # argon init first value is 7000 to be distinguishedfrom other lines
        # the list first value must be unique
        self.runtime = 100000
        self.argon40_weight =97.4/0.93
        self.argon36_weight = 2.5
        self.Energy_factor = 10 ** 3 * 1.602 * 10 ** (-19)
        self.ev=1.60218e-19
        self.c = 3 * 10 ** 8  # in m/s
        self.m_41 = 40.98 * 10 ** (-3) / (6.023 * 10 ** 23)  # Ar40 mass in kg
        self.m_40 = 39.98 * 10 ** (-3) / (6.023 * 10 ** 23)  # Ar40 mass in kg
        self.m_37 = 36.97 * 10 ** (-3) / (6.023 * 10 ** 23)  # Ar40 mass in kg
        self.time_factor = 10**(-3)  # time factor from ps to ns
        # self.address = "/data/runzezhang/result/SRIM_MC/MC_argon_el_full_20231107"
        # self.address = "/data/runzezhang/result/SRIM_MC/MC_argon_full_20231129_6299_-01"
        self.address = "/data/runzezhang/result/SRIM_MC/MC_argon_full_20231204_full"


        # self.argon_init = [10000, 0*self.time_factor ,{6098.9:[93.57*self.argon40_weight,0],3732:[0.121*self.argon40_weight,0],3702.9:[0.474*self.argon40_weight,0],3573:[0.0744*self.argon40_weight,0],
        #                               3564.9:[0.121*self.argon40_weight,0],3278.7:[0.372*self.argon40_weight,0],
        # #                               3111.4:[0.372*self.argon40_weight,0], 8791.2:[100.6*self.argon36_weight, 0]},self.m_41]
        self.argon_init = [10000, 0 * self.time_factor,
                           {6098.9: [93.57 * self.argon40_weight, 0], 3732: [0.121 * self.argon40_weight, 0],
                            3702.9: [0.474 * self.argon40_weight, 0], 3573: [0.0744 * self.argon40_weight, 0],
                             3278.7: [0.372 * self.argon40_weight, 0], 8791.2: [100.6 * self.argon36_weight, 0]},
                           self.m_41]
        # self.argon_init = [10000, 0 * self.time_factor,
        #                    {6098.9: [93.57 * self.argon40_weight, 0], 3732: [0.121 * self.argon40_weight, 0],
        #                     3702.9: [0.474 * self.argon40_weight, 0], 3573: [0.0744 * self.argon40_weight, 0],
        #                      3278.7: [0.372 * self.argon40_weight, 0], 8791.2: [100000000.6 * self.argon36_weight, 0]},
        #                    self.m_41]
        self.level60989 = [6098.9, 0 * self.time_factor,
                           {516.1: [10.8, 5582.0], 1034.7: [0.242, 5063.7], 1353.9: [51.2, 4745.0],
                            2398.1: [9.11, 3700.4], 2693: [0.0744, 3405.3], 2733.4: [3.91, 3365.5],
                            2948.7: [3.72, 3150.2],
                            3009.6: [1.02, 3089.4], 3326.8: [8, 2771.8], 3430.7: [0.474, 2668.1],
                            3968.2: [4.09, 2130.7], 4170.0: [0.93, 1828.8]}, self.m_41]
        # self.level60989 = [6098.9, 0 * self.time_factor,
        #                    {516.1: [10.8, 5582.0], 1034.7: [0.242, 5063.7], 1353.9: [51.2, 4745.0],
        #                     2398.1: [9.11, 3700.4], 2733.4: [3.91, 3365.5],
        #                     2948.7: [3.72, 3150.2],
        #                     3009.6: [1.02, 3089.4], 3326.8: [8, 2771.8], 3430.7: [0.474, 2668.1],
        #                     3968.2: [4.09, 2130.7], 4170.0: [0.93, 1828.8]}, self.m_41]
        self.level42700 = [4270, 0.021*self.time_factor ,{167.3:[0.279,4102.5]},self.m_41]
        # no data of 42700
        self.level39682 = [3968.2, 0.021*self.time_factor ,{516.1:[1.86,3451.8],1353.9:[2.7,2614.3]},self.m_41]
        self.level37320 = [3732.0, 0 *self.time_factor,{167.3:[0.121,3564.5]},self.m_41]
        self.level37029 = [3702.9, 0 *self.time_factor,{1034.7:[0.474, 2668.1]},self.m_41]
        self.level35730 = [3573.0, 0*self.time_factor, {167.3: [0.0744, 3405.5]},self.m_41]
        self.level35649 = [3564.9, 0*self.time_factor, {0: [0.121, 3564.7]},self.m_41]
        self.level33268 = [3326.8, 0.017*self.time_factor ,{516.1:[5.49,2810.5],1034.7:[0.186,2291.6],1353.9:[0.502,1972.6]},self.m_41]
        self.level32787 = [3278.7, 0*self.time_factor, {167.3: [0.372, 3111.3]},self.m_41]
        self.level31114 = [3111.4, 0*self.time_factor, {0: [0.372, 3111.3]},self.m_41]
        self.level30096 = [3009.6, 0.111*self.time_factor ,{167.3:[0.818,2842.5]},self.m_41]
        self.level29487 = [2948.7, 62*0.001*self.time_factor, {167.3: [1.58, 2781.8], 516.1:[0.781,2432.5]},self.m_41]
        self.level27334 = [2733.4, 31*0.001*self.time_factor, {167.3: [2.6, 2566.1]},self.m_41]
        self.level26930 = [2693, 0*self.time_factor, {0: [1, 0]},self.m_41]
        self.level23981 = [2398.1, 0.12*self.time_factor, {167.3: [0.27,2229.5], 516.1:[1.3,1881.5], 1353.9: [5.58, 1044.3]},self.m_41]
        self.level13539 = [1353.9, 0.40*self.time_factor, {0: [2.14,1354.0], 167.3:[48.5,1186.8], 516.1: [8.93, 837.7]},self.m_41]
        self.level10347 = [1034.7, 5*self.time_factor, {167.3 :[1.02,867.3]},self.m_41]
        self.level5161 = [516.1, 260*self.time_factor, {0 :[23.5,516], 167.3:[6.14,348.7]},self.m_41]
        self.level1673 = [167.3, 315*self.time_factor, {0 :[74, 167.3]},self.m_41]
        # argon 36
        self.level87912 = [8791.2, 0*self.time_factor ,{0:[10.9,8790.4],2490.9:[37.5,6299.7],3518.0:[25,5272.6],3938.5:[0.9,4851.8],3981.1:[1.7,4810.3],4448.6:[3.3,4342.3],
                                      4578.7:[0.9,4211.6],4637.6:[1.9,4153],5090.5:[13.4,3700.2],6583.7:[2.7,2207.6], 6826.2:[2.4,1966.7]},self.m_37]

        # self.level87912 = [8791.2, 0 * self.time_factor,
        #                    {0: [10.9, 8790.4], 2490.9: [1000000, 6299.7], 3518.0: [25, 5272.6], 3938.5: [0.9, 4851.8],
        #                     3981.1: [1.7, 4810.3], 4448.6: [3.3, 4342.3],
        #                     4578.7: [0.9, 4211.6], 4637.6: [1.9, 4153], 5090.5: [13.4, 3700.2], 6583.7: [2.7, 2207.6],
        #                     6826.2: [2.4, 1966.7]}, self.m_37]

        self.level68262 = [6826.2, 0.001*self.time_factor, {4578.7 :[5, 2247.9]},self.m_37]
        # nodata 6826
        self.level65837 = [6583.7, 0.001*self.time_factor, {4448.6: [0.2, 2135.3]},self.m_37]
        # nodata 6583
        self.level50905 = [5090.5, 0.01*self.time_factor, {1410.6: [7.9, 3679.3], 2490.9:[3.1,2599.6]},self.m_37]
        self.level46376 = [4637.6, 0.021*self.time_factor, {1410.6: [1, 3226.9],2490.9:[0.2,2145.2]},self.m_37]
        # only 4634
        self.level45787 = [4578.7, 0.014*self.time_factor, {2490.9:[0.3,2087.3]} ,self.m_37]
        # only4573
        self.level44486 = [4578.7, 0.014*self.time_factor, {2490.9: [2, 1957.3]},self.m_37]
        self.level39811 = [3981.1, 0.028*self.time_factor, {0: [1, 3981.1]},self.m_37]
        self.level39385 = [3938.5, 0.017*self.time_factor, {0: [1, 3938]},self.m_37]
        self.level35180 = [3518.0, 0.041*self.time_factor, {1410.6:[23.7,2107.5],2490.9: [2.7, 1026.7]},self.m_37]
        self.level24909 = [2490.9, 0.462*self.time_factor, {0: [57, 2490.6],1611.9:[0.5, 878.5]},self.m_37]
        # self.level24909 = [2490.9, 0.1*0.462 * self.time_factor, {0: [57, 2490.6], 1611.9: [0.5, 878.5]}, self.m_37]
        self.level16119 = [1611.9, 4.38*1000*self.time_factor, {0: [3.4, 1611.7]},self.m_37]
        self.level14106 = [1410.6, 0.59*self.time_factor, {0: [33, 1410.3]},self.m_37]

        self.level0 = [0, 10**6, {0 :[1, 0]}]
        self.argon_list = [self.argon_init,self.level60989, self.level42700, self.level39682, self.level37320 , self.level37029, self.level35730,self.level35649 ,self.level33268 ,self.level32787 ,self.level31114,self.level30096,self.level29487,
                           self.level27334,self.level26930 ,self.level23981 ,self.level13539 ,self.level10347 ,self.level5161  ,self.level1673 ,self.level87912,self.level68262,self.level65837,self.level50905,self.level46376,
                           self.level45787,self.level44486,self.level39811,self.level39385,self.level35180,self.level24909,self.level16119,self.level14106, self.level0]

        self.data_preparation()
        self.gamma_emission_list_1d = []
        self.gamma_emission_list_2d = []
        # self.gamma_sim(10000)
        # self.MC_sim(self.runtime)
        self.data_analysis(self.address)
        # self.plot_pile_up()

    def data_preparation(self):
        for i in range(len(self.argon_list)):# for each chain
            total_BR = 0
            for key in self.argon_list[i][2]:
                # calculate total branch ratio
                total_BR += self.argon_list[i][2][key][0]
                # nomalize the brach ratio
            for key in self.argon_list[i][2]:
                self.argon_list[i][2][key][0]=self.argon_list[i][2][key][0]/total_BR
            # calcualte possiblity boundary
            pointer = 0
            for key in self.argon_list[i][2]:
                pointer += self.argon_list[i][2][key][0]
                self.argon_list[i][2][key].append(pointer)
        print(self.argon_list)
        print("init",self.argon_init)
    def gamma_sim(self, N):
        max_step = 10
        for i in range(N):
            # print(i)
            state = self.argon_init
            # state = self.level60989
            temp_gamma_list = []
            for j in range(max_step):
                if state == self.level0:# # jump out of the loop and go to next event if gamma touch ground state
                    break
                else:# other wise continue the state
                    chance = random.uniform(0, 1)
                    for key in state[2]:
                        if chance <= state[2][key][2] :
                            # success and go to next event
                            if state[2][key][1] > 0:
                                self.gamma_emission_list_1d.append(state[2][key][1])
                                temp_gamma_list.append(state[2][key][1])
                            # search which level match the previous key and change to that state
                            for k in range(len(self.argon_list)):
                                if self.argon_list[k][0] ==key:
                                    state = self.argon_list[k]
                            break
                        else:
                            # move to next
                            continue

            self.gamma_emission_list_2d.append(temp_gamma_list)



        # print("1d",self.gamma_emission_list_1d)
        # print("2d",self.gamma_emission_list_2d)
        plt.hist(self.gamma_emission_list_1d, bins= 500, density = True, color = "red")
        # plt.xlim(0,6000)
        plt.show()

    def MC_sim(self, N):
        max_step = 10
        solve_tool = Eq_sol.E_loss_solve()
        self.E_deposit_1d = []

        for i in range(N):

            print(i)
            state = self.argon_init
            # state = self.level60989
            temp_gamma_list = []
            vx = 0
            vy = 0
            vz = 0
            vx_list = []
            vy_list = []
            vz_list = []
            E_deposit_list = []
            for j in range(max_step):
                if j==1:
                    print("event",i,"leading",state[0])


                if state != self.level0:# otherwise continue the state
                    chance = random.uniform(0, 1)
                    for key in state[2]:
                        # judge which mass it is
                        self.mass = state[3]
                        if chance <= state[2][key][2] :
                            # success and go to next event
                            gamma_energy = state[2][key][1]
                            # print("gamma",gamma_energy)
                            temp_energy = 0.5 * self.mass * (vx ** 2 + vy ** 2 + vz ** 2) /self.ev # energy in ev
                            # print("pre kenit", temp_energy, "t in ns", state[1])
                            # E_final = solve_tool.E_el_loss_result(temp_energy,state[1])
                            E_final = solve_tool.E_loss_result(temp_energy, state[1])
                            # print("post knit", E_final)
                            E_deposit = temp_energy- E_final
                            E_deposit_list.append(E_deposit)
                            if (vx**2+vy**2+vz**2) !=0:
                                velocity = float(np.sqrt(2*E_final*self.ev/self.mass))
                                # print(velocity)
                                # change sympy float to python float
                                vx = float(vx)
                                vy = float(vy)
                                vz = float(vz)
                                # print(np.sqrt(vx**2+vy**2+vz**2))
                                vx = vx * velocity/np.sqrt(vx**2+vy**2+vz**2)
                                vy = vy * velocity/np.sqrt(vx**2+vy**2+vz**2)
                                vz = vz * velocity/np.sqrt(vx**2+vy**2+vz**2)
                            else:
                                vx = 0
                                vy = 0
                                vz = 0
                            random_angle = self.generate_gamma_vec()
                            result = self.two_body_collision_xyz(self.mass, random_angle[0], random_angle[1], vx, vy, vz, gamma_energy*self.Energy_factor)


                            vx = float(result[0])
                            vy = float(result[1])
                            vz = float(result[2])
                            vx_list.append(vx)
                            vy_list.append(vy)
                            vz_list.append(vz)
                            # print("v in m/s", vx, vy, vz)
                            for k in range(len(self.argon_list)): # search which level match the previous key and change to that state
                                if self.argon_list[k][0] == key:
                                    state = self.argon_list[k]
                            break # stopping looping in level dictionary to find which level it goes to next.
                        else:
                            # move to next
                            continue

                elif state == self.level0: #if state is ground, jump to next event
                    E_last = 0.5 * self.mass * (vx ** 2 + vy ** 2 + vz ** 2)/self.ev
                    # print(vx_list)
                    # print(vy_list)
                    # print(vz_list)
                    # print(E_deposit_list)
                    E_deposit_sum = E_last  + sum(E_deposit_list)  # in ev
                    self.E_deposit_1d.append(E_deposit_sum)
                    # print("depositlist", E_deposit_list)
                    # print("deposit", E_last , E_deposit_sum)
                    break

        # print(self.E_deposit_1d)
        with open(self.address, "wb") as fp:  # Pickling
            pickle.dump(self.E_deposit_1d, fp)







    def two_body_collision_xyz(self, m, theta, phi, v_xi, v_yi, v_zi, E):
        # theta phi are the gamma's vectors and v_xi, v_yi, v_zi are initial state of LAr
        v_xf = Symbol('v_xf')
        v_yf = Symbol('v_yf')
        v_zf = Symbol('v_zf')
        f1 = m * v_zi - (m * v_zf + E * cos(theta) / self.c)
        f2 = m * v_xi - (m * v_xf + E * sin(theta) * cos(phi) / self.c)
        f3 = m * v_yi - (m * v_yf + E * sin(theta) * sin(phi) / self.c)
        result = nsolve((f1, f2, f3), (v_xf, v_yf, v_zf), (0, 0, 0))
        # print(result)
        return result

    def generate_gamma_vec(self):
        theta = random.uniform(0, PI)
        phi = random.uniform(0, 2 * PI)
        return (theta, phi)


    def run_generator(self, N, address, E_chain, t_chain, m):
        # given the address of the file, number of events, E chain and t chain and mass
        E_deposit_list = []
        for i in range(N):
            print(E_chain[0] / self.Energy_factor, "keV event:", i)
            E_deposit_list.append(self.event_generator(E_chain, t_chain, m))
        # print(E_deposit_list[:50])
        with open(address, "wb") as fp:  # Pickling
            pickle.dump(E_deposit_list, fp)

        return E_deposit_list

    def data_analysis(self, address):
        start = 0
        end = 1200
        x_bins = []
        with open(self.address, "rb") as fp:  # Unpickling
            MC_full = pickle.load(fp)
            print("read",MC_full)
        bin_n =500

        hist_result = plt.hist(MC_full, bins =bin_n, range=(start, end) ,density = True)
        plt.clf()
        for i in range(len(hist_result[1]) - 1):
            x_bins.append((hist_result[1][i] + hist_result[1][i + 1]) / 2)

        #plot the previous
        total_8spectrum_address = "/data/runzezhang/result/SRIM_MC/MC_argon_8cascades_20231107"
        with open(total_8spectrum_address, "rb") as fp:  # Unpickling
            MC_8 = pickle.load(fp)
            # thedata is stored as [[xbins],[y value]]
            print("read",MC_8)

        plt.plot(x_bins, hist_result[0], color="blue",label= "full chains")
        plt.plot(MC_8[0], MC_8[1], color="orange", label = "8 main chain")
        plt.xlabel("energy/eV")
        plt.ylabel("P")
        plt.yscale("log")
        plt.yticks(fontsize=18)
        plt.xticks(fontsize=18)
        plt.xlim([0, 1200])
        plt.ylim([1E-5,0.1])
        plt.legend()
        plt.show()



    def plot_pile_up(self):
        address1 = "/data/runzezhang/result/SRIM_MC/MC_argon_full_20231129_6299_-01"
        address2 = "/data/runzezhang/result/SRIM_MC/MC_argon_full_20231129_6299_-05"
        address3 = "/data/runzezhang/result/SRIM_MC/MC_argon_full_20231129_6299_00"
        start = 0
        end = 1200
        x1_bins = []
        x2_bins = []
        x3_bins = []
        with open(address1, "rb") as fp:  # Unpickling
            pile1 = pickle.load(fp)
            print("read", pile1)
        with open(address2, "rb") as fp:  # Unpickling
            pile2 = pickle.load(fp)
            print("read", pile2)
        with open(address3, "rb") as fp:  # Unpickling
            pile3 = pickle.load(fp)
            print("read", pile3)
        bin_n = 500

        hist_result1 = plt.hist(pile1, bins=bin_n, range=(start, end), density=True)
        plt.clf()
        for i in range(len(hist_result1[1]) - 1):
            x1_bins.append((hist_result1[1][i] + hist_result1[1][i + 1]) / 2)

        hist_result2 = plt.hist(pile2, bins=bin_n, range=(start, end), density=True)
        plt.clf()
        for i in range(len(hist_result2[1]) - 1):
            x2_bins.append((hist_result2[1][i] + hist_result2[1][i + 1]) / 2)

        hist_result3 = plt.hist(pile3, bins=bin_n, range=(start, end), density=True)
        plt.clf()
        for i in range(len(hist_result3[1]) - 1):
            x3_bins.append((hist_result3[1][i] + hist_result3[1][i + 1]) / 2)

        plt.plot(x1_bins, hist_result1[0], color="blue", label="0.1 lifetime")
        plt.plot(x2_bins, hist_result2[0], color="orange", label="0.5 lifetime ")
        plt.plot(x3_bins, hist_result3[0], color="red", label="origin lifetime")
        plt.xlabel("energy/eV")
        plt.ylabel("P")
        plt.yscale("log")
        plt.yticks(fontsize=18)
        plt.xticks(fontsize=18)
        plt.xlim([0, 1200])
        # plt.ylim([1E-5,0.1])
        plt.legend()
        plt.show()




if __name__ =="__main__":

    #data analysis
    # full_spectrum_run(address_full_list, run_number)
    sim = MC_sim_full_argon()