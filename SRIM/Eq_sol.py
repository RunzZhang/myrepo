import numpy as np
import matplotlib.pyplot as plt
from sympy import Function, dsolve, Derivative, checkodesol
from sympy.abc import x
from scipy.integrate import solve_ivp
# y = Function('y')
# # Solve the ODE
# result = dsolve(Derivative(y(x), x) + 9*y(x), y(x))
# print(result)
# # Check that the solution is correct
# print(checkodesol(Derivative(y(x), x) + 9*y(x), result))

PI = np.pi


def exponential_decay(t, y): return -0.5 * y

class E_loss_solve():
    def __init__(self):
        self.a = 1.1383
        self.b = 0.01321
        self.c = 0.21226
        self.d = 0.19593
        self.Z_tp = 18
        self.C_tf = (9 * PI ** 2 / (2 ** 7)) ** (1 / 3)
        self.A_tp = 41
        self.Z = (2 * self.Z_tp ** 0.23) ** 2
        self.EB = 13.6  # in eV
        self.a0 = 5.291 * 10 ** (-11)  # in m
        # self.a0 = 0.5291  # in A
        # self.energy = range(self.min, self.max, self.step)
        self.F1 = self.C_tf * 0.5 * (0.5 / 13.6) / (self.Z_tp ** 2 * self.Z ** (0.5))
        print("F1", self.F1)

        self.au = 0.8853 * self.a0 / (2 * self.Z_tp ** 2)  # in m
        self.au = self.au * 100  # in cm
        print("au", self.au)
        self.gam = 4 / 2 ** 2
        self.ev = 1.60218e-19
        self.Tar_Den = 2.1361E+22  # atoms/cm3
        self.mean_ratio =  2.1246181979237e-12 # adjust the formula unit to ev/A
        self.factor = PI * self.au ** 2 * self.gam * self.mean_ratio /(self.F1)
        self.k = 1.541e-05
        self.total_k = 0.0338721
        self.mass = 6.63551406835257e-26 # argon 40 in kg
        self.T_factor = 0.1 # change m/s to A/ns
        self.Tar_Den = 2.1361E+22 # atoms/cm3

        # self.main_fun()
        # self.test()

    def E_loss_t_fun_ODE(self, t, y):
        ep = self.C_tf*0.5* (y*0.5/13.6)/(self.Z_tp**2*self.Z**(0.5))
        part_a = - np.log(1+self.a * ep) / (2 * (ep + self.b * (ep) ** self.c) + self.d * (ep) ** 0.5)* self.factor *self.Tar_Den*10**8
        # part_a = - self.a * ep/ (
        #             2 * (ep + self.b * (ep) ** self.c) + self.d * (ep) ** 0.5) * self.factor * self.Tar_Den * 10 ** 8
        part_b = - (y) ** 0.5 * self.total_k
        part_c = (y) ** 0.5 * np.sqrt(2 * self.ev / self.mass) * (1/self.T_factor) # change dE/dx to dE/dt
        return (part_b+part_a)*part_c

    def E_loss_N_x_fun_ODE(self, t, y):
        #test part a
        ep = self.C_tf*0.5* (y*0.5/13.6)/(self.Z_tp**2*self.Z**(0.5))
        part_a = - np.log(1+self.a * ep) / (2 * (ep + self.b * (ep) ** self.c) + self.d * (ep) ** 0.5)* self.factor *self.Tar_Den*10**8
        # part_a = - self.a * ep/ (
        #             2 * (ep + self.b * (ep) ** self.c) + self.d * (ep) ** 0.5) * self.factor * self.Tar_Den * 10 ** 8

        return part_a

    def E_loss_el_t_fun_ODE(self, t, y):
        # test part c

        part_b = - (y) ** 0.5 * self.total_k
        part_c = (y) ** 0.5 * np.sqrt(2 * self.ev / self.mass) * (1/self.T_factor)  # change dE/dx to dE/dt
        return part_c*part_b

    def E_loss_el_x_fun_ODE(self, t, y):
        # test part b
        # part_b = - (y) ** 0.5 * np.sqrt(2 * 1000*self.ev / self.mass) * self.k
        part_b = - (y) ** 0.5 * self.total_k

        return part_b

    def main_fun(self):
        print("begin")
        self.E_loss_total_t_fun_ODE_posttest()
        # self.E_loss_N_x_fun_ODE_pretest()
        # self.E_loss_el_x_fun_ODE_pretest()
        # self.E_loss_el_t_fun_ODE_posttest()

    def E_loss_result(self, init_E, t):
        #given t in ns and E in ev, return the final energy
        self.last_t = t*10
        self.ini_E = init_E
        if t != 0:
            solve = solve_ivp(self.E_loss_el_t_fun_ODE, [0, self.last_t], [self.ini_E], t_eval=[t]) # check one point's value

            array = solve.y
            sol_y = array[0][0]
        else:
            sol_y = init_E
        return sol_y # return the E value after travels t


    def E_loss_total_t_fun_ODE_posttest(self):

        self.last_t = 0.005
        self.ini_E = 500
        self.threshold_v = 0.1
        self.t_crit = 1343 * 10 ** (-6) * np.log(self.ini_E / 0.1)
        print(self.t_crit)
        # solve = solve_ivp(self.E_loss_el_t_fun_ODE, [0, self.last_t], [self.ini_E], t_eval=[self.t_crit]) # check one point's value
        solve = solve_ivp(self.E_loss_t_fun_ODE, [0, self.last_t], [self.ini_E])
        array = solve.y
        sol_y = array[0]
        print("y", sol_y)
        print("t", solve.t)
        # print("threshold t", np.interp(self.threshold_v, sol_y, solve.t))
        plt.plot(solve.t, sol_y)
        Dy_list = []
        for i in range(1, len(sol_y)):  # double check with solution result
            value = (sol_y[i] - sol_y[i - 1]) / (solve.t[i] - solve.t[i - 1])
            Dy_list.append(value)
        # plt.plot(sol_y[1:], Dy_list)
        # print("Dydt", Dy_list)
        plt.show()

    def E_loss_N_x_fun_ODE_pretest(self):
        # part a works
        E_list = range(0, 500)
        DEDX_list = []
        for i in range(len(E_list)):
            DEDX_list.append(-self.E_loss_N_x_fun_ODE(0,E_list[i]))
        plt.plot(E_list, DEDX_list)
        plt.show()

    def E_loss_N_x_fun_ODE_posttest(self):
        self.last_t = 4000
        self.ini_E = 500
        solve = solve_ivp(self.E_loss_N_x_fun_ODE, [0, self.last_t], [self.ini_E], t_eval=range(0, self.last_t))
        array = solve.y
        sol_y = array[0]
        print("y", sol_y)
        print("t", solve.t)
        # plt.plot(solve.t, sol_y)
        Dy_list = []
        for i in range(1, len(sol_y)):  # double check with solution result
            value = (sol_y[i] - sol_y[i - 1]) / (solve.t[i] - solve.t[i - 1])
            Dy_list.append(value)
        plt.plot(sol_y[1:], Dy_list)
        print("Dydt", Dy_list)
        plt.show()


    def E_loss_el_x_fun_ODE_pretest(self):
        E_list = range(0, 500)
        DEDX_list = []
        for i in range(len(E_list)):
            DEDX_list.append(-self.E_loss_el_x_fun_ODE(0,E_list[i]))
        plt.plot(E_list, DEDX_list)
        plt.show()

    def E_loss_el_x_fun_ODE_posttest(self):
        self.last_t = 4000
        self.ini_E = 500
        solve = solve_ivp(self.E_loss_el_x_fun_ODE, [0, self.last_t], [self.ini_E], t_eval=range(0, self.last_t))
        array = solve.y
        sol_y = array[0]
        print("y", sol_y)
        print("t", solve.t)
        # plt.plot(solve.t, sol_y)
        Dy_list = []
        for i in range(1, len(sol_y)):  # double check with solution result
            value = (sol_y[i] - sol_y[i - 1]) / (solve.t[i] - solve.t[i - 1])
            Dy_list.append(value)
        plt.plot(sol_y[1:], Dy_list)
        print("Dydt", Dy_list)
        plt.show()

    def E_loss_el_t_fun_ODE_posttest(self):
        self.last_t = 0.01
        self.ini_E = 500
        self.t_crit = 1343 * 10 ** (-6) * np.log(self.ini_E / 0.1)
        print(self.t_crit)
        # solve = solve_ivp(self.E_loss_el_t_fun_ODE, [0, self.last_t], [self.ini_E], t_eval=[self.t_crit]) # check one point's value
        solve = solve_ivp(self.E_loss_el_t_fun_ODE, [0, self.last_t], [self.ini_E])
        array = solve.y
        sol_y = array[0]
        print("y", sol_y)
        print("t", solve.t)
        plt.plot(solve.t, sol_y)
        Dy_list = []
        for i in range(1, len(sol_y)):  # double check with solution result
            value = (sol_y[i] - sol_y[i - 1]) / (solve.t[i] - solve.t[i - 1])
            Dy_list.append(value)
        # plt.plot(sol_y[1:], Dy_list)
        # print("Dydt", Dy_list)
        plt.show()


    def test(self):
        sol = solve_ivp(exponential_decay, [0, 10], [2, 4, 8])
        print(sol.t)

        print(sol.y)




if __name__=="__main__":
    solve = E_loss_solve()
