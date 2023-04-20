import os
class momentum_check():
    def __init__(self, parent=None, join='Q_sigma_Info.txt'):
        # super().__init__(parent)
        print(os.getcwd())
        self.base = os.getcwd()
        #gamma info, [[initialxyz],[second xyz], energy]
        #all in Mev
        self.gamma1_Info=[[],[],]
        self.gamma2_Info=[[],[],]
        self.gamma3_Info = [[], [], ]
        self.gamma4_Info = [[0,0,0], [0,0,0],0 ]
        self.neutron_Info = [[0,0,0],[0,0,0],0]
        self.gamma1_p = [0,0,0,0]
        self.gamma2_p = [0, 0, 0, 0]
        self.gamma3_p = [0, 0, 0, 0]
        self.gamma4_p = [0, 0, 0, 0]
        self.gamma_list = [self.gamma1_p,self.gamma2_p,self.gamma3_p,self.gamma4_p]
        self.neutron_p = [0, 0, 0, 0]
    def get_momentum(self):
        for i in range(0,3):
            self.gamma1_p[i]=self.gamma1_Info[1][i]-self.gamma1_Info[0][i]
            self.gamma2_p[i] = self.gamma2_Info[1][i] - self.gamma2_Info[0][i]
            self.gamma3_p[i] = self.gamma3_Info[1][i] - self.gamma3_Info[0][i]
            self.gamma4_p[i] = self.gamma4_Info[1][i] - self.gamma4_Info[0][i]
            self.neutron_p[i] = self.neutron_Info[1][i] - self.neutron_Info[0][i]
        self.gamma1_p[3]=self.gamma1_Info[2]
        self.gamma2_p[3] = self.gamma2_Info[2]
        self.gamma3_p[3] = self.gamma3_Info[2]
        self.gamma4_p[3] = self.gamma4_Info[2]
        self.neutron_p[3] = self.neutron_Info[2]
    def check_momentum(self):
        neutron_x_s = (2*939*self.neutron_p[3])**0.5*self.neutron_p[0]/(self.neutron_p[0]**2+self.neutron_p[1]**2+self.neutron_p[2]**2)**0.5
        neutron_y_s = (2 * 939 * self.neutron_p[3]) ** 0.5 * self.neutron_p[1] / (
                    self.neutron_p[0] ** 2 + self.neutron_p[1] ** 2 + self.neutron_p[2] ** 2) ** 0.5
        neutron_z_s = (2 * 939 * self.neutron_p[3]) ** 0.5 * self.neutron_p[2] / (
                    self.neutron_p[0] ** 2 + self.neutron_p[1] ** 2 + self.neutron_p[2] ** 2) ** 0.5
        gamma_x_s = 0
        gamma_y_s = 0
        gamma_z_s = 0

        for i in range(0,4):
            gamma_x_s += self.gamma_list[i][0]*self.gamma_list[i][3]/(self.gamma_list[i][0]**2+self.gamma_list[i][1]**2+self.gamma_list[i][2]**2)**0.5
            gamma_y_s += self.gamma_list[i][1] * self.gamma_list[i][3] / (
                        self.gamma_list[i][0] ** 2 + self.gamma_list[i][1] ** 2 + self.gamma_list[i][2] ** 2) ** 0.5
            gamma_z_s += self.gamma_list[i][2] * self.gamma_list[i][3] / (
                        self.gamma_list[i][0] ** 2 + self.gamma_list[i][1] ** 2 + self.gamma_list[i][2] ** 2) ** 0.5

        print("x check", neutron_x_s, gamma_x_s)
        print("y check", neutron_y_s, gamma_y_s)
        print("z check", neutron_z_s, gamma_z_s)



