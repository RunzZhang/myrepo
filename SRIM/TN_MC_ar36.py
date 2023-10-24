# to run a new event, please check the event number, output file name, energy chain and time chain
import scipy,random, pickle
import sympy
from sympy import cos, sin, nsolve, Symbol
import numpy as np
from scipy.optimize import fsolve
#variable declaration
#energy in MeV, MeV to SI
Energy_factor = 10**6*1.602*10**(-19)
c = 3*10**8 #in m/s
m_41 = 40.98*10**(-3)/(6.023*10**23) # Ar40 mass in kg
m_37 = 36.97*10**(-3)/(6.023*10**23) # Ar40 mass in kg
tau = 2.6864*10**(-12) # time constant in s of LAr attenuation
time_factor = 10**(-12) # time factor in ps
# energy in Mev and change it into J

E_8790_chain = [8.790 * Energy_factor, 0 * Energy_factor,0 * Energy_factor, 0 * Energy_factor]
t_8790_chain = [0* time_factor, 0 * time_factor, 0 * time_factor, 0 * time_factor]
E_6299_chain = [6.300 * Energy_factor, 2.490 * Energy_factor,0.167 * Energy_factor, 0 * Energy_factor]
t_6299_chain = [0* time_factor, 0.462 * time_factor, 315 * time_factor, 0 * time_factor]
E_5272_chain = [5.272 * Energy_factor, 2.107 * Energy_factor ,1.410 * Energy_factor, 0 * Energy_factor]
t_5272_chain = [0* time_factor, 0.041 * time_factor, 0.59 * time_factor, 0 * time_factor]
E_3700_chain = [3.700 * Energy_factor, 3.679 * Energy_factor,1.1410 * Energy_factor, 0 * Energy_factor]
t_3700_chain = [0* time_factor, 0.010 * time_factor, 0.462 * time_factor, 0 * time_factor]
E_full_chain = [E_8790_chain, E_6299_chain, E_5272_chain, E_3700_chain]
t_full_chain = [t_8790_chain, t_6299_chain, t_5272_chain, t_3700_chain]
run_number = 100000
PI = scipy.pi
print(PI)
# set seed
random.seed(10)

address_8790="/data/runzezhang/result/SRIM_MC/MC_argon36_20231024_8790"
address_6299="/data/runzezhang/result/SRIM_MC/MC_argon36_20231024_6299"
address_5272="/data/runzezhang/result/SRIM_MC/MC_argon36_20231024_5272"
address_3700="/data/runzezhang/result/SRIM_MC/MC_argon36_20231024_3700"

address_full_list =[address_8790,address_6299,address_5272,address_3700]



def two_body_collision_xyz(m,theta, phi, v_xi, v_yi, v_zi, E):
    # theta phi are the gamma's vectors and v_xi, v_yi, v_zi are initial state of LAr
    v_xf=Symbol('v_xf')
    v_yf = Symbol('v_yf')
    v_zf = Symbol('v_zf')
    f1 = m*v_zi - (m * v_zf + E * cos(theta)/c)
    f2 = m*v_xi - (m * v_xf + E * sin(theta) * cos(phi)/c)
    f3 = m*v_yi - (m * v_yf + E * sin(theta) * sin(phi)/c)
    result = nsolve((f1, f2, f3), (v_xf, v_yf, v_zf),(0,0,0))
    # print(result)
    return result

def generate_gamma_vec():
    theta=random.uniform(0,PI)
    phi = random.uniform(0,2*PI)
    return (theta, phi)

def event_generator(E_chain, t_chain,m):
    #given the mass, E energy chain and its time chain, return the recoiled for 1 event
    #generate 1 event
    vx = 0
    vy = 0
    vz = 0
    vx_list = []
    vy_list = []
    vz_list = []
    E_deposit_list =[]
    for i in range(len(E_chain)):
        E_deposit =0.5*m*(vx**2+vy**2+vz**2)*(1-np.exp(-2*t_chain[i]/tau))
        E_deposit_list.append(E_deposit)
        vx=vx * np.exp(-t_chain[i]/tau)
        vy=vy * np.exp(-t_chain[i]/tau)
        vz=vz * np.exp(-t_chain[i] / tau)
        random_angle = generate_gamma_vec()
        result= two_body_collision_xyz(m,random_angle[0],random_angle[1],vx,vy,vz,E_chain[i])
        vx = result[0]
        vy = result[1]
        vz = result[2]
        vx_list.append(vx)
        vy_list.append(vy)
        vz_list.append(vz)
    E_last = 0.5*m*(vx**2+vy**2+vz**2)
    # print(vx_list)
    # print(vy_list)
    # print(vz_list)
    # print(E_deposit_list)
    E_deposit_sum = E_last+sum(E_deposit_list)
    # print("deposit", E_last, E_deposit_sum)

    return E_deposit_sum

def run_generator(N, address, E_chain, t_chain,m):
    #given the address of the file, number of events, E chain and t chain and mass
    E_deposit_list = []
    for i in range(N):
        print(E_chain[0]/Energy_factor,"keV event:",i)
        E_deposit_list.append(event_generator(E_chain, t_chain,m))
    # print(E_deposit_list[:50])
    with open(address, "wb") as fp:  # Pickling
        pickle.dump(E_deposit_list, fp)

    return E_deposit_list

def data_analysis(address):
    with open(address, "rb") as fp:  # Unpickling
        b = pickle.load(fp)

    return b

def full_spectrum_run(address_list ,N):
    for i in range(len(address_list)):
        run_generator(N, address_list[i],E_full_chain[i], t_full_chain[i], m_37)
    print("finished full run")





if __name__ =="__main__":
    # print("E", E1)
    # two_body_collision_xyz(0, 0, 0, 0, 0, 0)

    # event_generator()

    # if data_analysis() ==result:
    #     print("success!")

    #data analysis
    full_spectrum_run(address_full_list, run_number)
