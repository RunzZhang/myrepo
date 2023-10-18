import scipy,random, pickle
import sympy
from sympy import cos, sin, nsolve, Symbol
import numpy as np
from scipy.optimize import fsolve
#variable declaration
#energy in MeV, MeV to SI
Energy_factor = 10**6*1.602*10**(-19)
c = 3*10**8 #in m/s
m = 39.948*10**(-3)/(6.023*10**23) # Ar40 mass in kg
tau = 2.6864*10**(-12) # time constant in s of LAr attenuation
time_factor = 10**(-12) # time factor in fs
# energy in Mev and change it into J
E1 = 1 * Energy_factor
E2 = 1 * Energy_factor
E3 = 1 * Energy_factor
E4 = 1 * Energy_factor
t1=  1* time_factor
t2 = 1* time_factor
t3 = 1* time_factor
t4 = 1* time_factor
E_chain = [E1, E2, E3, E4]
t_chain = [t1, t2, t3, t4]
run_number = 1000
PI = scipy.pi
print(PI)
# set seed
random.seed(10)
address ="/data/runzezhang/result/SRIM_MC/MC_20231018"

def two_body_collision_xyz(theta, phi, v_xi, v_yi, v_zi, E):
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

def event_generator():
    #generate 1 event
    vx = 0
    vy = 0
    vz = 0
    vx_list = []
    vy_list = []
    vz_list = []
    E_deposit_list =[]
    for i in range(len(E_chain)):
        E_deposit =0.5*m*(vx**2+vy**2+vz**2)*(1-np.exp(-t_chain[i]/tau))
        E_deposit_list.append(E_deposit)
        vx=vx * np.exp(-t_chain[i]/tau)
        vy=vy * np.exp(-t_chain[i]/tau)
        vz=vz * np.exp(-t_chain[i] / tau)
        random_angle = generate_gamma_vec()
        result= two_body_collision_xyz(random_angle[0],random_angle[1],vx,vy,vz,E_chain[i])
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

def run_generator(N):
    E_deposit_list = []
    for i in range(N):
        print(i)
        E_deposit_list.append(event_generator())
    # print(E_deposit_list[:50])
    with open("test", "wb") as fp:  # Pickling
        pickle.dump(E_deposit_list, fp)

    return E_deposit_list

def data_analysis():
    with open("test", "rb") as fp:  # Unpickling
        b = pickle.load(fp)



if __name__ =="__main__":
    # print("E", E1)
    # two_body_collision_xyz(0, 0, 0, 0, 0, 0)

    # event_generator()
    Number=60
    result = run_generator(Number)
    result2 = data_analysis()
    if result == result2:
        print("save successfully!")