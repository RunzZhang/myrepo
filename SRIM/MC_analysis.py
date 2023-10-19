import scipy,random, pickle
import sympy
from sympy import cos, sin, nsolve, Symbol
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter

from scipy.optimize import fsolve
e = 1.602*10**(-19)
address1 ="/data/runzezhang/result/SRIM_MC/MC_20231018_5582"
address2 ="/data/runzezhang/result/SRIM_MC/MC_20231018_4745"
address3 ="/data/runzezhang/result/SRIM_MC/MC_20231018_3700"
weight1 = 10.8
weight2 = 51.25
weight3 = 9.11
address_list = [address1,address2,address3]
weight_list = [weight1, weight2, weight3]


def data_pick(address):
    with open(address, "rb") as fp:  # Unpickling
        b = pickle.load(fp)

    return b

def plot_chains():
    raw_data = []
    bin_n = 500
    raw_data_ev = []
    for address in address_list:
        raw_data.append(data_pick(address))
    # raw data is a data list
    #change data value from J into eV
    for data in raw_data:
        for data_ele in data:
            data_ele = data_ele/e
        raw_data_ev.append(data)

    for i in range(len(raw_data_ev)):
        print(raw_data_ev)
        plt.hist(raw_data_ev[i],bins= bin_n,label=address_list[i][-5:])
    plt.legend()
    plt.xlabel("energy/eV")
    plt.ylabel("N/bin")
    plt.show()

def plot_chain_sum():
    raw_data = []
    bin_n = 500
    raw_data_ev = []
    consistant_data = []
    consistant_weight =[]
    for address in address_list:
        raw_data.append(data_pick(address))
    # raw data is a data list
    # change data value from J into eV
    for i in range(len(raw_data)):
        for j in range(len(raw_data[i])):
            # add eV energy to a 1D list
            #together with its weights
            consistant_data.append(raw_data[i][j]/e)
            consistant_weight.append(weight_list[i])
    plt.hist(consistant_data,bins = bin_n, density = True, weight = consistant_weight)
    plt.xlabel("energy/eV")
    plt.ylabel("P")
    plt.show()

if __name__ =="__main__":
    plot_chains()
    # plot_chain_sum()







