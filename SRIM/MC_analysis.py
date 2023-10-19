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

test_address ="/data/runzezhang/result/SRIM_MC/MC_20231018_4745_2"


def data_pick(address):
    with open(address, "rb") as fp:  # Unpickling
        b = pickle.load(fp)

    return b

def plot_test():
    raw_data = []
    bin_n = 500
    raw_data_ev = []
    start = 0
    end = 600

    raw_data.append(data_pick(test_address))
    # print("type", type(raw_data[0][1]))
    # raw data is a data list
    #change data value from J into eV
    data_buffer = []
    for data_ele in raw_data:
        # change sympy float into float, otherwise the data cannot been plot by matplotlib
        data_buffer.append(float(round(data_ele/e,3)))
        # print(round(data_ele/e,3))
    print("finish one run")


    plt.hist(data_buffer,bins= bin_n,range=(start, end),label="test_run")
    plt.legend()
    plt.xlabel("energy/eV")
    plt.ylabel("N/bin")
    plt.show()
def plot_chains():
    raw_data = []
    bin_n = 500
    raw_data_ev = []
    start = 0
    end = 600

    for address in address_list:
        raw_data.append(data_pick(address))
    # print("type", type(raw_data[0][1]))
    # raw data is a data list
    #change data value from J into eV
    for data in raw_data:
        data_buffer = []
        for data_ele in data:
            # change sympy float into float, otherwise the data cannot been plot by matplotlib
            data_buffer.append(float(round(data_ele/e,3)))
            # print(round(data_ele/e,3))
        raw_data_ev.append(data_buffer)
        print("finish one run")
    # print(raw_data_ev[0])
    # print(address_list[0][-4:])
    # print("min and max", min(raw_data_ev[0]), max(raw_data_ev[0]))
    # plt.hist(np.array(raw_data_ev[0][:500]), bins=50, range=(min(raw_data_ev[0]),max(raw_data_ev[0])),label=address_list[0][-4:])
    figure, axis = plt.subplots(3)

    for i in range(len(raw_data_ev)):
        print(len(raw_data_ev))
        axis[i].hist(raw_data_ev[i],bins= bin_n,range=(start, end),label=address_list[i][-4:])
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
            consistant_data.append(round(float(raw_data[i][j]/e),3))
            consistant_weight.append(weight_list[i])
    plt.hist(consistant_data,bins = bin_n, density = True, weights = consistant_weight)
    plt.xlabel("energy/eV")
    plt.ylabel("P")
    plt.show()

if __name__ =="__main__":

    plot_test()
    # plot_chains()
    # plot_chain_sum()







