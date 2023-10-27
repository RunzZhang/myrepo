import scipy,random, pickle
import sympy
from sympy import cos, sin, nsolve, Symbol
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter

from scipy.optimize import fsolve
e = 1.602*10**(-19)

address_Ar40_5582="/data/runzezhang/result/SRIM_MC/MC_argon40_20231024_5582"
address_Ar40_4745="/data/runzezhang/result/SRIM_MC/MC_argon40_20231024_4745"
address_Ar40_3700="/data/runzezhang/result/SRIM_MC/MC_argon40_20231024_3700"
address_Ar40_2771="/data/runzezhang/result/SRIM_MC/MC_argon40_20231024_2771"

address_Ar36_8790="/data/runzezhang/result/SRIM_MC/MC_argon36_20231024_8790"
address_Ar36_6299="/data/runzezhang/result/SRIM_MC/MC_argon36_20231024_6299"
address_Ar36_5272="/data/runzezhang/result/SRIM_MC/MC_argon36_20231024_5272"
address_Ar36_3700="/data/runzezhang/result/SRIM_MC/MC_argon36_20231024_3700"
# 93 is norm factor of argon 40
ar_40_percent =0.974
ar_36_percent = 0.025
weight_Ar40_5582 = ar_40_percent * 10.8/93
weight_Ar40_4745 = ar_40_percent * 51.2/93
weight_Ar40_3700 = ar_40_percent * 9.11/93
weight_Ar40_2771 = ar_40_percent * 8/93

weight_Ar36_8790 = ar_36_percent * 10.9/100
weight_Ar36_6299 = ar_36_percent * 37.5/100
weight_Ar36_5272 = ar_36_percent * 25/100
weight_Ar36_3700 = ar_36_percent * 8/100


# address_list = [address_Ar40_5582,address_Ar40_4745,address_Ar40_3700,address_Ar40_2771,address_Ar36_8790,address_Ar36_6299,address_Ar36_5272,address_Ar36_3700]
# weight_list = [weight_Ar40_5582,weight_Ar40_4745,weight_Ar40_3700,weight_Ar40_2771,weight_Ar36_8790,weight_Ar36_6299,weight_Ar36_5272,weight_Ar36_3700]

address_list = [address_Ar40_3700,address_Ar40_2771,address_Ar36_8790,address_Ar36_5272,address_Ar36_3700]
weight_list = [weight_Ar40_3700,weight_Ar40_2771,weight_Ar36_8790,weight_Ar36_5272,weight_Ar36_3700]
address_list_40 =[address_Ar40_5582,address_Ar40_4745,address_Ar40_3700,address_Ar40_2771]
weight_list_40 = [weight_Ar40_5582,weight_Ar40_4745,weight_Ar40_3700,weight_Ar40_2771]
address_list_36 = [address_Ar36_8790,address_Ar36_6299,address_Ar36_5272,address_Ar36_3700]
weight_list_36 = [weight_Ar36_8790,weight_Ar36_6299,weight_Ar36_5272,weight_Ar36_3700]
color_list = ["red", "orange", "blue", "green"]
PI = np.pi

test_address ="/data/runzezhang/result/SRIM_MC/MC_argon36_20231024_8790"


def data_pick(address):
    with open(address, "rb") as fp:  # Unpickling
        b = pickle.load(fp)

    return b

def plot_test():
    raw_data = []
    bin_n = 500
    raw_data_ev = []
    start = 0
    end = 1200
    test_N = 100

    raw_data.append(data_pick(test_address))
    # print("type", type(raw_data[0][1]))
    # raw data is a data list
    #change data value from J into eV
    data_buffer = []
    for data_ele in raw_data[0][:test_N]:
        # change sympy float into float, otherwise the data cannot been plot by matplotlib
        data_buffer.append(float(round(data_ele/e,3)))
        print(data_ele)
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
    end = 1200

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
    figure, axis = plt.subplots(len(raw_data_ev))

    for i in range(len(raw_data_ev)):
        print(len(raw_data_ev))
        axis[i].hist(raw_data_ev[i],bins= bin_n,range=(start, end),label=address_list[i][-4:])
    plt.legend()
    plt.xlabel("energy/eV")
    plt.ylabel("N/bin")
    plt.show()
def plot_chains_ingroup():
    raw_data = []
    bin_n = 500
    raw_data_ev = []
    start = 0
    end = 1200
    data_result = []
    energy_x = []
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
    figure, axis = plt.subplots(len(raw_data_ev))

    for i in range(len(raw_data_ev)):
        print(len(raw_data_ev))
        data_result.append(axis[i].hist(raw_data_ev[i],bins= bin_n,range=(start, end),label=address_list[i][-4:]))
    plt.clf()

    figure_sum, axis_sum = plt.subplots(2)
    for i in range(len(address_list_40)):
        axis_sum[0].plot(data_result[i][1][:-1],data_result[i][0],label = address_list_40[i][-4:], color = color_list[i])
    for i in range(len(address_list_36)):
        axis_sum[1].plot(data_result[i+4][1][:-1],data_result[i+4][0],label = address_list_36[i][-4:], color = color_list[i])

    plt.yticks(fontsize=18)
    plt.xticks(fontsize=18)
    plt.legend(prop={'size': 10})
    plt.xlabel("energy/eV", fontsize = 18)
    plt.ylabel("N/bin", fontsize = 18)
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
    hist_result = plt.hist(consistant_data,bins = bin_n, density = True, weights = consistant_weight)
    plt.clf()
    # replot the hist as lineplot
    x_bins= []
    for i in range(len(hist_result[1])-1):
        x_bins.append((hist_result[1][i]+hist_result[1][i+1])/2)
    plt.plot(x_bins,hist_result[0])

    plt.xlabel("energy/eV")
    plt.ylabel("P")
    plt.yscale("log")
    plt.ylim([10**(-6),0.1])
    plt.xlim([0,1200])
    plt.show()

def plot_chain_separate():
    raw_data1 = []
    bin_n = 500
    raw_data_ev1 = []
    consistant_data1 = []
    consistant_weight1 =[]
    raw_data2 = []
    raw_data_ev2 = []
    consistant_data2 = []
    consistant_weight2 = []
    for address in address_list_40:
        raw_data1.append(data_pick(address))
    print(raw_data1)
    # raw data is a data list
    # change data value from J into eV
    for i in range(len(raw_data1)):
        for j in range(len(raw_data1[i])):
            # add eV energy to a 1D list
            #together with its weights
            consistant_data1.append(round(float(raw_data1[i][j]/e),3))
            consistant_weight1.append(weight_list_40[i])
    for address in address_list_36:
        raw_data2.append(data_pick(address))
    # raw data is a data list
    # change data value from J into eV
    for i in range(len(raw_data2)):
        for j in range(len(raw_data2[i])):
            # add eV energy to a 1D list
            #together with its weights
            consistant_data2.append(round(float(raw_data2[i][j]/e),3))
            consistant_weight2.append(weight_list_36[i])
    hist_result1 = plt.hist(consistant_data1,bins = bin_n, density = False, weights = consistant_weight1)
    hist_result2 = plt.hist(consistant_data2, bins=bin_n, density=False, weights=consistant_weight2)
    plt.clf()
    # replot the hist as lineplot
    x_bins1= []
    for i in range(len(hist_result1[1])-1):
        x_bins1.append((hist_result1[1][i]+hist_result1[1][i+1])/2)
    # make plot look better, the last bin's value should be a zero
    hist_result1_modified = np.append(hist_result1[0], 0)
    x_bins1.append(hist_result1[1][-1])
    x_bins2 = []
    for i in range(len(hist_result2[1]) - 1):
        x_bins2.append((hist_result2[1][i] + hist_result2[1][i + 1]) / 2)
    plt.plot(x_bins1,hist_result1_modified,color = "blue", label = "Argon40")
    plt.plot(x_bins2, hist_result2[0], color = "orange",label = "Argon36")

    plt.xlabel("energy/eV" , fontsize=18)
    plt.ylabel("evnets", fontsize=18)
    plt.yscale("log")
    plt.legend(prop={'size': 18})
    plt.yticks(fontsize = 18)
    plt.xticks(fontsize=18)
    # plt.ylim([10**(-6),0.1])
    plt.xlim([0,1200])
    plt.show()

def two_body_E_spectrum_47_func(x):
    a = 300
    b = 50
    deno = np.sqrt(1-(1/(2*a**2)+1/(2*b**2)-x/(a*b))**2)

    f = 1/deno
    return f

def Energy_vs_theta_47(theta):
    a = 300
    b = 50
    return a**2+b**2-2*a*b*np.cos(theta)

def plot_two_body_E_spectrum_47():
    angle = []
    E_list = []
    N = 100000
    for i in range(N):
        theta = random.uniform(0, 2 * PI)
        E = Energy_vs_theta_47(theta)
        angle.append(theta)
        E_list.append(E)

    plt.hist(E_list, bins = 500)
    plt.show()
if __name__ =="__main__":

    # plot_two_body_E_spectrum_47()
    # plot_test()
    # plot_chains()
    # plot_chain_sum()
    plot_chain_separate()






