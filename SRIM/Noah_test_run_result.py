import matplotlib.pyplot as plt
import numpy as np
import math
"""sum node [2183.504675142069, 6513.627576376943, 7539.910017746203, 8609.913301834209]
mean [272.93808439275864, 814.2034470471178, 942.4887522182754, 1076.2391627292761]
max [414.41989487614205, 1605.0256759410065, 1831.1480933101277, 1948.0116150885958]
min [119.60925946529423, 439.33535699254696, 455.47914544732305, 589.6382581080376]
T = 400 sigmahigh = sigmalow =30


sum node [449.6710770124641, 515.893047509419, 585.5963170420065, 716.2230694498137]
mean [64.23872528749487, 73.69900678705986, 83.65661672028664, 102.31758134997338]
max [70.56409963178898, 79.26294070015868, 90.04589742140719, 154.97070684283358]
min [55.30047822106481, 69.84327524360415, 79.88133999206757, 81.42715198180119]

T = 80 sigmahigh = low =10
guessfloor = 50
guessceil = 160
node [0, 0.3, 0.7, 1]"""


def NucleationEfficiencyTrue(r, T, sigLow, sigUp):
    if r < T:
        R = 1 / 2 * (1 + math.erf((r - T) / (sigLow * 2 ** (1 / 2))))
    else:
        # 1-R?
        R = 1 / 2 * (1 + math.erf((r - T) / (sigUp * 2 ** (1 / 2))))
    return R

def plot_error():
    x = range(0, 200, 1)
    y = []
    for i in x:
        y.append(NucleationEfficiencyTrue(i,T=80, sigLow= 10, sigUp=10))
    plt.plot(x, y, color='green', label ='erf curve')
    data_y = [0, 0.3, 0.7, 1]
    data_x = [64.23872528749487, 73.69900678705986, 83.65661672028664, 102.31758134997338]
    plt.plot(data_x, data_y,color = 'blue', label="data")
    max_y = [0, 0.3, 0.7, 1]
    max_x = [70.56409963178898, 79.26294070015868, 90.04589742140719, 154.97070684283358]
    plt.plot(max_x, max_y, color='red', label='max')
    min_y = [0, 0.3, 0.7, 1]
    min_x = [55.30047822106481, 69.84327524360415, 79.88133999206757, 81.42715198180119]
    plt.plot(min_x, min_y, color='red', label ='min')
    plt.legend()


    plt.show()


if __name__ =="__main__":
    plot_error()