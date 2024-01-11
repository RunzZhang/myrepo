import matplotlib.pyplot as plt
import numpy as np
import math
"""sum node [2183.504675142069, 6513.627576376943, 7539.910017746203, 8609.913301834209]
mean [272.93808439275864, 814.2034470471178, 942.4887522182754, 1076.2391627292761]
max [414.41989487614205, 1605.0256759410065, 1831.1480933101277, 1948.0116150885958]
min [119.60925946529423, 439.33535699254696, 455.47914544732305, 589.6382581080376]
T = 400 sigmahigh = sigmalow =30
node [0, 0.3, 0.7, 1]"""


def NucleationEfficiencyTrue(r, T, sigLow, sigUp):
    if r < T:
        R = 1 / 2 * (1 + math.erf((r - T) / (sigLow * 2 ** (1 / 2))))
    else:
        # 1-R?
        R = 1 / 2 * (1 + math.erf((r - T) / (sigUp * 2 ** (1 / 2))))
    return R

def plot_error():
    x = range(0, 2000, 10)
    y = []
    for i in x:
        y.append(NucleationEfficiencyTrue(i,T=400, sigLow= 30, sigUp=30))
    plt.plot(x, y, color='green')
    data_y = [0, 0.3, 0.7, 1]
    data_x = [272.93808439275864, 814.2034470471178, 942.4887522182754, 1076.2391627292761]
    plt.plot(data_x, data_y,color = 'blue')
    max_y = [0, 0.3, 0.7, 1]
    max_x = [414.41989487614205, 1605.0256759410065, 1831.1480933101277, 1948.0116150885958]
    plt.plot(max_x, max_y, color='red')
    min_y = [0, 0.3, 0.7, 1]
    min_x = [119.60925946529423, 439.33535699254696, 455.47914544732305, 589.6382581080376]
    plt.plot(min_x, min_y, color='red')


    plt.show()


if __name__ =="__main__":
    plot_error()