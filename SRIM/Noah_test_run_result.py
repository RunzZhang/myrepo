import matplotlib.pyplot as plt
import numpy as np
import math
"""
2d nodes [array([681.49118484, 707.6934431 , 733.94152889, 767.15565892]), array([697.0097815 , 781.03986393, 826.74229317, 935.19737187]), array([640.75658371, 766.40864333, 828.52539069, 902.06823053]), array([670.72035467, 713.46684167, 759.46686377, 808.52608359]), array([692.5781226 , 732.16897264, 782.4619356 , 838.96359105]), array([711.46531127, 719.82470001, 736.32741965, 826.38246438]), array([664.63396039, 775.61054676, 856.72485051, 901.85587095]), array([677.69865661, 713.81705671, 850.04811099, 873.73340577]), array([648.90995569, 727.3742491 , 812.71946701, 861.48759327]), array([693.68564013, 737.37693997, 859.36490041, 896.05904582]), array([656.53446179, 727.79645351, 829.49684869, 876.48850571]), array([659.06108953, 743.298589  , 751.51592437, 772.25332301]), array([671.77010017, 744.79789671, 899.66409945, 905.89923311]), array([749.40798001, 767.05925497, 795.87184971, 804.60179529]), array([698.38452769, 734.53289178, 761.26357846, 831.80524509]), array([674.6293317 , 707.93695444, 811.65138024, 924.41762921]), array([733.95169748, 763.29606258, 795.3986018 , 811.57336024])]
sum node [11622.688739788275, 12563.499360211716, 13691.185043418734, 14538.468407799923]
mean [683.687572928722, 739.0293741301009, 805.363826083455, 855.2040239882308]
max [749.4079800131791, 781.0398639257841, 899.6640994495501, 935.1973718652577]
min [640.7565837077811, 707.6934431046741, 733.9415288859597, 767.1556589167683]

T = 760+-40

sum node [8843.502398195737, 9228.589350577864, 10227.960865878678, 12819.459994011153]
mean [520.2060234232786, 542.8581970928155, 601.6447568163928, 754.085882000656]
max [627.8005310696678, 630.9385321013965, 717.1263923101665, 1008.5716150284896]
min [416.9465269897052, 425.19543847251936, 537.3969397574829, 581.7533781645059]
T = 600 +-40

sum node [2183.504675142069, 6513.627576376943, 7539.910017746203, 8609.913301834209]
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
    x = range(0, 1000, 1)
    y = []
    for i in x:
        y.append(NucleationEfficiencyTrue(i,T=760, sigLow= 40, sigUp=40))
    plt.plot(x, y, color='green', label ='erf curve')
    data_y = [0, 0.3, 0.7, 1]
    data_x = [683.687572928722, 739.0293741301009, 805.363826083455, 855.2040239882308]
    plt.plot(data_x, data_y,color = 'blue', label="data")
    max_y = [0, 0.3, 0.7, 1]
    max_x = [627.8005310696678, 630.9385321013965, 717.1263923101665, 1008.5716150284896]
    plt.plot(max_x, max_y, color='red', label='max')
    min_y = [0, 0.3, 0.7, 1]
    min_x =[640.7565837077811, 707.6934431046741, 733.9415288859597, 767.1556589167683]
    plt.plot(min_x, min_y, color='red', label ='min')
    plt.legend()


    plt.show()


if __name__ =="__main__":
    plot_error()