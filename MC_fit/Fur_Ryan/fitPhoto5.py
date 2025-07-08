import matplotlib.pyplot as plt
import sys
import numpy as np
import math
import scipy.optimize
import random

################################################
#flist = ["../Co60/JAEA.txt","../Sb124/JAEA.txt","../Th228/JAEA.txt"]
#flist = ["../JAEA_Co60.txt","../JAEA_Sb124.txt","../JAEA_Th228.txt"]
#flist = ["../JAEA_Cal/Co60_JAEA.txt","../JAEA_Cal/Sb124_JAEA.txt","../Th228/JAEA.txt"] 
#flist = ["../Compton_Th228/Compt.txt"]
flist = ["../Compton_Co60/Compt.txt","../Compton_Sb124/Compt.txt","../Compton_Th228/Compt.txt"]
#flist = ["../Eu152/JAEA.txt","../Bi207/JAEA.txt","../Y88/JAEA.txt"] 
################################################
#pnlist = ["../Bi"]
pnlist = ["../Bi","../Sb"]
################################################
#alist = [10]
alist = [50,20,10]
################################################
aplist = [8]  # this is not yet implimented

try:
    fileprefix = sys.argv[3]
except:
    fileprefix = "Test_Dump/test2"
""" 
def NucleationEfficiency(r,T,sigma):
    #A=1/2
    #denom=sigma*2**(1/2)
    #num=r-T
    #print ("NucleationEfficiency: ", r,T,sigma)
    R=1/2*(1+math.erf((r-T)/(sigma*2**(1/2))))
    #print (R)
    return R
"""


def NucleationEfficiencyTrue(r, T, sigLow, sigUp):
    if r < T:
        R = 1 / 2 * (1 + math.erf((r - T) / (sigLow * 2 ** (1 / 2))))
    else:
        R = 1 / 2 * (1 + math.erf((r - T) / (sigUp * 2 ** (1 / 2))))
    return R


def NucleationEfficiency(r, energies, efficiencies="auto"):
    n = len(energies)
    if efficiencies == "auto":
        print("Automation")
        efficiencies = np.zeros(n)
        for i in range(n):
            efficiencies[i] = i / (n - 1)
    # print (efficiencies)
    for i in range(n):
        if r < energies[i]:
            if i == 0:
                return 0
            else:
                return (efficiencies[i] - efficiencies[i - 1]) * (
                    r - energies[i - 1]
                ) / (energies[i] - energies[i - 1]) + efficiencies[i - 1]
    return 1


def VectorNucleationEfficiency(r, energies, efficiencies):
    n = len(energies)
    en0 = energies[0]
    en20 = energies[1]
    en50 = energies[2]
    en80 = energies[3]
    en100 = energies[4]
    return np.piecewise(
        r,
        [
            r <= en0,
            (r > en0) & (r <= en20),
            (r > en20) & (r <= en50),
            (r > en50) & (r <= en80),
            (r > en80) & (r < en100),
            r >= en100,
        ],
        [
            lambda r: 0,
            lambda r: (0.2) / (en20 - en0) * (r - en0),
            lambda r: (0.3) / (en50 - en20) * (r - en20) + 0.2,
            lambda r: (0.3) / (en80 - en50) * (r - en50) + 0.5,
            lambda r: (0.2) / (en100 - en80) * (r - en80) + 0.8,
            lambda r: 1,
        ],
    )


def rateFinder(recoil, energies, efficiencies, t, weight):
    # print ("rateFinder: ",data,energy,T,sigma,t,weight)
    W = np.ones(len(recoil))
    W = VectorNucleationEfficiency(recoil, energies, efficiencies)
    W *= weight
    W /= t
    #for i in range(len(recoil)):
        # print (np.shape(data))
        # print (data[i])
        # print (np.shape(data[i]))
        #W[i] = (
            #NucleationEfficiency(recoil[i], energies, efficiencies=efficiencies)
            #* weight[i]
            #/ t
        #)
    return sum(W)



def specrafy(recoils, weights, binsize=1):
    minr = np.floor(np.min(recoils))
    maxr = np.ceil(np.max(recoils))
    # if maxr>600:
    # maxr=600
    spread = int(maxr - minr)

    recoilS = np.arange(minr, maxr, binsize)
    lenS = len(recoilS)
    weightS = np.zeros(lenS)

    for i in range(len(recoils)):
        for j in range(lenS):
            if (
                recoils[i] > recoilS[j] - binsize / 2
                and recoils[i] <= recoilS[j] + binsize / 2
            ):
                weightS[j] += weights[i]
    # print (recoilS,weightS)
    return recoilS, weightS

def neutron_rates(spec1, T, sigLow, sigUp, t):
    recoil = spec1[:, 0]
    weight = spec1[:, 1]
    rate = rateFinderTrue(recoil, T, sigLow, sigUp, t, weight)
    return rate


def rateFinderTrue(recoil, T, sigLow, sigUp, t, weight):
    # print ("rateFinder: ",data,energy,T,sigma,t,weight)
    weight = weight
    length=len(recoil)
    W = np.ones(len(recoil))
    for i in range(length):
        # print (np.shape(data))
        # print (data[i])
        # print (np.shape(data[i]))
        W[i] = NucleationEfficiencyTrue(recoil[i], T, sigLow, sigUp) * weight[i] / t
    return sum(W)


def test(
    recoil,
    rate,
    energies,
    efficiencies,
    r_nuis,
    t,
    weight,
    mode=0,
    background=500,
    time=100,
):
    r_nuis=logErr(r_nuis)
    mode=logErr(mode)
    rateT = rateFinder(recoil, energies, efficiencies, t, weight) * (r_nuis) * (mode) + (background / time)
    # print (rateT)
    # print (rate)
    dif = rateT-rate
    chi = ((rateT - rate) ** 2 / (rate)) * time
    return chi,dif,rateT,rate


def testLoud(
    recoil, rate, energies, efficiencies, r_nuis, t, weight, background=500, time=100
):
    rateT = rateFinder(recoil, energies, efficiencies, t, weight) * (r_nuis) + (background / time)
    print("Model Rate (No Background): ", rateT)
    rateT += background / time
    print("Model Rate: ", rateT)
    chi = (rateT - rate) ** 2 / rate
    print("Data Rate: ", rate)
    return chi

def logErr(err):
    #if err==1:
        #print ("Pre Log Error: ",err)
    error=np.exp(err)
    #print ("Post Log Error: ",error)
    return error

def rateJitter(
    rate, count, commonMode, time=100, sourceErr=0.03, background=500, s_err="False"
):
    if s_err == "False":
        s_err = np.squeeze(np.random.normal(0, sourceErr, 1))
        #s_err = np.random.normal(0, sourceErr, 1)
    # s_err=1+random.uniform(-1,1)*sourceErr
    s_err=logErr(s_err)
    commonMode=logErr(commonMode)
    count = count * s_err * commonMode
    count = np.random.poisson(count)
    backgroundexcess = np.random.poisson(background)
    count += backgroundexcess
    rate = count / time
    b_err = backgroundexcess / background - 1
    print("Background: ", backgroundexcess, background)
    print("Source Error: ", s_err)
    return rate, count, backgroundexcess, s_err
    # (T,sigLow,sigUp,energiesb,efficiencies0,start=T-4*sigLow,end=T-4*sigUp)


def fittest(T, sigLow, sigUp, energiesM, efficienciesM, M=2000, start=0, end=200):
    rt = 0
    rg = 0
    dif = 0
    tot = 0
    calib = 0
    fifty = 0
    fiftydif = 0.5
    stepwidth=end/M
    # print (T,sigLow,sigUp)
    # print (energiesM,efficienciesM)
    for i in range(M):
        #stepwidth = end / M
        r = i * stepwidth
        RT = NucleationEfficiencyTrue(r, T, sigLow, sigUp)
        # print ("Real Efficiency: ",RT)
        RG = NucleationEfficiency(r, energiesM, efficiencies=efficienciesM)
        # print ("Model Efficiency: ",RG)
        # print ("Difference: ",RT-RG)
        dif += abs(RT - RG) * stepwidth
        tot += RT * stepwidth
        if fiftydif >= abs(0.5 - RG):
            fiftydif = 0.5 - RG
            fifty = r
        RT *= stepwidth
        RG *= stepwidth
        rt += RT
        rg += RG
        # print ("Absolute Difference: ",dif)
    return dif, fifty, rt, rg, tot


t = 0.03754 / 100  # hours to run the whole thing,

T = 105
sigma = 15


def analyze(file, T, sigLow, sigUp, N=2 * 10**9, Activity=100, time=100):
    Data = np.loadtxt(file)
    SourceRate = (3.7 * 10**6) * Activity / 100
    Energy = Data[:, 2]
    Recoils = Data[:, 4]
    Weights = Data[:, 5]

    t = N / SourceRate  # live time in seconds
    t /= 3600  # live time in hours

    Rate = rateFinderTrue(Recoils, T, sigLow, sigUp, t, Weights)

    Count = Rate * time

    return Recoils, Weights, Rate, Count, t, time

def analyze_Compton(file, T, sigLow, sigUp, N=5.0 * 10**7, Activity=100, time=100):
    Data = np.loadtxt(file)
    SourceRate = (3.7 * 10**6) * Activity / 100
    Energy = Data[:, 0]
    Recoils = Data[:, 1]
    Weights = Data[:, 2]

    t = N / SourceRate  # live time in seconds                                                                                                                                                              
    t /= 3600  # live time in hours                                                                                                                                                                         

    Rate = rateFinderTrue(Recoils, T, sigLow, sigUp, t, Weights)

    Count = Rate * time

    return Recoils, Weights, Rate, Count, t, time

def stepper(energies, nuisance, mode, step, n_step, m_step, min=40, max=2000):
    energiesMem = np.zeros(len(energies))
    nuisanceMem = np.zeros(len(nuisance))
    modeMem = np.zeros(len(mode))
    energiesMem[:] = energies[:]
    nuisanceMem[:] = nuisance[:]
    modeMem[:] = mode[:]
    chosen = np.random.randint(0, len(energies))
    for i in range(len(energiesMem)):
        energiesMem[i] = np.squeeze(np.random.normal(0, step, 1)) + energiesMem[i]
        #energiesMem[i] = np.random.normal(0, step, 1) + energiesMem[i]
        if i == 0:
            if energiesMem[0] < min:
                energiesMem[0] = min + 5
            if energiesMem[1] < min:
                energiesMem[1] = min + 5
            if energiesMem[0] > energiesMem[1]:
                low = energiesMem[1]
                energiesMem[0] = energiesMem[1]
                energiesMem[1] = low
        elif i == len(energiesMem) - 1:
            if energiesMem[i] < energiesMem[i - 1]:
                energiesMem[i] = energiesMem[i - 1]
            if energiesMem[i] < min:
                energiesMem[i] = min + 5
            if energiesMem[i] > max:
                energiesMem[i] = max - 10
            if energiesMem[i] < energiesMem[i - 1]:
                low = energies[i]
                energiesMem[i] = energiesMem[i - 1]
                energiesMem[i - 1] = low
        else:
            if energiesMem[i] > max:
                energiesMem[i] = max - 10
            if energiesMem[i] < min:
                energiesMem[i] = min + 5
            if energiesMem[i] < energiesMem[i - 1]:
                low = energiesMem[i]
                energiesMem[i] = energiesMem[i - 1]
    chosen = np.random.randint(0, len(nuisance))
    for i in range(len(nuisance)):
        # if i==chosen:
        #nuisanceMem[i]=logErr(nuisanceMem[i])
        if True:
            nuisanceMem[i] += np.squeeze(np.random.normal(0, n_step, 1))
            #nuisanceMem[i] += np.random.normal(0, n_step, 1)
        if nuisanceMem[i] > 0.2:
            nuisanceMem[i] = 0.2
        if nuisanceMem[i] < -0.2:
            nuisanceMem[i] = -0.2
    for i in range(len(mode)):
        modeMem[i] += np.squeeze(np.random.normal(0, m_step, 1))
        #modeMem[i] += np.random.normal(0, m_step, 1)
        if modeMem[i] > 1:
            modeMem[i] = 1
        if modeMem[i] < -1:
            modeMem[i] = -1
    return energiesMem[:], nuisanceMem[:], modeMem[:]


def photoread(pnc_filepath):
    r_list = []
    m_list = []
    data = np.loadtxt(pnc_filepath, skiprows=1, dtype="str")
    m = 0
    for i in range(len(data)):
        if i == 0:
            m = 0
        elif data[i, 0] == data[i - 1, 0]:
            m = m
        else:
            m += 1
        r = float(data[i, 1])
        u = data[i, 2]
        if u == "eV":
            r = r
        elif u == "keV":
            r = 1000 * r
        elif u == "MeV":
            r = 1000000 * r
        r = int(np.round(r, 0))
        r_list += [r]
        m_list += [m]
    pnc_array = np.empty([len(r_list), 2])
    pnc_array[:, 0] = m_list
    pnc_array[:, 1] = r_list
    print("crosscheck: ", len(data), m)
    return pnc_array

def photo_speed_read(pnc_filepath):
    r_list = []
    m_list = []
    data = np.loadtxt(pnc_filepath)
    m = 0
    i=0
    for line in data:
        if i == 0:
            m = 0
        elif line[0] == line_old[0]:
            m = m
        else:
            m += 1
        r = line[1]
        r_list.append(r)
        m_list.append(m)
        line_old=line
        i+=1
    print (r_list[0])
    print (m_list[0])
    print (len(r_list))
    pnc_array = np.empty([len(r_list), 2])
    print (np.shape(pnc_array))
    pnc_array[:, 0] = m_list
    pnc_array[:, 1] = r_list
    print("crosscheck: ", len(data), m)
    return pnc_array

def photo_spectra_read(pnc_prefix):
    highest=np.loadtxt(pnc_prefix+"highestrecoils.txt")
    all=np.loadtxt(pnc_prefix+"allrecoils.txt")
    r_list = []
    m_list = []
    data = np.loadtxt(pnc_filepath)
    m = 0
    i=0
    for line in data:
        if i == 0:
            m = 0
        elif line[0] == line_old[0]:
            m = m
        else:
            m += 1
        r = line[1]
        r_list.append(r)
        m_list.append(m)
        line_old=line
        i+=1
    print (r_list[0])
    print (m_list[0])
    print (len(r_list))
    pnc_array = np.empty([len(r_list), 2])
    print (np.shape(pnc_array))
    pnc_array[:, 0] = m_list
    pnc_array[:, 1] = r_list
    print("crosscheck: ", len(data), m)
    return pnc_array

def photoEval(pn, pnlist, photoArray, photoArrayTrue):
    for p in range(pn):
        photoArray[p, :] = phototestmean(
            pnlist[p] + "_fast.txt", 500, energies0, efficiencies0
        )
        for mb in range(len(photoArray)):
            if photoArrayTrue[p, mb] > 0:
                chipn[p] += (
                    photoArray[p, mb] - photoArrayTrue[p, mb]
                ) ** 2 / photoArrayTrue[p, mb]


def photoEval2(
    pn,
    photoneutrondata,
    livetime,
    energies,
    efficiencies,
    photoArrayTrue,
    photoBackMean,
    loud=False,
    pnuisance=[1, 1],
    mode=1,
):
    chi = np.zeros(pn)
    difTotal= np.zeros(pn)
    difCount= np.zeros(pn)
    totalModel= np.zeros(pn)
    totalTrue= np.zeros(pn)
    countModel = np.zeros(pn)
    countTrue = np.zeros(pn)
    for p in range(pn):
        #trueCount = photoArrayTrue[p, 1]
        trueMean = 0
        trueTot = 0
        for i in range(20):
            if i != 1:
                trueCount += photoArrayTrue[p, i] * (i)
                trueTot += photoArrayTrue[p, i] * (i)
        # trueMean=np.mean(photoArrayTrue[p,:])
        # print ("photoEval2")
        # print("photoneutrondata",photoneutrondata)
        # print ("photoneutrondata[p]", photoneutrondata[p])
        total, count, mean, czero = phototest3(
            photoneutrondata[p],
            livetime,
            energies,
            efficiencies,
            photoBackMean,
            sourceStrength=(logErr(pnuisance[p]) * logErr(mode)))
        # modelCount=count-czero
        totalTrue[p]=trueTot
        countTrue[p]=trueCount
        totalModel[p]=total
        countModel[p]=count
        difTotal[p]=total-trueTot
        difCount[p]=count-trueCount
        chiTotal = ((total - trueTot) ** 2 / trueTot) * 100
        chiCount = ((count - trueCount) ** 2 / trueCount) * 100
        chi[p] = chiTotal + chiCount
        if loud:
            print("True Mults: ", photoArrayTrue[p, :])
            print("True Count: ", trueCount)
            print("True Zeros: ", photoArrayTrue[p, 0])
            print("Model Zeros: ", czero)
            print(">0 True Count: ", sum(photoArrayTrue[p, 1:]))
            print("Model Count: ", count)
            print("Total/Mean True: ", trueTot / trueCount)
            print("Total/Mean Model: ", total / count)
            print("True Mean: ", trueMean)
            print("Model Mean: ", mean)
            print("True Sum: ", trueTot)
            print("Model Sum: ", total)
            print("Chi total: ", chiTotal)
            print("Chi count: ", chiCount)
            print(chi[p])
    return chi,difTotal,difCount,totalTrue,totalModel,countTrue,countModel
#chipn,dtotali[n:],dcounti[n:],totalTrue[n:],totalModeli[n:],countTrue[n:],countModeli[n:]

def photoEval3(
    pn,
    photoneutrondata,
    livetime,
    energies,
    efficiencies,
    photoArrayTrue,
    meanBackBubble,
    meanBackEvent,
    loud=False,
    pnuisance=[1, 1],
    mode=1,
):
    chi = np.zeros(pn)
    difTotal= np.zeros(pn)
    difCount= np.zeros(pn)
    totalModel= np.zeros(pn)
    totalTrue= np.zeros(pn)
    countModel = np.zeros(pn)
    countTrue = np.zeros(pn)
    for p in range(pn):
        trueCount = 0
        trueMean = 0
        trueTot = 0
        for i in range(20):
            if i>=1:
                trueCount += photoArrayTrue[p, i]
                trueTot += photoArrayTrue[p, i] * (i)
        # trueMean=np.mean(photoArrayTrue[p,:])    
        # print ("photoEval2")
        # print("photoneutrondata",photoneutrondata)
        # print ("photoneutrondata[p]", photoneutrondata[p]) 
        total,count = phototest3(
            photoneutrondata[p],
            energies,
            efficiencies,
            sourceStrength=(logErr(pnuisance[p]) * logErr(mode)))
        # modelCount=count-czero
        total+=meanBackBubble
        count+=meanBackEvent
        totalTrue[p]=trueTot
        countTrue[p]=trueCount
        totalModel[p]=total
        countModel[p]=count
        meanModel=total/count
        meanTrue=trueTot/trueCount
        difTotal[p]=meanModel - meanTrue
        difCount[p]=count - trueCount
        chiTotal = ((meanModel - meanTrue) ** 2 / meanTrue) * 100 * trueCount
        chiCount = ((count - trueCount) ** 2 / trueCount) * 100
        chi[p] = chiTotal + chiCount
        if loud:
            print("True Mults: ", photoArrayTrue[p, :])
            print("Mean Model: ",meanModel)
            print("Mean True: ",meanTrue)
            print("True Count: ", trueCount)
            print("Model Count: ", count)
            print(">0 True Count: ", sum(photoArrayTrue[p, 1:]))
            print("True Sum: ", trueTot)
            print("Model Sum: ", total)
            #print("True Zeros: ", photoArrayTrue[p, 0])
            #print("Model Zeros: ", czero)
            print("Total/Mean True: ", trueTot / trueCount)
            print("Total/Mean Model: ", total / count)
            #print("True Mean: ", trueMean)
            #print("Model Mean: ", mean)
            print("Chi total: ", chiTotal)
            print("Chi count: ", chiCount)
            print(chi[p])
    return chi,difTotal,difCount,totalTrue,totalModel,countTrue,countModel


def photoJitter(photoSourceTrue, photoBackMean, sourceError,pnCommonMode, time=100):
    # poisson uncertainty in background + events
    # source strength uncertainty
    photoBack = np.zeros(len(photoBackMean))
    pnCommonMode=logErr(pnCommonMode)
    sourceError=logErr(sourceError)
    photoSourceTrue *= sourceError
    photoSourceTrue *= pnCommonMode
    photoSourceTrue *= time
    photoBack = time * photoBackMean
    for i in range(len(photoSourceTrue)):
        photoSourceTrue[i] = np.random.poisson(photoSourceTrue[i])
        photoBack[i] = np.random.poisson(photoBack[i]) 
        #photoSourceTrue[i] = photoSourceTrue[i]
        #photoBack[i] = photoBack[i]
    photoSourceTrue = (photoSourceTrue + photoBack) / time
    return photoSourceTrue


def phototestTrue(pnc_filepath, livetime, T, sigLow, sigUp):
    # data1=np.loadtxt("../photoNC/Informacion_Sb124_high1.txt",skiprows=1,dtype="str")
    data = np.loadtxt(pnc_filepath, dtype="int")
    # data=data[0:2000]
    # print("Neutron output shape: ",np.shape(data))
    m = 0
    M = np.zeros(20)
    bsum = 0
    i=0
    # for i in range(500):
    for line in data:
        if i == 0:
            m = 0
        elif line[0] == line_old[0]:
            # print ("Same: ", i, data[i,0])
            m = m
        else:
            if m == 0:
                m = 0
                M[0] += 1
            else:
                M[m] += 1
                # print ("Recorded: ",m,i, M[m])
            m = 0
        r = line[1]
        t = 1
        weight = 1
        popper = np.random.uniform(0, 1)
        blower = NucleationEfficiencyTrue(r, T, sigLow, sigUp)
        bsum += blower
        # print ("blower: ",blower)
        # print ("popper: ",popper, i, r)
        if popper < blower:
            m += 1
            # print ("Tracker: ",i,m)
        if i == len(data) - 1:
            if m == 0:
                m = 0
                M[0] += 1
            else:
                M[m] += 1
        line_old=line
        i+=1
                # print ("Recorded: ",m,i, M[m])
    # print ("M :",M/np.sum(M))
    average = 0
    for i in range(len(M)):
        average += M[i] * i / np.sum(M)
    print("True Sum: ", bsum)
    print("Mult Sum: ", np.sum(M) / average)
    print("True Mean: ", bsum / average)
    print("Mult Sum: ", average)
    # print ("Average m: ",average)
    return M / livetime


def phototest(pnc_filepath, livetime, energies, efficiencies):
    # data1=np.loadtxt("../photoNC/Informacion_Sb124_high1.txt",skiprows=1,dtype="str")
    data = np.loadtxt(pnc_filepath, skiprows=1, dtype="str")
    # data=data[0:2000]
    # print("Neutron output shape: ",np.shape(data))
    m = 0
    M = np.zeros(20)
    # for i in range(500):
    for i in range(len(data)):
        if i == 0:
            m = 0
        elif data[i, 0] == data[i - 1, 0]:
            # print ("Same: ", i, data[i,0])
            m = m
        else:
            if m == 0:
                m = 0
            else:
                M[m] += 1
                # print ("Recorded: ",m,i, M[m])
            m = 0
        r = float(data[i, 1])
        u = data[i, 2]
        if u == "eV":
            r = r
        elif u == "keV":
            r = 1000 * r
        elif u == "MeV":
            r = 1000000 * r
        else:
            # print ("Help!",i+1,u)
            fiftyfive = 55
        t = 1
        weight = 1
        popper = np.random.uniform(0, 1)
        blower = NucleationEfficiency(r, energies, efficiencies)
        # print ("blower: ",blower)
        # print ("popper: ",popper, i, r)
        if popper < blower:
            m += 1
            # print ("Tracker: ",i,m)
        if i == len(data) - 1:
            if m == 0:
                m = 0
            else:
                M[m] += 1
                # print ("Recorded: ",m,i, M[m])
    # print ("M :",M/np.sum(M))
    average = 0
    for i in range(len(M)):
        average += M[i] * i / np.sum(M)
    # print ("Average m: ",average)
    return M[1:] / livetime

def phototest3(data, energies, efficiencies, sourceStrength=1):
    m = 0
    M = []
    PZ = 0
    lineO = [0, 0]
    pz = 0
    m = 0
    p1=0
    # print (data)                                                                                                                                                                                          
    # print (data[:,0])                                                                                                                                                                                     
    # print (data[:,1])                                                                                                                                                                                     
    # datOrg=data[:,1]
    #print (data)
    #print (np.shape(data))
    dataMem = np.zeros(np.shape(data))
    dataMem[:, 1] = data[:, 1]
    dataMem[:, 2] = data[:, 2]

    #T=80                                                                                                                                                                                                   
    #sigLow=10                                                                                                                                                                                              
    #sigUp=sigLow                                                                                                                                                                                           
    dataMem[:, 0] = VectorNucleationEfficiency(data[:, 0], energies, efficiencies)
    bubble_spectra=dataMem[:, 0] * dataMem[:, 1]
    event_spectra=dataMem[:, 0] * dataMem[:, 2]

    bubble_rate=np.sum(bubble_spectra)*sourceStrength
    event_rate=np.sum(event_spectra)*sourceStrength
    return bubble_rate,event_rate
    
def phototest2(data, livetime, energies, efficiencies, photoBackMean, sourceStrength=1):
    m = 0
    M = []
    PZ = 0
    lineO = [0, 0]
    pz = 0
    m = 0
    p1=0
    # print (data)
    # print (data[:,0])
    # print (data[:,1])
    # datOrg=data[:,1]
    dataMem = np.zeros(np.shape(data))
    dataMem[:, 0] = data[:, 0]
    
    #T=80
    #sigLow=10
    #sigUp=sigLow
    dataMem[:, 1] = VectorNucleationEfficiency(data[:, 1], energies, efficiencies)
    #print (dataMem)
    # for i in range(len(data)):
    # r=data[i,1]
    # dataMem[i]=NucleationEfficiencyTrue(r,T,sigLow,sigUp)
    pocket=[]
    for line in dataMem:
        #print ("Line: ", line)
        if line[0] == lineO[0]:
            pocket.append(lineO[1])
        if line[0] != lineO[0]:
            pocket.append(lineO[1])
            PZ += pz
            pz = 1
            m += 1
            p1p=0
            #print ("Pocket Full: ", pocket)
            plen=len(pocket)
            for i in range(plen):
                p1i=1
                for j in range(plen):
                    if i==j:
                        p1j=pocket[j]
                        #print ("p1j: ",p1j)
                    else:
                        p1j=1-pocket[j]
                        #print ("p1j: ",p1j)
                    p1i*=p1j
                #print ("p1i: ",p1i)
                p1p+=p1i
            p1+=p1p
            #print ("p1 pocket: ", p1p)
            pocket=[]
        pzi = 1 - line[1]
        pz *= pzi
        lineO = line
    PZ += pz
    m += 1
    backTot = 0
    for i in range(len(photoBackMean)):
        if i!=1:
            backTot += photoBackMean[i] * i
    datasum = sum(dataMem[:, 1])
    backCount = sum(photoBackMean)
    # print ("############")
    # print (dataMem)
    # print (len(dataMem))
    # print(energies,efficiencies)
    # print("Data Sum: ",datasum)
    # print ("Back Sum: ",backCount)
    # print ("############")
    #total = (datasum * sourceStrength) / livetime + backTot
    #count = (m * sourceStrength - PZ) / livetime + backCount
    total = ( (datasum - p1) * sourceStrength) / livetime + backTot
    count = (p1 * sourceStrength) / livetime + photoBackMean[1]   
    mean = total / count
    czero = PZ / livetime
    return total, count, mean, czero
    """    
        if i==0:
            m=0
            pz=1
        elif data[0,i]==data[0,i-1]:
            m=m
        else:
            if m==0:
                m=0
                M+=[m]
                PZ+=pz
            else:
                M+=[m]
                PZ+=pz
            m=0
            pz=1
        #print ("i: ",i)
        #print ("data: ",data)
        r=data[1,i]
        r=float(r)
        t=1
        weight=1
        mi=NucleationEfficiency(r,energies,efficiencies)
        pzi=(1-mi)
        pz*=pzi
        m+=mi
        if i==len(data[0,:])-1:
            if m==0:
                m=0
            else:
                M+=[m]
                PZ+=pz
        backTot=0
        for i in range(len(photoBackMean)):
            backTot+=photoBackMean[i]*i
        backCount=sum(photoBackMean)
        total=(sum(M)*sourceStrength)/livetime+backTot
        count=(len(M)*sourceStrength-PZ)/livetime+backCount
        mean=total/count
        czero=PZ/livetime
    return total,count,mean,czero
"""


def main(flist, alist, pnlist, aplist):
    # M=phototest("../Sb_fast.txt",500)
    # filename=sys.argv[1]

    Thomson = True
    Photoneutron = True
    Nuisance = True

    T = int(sys.argv[1])
    sigLow = int(sys.argv[2])
    sigUp = sigLow
    binsize = 1

    sourceErr = 0.05
    modeErrT = 0.1
    modeErrPN = 0.2
    backscale=1
    #backscale=1/3.473
    #backscale=1/1000 #reduce background
    
    background = 500*backscale
    backErr = np.round(background ** (1 / 2))
    # energies=[75,100,115,120,140]
    # efficiencies=[0,.2,.50,.8,1]
    # energies=[100,120]
    # efficiencies=[0,1]
    n = len(flist)
    RecoilList = []
    WeightList = []
    nuisanceT = np.zeros(n)
    TrueArray = np.zeros([3, n])
    InArray = np.zeros([4, n])
    UsePhotoError = []
    print(pnlist)
    pn = len(pnlist)
    #pnrateList = []
    #pnrecoilList = []
    #pnweightList = []
    time = 100
    # photoBackMean=np.zeros(20)
    # photoBackMean=np.ones(20)*50
    # photoBackMean[1]=5
    photoBackMean = np.array(
        [0, 5, 1.59, 0.69, 0.37, 0.16, 0.11, 0.08, 0.021, 0.016, 0.019, 0.022, 0, 0, 0, 0, 0, 0, 0, 0]
    )*backscale
    meanBackBubble=0
    meanBackEvent=np.sum(photoBackMean)
    for i in range(len(photoBackMean)):
        meanBackBubble+=photoBackMean[i]*i
    print("meanBackEvent: ", meanBackEvent)
    print ("meanBackBubble: ", meanBackBubble)
    pnlivetime = 20
    print("Photo-neutron livetime: ", pnlivetime)

    pn_sourceError = np.empty(pn)

    photoArrayTrue = np.zeros([pn, 20])

    #eventIDTestData=np.array([0,1,1,3,3,3,4,5,6,7,8,8,10,11])
    #recoilsTestData=np.array([50,50,150,50,150,150,50,150,150,50,50,150,150,50])
    #data=np.empty([14,2])
    #data[:,0]=eventIDTestData
    #data[:,1]=recoilsTestData
    #livetime=1
    #energies=[96,97,98,99,100]
    #efficiencies=[0, 0.2, 0.5, 0.8, 1]
    #print (data)
    #print (np.shape(data))
    #print (data[:,0])
    #total, count, mean, czero = phototest2(data, livetime, energies, efficiencies, photoBackMean, sourceStrength=1)
    #print ("#################################")
    #print ("Adjusted total/mean test: ", count,mean)
    #print (Chi_squared_11)
    # photoArray0=np.zeros([2,20])
    # photoArray=np.zeros([pn,20])
    # photoArrayC=np.zeros([2,20])
    # photoArrayB=np.zeros([2,20])
    photoneutrondata = []
    pnCommonMode = 0
    #while pnCommonMode<0.1:
    pnCommonMode = np.squeeze(np.random.normal(0, modeErrPN, 1))
    #pnCommonMode = np.random.normal(0, modeErrPN, 1)    
    for p in range(pn):
        print(p, pn)
        print(photoBackMean)
        pnc_array = np.loadtxt(pnlist[p]+"_recoils.txt")
        photoneutrondata += [pnc_array]
        print("")
        print("")
        print("********************")
        print("Photoneutron Data p: ")
        print(photoneutrondata[p])
        print("********************")
        print("")
        print("")
        photoArrayTrue[p, :] = phototestTrue(
            pnlist[p] + "_slow.txt", 2000, T, sigLow, sigUp
        )
        pn_sourceError[p] = np.squeeze(np.random.normal(0, sourceErr, 1))
        #pn_sourceError[p] = np.random.normal(0, sourceErr, 1)
        photoArrayTrue[p, :] = photoJitter(
            photoArrayTrue[p, :], photoBackMean, pn_sourceError[p],pnCommonMode, time=time
        )
        # print ("M:",M)
        print(photoBackMean)
    Mspecx = [1, 2, 3, 4]
    Mx = np.arange(1, 11, 1)
    # n4plus=1-(M[1]+M[2]+M[3])/np.sum(M)
    # print ("M 4+: ",n4plus)
    # n4plus=(M[4]+M[5]+M[6]+M[7]+M[8]+M[9])/np.sum(M)
    # print ("M 4+: ",n4plus)
    # plt.plot(Mspecx,Mspec/np.sum(Mspec))
    # plt.plot(Mx,M[1:11]/np.sum(M))
    # plt.show()
    # plt.savefig("photomult_bi207.png")
    # plt.clf()
    tCommonMode = 0
    #while tCommonMode<0.1:
    tCommonMode = np.squeeze(np.random.normal(0, modeErrT, 1))
    #tCommonMode = np.random.normal(0, modeErrT, 1)
    for i in range(n):
        print(i, flist[i])
        #Recoils, Weights, Rate, Count, t, time = analyze( flist[i], T, sigLow, sigUp, Activity=alist[i])
        Recoils, Weights, Rate, Count, t, time = analyze_Compton( flist[i], T, sigLow, sigUp, Activity=alist[i])  
        Recoils, Weights = specrafy(Recoils, Weights, binsize=binsize)
        print ("#######################")
        #print ("Model Energies: ", Recoils)
        print ("Max: ", np.max(Recoils), "Min: ", np.min(Recoils), "Shape: ",  np.shape(Recoils))
        #print ("Model Rates: ",Weights)
        print ("Sum: ", sum(Weights), "Mean: ", np.mean(Weights))
        print ("#######################")
        Rate2 = rateFinderTrue(Recoils, T, sigLow, sigUp, t, Weights)
        print("Recoils: ", Rate, Rate2, Rate - Rate2)
        Rate = np.squeeze(Rate2)
        Count = Rate * time
        TrueArray[0, i] = Rate
        TrueArray[1, i] = Count
        TrueArray[2, i] = t
        print ("#######################")
        print ("True Rates: ",Rate,Rate2,Count)
        print ("#######################")
        rate, count, backgroundexcess, nuisanceT[i] = rateJitter(
            Rate,
            Count,
            tCommonMode,
            time=time,
            sourceErr=sourceErr,
            background=background,
        )
        print(nuisanceT[i])
        # rate=Rate
        # count=Count
        InArray[0, i] = rate
        InArray[1, i] = count
        InArray[2, i] = t
        InArray[3, i] = backgroundexcess
        RecoilList += [Recoils]
        WeightList += [Weights]
        print ("Average Rate: ", Rate)
        print ("Count: ",count)
        print ("Background: ",backgroundexcess)
        print ("Source Error: ", nuisanceT[i])
        print ("Common Error: ", logErr(tCommonMode) )
        print ("Jittered Rate: ", rate)
        print ("Jittered Rate Check: ", backgroundexcess/time + count/time*logErr(tCommonMode)*nuisanceT[i])
        print(TrueArray[:, i])
        print(InArray[:, i])
    np.savetxt(fileprefix + "start.txt", InArray)
    print("THIS HERE")
    print (p)
    #print (np.shape(pnweightList))
    #print (pnweightList)
    #print(p, len(pnweightList[p]), len(pnrecoilList), len(pnweightList))
    #NitersRough = 5000
    #Niters = 1000
    
    NitersRough = 20000
    Niters = 4000000

    #NitersRough=10000
    #Niters=2000000
    
    #-3.426022520971181567e-02 -3.076520190371901178e-02
    #4.939968085322321567e-03 6.562545177614874381e-02
    #-2.400534680249160474e-02 2.250114429124067811e-02
    #1.403505933079687296e-01 8.379251897369664748e-02
    #1.922411485824054300e-01 1.268766240830859482e-01
    #1.431957932366064445e-01 4.553706020688447209e-02
    #-2.290043294275994137e-02 5.446431777108544890e-03
    
    #Niters = 150000
    Nwalkers = 4
    printstep = 100000
    energies0 = [0, 0, 0, 0 , 0]  # definition
    nuisance0 = np.zeros(n + pn)
    #nuisance0=[-3.426022520971181567e-02,4.939968085322321567e-03,-2.400534680249160474e-02,1.403505933079687296e-01,1.922411485824054300e-01]
    mode0 = [0, 0]
    #mode0=[1.431957932366064445e-01,-2.290043294275994137e-02]
    
    sample_n=10000
    
    samplestep=(Niters+1)//(sample_n)
    if samplestep==0:
        sample_n=Niters
        samplestep=1
    
    plater=np.zeros(len(energies0)+len(nuisance0)+len(mode0)+1)
    costco=np.zeros([sample_n,len(plater)])
    
    # pnuisance0=np.zeros(pn)
    bestChiList=[]

    efficiencies0 = [0, 0.2, 0.5, 0.8, 1]
    paracount = len(efficiencies0)
    energies0 = np.zeros(paracount)
    buffer = 5
    guessfloor = 50
    guessceil = 500
    bound = (guessceil - (paracount - 1) * buffer - guessfloor) / (paracount + 2)
    # print (paracount)
    print("Random Threshold Step: ", bound)
    energies0 = np.zeros(paracount)
    step = 1
    #n_step = sourceErr/10
    #m_step = modeErr/10
    #n_step = 0.01
    #m_step = 0.02
    #step = 1
    n_step = 0.003
    m_step = 0.001
    RoughFact = 2
    roughChiPenalty = 0
    X = 20
    step_shift=1

    #max_shift=1
    #min_shift=.1
    #baseChi=1
    
    #max_shift2=max_shift**2
    #min_shift2=min_shift**2
    # Chi0/T0/sigma0 is the initial guess
    # Chic/Tc/sigmac is the current basis, the iteration that the model will default back on if it does not select the new guess
    # Chii/Ti/sigmai is the most recent iteration
    # Chib/Tb/sigmab is the best fit so far, being stored to compare to future models
    OutArray = np.zeros([len(energies0) + 4, Nwalkers])
    chimode = np.zeros(2)
    chi = np.zeros(n)
    chiNuis = np.zeros(n+pn)
    chipn = np.zeros(pn)

    fullnuisance = np.zeros([n+pn+len(mode0),Nwalkers])
    DTotalB=np.zeros([n + pn,Nwalkers])
    DCountB=np.zeros([n + pn,Nwalkers])
    for nw in range(Nwalkers):
        energies0 = np.zeros(paracount)
        for j in range(paracount):
            if j == 0:
                energies0[j] = guessfloor + 3 * bound * np.random.rand()
            else:
                energies0[j] = energies0[j - 1] + bound * np.random.rand() + buffer
        #energies0 = [6.252983462309628493e+01,7.599864648098314035e+01,7.619371624296209689e+01,8.055708179978147143e+01,1.272788497247918968e+02] 
        #energies0 = [60,72,80,88,100]  # fixed initial guess
        #energies0 = [120,144,160,176,200]  # fixed initial guess
        #energies0 = [240,288,320,352,400]
        #energies0=[61.07507532,69.51169466,81.90851564,95.18992111,1037.24945314]
        #energies0=[70.03499573, 74.51155923, 75.96067757, 82.42290521, 1425.13382012]
        #energies0=[57.10434916, 77.46442932, 80.46113462, 82.23119383, 1639.12546188]
        print("Initial Guess: ", energies0)
        np.savetxt(fileprefix + "guess_" + str(nw + 1) + ".txt", energies0)
        totalTrue=np.zeros(n+pn)
        countTrue=np.zeros(n+pn)

        dtotalb = np.zeros(n + pn)
        dcountb = np.zeros(n + pn)
        dtotalc = np.zeros(n + pn)
        dcountc = np.zeros(n + pn)
        dtotali = np.zeros(n + pn)
        dcounti = np.zeros(n + pn)
        
        totalModeli = np.zeros(n + pn)
        countModeli = np.zeros(n + pn)
        totalModelc = np.zeros(n + pn)
        countModelc = np.zeros(n + pn)
        totalModelb = np.zeros(n + pn)
        countModelb = np.zeros(n + pn)
        for i in range(n):
            chi[i],dtotali[i],totalModeli[i],totalTrue[i] = test(
                RecoilList[i],
                InArray[0, i],
                energies0,
                efficiencies0,
                nuisance0[i],
                InArray[2, i],
                WeightList[i],
                background=background,
                time=time,
                mode=mode0[0],
            )
            #chi,dif,rateT,rate
            dcounti[i],countModeli[i],countTrue[i]=dtotali[i],totalModeli[i],totalTrue[i]
        chipn,dtotali[n:],dcounti[n:],totalTrue[n:],totalModeli[n:],countTrue[n:],countModeli[n:] = photoEval3(
            pn,
            photoneutrondata,
            pnlivetime,
            energies0,
            efficiencies0,
            photoArrayTrue,
            meanBackBubble,
            meanBackEvent,
            loud=True,
            pnuisance=nuisance0[i + 1 :],
            mode=mode0[1],
        )
        #chi,difTotal,difCount,totalTrue,totalModel,countTrue,countModel
        #chi,difTotal,difCount,total,trueTot,count,trueCount
        """
        for p in range(pn):
            #test(recoil,rate,energies,efficiencies,r_nuis,t,weight,background=500,time=100)
            #print (len(pnrecoilList[i]),len(RecoilList[i]))
            #chipn[p]=test(pnrecoilList[p],pnrateList[p],energies0,efficiencies0,pnuisance0[p],1,pnweightList[p],background=background,time=time)
            photoArray[p,:]=phototest(pnlist[p]+"_ultrafast.txt",10,energies0,efficiencies0)
            #photoArray[p,:]=phototest(pnlist[p]+"_ultrafast.txt",500,energiesi,efficiencies0)
            for mb in range(len(photoArray)):
                if photoArrayTrue[p,mb]>0:
                    chipn[p]+=((photoArray[p,mb]-photoArrayTrue[p,mb])**2/photoArrayTrue[p,mb])
            #for mb in range(len(photoArray)):
                #if photoArrayTrue[p,mb]>0:
                    #chipn[p]+=((photoArray[p,mb]-photoArrayTrue[p,mb])**2/photoArrayTrue[p,mb])
        """
        print("Chi Thomson: ", chi)
        print("Chi Photo-Neutron: ", np.sum(chipn))
        print("Chi photo: ", chipn)
        for i in range(n+pn):
            chiNuis[i] = ((nuisance0[i]) / sourceErr) ** 2
        Chi0 = 0
        if Thomson == True:
            Chi0 += sum(chi)
        if Photoneutron == True:
            Chi0 += sum(chipn)
        if Nuisance == True:
            Chi0 += sum(chiNuis)
        # Chi0=sum(chi)+sum(chipn) #both
        # Chi0=sum(chi) #thomson only
        # Chi0=sum(chipn) #photoneutron only
        Chib = Chi0
        Chic = Chi0
        MM = len(energies0)
        modei= np.empty(len(mode0))
        modec=np.empty(len(mode0))
        modeb=np.empty(len(mode0))
        modei = mode0[:]
        modec = mode0[:]
        modeb = mode0[:]
        energiesi = np.empty(MM)
        energiesb = np.empty(MM)
        energiesc = np.empty(MM)
        nuisancei = np.zeros(n + pn)
        nuisanceb = np.zeros(n + pn)
        nuisancec = np.zeros(n + pn)
        nuisancei[:]=nuisance0[:]
        nuisancec[:]=nuisance0[:]
        nuisanceb[:]=nuisance0[:]
        # pnuisancei=np.zeros(pn)
        # pnuisanceb=np.zeros(pn)
        # pnuisancec=np.zeros(pn)
        bestChiThom=np.empty(len(chi))
        bestChiPhot=np.empty(len(chipn))
        bestChiNuis=np.empty(len(chiNuis))
        bestChiMode=np.empty(len(chimode))
        
        for i in range(MM):
            value = energies0[i]
            energiesi[i] = value
            energiesb[i] = value
            energiesc[i] = value
        # print ("Defined")
        # print (energiesb)
        # print (energiesc)
        # print (Chi0,energies0)
        end = T + X * sigUp
        if end < 500:
            end=500
        grade, fifty, RT, RG, TOT = fittest(
            T, sigLow, sigUp, energiesb, efficiencies0, start=0, end=end
        )
        print("")
        print("")
        print("Guess Score: ", grade)
        print("Total Value: ", TOT)
        print("Guess Chi: ", Chi0)
        print(photoArrayTrue)
        # print (photoArray)
        # print ("Chi photo: ",chipn)
        print("50% Efficiency Point: ", fifty, T)
        print("Guess Parameters: ", energies0, efficiencies0)
        print("")
        print("")
        for ni in range(NitersRough):
            # print ("*****")
            # print ("Current Parameters: ", energiesc)
            # print ("*****")
            # energiesi[:],nuisancei[:]=stepper(energiesc,nuisancec,step*RoughFact,n_step)
            #step_shift=Chic/baseChi
            #if step_shift>max_shift2:
                #step_shift=max_shift
            #elif step_shift<min_shift2:
                #step_shift=min_shift
            #else:
                #step_shift=(step_shift)**(1/2)
            #stepi = step
            # energiesi[:],nuisancei[:]=stepper(energiesc,nuisancec,stepi,n_step)                                                                                                                                                                                               
            energiesi[:], nuisancei[:], modei[:] = stepper(
                energiesc, nuisancec, modec, step*step_shift*RoughFact, n_step*step_shift, m_step*step_shift, max=2000
            )
            if Nuisance == False:
                nuisancei[:] = nuisance0[:]  # turn nuisance parameters off for the whole fit 
                modei[:] = mode0[:] # whole fit ignores common mode uncertainty 
            nuisancei[:] = nuisance0[:]  # rough fit ignores nuisance parameters
            #modei[:] = mode0[:] # rough fit ignores common mode uncertainty
            # pnuisancei[:]=pnuisance0[:]
            for i in range(n):
                chi[i],dtotali[i],totalModeli[i],totalTrue[i] = test(
                    RecoilList[i],
                    InArray[0, i],
                    energiesi,
                    efficiencies0,
                    nuisancei[i],
                    InArray[2, i],
                    WeightList[i],
                    background=background,
                    time=time,
                    mode=modei[0],
                )
                dcounti[i],countModeli[i],countTrue[i]=dtotali[i],totalModeli[i],totalTrue[i]
                #dcountb[i],countModelb[i],dtotalb[i],totalModelb[i]=dcounti[i],countModeli[i],dtotali[i],totalModeli[i]
                #dcountc[i],countModelc[i],dtotalc[i],totalModelb[i]=dcounti[i],countModeli[i],dtotali[i],totalModeli[i]
                #chi,dif,rateT,rate
                if chi[i] < 0:
                    print("chi i Negative!: ", i)
            for i in range(n+pn):
                chiNuis[i] = ((nuisancei[i]) / sourceErr) ** 2
            if ni % printstep == 0:
                pnverb = True
            else:
                pnverb = False
            chipn,dtotali[n:],dcounti[n:],totalTrue[n:],totalModeli[n:],countTrue[n:],countModeli[n:]= photoEval3(
                pn,
                photoneutrondata,
                pnlivetime,
                energiesi,
                efficiencies0,
                photoArrayTrue,
                meanBackBubble,
                meanBackEvent,
                loud=pnverb,
                pnuisance=nuisancei[n:],
                mode=modei[1],
            )
            #chipn,dtotali[n:],dcounti[n:],totalTrue[n:],totalModeli[n:],countTrue[n:],countModeli[n:]
            dcountb[:],countModelb[:],dtotalb[:],totalModelb[:]=dcounti[:],countModeli[:],dtotali[:],totalModeli[:]
            dcountc[:],countModelc[:],dtotalc[:],totalModelb[:]=dcounti[:],countModeli[:],dtotali[:],totalModeli[:]
            chimode[0] = ( (modei[0]) /modeErrT) ** 2
            chimode[1] = ( (modei[1]) /modeErrPN) ** 2
            """
            for p in range(pn):
                #test(recoil,rate,energies,efficiencies,r_nuis,t,weight,background=500,time=100)
                #chipn[p]=test(pnrecoilList[p],pnrateList[p],energiesi,efficiencies0,pnuisancei[p],1,pnweightList[p],background=background,time=time)
                photoArray[p,:]=phototest(pnlist[p]+"_ultrafast.txt",10,energiesi,efficiencies0)
                for mb in range(len(photoArray)):
                    if photoArrayTrue[p,mb]>0:
                        chipn[p]+=((photoArray[p,mb]-photoArrayTrue[p,mb])**2/photoArrayTrue[p,mb])
            """
            Chii = 0
            if Thomson:
                Chii += sum(chi)
            if Photoneutron:
                Chii += sum(chipn)
            if Nuisance:
                if Thomson and Photoneutron:
                    Chii += sum(chiNuis)
                    Chii += sum(chimode)
                elif Thomson:
                    Chii += sum(chiNuis[:n])
                    Chii += chimode[0]
                elif Photoneutron:
                    Chii += sum(chiNuis[n:])
                    Chii += chimode[1]
                else:
                    print ("EMPTY RUN EMPTY RUN")
                    print ("no calibration methods")
                    x=redapple
            # Chii=sum(chi)+sum(chipn)+sum(chiNuis) #both
            # Chii=sum(chi)+roughChiPenalty #thomson only
            # Chii=sum(chipn)+sum(chiNuis) #photoneutron only
            #print (Chii)
            if Chii < 0:
                print("Negative!")
                Chii = 10**6
            if Chii == 0:
                print("Zero!")
                Chii == 10 ** (-1)
            bar = Chii / (Chic + Chii)
            #bar = Chic/Chii
            judge = np.random.uniform(0, 1)
            if ni % printstep == 0:
                print("*****")
                print("*****")
                print("Iteration: ", ni)
                print("Chi Thomson: ", sum(chi))
                print("Chi photo-n: ", sum(chipn))
                print("Chi photo: ", chipn)
                print("Chi nuisance: ", sum(chiNuis))
                print("Chi mode: ",sum(chimode))
                print("Chi Iteration: ", Chii)
                print("Current Chi: ", Chic)
                print("Best Chi: ", Chib)
                print("*****")
                print("Iteration Parameters: ", energiesi)
                print("Current Parameters: ", energiesc)
                print("Best Parameters: ", energiesb)
                print("*****")
                print("Iteration Nuisance: ", nuisancei)
                print("Current Nuisance: ", nuisancec)
                print("Best Nuisance: ", nuisanceb)
                print("Bar/Judge: ", bar, "/", judge)
                print("Iteration Mode Nuiance: ", modei)
                print("Current Mode Nuiance: ", modec)
                print("Best Mode Nuisance: ", modeb)
                print("Mode Thomson Error (real/model): ", tCommonMode, " / ", modeb[0])
                print("Mode Photoneutron Error (real/model): ", pnCommonMode, " / ", modeb[1])
                bestChiList+=[Chib]
                print("*****")
                print("*****")
            if Chii < Chib:
                Chib = Chii
                Chic = Chii
                energiesb[:] = energiesi[:]
                energiesc[:] = energiesi[:]
                nuisanceb[:] = nuisancei[:]
                nuisancec[:] = nuisancei[:]
                modeb[:] = modei[:]
                modec[:] = modei[:]
                dtotalc[:] = dtotali[:]
                dcountc[:] = dcounti[:]
                totalModelc[:] = totalModeli[:]
                countModelc[:] = countModeli[:]
                dtotalb[:] = dtotali[:]
                dcountb[:] = dcounti[:]
                totalModelb[:] = totalModeli[:]
                countModelb[:] = countModeli[:]
                bestChiThom[:]=chi[:]
                bestChiPhot[:]=chipn[:]
                bestChiNuis[:]=chiNuis[:]
                bestChiMode[:]=chimode[:]
            elif Chii < Chic:
                Chic = Chii
                energiesc[:] = energiesi[:]
                nuisancec[:] = nuisancei[:]
                modec[:] = modei[:]
                dtotalc[:] = dtotali[:]
                dcountc[:] = dcounti[:]
                totalModelc[:] = totalModeli[:]
                countModelc[:] = countModeli[:]
            elif judge > bar:
                Chic = Chii
                energiesc[:] = energiesi[:]
                nuisancec[:] = nuisancei[:]
                modec[:] = modei[:]
                dtotalc[:] = dtotali[:]
                dcountc[:] = dcounti[:]
                totalModelc[:] = totalModeli[:]
                countModelc[:] = countModeli[:]
            else:
                continue
        grade, fifty, RT, RG, TOT = fittest(
            T, sigLow, sigUp, energiesb, efficiencies0, start=0, end=end
        )
        print("")
        print("")
        print("Rough Score: ", grade)
        print("Total Value: ", TOT)
        print("Rough Chii :", Chib)
        print("Rough Parameters: ", energiesb)
        print("50% Efficiency Point: ", fifty, T)
        print("Model Effective Threshold: ", RG)
        print("True Effective Threshold: ", RT)
        print("Mode Photoneutron Error: ", pnCommonMode)
        print("Mode Thomson Error (real/model): ", tCommonMode, " / ", modeb[0])
        print("Mode Photoneutron Error (real/model): ", pnCommonMode, " / ", modeb[1])
        print("")
        print("")
        Chic = Chib
        energiesc[:] = energiesb[:]
        nuisancec[:] = nuisanceb[:]
        sample_i=0
        for ni in range(Niters):
            # stepi=step*(Niters-ni)/Niters
            #step_shift=Chic/baseChi
            #if step_shift>max_shift2:
                #step_shift=max_shift
            #elif step_shift<min_shift2:
                #step_shift=min_shift
            #else:
                #step_shift=(step_shift)**(1/2)
            #stepi = step
            # energiesi[:],nuisancei[:]=stepper(energiesc,nuisancec,stepi,n_step)
            energiesi[:], nuisancei[:], modei[:] = stepper(
                energiesc, nuisancec, modec, step*step_shift, n_step*step_shift, m_step*step_shift, max=2000
            )
            #modei[:] = mode0[:] # whole fit ignores common mode uncertainty 
            if Nuisance == False:
                nuisancei[:] = nuisance0[:]  # turn nuisance parameters off for the whole fit
                modei[:] = mode0[:] # whole fit ignores common mode uncertainty   
            for i in range(n):
                chi[i],dtotali[i],totalModeli[i],totalTrue[i] = test(
                    RecoilList[i],
                    InArray[0, i],
                    energiesi,
                    efficiencies0,
                    nuisancei[i],
                    InArray[2, i],
                    WeightList[i],
                    background=background,
                    time=time,
                    mode=modei[0],
                )
                dcounti[i],countModeli[i],countTrue[i]=dtotali[i],totalModeli[i],totalTrue[i]
                if chi[i] < 0:
                    print("chi i Negative!: ", i)
                #print (nuisancei[i],sourceErr)
            for i in range(n+pn):
                chiNuis[i] = ((nuisancei[i]) / sourceErr) ** 2
                #print (chiNuis[i])
                if chiNuis[i] < 0:
                    print("chi nuisance Negative!: ", i)
                    x = red
            if ni % printstep == 0:
                pnverb = True
            else:
                pnverb = False
            chipn,dtotali[n:],dcounti[n:],totalTrue[n:],totalModeli[n:],countTrue[n:],countModeli[n:] = photoEval3(
                pn,
                photoneutrondata,
                pnlivetime,
                energiesi,
                efficiencies0,
                photoArrayTrue,
                meanBackBubble,
                meanBackEvent,
                loud=pnverb,
                pnuisance=nuisancei[n:],
                mode=modei[1],
            )
            chimode[0] = ( (modei[0]) /modeErrT) ** 2
            chimode[1] = ( (modei[1]) /modeErrPN) ** 2
            # for p in range(pn):
            """
            for p in range(pn):
                #test(recoil,rate,energies,efficiencies,r_nuis,t,weight,background=500,time=100)
                #chipn[p]=testLoud(pnrecoilList[p],pnrateList[p],energiesi,efficiencies0,pnuisance0[p],1,pnweightList[p],background=background,time=time)
                #chipn[p]=test(pnrecoilList[p],pnrateList[p],energiesi,efficiencies0,pnuisancei[p],1,pnweightList[p],background=background,time=time)
                photoArray[p,:]=phototest(pnlist[p]+"_ultrafast.txt",10,energiesi,efficiencies0)
                #for mb in range(len(photoArray)):
                    #if photoArrayTrue[p,mb]>0:
                        #chipn[p]+=((photoArray[p,mb]-photoArrayTrue[p,mb])**2/photoArrayTrue[p,mb])
                for mb in range(len(photoArray)):
                    if photoArrayTrue[p,mb]>0:
                        chipn[p]+=((photoArray[p,mb]-photoArrayTrue[p,mb])**2/photoArrayTrue[p,mb])
            """
            Chii = 0
            if Thomson:
                Chii += sum(chi)
            if Photoneutron:
                Chii += sum(chipn)
            if Nuisance:
                if Thomson and Photoneutron:
                    Chii += sum(chiNuis)
                    Chii += sum(chimode)
                elif Thomson:
                    Chii += sum(chiNuis[:n])
                    Chii += chimode[0]
                elif Photoneutron:
                    Chii += sum(chiNuis[n:])
                    Chii += chimode[1]
                else:
                    print ("EMPTY RUN EMPTY RUN")
                    print ("no calibration methods")
                    x=redapple
            # Chii=sum(chi)+sum(chiNuis)+sum(chipn) #both
            # Chii=sum(chi)+sum(chiNuis) #no nuisance both
            # Chii=sum(chi)+sum(chiNuis) #thomson only
            # Chii=sum(chiNuis)+sum(chipn) #photoneutron only
            #print (Chii)
            if Chii < 0:
                Chii = 10**6
                print("Negative!")
            if Chii == 0:
                print("Zero!")
                Chii = 10 ** (-3)
            bar = Chii / (Chic + Chii)
            judge = np.random.uniform(0, 1)
            if ni % samplestep==0:
                sample_i+=1
                plater[0]=Chii
                plater[1:MM+1]=energiesi
                plater[MM+1:MM+n+pn+1]=nuisancei
                plater[MM+n+pn+1:]=modei
                #sample_i=int(ni/samplestep)
                if sample_i<sample_n:
                    costco[sample_i,:]=plater
            if ni % printstep == 0:
                print("*****")
                print("*****")
                print("Iteration: ", ni)
                print("Chi Thomson: ", sum(chi))
                print("Chi photo-n: ", sum(chipn))
                print("Chi photo: ", chipn)
                print("Chi nuisance: ", sum(chiNuis))
                print("Chi mode: ",sum(chimode))
                print("Chi Iteration: ", Chii)
                print("Current Chi: ", Chic)
                print("Best Chi: ", Chib)
                # print (photoArrayTrue)
                # print (photoArray)
                print("*****")
                print("Iteration Parameters: ", energiesi)
                print("Current Parameters: ", energiesc)
                print("Best Parameters: ", energiesb)
                print("*****")
                print("Iteration Nuisance: ", nuisancei)
                print("Current Nuisance: ", nuisancec)
                print("Best Nuisance: ", nuisanceb)
                print("Bar/Judge: ", bar, "/", judge)
                print("Iteration Mode Nuiance: ", modei)
                print("Current Mode Nuiance: ", modec)
                print("Best Mode Nuisance: ", modeb)
                print("Mode Thomson Error (real/model): ", tCommonMode, " / ", modeb[0])
                print("Mode Photoneutron Error (real/model): ", pnCommonMode, " / ", modeb[1])
                bestChiList+=[Chib]
                print("*****")
                print("*****")
            if Chii < Chib:
                Chib = Chii
                energiesb[:] = energiesi[:]
                energiesc[:] = energiesi[:]
                nuisanceb[:] = nuisancei[:]
                nuisancec[:] = nuisancei[:]
                modeb[:] = modei[:]
                modec[:] = modei[:]
                
                dtotalc[:] = dtotali[:]
                dcountc[:] = dcounti[:]
                totalModelc[:] = totalModeli[:]
                countModelc[:] = countModeli[:]
                dtotalb[:] = dtotali[:]
                dcountb[:] = dcounti[:]
                totalModelb[:] = totalModeli[:]
                countModelb[:] = countModeli[:]
                bestChiThom[:]=chi[:]
                bestChiPhot[:]=chipn[:]
                bestChiNuis[:]=chiNuis[:]
                bestChiMode[:]=chimode[:]
            elif Chii < Chic:
                Chic = Chii
                energiesc[:] = energiesi[:]
                nuisancec[:] = nuisancei[:]
                modec[:] = modei[:]
                dtotalc[:] = dtotali[:]
                dcountc[:] = dcounti[:]
                totalModelc[:] = totalModeli[:]
                countModelc[:] = countModeli[:]
                #dtotalb[:] = dtotali[:]
                #dcountb[:] = dcounti[:]
                #totalModelb[:] = totalModeli[:]
                #countModelb[:] = countModeli[:]
            elif judge > bar:
                Chic = Chii
                energiesc[:] = energiesi[:]
                nuisancec[:] = nuisancei[:]
                modec[:] = modei[:]
                dtotalc[:] = dtotali[:]
                dcountc[:] = dcounti[:]
                totalModelc[:] = totalModeli[:]
                countModelc[:] = countModeli[:]
            else:
                continue
        grade, fifty, RT, RG, TOT = fittest(
            T, sigLow, sigUp, energiesb, efficiencies0, start=0, end=end
        )
        teff = end - RT
        teffG = end - RG
        print("Real Effective Threshold:  ", teff)
        print("Model Effective Threshold: ", teffG)
        teffdiff = teffG / teff - 1
        OutArray[0, nw] = grade
        OutArray[1, nw] = Chib
        OutArray[2, nw] = fifty / T - 1
        OutArray[3, nw] = teffdiff
        OutArray[4:, nw] = energiesb[:]
        print("")
        print("")
        print("Post Run Score: ", grade)
        print("Total Value: ", TOT)
        print("Final Chii :", Chib)
        print("50% Efficiency Point: ", fifty, T)
        print("Final Parameters: ", energiesb)
        print("Mode Thomson Error (real/model): ", tCommonMode, " / ", modeb[0])
        print("Mode Photoneutron Error (real/model): ", pnCommonMode, " / ", modeb[1])
        print(nuisanceT)
        print(nuisanceb)
        print("Mode  Error: ", pnCommonMode)
        print("Mode Thomson Error: ", tCommonMode)
        print("Mode Nuisance: ", modeb)
        for i in range(n):
            print(i, flist[i])
            print("Source Strength Uncertainty")
            print("True: ", (nuisanceT[i] - 1) * 100, "%")
            print("Model: ", nuisanceb[i] * 100, "%")
        print("Total Nuisance Squared: ", sum(nuisanceT**2))
        # print ("Total Unfit Nuisance: ",sum((nuisanceT-nuisanceb)**2))
        print("")
        print("")
        # postfitplot(Tb,sigmab,zero,twenty,fifty,eighty,onehundred,start=zero-10,end=onehundred+10)

        fullnuisance[:n+pn,nw]=nuisanceb[:]
        fullnuisance[n+pn,nw]=modeb[0]
        fullnuisance[n+pn+1,nw]=modeb[1]

        for ex in range(len(dcountb)):
            print ("Source #",ex+1)
            print ("Total - Model - True - Excess: ",countModelb[ex],countTrue[ex],dcountb[ex])
            print ("Total - Model - True - Excess: ",totalModelb[ex],totalTrue[ex],dtotalb[ex])
            print ("#####################################################################")

        DTotalB[:,nw]=dtotalb[:]
        DCountB[:,nw]=dcountb[:]
        np.savetxt(fileprefix + "fullnuisance.txt", fullnuisance)
        
        np.savetxt(fileprefix + "count_difference.txt", DCountB)
        np.savetxt(fileprefix + "total_difference.txt", DTotalB)

        np.savetxt(fileprefix + "photoArrayTrue.txt", photoArrayTrue)
        np.savetxt(fileprefix + "InArray.txt", InArray)
        np.savetxt(fileprefix + "fit.txt", OutArray)
        np.savetxt(fileprefix + "convergence.txt",bestChiList)

        np.savetxt(fileprefix + "samples_"+str(nw)+".txt",costco)
        
        print ("Best Chi Thomson: ", sum(bestChiThom),bestChiThom)
        print ("Best Chi PN: ", sum(bestChiPhot),bestChiPhot)
        print ("Best Chi Nuisance: ", sum(bestChiNuis),bestChiNuis)
        print ("Best Chi Mode: ", sum(bestChiMode),bestChiMode)
    print(OutArray)


# flist = ["NewRuns/Eu152JAEA.txt","NewRuns/Bi207JAEA.txt","NewRuns/Y88JAEA.txt","NewRuns/Th228JAEA.txt"]
# flist = ["Eu152/JAEA.txt", "Bi207/JAEA.txt", "Sb124/JAEA.txt", "Th228/JAEA.txt", "Y88/JAEA.txt"]
# flist = ["Sb124/JAEA.txt","Y88/JAEA.txt","Th228/JAEA.txt"]
# "Eu152/JAEA.txt", "Bi207/JAEA.txt", "Sb124/JAEA.txt", "Th228/JAEA.txt", "Y88/JAEA.txt"
# flist = ["Bi207/JAEA.txt","Sb124/JAEA.txt","Eu152/JAEA.txt","Y88/JAEA.txt"]
# flist=["Bi207/JAEA.txt", "Sb124/JAEA.txt", "Th228/JAEA.txt", "Y88/JAEA.txt"]
# flist=["Bi207/JAEA.txt","Th228/JAEA.txt"]
# flist=["../Th228/JAEA.txt","../Bi207/JAEA.txt"]

# phototest()

main(flist, alist, pnlist, aplist)
