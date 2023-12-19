import matplotlib.pyplot as plt
import sys
import numpy as np
import math
import scipy.optimize
import random

################################################
flist=["./Sb124JAEA.txt","./Co60JAEA.txt","./Th228JAEA.txt"]
################################################
pnlist=["./Sb","./Bi"]
################################################
alist=[20,100,10]
################################################
aplist=[1,8] #this is not yet implimented

try:
    fileprefix=sys.argv[1]
except:
    fileprefix="./Test_Dump/test2"
save_path = "/data/runzezhang/result/chi2_test/"
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

def NucleationEfficiencyTrue(r,T,sigLow,sigUp):
    if r<T:
        R=1/2*(1+math.erf((r-T)/(sigLow*2**(1/2))))
    else:
        # 1-R?
        R=1/2*(1+math.erf((r-T)/(sigUp*2**(1/2))))
    return R

def NucleationEfficiency(r,energies,efficiencies="auto"):
    # unknow energy efficiency 2-d array, get r's efficiency value by interception
    n=len(energies)
    if efficiencies=="auto":
        print ("Automation")
        efficiencies=np.zeros(n)
        for i in range(n):
            efficiencies[i]=i/(n-1)
    #print (efficiencies)
    for i in range(n):
        if r<energies[i]:
            if i==0:
                return 0
            else:
                return (efficiencies[i]-efficiencies[i-1])*(r-energies[i-1])/(energies[i]-energies[i-1])+efficiencies[i-1]
    return 1
def VectorNucleationEfficiency(r,energies,efficiencies):
    n=len(energies)
    en0=energies[0]
    en30=energies[1]
    en70=energies[2]
    en100=energies[3]
    return np.piecewise(r,[r<=en0,(r>en0) & (r<=en30),(r>en30) & (r<=en70),(r>en70) & (r<en100),r>=en100],[lambda r: 0,lambda r: (0.3)/(en30-en0)*(r-en0),lambda r: (0.4)/(en70-en30)*(r-en30)+0.3,lambda r: (0.3)/(en100-en70)*(r-en70)+0.7,lambda r: 1])
    
def rateFinder(recoil,energies,efficiencies,t,weight):
    #print ("rateFinder: ",data,energy,T,sigma,t,weight)
    W=np.ones(len(recoil))
    for i in range(len(recoil)):
        #print (np.shape(data))
        #print (data[i])
        #print (np.shape(data[i]))
        W[i]=NucleationEfficiency(recoil[i],energies,efficiencies=efficiencies)*weight[i]/t
    return sum(W)

def specrafy(recoils,weights,binsize=1):
    # print("specrafy")
    minr=np.floor(np.min(recoils))
    maxr=np.ceil(np.max(recoils))
    #if maxr>600:
        #maxr=600
    spread=int(maxr-minr)
    
    recoilS=np.arange(minr,maxr,binsize)
    lenS=len(recoilS)
    weightS=np.zeros(lenS)
    
    for i in range(len(recoils)):
        for j in range(lenS):
            if recoils[i]>recoilS[j]-binsize/2 and recoils[i]<=recoilS[j]+binsize/2:
                weightS[j]+=weights[i]
    print (recoilS,weightS)
    return recoilS,weightS
# weightS is finally recoils counts by a factor of its weights
    
def load_photN(name):
    spec1=np.loadtxt(name+"_spec1.txt")
    spec2=np.loadtxt(name+"_spec2.txt")
    spec3=np.loadtxt(name+"_spec3.txt")
    specn=np.loadtxt(name+"_specn.txt")
    return spec1,spec2,spec3,specn

def load_photN2(name):
    spec11=np.loadtxt(name+"_spec11.txt")
    spec21=np.loadtxt(name+"_spec21.txt")
    spec22=np.loadtxt(name+"_spec22.txt")
    spec31=np.loadtxt(name+"_spec31.txt")
    spec32=np.loadtxt(name+"_spec32.txt")
    spec33=np.loadtxt(name+"_spec33.txt")
    specn1=np.loadtxt(name+"_specn1.txt")
    specn2=np.loadtxt(name+"_specn2.txt")
    specn3=np.loadtxt(name+"_specn3.txt")
    specnn=np.loadtxt(name+"_specnn.txt")
    return spec11,spec21,spec22,spec31,spec32,spec33,specn1,specn2,specn3,specnn

def neutron_rates(spec1,T,sigLow,sigUp,t):
    recoil=spec1[:,0]
    weight=spec1[:,1]
    rate=rateFinderTrue(recoil,T,sigLow,sigUp,t,weight)
    return rate

def multiplicer(spec1,spec2,spec3,specn,T,sigLow,sigUp,t):
    r=spec1[:,0]
    W1,W2,W3,Wn=spec1[:,1],spec2[:,1],spec3[:,1],specn[:,1]
    #print (spec2)
    R1=np.sum(W1)
    R2=np.sum(W2)
    R3=np.sum(W3)
    Rn=np.sum(Wn)

    print (R1,R2,R3,Rn)

    w1,w2,w3,wn=W1/R1,W2/R2,W3/R3,Wn/Rn

    
    r1=R1-R2-R3-Rn
    r2=R2-R3-Rn
    r3=R3-Rn
    rn=Rn
    
    B1=0
    B2=0
    B3=0
    Bn=0
    
    p1=rateFinderTrue(r,T,sigLow,sigUp,1,w1)
    p2=rateFinderTrue(r,T,sigLow,sigUp,1,w2)
    p3=rateFinderTrue(r,T,sigLow,sigUp,1,w3)
    pn=rateFinderTrue(r,T,sigLow,sigUp,1,wn)

    n1,n2,n3,nn=1-p1,1-p2,1-p3,1-pn
    print (p1,p2,p3,pn)
    print (n1,n2,n3,nn)
    
    B11=r1*p1
    B12=r2*(p1*n2+n1*p2)
    B13=r3*(p1*n2*n3 + n1*p2*n3 + n1*n2*p3)
    B1n=rn*(p1*n2*n3*nn + n1*p2*n3*nn + n1*n2*p3*nn + n1*n2*n3*pn)
    B1=(B11+B12+B13+B1n)/t
    
    B22=r2*p1*p2
    B23=r3*(p1*p2*n3 + p1*n2*p3 + n1*p2*p3)
    B2n=rn*(p1*p2*n3*nn + p1*n2*p3*nn + p1*n2*n3*pn + n1*p2*p3*nn + n1*p2*n3*pn + n1*n2*p3*pn)
    B2=(B22+B23+B2n)/t
    
    B33=r3*p1*p2*p3
    B3n=rn*(p1*p2*p3*nn + p1*p2*n3*pn + p1*n2*p3*pn + n1*p2*p3*pn)
    B3=(B33+B3n)/t
    
    Bnn=rn*p1*p2*p3*pn
    Bn=Bnn/t

    return B1,B2,B3,Bn

#np.savetxt(name+"spec1.txt",spec1)

def multiplicer2(spec11,spec21,spec22,spec31,spec32,spec33,specn1,specn2,specn3,specnn,T,sigLow,sigUp,t):
    r=spec11[:,0]
    W11,W21,W22,W31,W32,W33,Wn1,Wn2,Wn3,Wnn=spec11[:,1],spec21[:,1],spec22[:,1],spec31[:,1],spec32[:,1],spec33[:,1],specn1[:,1],specn2[:,1],specn3[:,1],specnn[:,1]
    #print (spec2)

    print ("Shape spec11: ",np.shape(spec11))
    print ("Shape W11: ",np.shape(W11))
    #print ("spec11[:,1] ",spec11[:,1])
    #print ("W11: ",W11)
    
    R1=np.sum(W11)
    R2=np.sum(W21)
    R3=np.sum(W31)
    Rn=np.sum(Wn1)
    Rnn=np.sum(Wnn)

    print (R1,R2,R3,Rn,Rnn)

    w11,w21,w22,w31,w32,w33,wn1,wn2,wn3,wnn=W11/R1,W21/R2,W22/R2,W31/R3,W32/R3,W33/R3,Wn1/Rn,Wn2/Rn,Wn3/Rn,Wnn/Rn #Wnn/Rn may be wrong, this might need to be over Rnn?

    r1=R1 #the rate of one recoil long events
    r2=R2 #the rate of two recoil long events
    r3=R3 #the rate of three recoil long events
    rn=Rn
    rnn=Rnn
    
    B1=0
    B2=0
    B3=0
    Bn=0

    #pij is the nucleation probability of the jth recoil in an event with i recoils
    p11=rateFinderTrue(r,T,sigLow,sigUp,1,w11) #the chance of the first (and only recoil) in a one recoil event of nucleating a bubble
    p21=rateFinderTrue(r,T,sigLow,sigUp,1,w21) #the chance of the first recoil in a two recoil event of nucleating a bubble
    p22=rateFinderTrue(r,T,sigLow,sigUp,1,w22) #the chance of the second recoil in a two recoil event of nucleating a bubble
    p31=rateFinderTrue(r,T,sigLow,sigUp,1,w31) #and so on
    p32=rateFinderTrue(r,T,sigLow,sigUp,1,w32)
    p33=rateFinderTrue(r,T,sigLow,sigUp,1,w33)
    pn1=rateFinderTrue(r,T,sigLow,sigUp,1,wn1)
    pn2=rateFinderTrue(r,T,sigLow,sigUp,1,wn2)
    pn3=rateFinderTrue(r,T,sigLow,sigUp,1,wn3)
    pnn=rateFinderTrue(r,T,sigLow,sigUp,1,wnn)

    #nij is the probability of the jth recoil in an i-long chain NOT forming a bubble, 1-pij
    n11,n21,n22,n31,n32,n33,nn1,nn2,nn3,nnn=1-p11,1-p21,1-p22,1-p31,1-p32,1-p33,1-pn1,1-pn2,1-pn3,1-pnn
    print ("Prob First Scatter: ",p11,p21,p31,pn1)
    
    B01=r1*n11
    B02=r2*n21*n22+n21*p22
    B03=r3*n31*n32*n33
    B0n=rn*nn1*nn2*nn3*nnn
    B0=(B01+B02+B03+B0n)/t

    B11=r1*p11
    B12=r2*(p21*n22+n21*p22)
    B13=r3*(p31*n32*n33 + n31*p32*n33 + n31*n32*p33)
    B1n=rn*(pn1*nn2*nn3*nnn + nn1*pn2*nn3*nnn + nn1*nn2*pn3*nnn + nn1*nn2*nn3*pnn)
    B1=(B11+B12+B13+B1n)/t

    B22=r2*p21*p22
    B23=r3*(p31*p32*n33 + p31*n32*p33 + n31*p32*p33)
    B2n=rn*(pn1*pn2*nn3*nnn + pn1*nn2*pn3*nnn + pn1*nn2*nn3*pnn + nn1*pn2*pn3*nnn + nn1*pn2*nn3*pnn + nn1*nn2*pn3*pnn)
    B2=(B22+B23+B2n)/t

    B33=r3*p31*p32*p33
    B3n=rn*(pn1*pn2*pn3*nnn + pn1*pn2*nn3*pnn + pn1*nn2*pn3*pnn + nn1*pn2*pn3*pnn)
    B3=(B33+B3n)/t

    Bnn=rnn*pn1*pn2*pn3*pnn
    Bn=Bnn/t
    BT=B0+B1+B2+B3+Bn
    rT=r1+r2+r3+rn
    print ("B vs R ", BT,rT)
    
    return B1,B2,B3,Bn

def rateFinderTrue(recoil,T,sigLow,sigUp,t,weight):
    #print ("rateFinder: ",data,energy,T,sigma,t,weight)
    weight=weight
    W=np.ones(len(recoil))
    for i in range(len(recoil)):
        #print (np.shape(data))
        #print (data[i])
        #print (np.shape(data[i]))
        W[i]=NucleationEfficiencyTrue(recoil[i],T,sigLow,sigUp)*weight[i]/t
    return sum(W)

def test(recoil,rate,energies,efficiencies,r_nuis,t,weight,background=500,time=100):
    rateT=rateFinder(recoil,energies,efficiencies,t,weight)*(1+r_nuis)+(background/time)
    #print (rateT)
    #print (rate)
    chi= ( (rateT-rate)**2/(rate) ) * time
    return chi

def testLoud(recoil,rate,energies,efficiencies,r_nuis,t,weight,background=500,time=100):
    rateT=rateFinder(recoil,energies,efficiencies,t,weight)*(1+r_nuis)+(background/time)
    print ("Model Rate (No Background): ",rateT)
    rateT+=background/time
    print ("Model Rate: ",rateT)  
    chi=(rateT-rate)**2/rate
    print ("Data Rate: ",rate)
    return chi

def rateJitter(rate,count, time=100, sourceErr=.03,background=500,s_err="False"):
    if s_err=="False":
        s_err=np.random.normal(1,sourceErr,1)
    #s_err=1+random.uniform(-1,1)*sourceErr
    count=count*s_err
    count=np.random.poisson(count)
    backgroundexcess = np.random.poisson(background)
    count+=backgroundexcess
    rate=count/time
    b_err=backgroundexcess/background - 1
    print ("Background: ",backgroundexcess,background)
    print ("Source Error: ",s_err)
    return rate, count,background, s_err
        #(T,sigLow,sigUp,energiesb,efficiencies0,start=T-4*sigLow,end=T-4*sigUp)

def fittest(T,sigLow,sigUp,energiesM,efficienciesM,M=10000,start=40,end=200):
    rt=0
    rg=0
    dif=0
    tot=0
    calib=0
    fifty=0
    fiftydif=.5
    #print (T,sigLow,sigUp)
    #print (energiesM,efficienciesM)
    for i in range(M):
        stepwidth=end/M
        r=i*stepwidth
        RT=NucleationEfficiencyTrue(r,T,sigLow,sigUp)
        #print ("Real Efficiency: ",RT)
        RG=NucleationEfficiency(r,energiesM,efficiencies=efficienciesM)
        #print ("Model Efficiency: ",RG)
        #print ("Difference: ",RT-RG)
        dif+=abs(RT-RG)*stepwidth
        tot+=RT*stepwidth
        if fiftydif>=abs(.5-RG):
            fiftydif=.5-RG
            fifty=r
        RT*=stepwidth
        RG*=stepwidth
        rt+=RT
        rg+=RG
        #print ("Absolute Difference: ",dif)
    return dif,fifty,rt,rg,tot
        


t=0.03754/100 #hours to run the whole thing,

T=105
sigma=15

def analyze(file,T,sigLow,sigUp, N=2*10**9, Activity=100,time=100):
    # N is number of events
    # time is calibration time
    # t is the simulated life time
    Data=np.loadtxt(file)
    # print("analysing...", Data[:10])
    SourceRate=(3.7*10**6)*Activity/100
    Energy=Data[:,2]
    Recoils=Data[:,4]
    Weights=Data[:,5]
    
    t = N/SourceRate #live time in seconds
    t /= 3600 #live time in hours
    
    Rate=rateFinderTrue(Recoils,T,sigLow,sigUp,t,Weights)
    
    Count = Rate * time
    
    return Recoils, Weights, Rate, Count,t, time

def stepper(energies,nuisance,step,n_step,min=20,max=2000):
    energiesMem=np.zeros(len(energies))
    nuisanceMem=np.zeros(len(nuisance))
    energiesMem[:]=energies[:]
    nuisanceMem[:]=nuisance[:]
    chosen=np.random.randint(0,len(energies))
    for i in range(len(energiesMem)):
        energiesMem[i]=np.random.normal(0,step,1)+energiesMem[i]
        if i==0:
            if energiesMem[0]<min:
                energiesMem[0]=min+5
            if energiesMem[1]<min:
                energiesMem[1]=min+5
            if energiesMem[0]>energiesMem[1]:
                low=energiesMem[1]
                energiesMem[0]=energiesMem[1]
                energiesMem[1]=low
        elif i==len(energiesMem)-1:
            if energiesMem[i]<energiesMem[i-1]:
                energiesMem[i]=energiesMem[i-1]
            if energiesMem[i]<min:
                energiesMem[i]=min+5
            if energiesMem[i]>max:
                energiesMem[i]=max-10
            if energiesMem[i]<energiesMem[i-1]:
                low=energies[i]
                energiesMem[i]=energiesMem[i-1]
                energiesMem[i-1]=low
        else:
            if energiesMem[i]>max:
                energiesMem[i]=max-10
            if energiesMem[i]<min:
                energiesMem[i]=min+5
            if energiesMem[i]<energiesMem[i-1]:
                low=energiesMem[i]
                energiesMem[i]=energiesMem[i-1]
    chosen=np.random.randint(0,len(nuisance))
    for i in range(len(nuisance)):
        #if i==chosen:
        if True:
            nuisance[i]+=np.random.normal(0,n_step,1)
        if nuisance[i]>0.1:
            nuisance[i]=.1
        if nuisance[i]<-.1:
            nuisance[i]=-.1
    return energiesMem[:],nuisanceMem[:]

def photoread(pnc_filepath):
    r_list=[]
    m_list=[]
    data=np.loadtxt(pnc_filepath,skiprows=1,dtype="str")
    m=0
    for i in range(len(data)):
        if i==0:
            m=0
        elif data[i,0]==data[i-1,0]:
            m=m
        else:
            m+=1
        r=float(data[i,1])
        u=data[i,2]
        if u=="eV":
            r=r
        elif u=="keV":
            r=1000*r
        elif u=="MeV":
            r=1000000*r
        r=int(np.round(r,0))
        r_list+=[r]
        m_list+=[m]
    pnc_array=np.empty([len(r_list),2])
    pnc_array[:,0]=m_list
    pnc_array[:,1]=r_list
    print ("crosscheck: ",len(data),m)
    return pnc_array

def photoEval(pn,pnlist,photoArray,photoArrayTrue):
    for p in range(pn):
        photoArray[p,:]=phototestmean(pnlist[p]+"_fast.txt",500,energies0,efficiencies0)
        for mb in range(len(photoArray)):
            if photoArrayTrue[p,mb]>0:
                chipn[p]+=((photoArray[p,mb]-photoArrayTrue[p,mb])**2/photoArrayTrue[p,mb])

def photoEval2(pn,photoneutrondata,livetime,energies,efficiencies,photoArrayTrue,photoBackMean,loud=False,pnuisance=[0,0]):
    chi=np.zeros(pn)
    for p in range(pn):
        trueCount=sum(photoArrayTrue[p,1:])
        trueMean=0
        trueTot=0
        for i in range(20):
            trueMean+=photoArrayTrue[p,i]*(i)/trueCount
            trueTot+=photoArrayTrue[p,i]*(i)
        #trueMean=np.mean(photoArrayTrue[p,:])
        #print ("photoEval2")
        #print("photoneutrondata",photoneutrondata)
        #print ("photoneutrondata[p]", photoneutrondata[p])
        total,count,mean,czero=phototest2(photoneutrondata[p],livetime,energies,efficiencies,photoBackMean,sourceStrength=pnuisance[p]+1)
        #modelCount=count-czero
        chiTotal=( (total-trueTot)**2/trueTot ) * 100
        chiCount=( (count - trueCount)**2/trueCount ) * 100
        chi[p]=chiTotal+chiCount
        if loud:
            print ("True Mults: ",photoArrayTrue[p,:])
            print ("True Count: ",trueCount)
            print ("True Zeros: ",photoArrayTrue[p,0])
            print ("Model Zeros: ",czero)
            print (">0 True Count: ",sum(photoArrayTrue[p,1:]))
            print ("Model Count: ",count)
            print ("Total/Mean True: ",trueTot/trueCount)
            print ("Total/Mean Model: ",total/count)
            print ("True Mean: ",trueMean)
            print ("Model Mean: ",mean)
            print ("True Sum: ", trueTot)
            print ("Model Sum: ",total)
            print ("Chi total: ", chiTotal)
            print ("Chi count: ", chiCount)
            print (chi[p])
    return chi
def photoJitter(photoSourceTrue,photoBackMean,sourceError,time=100):
    #poisson uncertainty in background + events
    #source strength uncertainty
    photoBack=np.zeros(len(photoBackMean))
    photoSourceTrue*=(sourceError)
    photoSourceTrue*=time
    photoBack=time*photoBackMean
    for i in range(len(photoSourceTrue)):
        photoSourceTrue[i]=np.random.poisson(photoSourceTrue[i])
        photoBack[i]=np.random.poisson(photoBack[i])
    photoSourceTrue= ( photoSourceTrue + photoBack ) / time
    return photoSourceTrue

def phototestTrue(pnc_filepath,livetime,T,sigLow,sigUp):
    #data1=np.loadtxt("../photoNC/Informacion_Sb124_high1.txt",skiprows=1,dtype="str")
    data=np.loadtxt(pnc_filepath,skiprows=1,dtype="str")
    # return a 1-d array with at i'th position, record number of events of multiplicity of i
    print(pnc_filepath, data[:10])
    #data=data[0:2000]
    #print("Neutron output shape: ",np.shape(data))
    m=0
    M=np.zeros(20)
    bsum=0
    #for i in range(500):
    for i in range(len(data)):
        if i==0:
            m=0
        elif data[i,0]==data[i-1,0]:
            #print ("Same: ", i, data[i,0])
            m=m
        else:
            if m==0:
                m=0
                M[0]+=1
            else:
                M[m]+=1
                #print ("Recorded: ",m,i, M[m])
            m=0
        r=float(data[i,1])
        u=data[i,2]
        if u=="eV":
            r=r
        elif u=="keV":
            r=1000*r
        elif u=="MeV":
            r=1000000*r
        else:
            #print ("Help!",i+1,u)
            fiftyfive=55
        t=1
        weight=1
        popper=np.random.uniform(0,1)
        blower=NucleationEfficiencyTrue(r,T,sigLow,sigUp)
        bsum+=blower
        #print ("blower: ",blower)
        #print ("popper: ",popper, i, r)
        if popper<blower:
            m+=1
            #print ("Tracker: ",i,m)
        if i==len(data)-1:
            if m==0:
                m=0
                M[0]+=1
            else:
                M[m]+=1
                #print ("Recorded: ",m,i, M[m])
    #print ("M :",M/np.sum(M))
    average=0
    for i in range(len(M)):
        average+=M[i]*i/np.sum(M)
    print ("True Sum: ",bsum)
    print ("Mult Sum: ",np.sum(M)/average)
    print ("True Mean: ", bsum/average)
    print ("Mult Sum: ", average)
    #print ("Average m: ",average)
    return M/livetime

def phototest(pnc_filepath,livetime,energies,efficiencies):
    #data1=np.loadtxt("../photoNC/Informacion_Sb124_high1.txt",skiprows=1,dtype="str")                                                                                                         
    data=np.loadtxt(pnc_filepath,skiprows=1,dtype="str")
    #data=data[0:2000]
    #print("Neutron output shape: ",np.shape(data))
    m=0
    M=np.zeros(20)
    #for i in range(500):                                                                                                                                                                                
    for i in range(len(data)):
        if i==0:
            m=0
        elif data[i,0]==data[i-1,0]:
            #print ("Same: ", i, data[i,0])                                                                                                                                         
            m=m
        else:
            if m==0:
                m=0
            else:
                M[m]+=1
                #print ("Recorded: ",m,i, M[m])                                                                                                                    
            m=0
        r=float(data[i,1])
        u=data[i,2]
        if u=="eV":
            r=r
        elif u=="keV":
            r=1000*r
        elif u=="MeV":
            r=1000000*r
        else:
            #print ("Help!",i+1,u)                                                                                                                                                                          
            fiftyfive=55
        t=1
        weight=1
        popper=np.random.uniform(0,1)
        blower=NucleationEfficiency(r,energies,efficiencies)
        #print ("blower: ",blower)                                                                                                                                                                          
        #print ("popper: ",popper, i, r)                                                                                                                                                                    
        if popper<blower:
            m+=1
            #print ("Tracker: ",i,m)        
        if i==len(data)-1:
            if m==0:
                m=0
            else:
                M[m]+=1
                #print ("Recorded: ",m,i, M[m]) 
    #print ("M :",M/np.sum(M))
    average=0
    for i in range(len(M)):
        average+=M[i]*i/np.sum(M)
    #print ("Average m: ",average)  
    return M[1:]/livetime

def phototest2(data,livetime,energies,efficiencies,photoBackMean,sourceStrength=1):
    m=0
    M=[]
    PZ=0
    lineO=[0,0]
    pz=0
    m=0
    #print (data)
    #print (data[:,0])
    #print (data[:,1])
    #datOrg=data[:,1]
    dataMem=np.zeros(np.shape(data))
    dataMem[:,0]=data[:,0]
    T=240
    sigLow=30
    sigUp=30
    dataMem[:,1]=VectorNucleationEfficiency(data[:,1],energies,efficiencies)
    #for i in range(len(data)):
        #r=data[i,1]
        #dataMem[i]=NucleationEfficiencyTrue(r,T,sigLow,sigUp)  
    for line in dataMem:
        if line[0] != lineO[0]:
            PZ+=pz
            pz=1
            m+=1
        pzi=(1-line[1])
        pz*=pzi
        lineO[:] = line[:]
    PZ+=pz
    m+=1
    backTot=0
    for i in range(len(photoBackMean)):
        backTot+=photoBackMean[i]*i
    datasum=sum(dataMem[:,1])
    backCount=sum(photoBackMean)
    #print ("############")
    #print (dataMem)
    #print (len(dataMem))
    #print(energies,efficiencies)
    #print("Data Sum: ",datasum)
    #print ("Back Sum: ",backCount)
    #print ("############")
    total=(datasum*sourceStrength)/livetime+backTot
    count=(m*sourceStrength-PZ)/livetime+backCount
    mean=total/count
    czero=PZ/livetime
    return total,count,mean,czero
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
def main(flist,alist,pnlist,aplist):
    #M=phototest("../Sb_fast.txt",500)
    #filename=sys.argv[1]
    
    Thomson=True
    Photoneutron=False
    Nuisance=True

    T=400
    sigLow=30
    sigUp=30
    binsize=.5
    sourceErr=.05
    background=500
    backErr=np.round(background**(1/2))
    #energies=[75,100,115,120,140]
    #efficiencies=[0,.2,.50,.8,1]
    #energies=[100,120]
    #efficiencies=[0,1]
    n = len(flist)
    RecoilList=[]
    WeightList=[]
    nuisanceT=np.zeros(n)
    TrueArray=np.zeros([3,n])
    InArray=np.zeros([4,n])
    UsePhotoError=[0]
    print (pnlist)
    pn=len(pnlist)
    pnrateList=[]
    pnrecoilList=[]
    pnweightList=[]
    time=100
    #photoBackMean=np.zeros(20)
    #photoBackMean=np.ones(20)*50
    #photoBackMean[1]=5
    photoBackMean=np.array([0,5,1.9,1.0,0.8,0.3,.2,.2,0,0,0,0,0,0,0,0,0,0,0,0])
    print ("photoBackMean length: ", len(photoBackMean))
    print ("!!!!!")
    pnlivetime=20
    print ("Photo-neutron livetime: ", pnlivetime)

    pn_sourceError=np.empty(pn)
    
    photoArrayTrue=np.zeros([pn,20])
    #photoArray0=np.zeros([2,20])
    #photoArray=np.zeros([pn,20])
    #photoArrayC=np.zeros([2,20])
    #photoArrayB=np.zeros([2,20])
    photoneutrondata=[]
    for p in range(pn):
        print (p,pn)
        print (photoBackMean)
        pnc_array=photoread(pnlist[p]+"_faster.txt")
        photoneutrondata+=[pnc_array]
        print ("")
        print ("")
        print ("********************")
        print ("Photoneutron Data p: ")
        print (photoneutrondata[p])
        print ("********************")
        print ("")
        print ("")
        photoArrayTrue[p,:]=phototestTrue(pnlist[p]+"_fast.txt",500,T,sigLow,sigUp)
        print("before",photoArrayTrue)
        pn_sourceError[p]=np.random.normal(1,sourceErr,1)
        photoArrayTrue[p,:]=photoJitter(photoArrayTrue[p,:],photoBackMean,pn_sourceError[p],time=time)
        #print ("M:",M)
        print("after",photoArrayTrue)
        print(photoBackMean)
        # comment spec1
        # spec11,spec21,spec22,spec31,spec32,spec33,specn1,specn2,specn3,specnn=load_photN2(pnlist[p])
        # B1,B2,B3,Bn=multiplicer2(spec11,spec21,spec22,spec31,spec32,spec33,specn1,specn2,specn3,specnn,T,sigLow,sigUp,1)
        # Mspec=np.array([B1,B2,B3,Bn])
        # print ("Multiplicity: ",Mspec)
        # Bsum=np.sum(Mspec)
        # print ("Mult Percent: ",Mspec/Bsum)
        # Bav=B1/Bsum+B2/Bsum*2+B3/Bsum*3+Bn/Bsum*4
        # print ("Mult Average: ", Bav)
        # rate=neutron_rates(spec11,T,sigLow,sigUp,1)
        # Count=rate*time
        # rate,count,background,nuisanceT[p]= rateJitter(rate,Count, time=time, sourceErr=sourceErr,background=background)
        # pnrateList+=[rate]
        # pnrecoilList+=[spec11[:,0]]
        # pnweightList+=[spec11[:,1]]
        # print (p,pnlist[p],rate)
    Mspecx=[1,2,3,4]
    Mx=np.arange(1,11,1)
    #n4plus=1-(M[1]+M[2]+M[3])/np.sum(M)
    #print ("M 4+: ",n4plus)
    #n4plus=(M[4]+M[5]+M[6]+M[7]+M[8]+M[9])/np.sum(M)
    #print ("M 4+: ",n4plus)
    #plt.plot(Mspecx,Mspec/np.sum(Mspec))
    #plt.plot(Mx,M[1:11]/np.sum(M))
    #plt.show()
    #plt.savefig("photomult_bi207.png")
    #plt.clf()
    for i in range(n):
        print (i,flist[i])
        Recoils, Weights, Rate, Count, t, time= analyze(flist[i],T,sigLow,sigUp,Activity=alist[i])
        Recoils,Weights=specrafy(Recoils,Weights,binsize=binsize)
        Rate2=rateFinderTrue(Recoils,T,sigLow,sigUp,t,Weights)
        print ("Recoils: ",Rate,Rate2,Rate-Rate2)
        Rate=Rate2
        Count=Rate*time
        TrueArray[0,i]=Rate
        TrueArray[1,i]=Count
        TrueArray[2,i]=t
        if i in UsePhotoError:
            rate ,count,background,nuisanceT[i]= rateJitter(Rate,Count, time=time, sourceErr=sourceErr,background=background,s_err=pn_sourceError[i])
        else:
            rate ,count,background,nuisanceT[i]=rateJitter(Rate,Count, time=time, sourceErr=sourceErr,background=background)
        print (nuisanceT[i])
        #rate=Rate
        #count=Count
        # InArray is jittered data
        InArray[0,i]=rate
        InArray[1,i]=count
        InArray[2,i]=t
        InArray[3,i]=background
        # original data
        RecoilList+=[Recoils]
        WeightList+=[Weights]

        print (TrueArray[:,i])
        print (InArray[:,i])
    print("recoillist", RecoilList)
    print("Weightlist", WeightList)
    # recoil and weight is modified spectrum
    np.savetxt(save_path+"start.txt",InArray)
    print ("THIS HERE")
    print("p",len(pnweightList),len(pnrecoilList),len(pnweightList))
    NitersRough = 25000
    Niters = 100000
    Nwalkers = 4
    printstep=5000    
    energies0=[0,0,0,0] #definition
    nuisance0=np.zeros(n+pn)
    #pnuisance0=np.zeros(pn)
    
    efficiencies0=[0,.3,.7,1]
    paracount=len(efficiencies0)
    energies0=np.zeros(paracount)
    buffer=5
    guessfloor=50
    guessceil=230
    bound=(guessceil-(paracount-1)*buffer-guessfloor)/(paracount+2)
    #print (paracount)
    print ("Random Threshold Step: ",bound)
    energies0=np.zeros(paracount)
    step=1
    n_step=.005
    RoughFact=2
    roughChiPenalty=0
    X=3
    
    #Chi0/T0/sigma0 is the initial guess
    #Chic/Tc/sigmac is the current basis, the iteration that the model will default back on if it does not select the new guess
    #Chii/Ti/sigmai is the most recent iteration
    #Chib/Tb/sigmab is the best fit so far, being stored to compare to future models
    OutArray=np.zeros([len(energies0)+4,Nwalkers])
    chi=np.zeros(n)
    chiNuis=np.zeros(n)
    chiNuisPN=np.zeros(pn)
    chipn=np.zeros(pn)
    for nw in range(Nwalkers):
        energies0=np.zeros(paracount)
        for j in range(paracount):
            if j==0:
                energies0[j]=guessfloor+3*bound*np.random.rand()
            else:
                energies0[j]=energies0[j-1]+bound*np.random.rand()+buffer
        #energies0=[180,210,270,300] #fixed intial guess
        print ("Initial Guess: ",energies0)
        np.savetxt(save_path+"guess_"+str(nw+1)+".txt",energies0)
        for i in range(n):
            # n is for flist index
            chi[i]=test(RecoilList[i],InArray[0,i],energies0,efficiencies0,nuisance0[i],InArray[2,i],WeightList[i],background=background,time=time)
        chipn=photoEval2(pn,photoneutrondata,pnlivetime,energies0,efficiencies0,photoArrayTrue,photoBackMean,loud=True,pnuisance=nuisance0[i+1:])
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
        print ("Chi Thomson: ",chi)
        print ("Chi Photo-Neutron: ",np.sum(chipn))
        print ("Chi photo: ",chipn)
        Chi0=0
        if Thomson==True:
            Chi0+=sum(chi)
        if Photoneutron==True:
            Chi0+=sum(chipn)
        if Nuisance==True:
            Chi0+=sum(chiNuis)
        #Chi0=sum(chi)+sum(chipn) #both
        #Chi0=sum(chi) #thomson only   
        #Chi0=sum(chipn) #photoneutron only
        Chib=Chi0
        Chic=Chi0
        MM=len(energies0)
        energiesi=np.empty(MM)
        energiesb=np.empty(MM)
        energiesc=np.empty(MM)
        nuisancei=np.zeros(n+pn)
        nuisanceb=np.zeros(n+pn)
        nuisancec=np.zeros(n+pn)
        #pnuisancei=np.zeros(pn)
        #pnuisanceb=np.zeros(pn)
        #pnuisancec=np.zeros(pn)
        for i in range(MM):
            value=energies0[i]
            energiesi[i]=value
            energiesb[i]=value
            energiesc[i]=value
        #print ("Defined")
        #print (energiesb)
        #print (energiesc)
        #print (Chi0,energies0)
        end=T+X*sigUp
        grade,fifty,RT,RG,TOT=fittest(T,sigLow,sigUp,energiesb,efficiencies0,start=0,end=end)
        print ("")
        print ("")
        print ("Guess Score: ",grade)
        print ("Total Value: ",TOT)
        print ("Guess Chi: ",Chi0)
        print (photoArrayTrue)
        #print (photoArray)
        #print ("Chi photo: ",chipn)
        print ("50% Efficiency Point: ",fifty,T)
        print ("Guess Parameters: ",energies0,efficiencies0)
        print ("")
        print ("")
        for ni in range(NitersRough):
            #print ("*****")
            #print ("Current Parameters: ", energiesc)
            #print ("*****")
            #energiesi[:],nuisancei[:]=stepper(energiesc,nuisancec,step*RoughFact,n_step)
            energiesi[:],nuisancei[:]=stepper(energiesc,nuisancec,step*RoughFact,n_step,max=2000)
            nuisancei[:]=nuisance0[:] #rough fit ignores nuisance parameters
            #pnuisancei[:]=pnuisance0[:]
            for i in range(n):
                chi[i]=test(RecoilList[i],InArray[0,i],energiesi,efficiencies0,nuisancei[i],InArray[2,i],WeightList[i],background=background,time=time)
                if chi[i]<0:
                    print("chi i Negative!: ",i)
                chiNuis[i]=(nuisancei[i]/sourceErr)**2
                if chiNuis[i]<0:
                    print("chi nuisance Negative!: ",i)
                    x=red
            if ni%printstep==0:
                pnverb=True
            else:
                pnverb=False
            chipn=photoEval2(pn,photoneutrondata,pnlivetime,energiesi,efficiencies0,photoArrayTrue,photoBackMean,loud=pnverb,pnuisance=nuisancei[i+1:])
            """
            for p in range(pn):
                #test(recoil,rate,energies,efficiencies,r_nuis,t,weight,background=500,time=100)
                #chipn[p]=test(pnrecoilList[p],pnrateList[p],energiesi,efficiencies0,pnuisancei[p],1,pnweightList[p],background=background,time=time)
                photoArray[p,:]=phototest(pnlist[p]+"_ultrafast.txt",10,energiesi,efficiencies0)
                for mb in range(len(photoArray)):
                    if photoArrayTrue[p,mb]>0:
                        chipn[p]+=((photoArray[p,mb]-photoArrayTrue[p,mb])**2/photoArrayTrue[p,mb])
            """
            Chii=0
            if Thomson==True:
                Chii+=sum(chi)
            if Photoneutron==True:
                Chii+=sum(chipn)
            if Nuisance==True:
                Chii+=sum(chiNuis)
            #Chii=sum(chi)+sum(chipn)+sum(chiNuis) #both
            #Chii=sum(chi)+roughChiPenalty #thomson only
            #Chii=sum(chipn)+sum(chiNuis) #photoneutron only
            if Chii<0:
                print ("Negative!")
                Chii=10**6
            if Chii==0:
                print ("Zero!")
                Chii==10**(-1)
            bar=Chii/(Chic+Chii)
            judge=np.random.uniform(0,1)
            if ni%printstep==0:
                print ("*****")
                print ("*****")
                print ("Iteration: ",ni)
                print ("Chi Thomson: ",sum(chi))
                print ("Chi photo-n: ",sum(chipn))
                print ("Chi photo: ",chipn)
                print ("Chi nuisance: ", sum(chiNuis))
                print ("Chi Iteration: ", Chii)
                print ("Current Chi: ", Chic)
                print ("Best Chi: ", Chib)
                #print (photoArrayTrue)
                #print (photoArray)
                print ("*****")
                print ("Iteration Parameters: ", energiesi)
                print ("Current Parameters: ", energiesc)
                print ("Best Parameters: ", energiesb)
                print ("*****")
                print ("Iteration Nuisance: ", nuisancei)
                print ("Current Nuisance: ", nuisancec)
                print ("Best Nuisance: ", nuisanceb)
                print ("Bar/Judge: ", bar,"/",judge)
                print ("*****")
                print ("*****")
            if Chii<Chib:
                Chib=Chii
                energiesb[:]=energiesi[:]
                energiesc[:]=energiesi[:]
                nuisanceb[:]=nuisancei[:]
                nuisancec[:]=nuisancei[:]
                #pnuisanceb[:]=pnuisancei[:]
                #pnuisancec[:]=pnuisancei[:]
            elif Chii<Chic:
                Chic=Chii
                energiesc[:]=energiesi[:]
                nuisancec[:]=nuisancei[:]
                #pnuisancec[:]=pnuisancei[:]
            elif judge>bar:
                Chic=Chii
                energiesc[:]=energiesi[:]
                nuisancec[:]=nuisancei[:]
                #pnuisancec[:]=pnuisancei[:]
            else:
                continue
        grade,fifty,RT,RG,TOT=fittest(T,sigLow,sigUp,energiesb,efficiencies0,start=0,end=end)
        print ("")
        print ("")
        print ("Rough Score: ",grade)
        print ("Total Value: ",TOT)
        print ("Rough Chii :",Chib)
        print ("Rough Parameters: ", energiesb)
        print ("50% Efficiency Point: ",fifty,T)
        print ("Model Effective Threshold: ", RG)
        print ("True Effective Threshold: ", RT)
        print ("")
        print ("")
        Chic=Chib
        energiesc[:]=energiesb[:]
        nuisancec[:]=nuisanceb[:]
        for ni in range(Niters):
            #stepi=step*(Niters-ni)/Niters
            stepi=step
            #energiesi[:],nuisancei[:]=stepper(energiesc,nuisancec,stepi,n_step)
            energiesi[:],nuisancei[:]=stepper(energiesc,nuisancec,step*RoughFact,n_step,max=2000)
            if Nuisance==False:
                nuisancei[:]=nuisance0[:] #turn nuisance parameters off for the whole fit  
            for i in range(n):
                chi[i]=test(RecoilList[i],InArray[0,i],energiesi,efficiencies0,nuisancei[i],InArray[2,i],WeightList[i],background=background,time=time)
                if chi[i]<0:
                    print("chi i Negative!: ",i)
                chiNuis[i]=(nuisancei[i]/sourceErr)**2
                if chiNuis[i]<0:
                    print("chi nuisance Negative!: ",i)
                    x=red
            if ni%printstep==0:
                pnverb=True
            else:
                pnverb=False
            chipn=photoEval2(pn,photoneutrondata,pnlivetime,energiesi,efficiencies0,photoArrayTrue,photoBackMean,loud=pnverb,pnuisance=nuisancei[i+1:])
            #for p in range(pn):
                #chiNuisPN[p]=(pnuisancei[p]/sourceErr)**2
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
            Chii=0
            if Thomson==True:
                Chii+=sum(chi)
            if Photoneutron==True:
                Chii+=sum(chipn)
            if Nuisance==True:
                Chii+=sum(chiNuis)
            #Chii=sum(chi)+sum(chiNuis)+sum(chipn) #both
            #Chii=sum(chi)+sum(chiNuis) #no nuisance both
            #Chii=sum(chi)+sum(chiNuis) #thomson only
            #Chii=sum(chiNuis)+sum(chipn) #photoneutron only
            if Chii<0:
                Chii=10**6
                print ("Negative!")
            if Chii==0:
                print ("Zero!")
                Chii=10**(-3)
            bar=Chii/(Chic+Chii)
            judge=np.random.uniform(0,1)
            if ni%printstep==0:                                                                                                                                                                            
                print ("*****")
                print ("*****")
                print ("Iteration: ",ni)
                print ("Chi Thomson: ",sum(chi))
                print ("Chi photo-n: ",sum(chipn))
                print ("Chi photo: ",chipn)
                print ("Chi nuisance: ", sum(chiNuis))
                print ("Chi Iteration: ", Chii)
                print ("Current Chi: ", Chic)
                print ("Best Chi: ", Chib)
                #print (photoArrayTrue)
                #print (photoArray)
                print ("*****")
                print ("Iteration Parameters: ", energiesi)
                print ("Current Parameters: ", energiesc)
                print ("Best Parameters: ", energiesb)
                print ("*****")
                print ("Iteration Nuisance: ", nuisancei)
                print ("Current Nuisance: ", nuisancec)
                print ("Best Nuisance: ", nuisanceb)
                print ("Bar/Judge: ", bar,"/",judge)
                print ("*****")
                print ("*****")                                                                                                                                                      
            if Chii<Chib:
                Chib=Chii
                energiesb[:]=energiesi[:]
                energiesc[:]=energiesi[:]
                nuisanceb[:]=nuisancei[:]
                nuisancec[:]=nuisancei[:]
            elif Chii<Chic:
                Chic=Chii
                energiesc[:]=energiesi[:]
                nuisancec[:]=nuisancei[:]
            elif judge>bar:
                Chic=Chii
                energiesc[:]=energiesi[:]
                nuisancec[:]=nuisancei[:]
            else:
                continue
        grade,fifty,RT,RG,TOT=fittest(T,sigLow,sigUp,energiesb,efficiencies0,start=0,end=end)
        teff=end-RT
        teffG=end-RG
        print ("Real Effective Threshold:  ",teff)
        print ("Model Effective Threshold: ",teffG)
        teffdiff=teffG/teff-1
        OutArray[0,nw]=grade
        OutArray[1,nw]=Chib
        OutArray[2,nw]=fifty/T-1
        OutArray[3,nw]=teffdiff
        OutArray[4:,nw]=energiesb[:]
        print ("")
        print ("")
        print ("Post Run Score: ",grade)
        print ("Total Value: ",TOT)
        print ("Final Chii :",Chib)
        print ("50% Efficiency Point: ",fifty,T)
        print ("Final Parameters: ", energiesb)
        print (nuisanceT)
        print (nuisanceb)
        for i in range(n):
            print (i,flist[i])
            print ("Source Strength Uncertainty")
            print ("True: ",(nuisanceT[i]-1)*100,"%")
            print ("Model: ",nuisanceb[i]*100,"%")
        print ("Total Nuisance Squared: ", sum(nuisanceT**2))
        #print ("Total Unfit Nuisance: ",sum((nuisanceT-nuisanceb)**2))
        print ("")
        print ("")
        #postfitplot(Tb,sigmab,zero,twenty,fifty,eighty,onehundred,start=zero-10,end=onehundred+10)
        np.savetxt(save_path+"fit.txt",OutArray)
    print (OutArray)
        
    
    
#flist = ["NewRuns/Eu152JAEA.txt","NewRuns/Bi207JAEA.txt","NewRuns/Y88JAEA.txt","NewRuns/Th228JAEA.txt"]
#flist = ["Eu152/JAEA.txt", "Bi207/JAEA.txt", "Sb124/JAEA.txt", "Th228/JAEA.txt", "Y88/JAEA.txt"]
#flist = ["Sb124/JAEA.txt","Y88/JAEA.txt","Th228/JAEA.txt"]
#"Eu152/JAEA.txt", "Bi207/JAEA.txt", "Sb124/JAEA.txt", "Th228/JAEA.txt", "Y88/JAEA.txt"
#flist = ["Bi207/JAEA.txt","Sb124/JAEA.txt","Eu152/JAEA.txt","Y88/JAEA.txt"]
#flist=["Bi207/JAEA.txt", "Sb124/JAEA.txt", "Th228/JAEA.txt", "Y88/JAEA.txt"]
#flist=["Bi207/JAEA.txt","Th228/JAEA.txt"]     
#flist=["../Th228/JAEA.txt","../Bi207/JAEA.txt"]

#phototest()

main(flist,alist,pnlist,aplist)
