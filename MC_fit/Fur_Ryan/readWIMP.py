import numpy as np
import sys
import math
import matplotlib.pyplot as plt
import scipy.stats

source=sys.argv[1]
m=0
n=50

try:
    n=int(sys.argv[4])
    print ("n=",n)
    m=int(sys.argv[5])
    print ("m=",m)
except:
    print ("Default values for n and m")

#n=1

print (m,n)
fname="fit.txt"

def fittest(T,sigLow,sigUp,energiesM,efficienciesM,M=10000,start=40,end=200,scorerange=90):
    rt=0
    rg=0
    dif=0
    tot=0
    calib=0
    fifty=0
    fiftydif=.5
    #print (T,sigLow,sigUp)                                                                                                                                                                                                                   
    print ("Energies and Efficiencies: ",energiesM,efficienciesM)                                                                                                                                                                                                          
    for i in range(M):
        stepwidth=(end-start)/M
        r=i*stepwidth+start
        RT=NucleationEfficiencyTrue(r,T,sigLow,sigUp)
        #print ("Real Efficiency: ",RT)                                                                                                                                                                                                       
        RG=NucleationEfficiency(r,energiesM,efficiencies=efficienciesM)
        #print ("Model Efficiency: ",RG)                                                                                                                                                                                                      
        #print ("Difference: ",RT-RG)
        #if r>T-scorerange and r<T+scorerange:
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

def fit_CEvNS(T,sigLow,sigUp,energiesM,efficienciesM,cevns,M=10000,start=50,end=1370,scorerange=90):
    rt=0
    rg=0
    dif=0
    cdif=0
    tot=0
    calib=0
    fifty=0
    fiftydif=.5
    #print (T,sigLow,sigUp)
    ci=0
    cr=cevns[ci,0]
    cwr=cevns[ci,1]
    ci_max=len(cevns[:,0])-1
    ctot=0
    print ("Energies and Efficiencies: ",energiesM,efficienciesM)
    for i in range(M):
        stepwidth=(end-start)/M
        cw=cwr*stepwidth
        r=i*stepwidth+start
        while r>cr:
            ci+=1
            cr=cevns[ci,0]
        if ci==0:
            cwr=cevns[ci,1]
            #print (cwr)
            cw=cwr*stepwidth
            #print (cw)
        elif r>cevns[ci-1,0]:
            x0=cevns[ci-1,0]
            x1=cevns[ci,0]
            #print ("x0: ",x0,"x1: ",x1)
            y0=cevns[ci-1,1]
            y1=cevns[ci,1]
            #print ("y0: ", y0,"y1: ",y1)
            slope = (y1-y0)/(x1-x0)
            #print ("Slope: ",slope)
            x = r-x0
            #print ("x: ",x)
            cwr = slope * x +y0
            #print ("slope * x + y0: ",cwr)
            cw = cwr*stepwidth
            #print (cw)
        else:
            pizza=cow
            cr0=cr
            ci0=ci
            while r<cr0:
                ci0 -= ci0
                cr0  = cevns[ci0,0]
            x0=cevns[ci0,0]
            x1=cevns[ci,0]
            y0=cevns[ci0,1]
            y1=cevns[ci,1]
            slope = (y1-y0)/(x1-x0)
            x =	r-x0
            cwr = slope	* x + y0
            cw = cwr*stepwidth
        """
        if r>cr:
            if ci==ci_max:
                ci=len(cevns[:,1])-1
                cr=cevns[ci,0]
                cwr=cevns[ci,1]
                cw=cwr*stepwidth
            else:
                if cevns[ci,0]>r>cevns[ci+1,0]:
                    x0=cevns[ci,0]
                    y0=cevns[ci,1]
                    x1=cevns[ci+1,0]
                    y1=cevns[ci+1,1]
                    slope=(y1-y0)/(x1-x0)
                    x=r-x0
                    cwr=slope*x+x0
                    cw=cwr*stepwidth
                else:
                    while r<cevns[ci+1,0] and ci<=ci_max:
                        print ("While", ci)
                        ci+=1
                    if ci==ci_max:
                        ci=len(cevns[:,1])-1
                        cr=cevns[ci,0]
                        cwr=cevns[ci,1]
                        cw=cwr*stepwidth
                    else:
                        x0=cevns[ci,0]
                        y0=cevns[ci,1]
                        x1=cevns[ci+1,0]
                        y1=cevns[ci+1,1]
                        slope=(y1-y0) / (x1-x0)
                        x=r-x0
                        cwr=slope*x+x0
                        cr=cwr*stepwidth
        """
        RT=NucleationEfficiencyTrue(r,T,sigLow,sigUp)
        """
        print ("*********************************")
        print ("Recoil Energy: ", r)
        print ("Recoil CEvNS: ", cr)
        print ("CEvNS Rate: ",cevns[ci,1])
        print ("CEvNS Raw: ",cwr)
        print ("Step Size (eV): ",stepwidth)
        print ("Recoil Weight: ",cw)
        print ("Real Efficiency: ",RT)
        print ("Bubble Weight: ",cw*RT)
        print ("Bubble Check: ",cevns[ci,1]*stepwidth*RT)
        """
        RG=NucleationEfficiency(r,energiesM,efficiencies=efficienciesM)
        #print ("Model Efficiency: ",RG)
        #print ("Difference: ",RT-RG)
        #if r>T-scorerange and r<T+scorerange:
        ctot+=RT*cw
        cdif+=(RT-RG)*cw
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
    return dif,cdif,fifty,rt,rg,tot,ctot

def NucleationEfficiencyTrue(r,T,sigLow,sigUp):
    if r<T:
        R=1/2*(1+math.erf((r-T)/(sigLow*2**(1/2))))
    else:
        R=1/2*(1+math.erf((r-T)/(sigUp*2**(1/2))))
    return R

def NucleationEfficiency(r,energies,efficiencies="auto"):
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

#cevns=np.loadtxt("dahlspectrum_10GeV_Ar.txt")
cevns=np.loadtxt("WIMP_recoil_Lamb_2GeV.txt")  

cevns[:,0]*=1000
cevns[:,1]/=1000

print("cevns: ", cevns)

f1=source+"/"+str(m+1)+"/"+fname
lines = np.loadtxt(f1, dtype="str",comments="#")
walkers=len(lines[0,:])
params=len(lines)-4
models=len(range(m,n))

print("Models: ",models)

score=np.zeros([models,walkers])
score2=np.zeros([models,walkers])
chi=np.zeros([models,walkers])
fifty=np.zeros([models,walkers])
teff=np.zeros([models,walkers])

tt=np.zeros([models,walkers])
tm=np.zeros([models,walkers])

T=float(sys.argv[2])
sigLow=float(sys.argv[3])
sigUp=sigLow

fiftyreal=80

parameters=np.zeros([models,params,walkers])

scuffmap=[]

blankmapi=[]
blankmapj=[]
for i in range(m,n):
    filename=source+"/"+str(i+1)+"/"+fname
    i=i-m
    print (filename)
    lines = np.loadtxt(filename, dtype="str",comments="#")
    #print (np.shape(lines))
    #print (len(lines))
    #print (len(lines[0,:]))
    #print (len(lines[:,0]))
    score[i,:]=lines[0,:]
    chi[i,:]=lines[1,:]
    fifty[i,:]=lines[2,:]
    teff[i,:]=lines[3,:]
    parameters[i,:,:]=lines[4:,:]
    paracheck=10
    paracheckHarsh=100
    if paracheck<0:
        print ("Run scuffed!")
        scuffmap+=[i]
        print ("Mean Parameter: ",paracheck)
    elif paracheckHarsh<0:
        print ("Run may be scuffed")
        print ("Mean Parameter: ",paracheck)
        print ("Lowest Quartile: ",paracheckHarsh)
    efficienciesM=[0,.2,.5,.8,1]
    X=10
    #end=T+X*sigUp
    #start=T-X*sigUp
    #if start<0:
        #start=0
    #start=0
    #if len(efficienciesM)==len(lines[5:,0]):
    start=50
    end=1370
    for j in range(walkers):
        grade,cgrade,fiftyV,RT,RG,TOT,CTOT=fit_CEvNS(T,sigLow,sigUp,parameters[i,:,j],efficienciesM,cevns,M=10000,scorerange=sigUp*10)
        score[i,j]=cgrade
        fifty[i,j]=(fiftyV-fiftyreal)/fiftyreal*100
        #print (RG,RT)
        tm[i,j]=RG
        tt[i,j]=RT
        adjustedT=end-RT
        adjustedG=end-RG
        teff[i,j]=(adjustedG-adjustedT)/adjustedT*100
    #else:
        #print ("Skipped the score recalculation")
#spread=end-start

blotcount=len(scuffmap)
scoreclean=np.zeros([models-blotcount,walkers])
chiclean=np.zeros([models-blotcount,walkers])
fiftyclean=np.zeros([models-blotcount,walkers])
teffclean=np.zeros([models-blotcount,walkers])

parametersclean=np.zeros([models-blotcount,params,walkers])

if blotcount>0:
    index=0
    for i in range(models-blotcount):
        if index in scuffmap:
            index+=1
        scoreclean[i,:]=score[index,:]
        chiclean[i,:]=chi[index,:]
        fiftyclean[i,:]=fifty[index,:]
        teffclean[i,:]=teff[index,:]
        parametersclean[i,:,:]=parameters[index,:,:]
        index+=1


scoremean=np.mean(score,1)
chimean=np.mean(chi,1)
fiftymean=np.mean(fifty,1)
teffmean=np.mean(teff,1)
teffmedian=np.median(teff,1)
paramean=np.mean(parameters,2)

#print (np.shape(paramean))

#print (np.shape(chi))
#print(chi)
arg_min=np.argmin(chi,axis=1)
#print ("Shape Argmin: ",np.shape(arg_min))
chi_min=chi[range(len(arg_min)),arg_min]
score_min=score[range(len(arg_min)),arg_min]
fifty_min=fifty[range(len(arg_min)),arg_min]
teff_min=teff[range(len(arg_min)),arg_min]
para_min=parameters[range(len(arg_min)),:,arg_min]
print ("Shape Para Min: ",np.shape(para_min))

score_med=np.median(score)

teff_mean=np.mean(teff)
teff_median=np.median(teff)
teff_std=np.std(teff)
teff_mad=scipy.stats.median_abs_deviation(teff)
teff_min_mean=np.mean(teff_min)
teff_min_median=np.median(teff_min)
teff_min_std=np.std(teff_min)
teff_min_mad=scipy.stats.median_abs_deviation(teff_min)

print ("Mean Effective Threshold (all): ", teff_mean, "Median Effective Threshold (all): ", teff_median)
print ("Standard Deviation (all): ", teff_std, "Median Absolute Deviation (all): ", teff_mad)
print ("Mean Effective Threshold (best): ", teff_min_mean, "Median Effective Threshold (best): ", teff_min_median)
print ("Standard Deviation (best): ", teff_min_std, "Median Absolute Deviation (best): ", teff_min_mad)

tmmean=np.mean(np.mean(tm,1))
ttmean=np.mean(np.mean(tt,1))

TMmean=end-tmmean
TTmean=end-ttmean

#print (ttmean/spread,1-(TTmean-start)/spread)

#print (tmmean/spread,1-(TMmean-start)/spread)

mean4=np.median(para_min,0)
#mean4=np.median(paramean,0)

scoremeanclean=np.mean(scoreclean,1)
chimeanclean=np.mean(chiclean,1)
fiftymeanclean=np.mean(fiftyclean,1)
teffmeanclean=np.mean(teffclean,1)
parameanclean=np.mean(parametersclean,2)

x=.16
efficiencies=[0,0.2,0.5,.8,1]

print (np.shape(parameters))

paralow=np.zeros(params)
parahigh=np.zeros(params)
for i in range(params):
    paralow[i]=np.quantile(parameters[:,i,:],x)
    parahigh[i]=np.quantile(parameters[:,i,:],1-x)
#print ("paralowi",paralowi)
#print (np.shape(paralow))
#paralow=np.quantile(parameters,x,0)
#parahigh=np.quantile(paramean,1-x,0)
print (paralow)
print (parahigh)
paralow2=np.quantile(paramean,.05,0)
parahigh2=np.quantile(paramean,.95,0)

paralowC=np.quantile(paramean,x,0)
parahighC=np.quantile(paramean,1-x,0)

paralow2C=np.quantile(paramean,.05,0)
parahigh2C=np.quantile(paramean,.95,0)

#print (np.shape(paralow))

score_median=np.quantile(score,.5)
score_median_min=np.quantile(score_min,.5)
score_median2=np.median(score)
print ("Median Test: ", score_median - np.median(score))
score_low=np.quantile(score,x)
score_high=np.quantile(score,1-x)

score_min_median=np.median(score_min)
score_mad=scipy.stats.median_abs_deviation(score_min)

teff_min_median=np.median(teff_min)
teff_mad=scipy.stats.median_abs_deviation(teff_min)

print (score_median,score_high,score_low)
print ( "score min quant .5", score_median_min, "score np.median", score_median2, "score min np.median", score_min_median, "old score median", score_med)
print ("score min quant .5", score_median_min/CTOT*100, "score np.median", score_median2/CTOT*100, "score min np.median", score_min_median/CTOT*100, "old score median", score_med/CTOT*100)

print ("WIMP median: ", score_median, " + ", score_high - score_median," - ",score_median - score_low)
print ("WIMP median/MAD: ", score_min_median, " +/- ", score_mad)
print ("WIMP median/MAD: ", score_min_median, " +/- ", scipy.stats.median_abs_deviation(score_min))
print ("Effective Threshold median/MAD: ", teff_min_median, " +/- ", teff_mad)
print ("Scuffmap: ",scuffmap)
print ("blotcount: ",blotcount)

print ("Efficiencies :",efficiencies)
print ("Parameters 2 Sigma High: ",parahigh2)
print ("Parameters 2 Sigma Low: ",paralow2)
print ("Parameters 1 Sigma High: ",parahigh)
print ("Parameters 1 Sigma Low: ",paralow)
print ("Mean Fit: ", mean4)
if blotcount>0:
    print ("Clean Parameters:")
    print ("Parameters 2 Sigma High: ",parahigh2C)
    print ("Parameters 2 Sigma Low: ",paralow2C)
    print ("Parameters 1 Sigma High: ",parahighC)
    print ("Parameters 1 Sigma Low: ",paralowC)

#print (np.shape(paralow))
#print (np.mean(parameters))
#print (np.mean(parameters,0))
#print (np.mean(parameters,1))
#print (np.mean(parameters,2))
M=10000

R=np.zeros(M)
RM=np.zeros(M)
RT=np.zeros(M)
RH=np.zeros(M)
RL=np.zeros(M)
Rmean4=np.zeros(M)

RH2=np.zeros(M)
RL2=np.zeros(M)

start=T-sigLow*10
end=T+10*sigUp

start=paralow[0]-20
end=parahigh[4]+20

if start<0:
    start = 0
print (np.shape(mean4))
print (np.shape(parahigh))

for i in range(M):
    r=i/M*end
    R[i]=r
    RT[i]=NucleationEfficiencyTrue(r,T,sigLow,sigUp)
    RM[i]=NucleationEfficiency(r,mean4,efficiencies)
    RH[i]=NucleationEfficiency(r,parahigh,efficiencies)
    RL[i]=NucleationEfficiency(r,paralow,efficiencies)
    RH2[i]=NucleationEfficiency(r,parahigh2,efficiencies)
    RL2[i]=NucleationEfficiency(r,paralow2,efficiencies)

#dif,cdif,fifty,rt,rg,tot,ctot
_,score_mean4,_,_,_,_,_ = fit_CEvNS(T,sigLow,sigUp,mean4,efficienciesM,cevns,M=10000,scorerange=sigUp*10)

RTadj=end-np.sum(RT)
RMadj=end-np.sum(RM)
RHadj=end-np.sum(RH)
RLadj=end-np.sum(RL)
RH2adj=end-np.sum(RH2)
RL2adj=end-np.sum(RL2)
Rmean4adj=end-np.sum(Rmean4)

TeffM=(RMadj-RTadj)/RTadj*100
TeffH=(RHadj-RTadj)/RTadj*100
TeffL=(RLadj-RTadj)/RTadj*100

#print ("Effective T: ", TeffM, "   1 sigma high: ", TeffL, "   1 sigma low: ", TeffH)

print (score_median , score_high ,score_low, CTOT )
print ("WIMP median: ", score_median/CTOT*100, " + ", (score_high - score_median)/CTOT*100," - ",(score_median - score_low)/CTOT*100)
print ("WIMP median: ", score_median/CTOT*100, " +/- ", ( (score_high - score_median)+(score_median - score_low) )/CTOT*50)

print ("CTOT: ", CTOT, 1/CTOT*100)
print ("WIMP median/MAD: ", score_min_median, " +/- ", score_mad)
print ("Effective Threshold median/MAD: ", teff_min_median, " +/- ", teff_mad)
print ("WIMP median/MAD: ", score_min_median/CTOT*100, " +/- ", score_mad/CTOT*100)
print ("Effective Threshold median/MAD: ", T+teff_min_median*T/100, " +/- ", teff_mad*T/100)

print ("WIMP plot: ", score_mean4/CTOT*100, " +/- ", ( (score_high - score_mean4)+(score_mean4 - score_low) )/CTOT*50)

print ("Mean Score: ",np.mean(score_min),np.std(score))
print ("Total Rate: ",CTOT,"events/kg/day")
print ("WIMP Systematic (lowest chi): ",np.mean(score_min)/CTOT*100,np.std(score)/CTOT*100)
print ("Mean Chi (lowest chi): ",np.mean(chi_min),np.std(chi_min))
print ("WIMP Systematic (all): ",np.mean(score)/CTOT*100,np.std(score)/CTOT*100)
print ("Mean Chi (all): ",np.mean(chi),np.std(chi))
print ("50% Crossing Point:",np.mean(fifty_min),np.std(fifty))
print ("Effective Threshold (%):",np.mean(teff_min),np.std(teff_min))
print ("Effective Threshold (eV): ",(np.mean(teff_min))*T/100+T,np.std(teff_min)*T/100)
print ("Lowest Chi Effective Threshold (%):",np.mean(teff),np.std(teff))
print ("All fits Effective Threshold (eV): ",(np.mean(teff))*T/100+T,np.std(teff)*T/100)

#effectiveT=115.98381738323732
effectiveT=80
effectiveG=80*(1+np.mean(teffmean)/100)
axis_font = {'size':'14'}
print (effectiveG,effectiveT)
#print (TMmean, TTmean)

plt.plot(R,RT,label="True Efficiency")
plt.plot(R,RM,color="red",label="Model Efficiency")
plt.plot(mean4,efficiencies,linestyle="",marker="o",color="red")
#plt.plot(cevns[:,0],cevns[:,1]/cevns[0,1],linestyle=":",marker="",color="black",label="WIMP Spectrum")
#plt.fill_between(R,RL,RH,alpha=.15,color="red",label="1-sigma error")
plt.xlim([50,T*2])
plt.ylabel("Nucleation Efficiency", **axis_font)
plt.xlabel("Energy (eV)", **axis_font)
plt.legend(fontsize=12)
plt.tick_params(labelsize=14)
plt.savefig("Sample_Efficiency.png")
plt.clf()

plt.plot(R,RT)
plt.fill_between(R,RL,RH,alpha=.4,color="red")
plt.fill_between(R,RL2,RH2,alpha=.2,color="red")
#plt.plot(R,RM)
plt.xlim([start,end])
plt.savefig(source+".png")
plt.clf()

print ("saved: ",source+".png")

#axis_font = {'size':'14'}

print ("Mean figure start/end: ",start,end)

plt.plot(R,RT,label="Monte Carlo Truth")
plt.plot(R,RM,color="red",label="Avg. Best-fit Model")
plt.plot(mean4,efficiencies,linestyle="",marker="o",color="red")
plt.plot(cevns[:,0],cevns[:,1]/cevns[0,1],linestyle=":",marker="",color="black",label="WIMP Spectrum")
plt.fill_between(R,RL,RH,alpha=.15,color="red",label="1-sigma error")
#plt.fill_between(R,RL2,RH2,alpha=.2,color="red") 
#plt.xlim([50,end])
plt.xlim([50,T*2]) 
plt.ylabel("Nucleation Efficiency", **axis_font)
plt.xlabel("Energy (eV)", **axis_font)
plt.legend(fontsize=12)
plt.tick_params(labelsize=14)
plt.savefig(source+"_WIMP.png")
plt.clf()

#plt.plot(R,RT,label="Monte Carlo Truth")
#plt.plot(R,RM,color="red",label="Avg. Best-fit Model")
#plt.plot(mean4,efficiencies,linestyle="",marker="o",color="red")
plt.plot(cevns[:,0]/1000,cevns[:,1],linestyle=":",marker="",color="green",label="WIMP Spectrum")
#plt.fill_between(R,RL,RH,alpha=.15,color="red",label="1-sigma error")
#plt.fill_between(R,RL2,RH2,alpha=.2,color="red")                                                                                                                                                                                             
#plt.xlim([start,end])
plt.ylabel("WIMP Rate (events/day/kg/keV)", **axis_font)
plt.xlabel("Energy (keV)", **axis_font)
plt.legend(fontsize=12)
plt.tick_params(labelsize=14)
plt.savefig("WIMP.png")
plt.clf()

print ("saved: ",source+"_WIMP.png")
if blotcount>0:
    for i in range(M):
        r=i/M*end
        R[i]=r
        RT[i]=NucleationEfficiencyTrue(r,T,sigLow,sigUp)
        RH[i]=NucleationEfficiency(r,parahighC,efficiencies=[0,.3,.7,1])
        RL[i]=NucleationEfficiency(r,paralowC,efficiencies=[0,.3,.7,1])
        RH2[i]=NucleationEfficiency(r,parahigh2C,efficiencies=[0,.3,.7,1])
        RL2[i]=NucleationEfficiency(r,paralow2C,efficiencies=[0,.3,.7,1])

    print ("Clean Score: ",np.mean(scoremeanclean),np.std(scoremeanclean))
    print ("Clean Chi: ",np.mean(chimeanclean),np.std(chimeanclean))
    print ("Clean 50% Crossing Point :",np.mean(fiftymeanclean),np.std(fiftymeanclean))
    print ("Clean Effective Threshold :",np.mean(teffmeanclean),np.std(teffmeanclean))

    effectiveT=80
    effectiveG=80*(np.mean(teffmeanclean))
    print (effectiveG,effectiveT)

    plt.plot(R,RT)
    plt.fill_between(R,RL,RH,alpha=.4,color="red")
    plt.fill_between(R,RL2,RH2,alpha=.2,color="red")
    plt.xlim([start,end])
    plt.savefig(source+"clean.png")
    plt.clf()
