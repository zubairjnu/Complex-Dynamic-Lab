#	
#		
#				 
#			
##########################################################################################

#from mpl_toolkits.mplot3d import Axes3D
#import pylab
#import sympy as sp
#from operator import itemgetter
#from scipy import linalg
#from pylab import genfromtxt
#import fit
#import networkx as nx
#from mayavi import mlab
#from multiprocessing import Pool, TimeoutError, Process, Queue, Pipe, Value, Array
import os, math, pdb, time, random as ra, numpy as np, matplotlib.pyplot as plt
from datetime import datetime
#import seaborn as sns; sns.set()
from scipy.integrate import odeint
#from scipy import signal
########################################################################################################
start_time = datetime.now()
print start_time, 'STARTING TIME OF THE PROGRAM'
print '####', '\n'

print 'D1 = P53A' 
print 'D2 = RNAn'
print 'D3 = RNAc'
print 'D4 = MDM2c'
print 'D5 = MDM2n'
print 'D6 = ARF'
print 'D7 = p53M'
print 'D8 = Oncogene', '\n'


def SystemC(Stress, KValue):
    # Old model parameters
    kp=0.50;                    k1=9.963*10**(-6);          dp=1.925*10**(-5)
    km=1.5*10**(-3);            k2=1.5*10**(-2);            kD=740.0;            k0=8*10**(-4)
    drc=1.444*10**(-4);         kT=1.66*10**(-2);           ki=9.0*10**(-4)
    dmm=1.66*10**(-7) ;         ka=0.5
    da=3.209*10**(-5) ;         k3=9.963*10**(-6) 

    n=3.0;                      Ks=3.0;                     Kg=10.0;
    Delta=3.5*10**(-6);         KK=50;                      h=4

    # New model parameters 
    gammap53A=9.963*10**(-7);   Deltap53A=5.963*10**(-4)
    Alphap53M=1.5*10**(-7);     Gammap53M=9.963*10**(-6);   Deltap53M=1.925*10**(-5)
    AlphaG=1.5*10**(-7);        BetaG=6.5*10**(-3);         GammaG=1.925*10**(-5); 

    alphaM=5.4*10**(-1);        KM=50.0; x=2.0;             betaM=8.4*10**(-3)
    gammaP=5.4*10**(-1);        deltaP=1.4*10**(-1)

    # Parameters for stress input
    T=3600*6;                   A=Stress;                      lemda=0.05/3600.0
    DeltaG=0.375*BetaG;         Karf=KValue

    def osci(Var, t):
        p53A, RNAn, RNAc, MDM2c, MDM2n, ARF, p53M, G = Var
        #SS = A + A*math.sin(2*math.pi*t/T)
        SS = A#*np.exp(-lemda*t)
        F1 = kp - k1*p53A*MDM2n - dp*p53A - gammap53A*p53A*p53M - Deltap53A*(G**h)/((KK**h)+(G**h))*p53A
        F2 = km + k2*(p53A**1.8)/((kD**1.8)+(p53A**1.8)) - k0*RNAn
        F3 = k0*RNAn - drc*RNAc
        F4 = kT*RNAc - ki*MDM2c
        F5 = ki*MDM2c - dmm*MDM2c**2 - k3*MDM2n*ARF
        F6 = ka - da*ARF - k3*MDM2n*ARF + Delta*(G**n)/((Kg**n)+(G**n))*ARF
        F7 = Alphap53M + Deltap53A*(G**h)/((KK**h)+(G**h))*p53A - Gammap53M*p53M*MDM2n - Deltap53M*p53M
        F8 = AlphaG + BetaG*(SS**n)/((Ks**n)+(SS**n)) - GammaG*G + DeltaG*(p53M**n)/((Karf**n)+(p53M**n))
        return [F1, F2, F3, F4, F5, F6, F7, F8]

    #Label = ['p53A', 'RNAn', 'RNAc', 'MDM2c', 'MDM2n', 'ARF', 'p53M', 'G']
    Initial = [5.0, 2.0, 3.0, 3.0, 4.0, 1.0, 0.0, 0.0]
    t = np.linspace(0, 6*10**5, 12*10**6)
    Sol = odeint(osci, Initial, t)

    Time=[]; D1=[]; D2=[]; D3=[]; D4=[]; D5=[]; D6=[]; D7=[]; D8=[]; D9=[]; D10=[]
    for index, time, d1, d5, d6, d7, d8 in zip(range(len(t)), t, Sol[:,0], Sol[:,4], Sol[:,5], Sol[:,6], Sol[:,7]):
        if index%100==0: Time.append(time); D1.append(d1); D5.append(d5); D6.append(d6); D7.append(d7); D8.append(d8)
    TIME = np.array(Time)/3600.0

    return TIME, D1, D7

def SystemS(Stress, KValue):
    # Old model parameters
    kp=0.50;                    k1=9.963*10**(-6);          dp=1.925*10**(-5)
    km=1.5*10**(-3);            k2=1.5*10**(-2);            kD=740.0;            k0=8*10**(-4)
    drc=1.444*10**(-4);         kT=1.66*10**(-2);           ki=9.0*10**(-4)
    dmm=1.66*10**(-7) ;         ka=0.5
    da=3.209*10**(-5) ;         k3=9.963*10**(-6)
    
    n=3.0;                      Ks=3.0;                     Kg=10.0;
    Delta=3.5*10**(-6);         KK=50;                      h=4
    
    # New model parameters
    gammap53A=9.963*10**(-7);   Deltap53A=5.963*10**(-4)
    Alphap53M=1.5*10**(-7);     Gammap53M=9.963*10**(-6);   Deltap53M=1.925*10**(-5)
    AlphaG=1.5*10**(-7);        BetaG=6.5*10**(-3);         GammaG=1.925*10**(-5);
    
    alphaM=5.4*10**(-1);        KM=50.0; x=2.0;             betaM=8.4*10**(-3)
    gammaP=5.4*10**(-1);        deltaP=1.4*10**(-1)
    
    # Parameters for stress input
    T=3600*6;                   A=Stress;                      lemda=0.05/3600.0
    DeltaG=0.375*BetaG;         Karf=KValue
    
    def osci(Var, t):
        p53A, RNAn, RNAc, MDM2c, MDM2n, ARF, p53M, G = Var
        SS = A + A*math.sin(2*math.pi*t/T)
        #SS = A#*np.exp(-lemda*t)
        F1 = kp - k1*p53A*MDM2n - dp*p53A - gammap53A*p53A*p53M - Deltap53A*(G**h)/((KK**h)+(G**h))*p53A
        F2 = km + k2*(p53A**1.8)/((kD**1.8)+(p53A**1.8)) - k0*RNAn
        F3 = k0*RNAn - drc*RNAc
        F4 = kT*RNAc - ki*MDM2c
        F5 = ki*MDM2c - dmm*MDM2c**2 - k3*MDM2n*ARF
        F6 = ka - da*ARF - k3*MDM2n*ARF + Delta*(G**n)/((Kg**n)+(G**n))*ARF
        F7 = Alphap53M + Deltap53A*(G**h)/((KK**h)+(G**h))*p53A - Gammap53M*p53M*MDM2n - Deltap53M*p53M
        F8 = AlphaG + BetaG*(SS**n)/((Ks**n)+(SS**n)) - GammaG*G + DeltaG*(p53M**n)/((Karf**n)+(p53M**n))
        return [F1, F2, F3, F4, F5, F6, F7, F8]
    
    #Label = ['p53A', 'RNAn', 'RNAc', 'MDM2c', 'MDM2n', 'ARF', 'p53M', 'G']
    Initial = [5.0, 2.0, 3.0, 3.0, 4.0, 1.0, 0.0, 0.0]
    t = np.linspace(0, 6*10**5, 12*10**6)
    Sol = odeint(osci, Initial, t)
    
    Time=[]; D1=[]; D2=[]; D3=[]; D4=[]; D5=[]; D6=[]; D7=[]; D8=[]; D9=[]; D10=[]
    for index, time, d1, d5, d6, d7, d8 in zip(range(len(t)), t, Sol[:,0], Sol[:,4], Sol[:,5], Sol[:,6], Sol[:,7]):
        if index%100==0: Time.append(time); D1.append(d1); D5.append(d5); D6.append(d6); D7.append(d7); D8.append(d8)
    TIME = np.array(Time)/3600.0
    
    return TIME, D1, D7

def SystemD(Stress, KValue):
    # Old model parameters
    kp=0.50;                    k1=9.963*10**(-6);          dp=1.925*10**(-5)
    km=1.5*10**(-3);            k2=1.5*10**(-2);            kD=740.0;            k0=8*10**(-4)
    drc=1.444*10**(-4);         kT=1.66*10**(-2);           ki=9.0*10**(-4)
    dmm=1.66*10**(-7) ;         ka=0.5
    da=3.209*10**(-5) ;         k3=9.963*10**(-6)
    
    n=3.0;                      Ks=3.0;                     Kg=10.0;
    Delta=3.5*10**(-6);         KK=50;                      h=4
    
    # New model parameters
    gammap53A=9.963*10**(-7);   Deltap53A=5.963*10**(-4)
    Alphap53M=1.5*10**(-7);     Gammap53M=9.963*10**(-6);   Deltap53M=1.925*10**(-5)
    AlphaG=1.5*10**(-7);        BetaG=6.5*10**(-3);         GammaG=1.925*10**(-5);
    
    alphaM=5.4*10**(-1);        KM=50.0; x=2.0;             betaM=8.4*10**(-3)
    gammaP=5.4*10**(-1);        deltaP=1.4*10**(-1)
    
    # Parameters for stress input
    T=3600*6;                   A=Stress;                      lemda=0.05/3600.0
    DeltaG=0.375*BetaG;         Karf=KValue
    
    def osci(Var, t):
        p53A, RNAn, RNAc, MDM2c, MDM2n, ARF, p53M, G = Var
        #SS = A + A*math.sin(2*math.pi*t/T)
        SS = A*np.exp(-lemda*t)
        F1 = kp - k1*p53A*MDM2n - dp*p53A - gammap53A*p53A*p53M - Deltap53A*(G**h)/((KK**h)+(G**h))*p53A
        F2 = km + k2*(p53A**1.8)/((kD**1.8)+(p53A**1.8)) - k0*RNAn
        F3 = k0*RNAn - drc*RNAc
        F4 = kT*RNAc - ki*MDM2c
        F5 = ki*MDM2c - dmm*MDM2c**2 - k3*MDM2n*ARF
        F6 = ka - da*ARF - k3*MDM2n*ARF + Delta*(G**n)/((Kg**n)+(G**n))*ARF
        F7 = Alphap53M + Deltap53A*(G**h)/((KK**h)+(G**h))*p53A - Gammap53M*p53M*MDM2n - Deltap53M*p53M
        F8 = AlphaG + BetaG*(SS**n)/((Ks**n)+(SS**n)) - GammaG*G + DeltaG*(p53M**n)/((Karf**n)+(p53M**n))
        return [F1, F2, F3, F4, F5, F6, F7, F8]
    
    #Label = ['p53A', 'RNAn', 'RNAc', 'MDM2c', 'MDM2n', 'ARF', 'p53M', 'G']
    Initial = [5.0, 2.0, 3.0, 3.0, 4.0, 1.0, 0.0, 0.0]
    t = np.linspace(0, 6*10**5, 12*10**6)
    Sol = odeint(osci, Initial, t)
    
    Time=[]; D1=[]; D2=[]; D3=[]; D4=[]; D5=[]; D6=[]; D7=[]; D8=[]; D9=[]; D10=[]
    for index, time, d1, d5, d6, d7, d8 in zip(range(len(t)), t, Sol[:,0], Sol[:,4], Sol[:,5], Sol[:,6], Sol[:,7]):
        if index%100==0: Time.append(time); D1.append(d1); D5.append(d5); D6.append(d6); D7.append(d7); D8.append(d8)
    TIME = np.array(Time)/3600.0
    
    return TIME, D1, D7

A = 1.0; T=3600*6; lemda=0.05/3600.0; Time=[time for index, time in zip(range(12*10**6), np.linspace(0, 6*10**5, 12*10**6)) if index%100==0]
def Stress(t): return [A, A + A*math.sin(2*math.pi*t/T), A*np.exp(-lemda*t)]


plt.figure(figsize=(17,14))
grid = plt.GridSpec(680, 500, left=0.025, right=0.99, top=0.98, bottom=0.048, hspace = 0.0, wspace=0.0)


lw=0.75; fs=18; Legfs=11; fs1=18; fs2=13
XtickS=[0,40,80,120,160]; YtickS=[0,600,1200,1800]
XliM=[-4.25,170]; YliM=[-50,2000]; Te0=0; Te1=1750


# Plot 1
###############
# Data ploting
Data1=SystemC(0.5,500.0); Data2=SystemC(1.6, 500.0)
Data3=SystemC(1.27, 100.0); Data4=SystemC(2.0, 200.0)

plt.subplot(grid[0:90, 300:395])
plt.plot(Data1[0],Data1[1],'-g' ,Data1[0],Data1[2],'-r', linewidth=lw)
plt.xlim(XliM); plt.ylim(YliM); plt.xticks([]); plt.yticks(YtickS, fontsize=fs-4); plt.text(Te0,Te1,r'  A$_1$', fontsize=fs)
plt.legend([r'p53$_{A}$','p53$_{M}$'], loc='best', prop={'size': Legfs},facecolor=None); plt.ylabel('Concentration', fontsize=fs1)

plt.subplot(grid[0:90, 405:500])
plt.plot(Data2[0],Data2[1],'-g' ,Data2[0],Data2[2],'-r', linewidth=lw)
plt.xlim(XliM); plt.ylim(YliM); plt.xticks([]); plt.yticks([]); plt.text(Te0,Te1,r'  A$_2$', fontsize=fs)
plt.legend([r'p53$_{A}$','p53$_{M}$'], loc='best', prop={'size': Legfs},facecolor=None)

plt.subplot(grid[110:200, 300:395])
plt.plot(Data3[0],Data3[1],'-g' ,Data3[0],Data3[2],'-r', linewidth=lw)
plt.xlim(XliM); plt.ylim(YliM); plt.xticks(XtickS, fontsize=fs-4); plt.yticks(YtickS, fontsize=fs-4); plt.ylabel('Concentration', fontsize=fs1)
plt.text(Te0,Te1,r'  A$_3$', fontsize=fs); plt.legend([r'p53$_{A}$','p53$_{M}$'], loc='best', prop={'size': Legfs},facecolor=None)
plt.xlabel('Time (hrs)', fontsize=fs)

plt.subplot(grid[110:200, 405:500])
plt.plot(Data4[0],Data4[1],'-g' ,Data4[0],Data4[2],'-r', linewidth=lw)
plt.xlim(XliM); plt.ylim(YliM); plt.xticks(XtickS, fontsize=fs-4); plt.yticks([]); plt.xlabel('Time (hrs)', fontsize=fs)
plt.text(Te0,Te1,r'  A$_4$', fontsize=fs); plt.legend([r'p53$_{A}$','p53$_{M}$'], loc='best', prop={'size': Legfs},facecolor=None)


# Plot 2
###############

# Data ploting
Data11=SystemS(0.5,500.0); Data21=SystemS(1.3, 500.0)
Data31=SystemS(0.96, 100.0); Data41=SystemS(2.0, 200.0)

plt.subplot(grid[240:330, 300:395])
plt.plot(Data11[0],Data11[1],'-g' ,Data11[0],Data11[2],'-r', linewidth=lw)
plt.xlim(XliM); plt.ylim(YliM); plt.xticks([]); plt.yticks(YtickS, fontsize=fs-4); plt.text(Te0,Te1,r'  B$_1$', fontsize=fs)
plt.legend([r'p53$_{A}$','p53$_{M}$'], loc='best', prop={'size': Legfs},facecolor=None); plt.ylabel('Concentration', fontsize=fs1)

plt.subplot(grid[240:330, 405:500])
plt.plot(Data21[0],Data21[1],'-g' ,Data21[0],Data21[2],'-r', linewidth=lw)
plt.xlim(XliM); plt.ylim(YliM); plt.xticks([]); plt.yticks([]); plt.text(Te0,Te1,r'  B$_2$', fontsize=fs)
plt.legend([r'p53$_{A}$','p53$_{M}$'], loc='best', prop={'size': Legfs},facecolor=None)

plt.subplot(grid[350:440, 300:395])
plt.plot(Data31[0],Data31[1],'-g' ,Data31[0],Data31[2],'-r', linewidth=lw)
plt.xlim(XliM); plt.ylim(YliM); plt.xticks(XtickS, fontsize=fs-4); plt.yticks(YtickS, fontsize=fs-4); plt.ylabel('Concentration', fontsize=fs1)
plt.text(Te0,Te1,r'  B$_3$', fontsize=fs); plt.legend([r'p53$_{A}$','p53$_{M}$'], loc='best', prop={'size': Legfs},facecolor=None)
plt.xlabel('Time (hrs)', fontsize=fs)

plt.subplot(grid[350:440, 405:500])
plt.plot(Data41[0],Data41[1],'-g' ,Data41[0],Data41[2],'-r', linewidth=lw)
plt.xlim(XliM); plt.ylim(YliM); plt.xticks(XtickS, fontsize=fs-4); plt.yticks([]); plt.xlabel('Time (hrs)', fontsize=fs)
plt.text(Te0,Te1,r'  B$_4$', fontsize=fs); plt.legend([r'p53$_{A}$','p53$_{M}$'], loc='best', prop={'size': Legfs},facecolor=None)


# Plot 3
###############
# Data ploting
Data12=SystemD(2.0,500.0); Data22=SystemD(4.0, 500.0)
Data32=SystemD(2.0, 100.0); Data42=SystemD(4.0, 100.0)

plt.subplot(grid[480:570, 300:395])
plt.plot(Data12[0],Data12[1],'-g' ,Data12[0],Data12[2],'-r', linewidth=lw)
plt.xlim(XliM); plt.ylim(YliM); plt.xticks([]); plt.yticks(YtickS, fontsize=fs-4); plt.text(Te0,Te1,r'  C$_1$', fontsize=fs)
plt.legend([r'p53$_{A}$','p53$_{M}$'], loc='best', prop={'size': Legfs},facecolor=None); plt.ylabel('Concentration', fontsize=fs1)

plt.subplot(grid[480:570, 405:500])
plt.plot(Data22[0],Data22[1],'-g' ,Data22[0],Data22[2],'-r', linewidth=lw)
plt.xlim(XliM); plt.ylim(YliM); plt.xticks([]); plt.yticks([]); plt.text(Te0,Te1,r'  C$_2$', fontsize=fs)
plt.legend([r'p53$_{A}$','p53$_{M}$'], loc='best', prop={'size': Legfs},facecolor=None)

plt.subplot(grid[590:680, 300:395])
plt.plot(Data32[0],Data32[1],'-g' ,Data32[0],Data32[2],'-r', linewidth=lw)
plt.xlim(XliM); plt.ylim(YliM); plt.xticks(XtickS, fontsize=fs-4); plt.yticks(YtickS, fontsize=fs-4); plt.ylabel('Concentration', fontsize=fs1)
plt.text(Te0,Te1,r'  C$_3$', fontsize=fs); plt.legend([r'p53$_{A}$','p53$_{M}$'], loc='best', prop={'size': Legfs},facecolor=None)
plt.xlabel('Time (hrs)', fontsize=fs)

plt.subplot(grid[590:680, 405:500])
plt.plot(Data42[0],Data42[1],'-g' ,Data42[0],Data42[2],'-r', linewidth=lw)
plt.xlim(XliM); plt.ylim(YliM); plt.xticks(XtickS, fontsize=fs-4); plt.yticks([]); plt.xlabel('Time (hrs)', fontsize=fs)
plt.text(Te0,Te1,r'  C$_4$', fontsize=fs); plt.legend([r'p53$_{A}$','p53$_{M}$'], loc='best', prop={'size': Legfs},facecolor=None)


lw=0.75; fs=18; Legfs=11; fs1=18; fs2=18
XtickS=[0,40,80,120,160]; YtickS=[0,500,1000,1500,2000]
XliM=[-4.25,170]; YliM=[-50,2000]; Te0=0; Te1=1750


# Stress ploting
###################
TIME = np.array(Time)/3600.0
plt.subplot(grid[12:62, 0:80])
plt.plot(TIME, [Stress(Ti)[0] for Ti in Time],'-b', linewidth=lw)
plt.xlim([0-max(TIME)*0.0125,max(TIME)*1.0125]); plt.ylim([-0.1, 1.1])
plt.yticks([]); plt.ylabel('Stress', fontsize=fs)
plt.text(0,1.25,'a', fontsize=fs); plt.xticks([0,40,80,120,160], fontsize=fs-4)
plt.title('Constant stress', fontsize=fs2)


plt.subplot(grid[240:290, 0:80])
plt.plot(TIME, [Stress(Ti)[1] for Ti in Time],'-b', linewidth=lw)
plt.xlim([0-max(TIME)*0.0125,max(TIME)*1.0125]); plt.ylim([-0.2, 2.2])
plt.yticks([]); plt.ylabel('Stress', fontsize=fs)
plt.text(0,2.50,'b', fontsize=fs); plt.xticks([0,40,80,120,160], fontsize=fs-4)
plt.title('Oscillatory stress', fontsize=fs2)


plt.subplot(grid[480:530, 0:80])
plt.plot(TIME, [Stress(Ti)[2] for Ti in Time],'-b', linewidth=lw)
plt.xlim([0-max(TIME)*0.0125,max(TIME)*1.0125]); plt.ylim([-0.1, 1.1])
plt.yticks([]); plt.xlabel('Time (hrs)', fontsize=fs); plt.ylabel('Stress', fontsize=fs)
plt.text(0,1.25,'c', fontsize=fs); plt.xticks([0,40,80,120,160], fontsize=fs-4)
plt.title('Decaying stress', fontsize=fs2)


# Heat plot
##################
import seaborn as sns; sns.set()
Xdat=np.linspace(0.5,2.0,200)
Ydat=np.linspace(100.0,500,200)[::-1]
Index=[0,66,132,199]
Xticks=[]; Yticks=[]
for i in range(len(Xdat)):
    if i in Index: Xticks.append(round(Xdat[i],3)); Yticks.append(round(Ydat[i],3))
    else: Xticks.append(''); Yticks.append('')

Data=np.genfromtxt('./data/Avgp53aC.dat')
plt.subplot(grid[0:200, 110:260])
ax=sns.heatmap(Data[:][:],linewidth=0.0,yticklabels=Yticks,xticklabels=Xticks,rasterized=True, cmap='RdYlGn')
cbar = ax.collections[0].colorbar
cbar.set_label(r'p53$_{A}$ concentration', labelpad=5, fontsize=fs1)
ax.set_xticks(list(np.linspace(3,200,200)))
plt.xticks(fontsize=fs-3); plt.yticks(fontsize=fs-3); plt.gcf().axes[-1].tick_params(labelsize=15)
plt.title('A', fontsize=fs, loc='left')
plt.xlabel('Magnitude of stress ($I$)', fontsize=fs);
plt.ylabel(r'$K_{3}$', fontsize=fs)


Xdat1=np.linspace(0.5,2.0,200)
Ydat1=np.linspace(100.0,500,200)[::-1]
Index1=[0,66,132,199]
Xticks1=[]; Yticks1=[]
for i in range(len(Xdat1)):
    if i in Index1: Xticks1.append(round(Xdat1[i],3)); Yticks1.append(round(Ydat1[i],3))
    else: Xticks1.append(''); Yticks1.append('')

Data1=np.genfromtxt('./data/Avgp53aS.dat')
plt.subplot(grid[240:440, 110:260])
ax=sns.heatmap(Data1[:][:],linewidth=0.0,yticklabels=Yticks1,xticklabels=Xticks1,rasterized=True, cmap='RdYlGn')
cbar = ax.collections[0].colorbar
cbar.set_label(r'p53$_{A}$ concentration', labelpad=5, fontsize=fs1)
ax.set_xticks(list(np.linspace(3,200,200)))
plt.xticks(fontsize=fs-3); plt.yticks(fontsize=fs-3); plt.gcf().axes[-1].tick_params(labelsize=15)
plt.title('B', fontsize=fs, loc='left')
plt.xlabel('Magnitude of stress ($I$)', fontsize=fs);
plt.ylabel(r'$K_{3}$', fontsize=fs)

Xdat2=np.linspace(2.0,4.0,200)
Ydat2=np.linspace(100.0,500,200)[::-1]
Index2=[0,66,132,199]
Xticks2=[]; Yticks2=[]
for i in range(len(Xdat2)):
    if i in Index2: Xticks2.append(round(Xdat2[i],3)); Yticks2.append(round(Ydat2[i],3))
    else: Xticks2.append(''); Yticks2.append('')

Data2=np.genfromtxt('./data/Avgp53aD.dat')
plt.subplot(grid[480:680, 110:260])
ax=sns.heatmap(Data2[:][:],linewidth=0.0,yticklabels=Yticks2,xticklabels=Xticks2,rasterized=True, cmap='RdYlGn')
cbar = ax.collections[0].colorbar
cbar.set_label(r'p53$_{A}$ concentration', labelpad=5, fontsize=fs1)
ax.set_xticks(list(np.linspace(3,200,200)))
plt.xticks(fontsize=fs-3); plt.yticks(fontsize=fs-3); plt.gcf().axes[-1].tick_params(labelsize=15)
plt.title('C', fontsize=fs, loc='left')
plt.xlabel(r'Magnitude of stress ($I$)', fontsize=fs); plt.ylabel(r'$K_{3}$', fontsize=fs)


plt.savefig('fig5.eps')
plt.savefig('fig5.pdf')
#plt.show()

end_time = datetime.now()
print 'MAIN PROGRAM IS COMPLETED||Duration||H:M:S||{}'.format(end_time - start_time), '\n'
###############################################################################################################
###############################################################################################################
