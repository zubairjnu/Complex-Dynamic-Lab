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
gammaP=5.4*10**(-1);        deltaP=1.4*10**(-1);        OMx = 9.963*10**(-8)

# Parameters for stress input
T=3600*6;                   A=2.50;                     lemda=0.05/3600.0
DeltaG=0.375*BetaG;         Karf=1000.0

SSS=A

def osci(Var, t):
    p53A, RNAn, RNAc, MDM2c, MDM2n, ARF, p53M, G = Var
    SS = A #+ A*math.sin(2*math.pi*t/T)
    #SS = A*np.exp(-lemda*t)
    F1 = kp - k1*p53A*MDM2n - dp*p53A - gammap53A*p53A*p53M - Deltap53A*((G**h)/((KK**h)+(G**h)))*p53A
    F2 = km + k2*(p53A**1.8)/((kD**1.8)+(p53A**1.8)) - k0*RNAn
    F3 = k0*RNAn - drc*RNAc
    F4 = kT*RNAc - ki*MDM2c
    F5 = ki*MDM2c - dmm*MDM2c**2 - k3*MDM2n*ARF
    F6 = ka - da*ARF - k3*MDM2n*ARF + Delta*(G**n)/((Kg**n)+(G**n))*ARF
    F7 = Alphap53M + Deltap53A*((G**h)/((KK**h)+(G**h)))*p53A - Gammap53M*p53M*MDM2n - Deltap53M*p53M
    F8 = AlphaG + BetaG*(SS**n)/((Ks**n)+(SS**n)) - GammaG*G + DeltaG*(p53M**n)/((Karf**n)+(p53M**n))
    return [F1, F2, F3, F4, F5, F6, F7, F8]

#Label = ['p53A', 'RNAn', 'RNAc', 'MDM2c', 'MDM2n', 'ARF', 'p53M', 'G']
#Initial = [5.0, 2.0, 3.0, 3.0, 4.0, 1.0, 0.0, 0.0]
Initial = [ra.uniform(0,10), ra.uniform(0,10), ra.uniform(0,10), ra.uniform(0,10), ra.uniform(0,10), ra.uniform(0,10), 0.0, 0.0]
t = np.linspace(0, 6*10**5, 6*10**7)
Sol = odeint(osci, Initial, t)

Time=[]; D1=[]; D2=[]; D3=[]; D4=[]; D5=[]; D6=[]; D7=[]; D8=[]; D9=[]; D10=[]
for index, time, d1, d5, d6, d7, d8 in zip(range(len(t)), t, Sol[:,0], Sol[:,4], Sol[:,5], Sol[:,6], Sol[:,7]):
    if index%100==0: Time.append(time); D1.append(d1); D5.append(d5); D6.append(d6); D7.append(d7); D8.append(d8)
TIME = np.array(Time)/3600.0


plt.figure(figsize=(4,2.0))
#plt.subplots_adjust(left=0.115, right=0.990, top=0.970, bottom=0.125, hspace = 0.0)
grid = plt.GridSpec(10, 1, left=0.115, right=0.990, top=0.970, bottom=0.125, hspace = 0.0, wspace=0.0)

LW=0.5;
#Tau = [A + A*math.sin(2*math.pi*tt/T) for tt in Time]
Tau = [A*np.exp(-lemda*tt) for tt in Time]

# Ploting
plt.subplot(grid[0:2, 0])
plt.plot(TIME,Tau,'-r',linewidth=LW)
plt.legend(['Stress'], loc='best', prop={'size': 6})
plt.yticks([SSS,SSS+11.5]); plt.ylim([0-SSS*0.05,SSS*1.05])
plt.xlim([0-max(TIME)*0.0125,max(TIME)*1.0125])

plt.subplot(grid[2:, 0])
plt.plot(TIME,np.array(D1),'-b',TIME,np.array(D7),'-r',linewidth=LW)
plt.legend([r'p53$_{A}$','p53$_{M}$','GENE'], loc='best', prop={'size': 6},facecolor=None)
plt.xlim([0-max(TIME)*0.0125,max(TIME)*1.0125])
#plt.savefig('p53TwoDecayS{}.pdf'.format(int(SSS*10)))
#plt.savefig('p53F.pdf')
plt.show()

end_time = datetime.now()
print 'MAIN PROGRAM IS COMPLETED||Duration||H:M:S||{}'.format(end_time - start_time), '\n'
###############################################################################################################
###############################################################################################################
