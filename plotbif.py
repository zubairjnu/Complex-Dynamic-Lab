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

Data=[np.genfromtxt('./data/Bif{}.dat'.format(Ind)) for Ind in range(1,10)]

A = 1.0; T=3600*6; lemda=0.05/3600.0; Time=[time for index, time in zip(range(12*10**6), np.linspace(0, 6*10**5, 12*10**6)) if index%100==0]
def Stress(t): return [A, A + A*math.sin(2*math.pi*t/T), A*np.exp(-lemda*t)]

plt.figure(figsize=(16,7.5))
grid = plt.GridSpec(300, 400, left=0.023, right=0.990, top=0.985, bottom=0.080, hspace = 0.0, wspace=0.0)

lw=1.25; fs=17; Legfs=11; fs1=17; ColStrength=0.25
XtickS=[0.5,1.0,1.5,2.0,2.5]; YtickS=[0,600,1200,1800]
XliM=[0.5,2.5]; YliM=[-50,2000]; Te0=1.05; Te1=1775

# Data ploting
plt.subplot(grid[0:87, 95:190])
DATAX=Data[0]
Ind1=0; Ind2=0; Ind3=0
for INDEX, data in zip(range(len(DATAX[:,1])),DATAX[:,1]):
    if data<(518*0.95): 
        Ind1=INDEX        
        print INDEX, DATAX[:,0][INDEX]; break

for INDEX, dataA,dataM in zip(range(len(DATAX[:,1])),DATAX[:,1],DATAX[:,2]):
    if dataM>dataA: 
        Ind2=INDEX            
        print INDEX, DATAX[:,0][INDEX]; break

for INDEX, dataA,dataM in zip(range(len(DATAX[:,1])),DATAX[:,1],DATAX[:,2]):
    if dataM>max(DATAX[:,1]):
        Ind3=INDEX  
        print INDEX, DATAX[:,0][INDEX], 'Ind3'; break

plt.plot(Data[0][:,0],Data[0][:,1],'-g' ,Data[0][:,0],Data[0][:,2],'-r', linewidth=lw)
plt.fill_between(Data[0][:Ind1,0], [-50 for AA in range(len(Data[0][:Ind1,0]))], 2000, facecolor='y', alpha=ColStrength)
plt.fill_between(Data[0][Ind1-1:Ind2,0], [-50 for AA in range(len(Data[0][Ind1-1:Ind2,0]))], 2000, facecolor='c', alpha=ColStrength)
plt.fill_between(Data[0][Ind2-1:Ind3,0], [-50 for AA in range(len(Data[0][Ind2-1:Ind3,0]))], 2000, facecolor='m', alpha=ColStrength)
plt.fill_between(Data[0][Ind3-1:,0], [-50 for AA in range(len(Data[0][Ind3-1:,0]))], 2000, facecolor='grey', alpha=ColStrength)
plt.xlim(XliM); plt.ylim(YliM); plt.xticks(XtickS, fontsize=fs1-4); plt.yticks(YtickS, fontsize=fs1-4); plt.ylabel('Concentration', fontsize=fs1)
plt.legend([r'p53$_{A}$','p53$_{M}$'], loc='best', prop={'size': Legfs},facecolor=None); plt.text(Te0,Te1,r'  A$_1$', fontsize=fs)
plt.text(0.525,1000,'Active state', fontsize=Legfs)
plt.text(1.2475,1000,'Apoptotic state', fontsize=Legfs)
plt.text(2.075,1800,'Pre-malignant state', fontsize=Legfs,rotation=90)
plt.text(2.4,1500,'Cancer state', fontsize=Legfs,rotation=90)

DATAX=Data[1]
Ind1=0; Ind2=0; Ind3=0
for INDEX, data in zip(range(len(DATAX[:,1])),DATAX[:,1]):
    if data<(518*0.95): 
        Ind1=INDEX        
        print INDEX, DATAX[:,0][INDEX]; break

for INDEX, dataA,dataM in zip(range(len(DATAX[:,1])),DATAX[:,1],DATAX[:,2]):
    if dataM>dataA: 
        Ind2=INDEX            
        print INDEX, DATAX[:,0][INDEX]; break

for INDEX, dataA,dataM in zip(range(len(DATAX[:,1])),DATAX[:,1],DATAX[:,2]):
    if dataM>max(DATAX[:,1]):
        Ind3=INDEX  
        print INDEX, DATAX[:,0][INDEX], 'Ind3'; break


plt.subplot(grid[0:87, 200:295])
plt.plot(Data[1][:,0],Data[1][:,1],'-g' ,Data[1][:,0],Data[1][:,2],'-r', linewidth=lw)
plt.fill_between(Data[1][:Ind1,0], [-50 for AA in range(len(Data[1][:Ind1,0]))], 2000, facecolor='y', alpha=ColStrength)
plt.fill_between(Data[1][Ind1-1:Ind2,0], [-50 for AA in range(len(Data[1][Ind1-1:Ind2,0]))], 2000, facecolor='c', alpha=ColStrength)
plt.fill_between(Data[1][Ind2-1:Ind3,0], [-50 for AA in range(len(Data[1][Ind2-1:Ind3,0]))], 2000, facecolor='m', alpha=ColStrength)
plt.fill_between(Data[1][Ind3-1:,0], [-50 for AA in range(len(Data[1][Ind3-1:,0]))], 2000, facecolor='grey', alpha=ColStrength)
plt.xlim(XliM); plt.ylim(YliM); plt.xticks(XtickS, fontsize=fs1-4); plt.yticks([]); plt.text(Te0,Te1,r'  A$_2$', fontsize=fs)
plt.legend([r'p53$_{A}$','p53$_{M}$'], loc='best', prop={'size': Legfs},facecolor=None)
plt.text(0.525,1000,'Active state', fontsize=Legfs)
plt.text(1.5,1600,'Apoptotic state', fontsize=Legfs,rotation=90)
plt.text(1.78,1800,'Pre-malignant state', fontsize=Legfs,rotation=90)
plt.text(1.9,1000,'Cancer state', fontsize=Legfs,rotation=0)

DATAX=Data[2]
Ind2=0
for INDEX, dataA,dataM in zip(range(len(DATAX[:,1])),DATAX[:,1],DATAX[:,2]):
    if dataM>dataA: 
        Ind2=INDEX            
        print INDEX, DATAX[:,0][INDEX]; break

plt.subplot(grid[0:87, 305:400])
plt.plot(Data[2][:,0],Data[2][:,1],'-g' ,Data[2][:,0],Data[2][:,2],'-r', linewidth=lw)
plt.fill_between(Data[2][:Ind2,0], [-50 for AA in range(len(Data[2][:Ind2,0]))], 2000, facecolor='y', alpha=ColStrength)
plt.fill_between(Data[2][Ind2-1:,0], [-50 for AA in range(len(Data[2][Ind2-1:,0]))], 2000, facecolor='grey', alpha=ColStrength)
plt.xlim(XliM); plt.ylim(YliM); plt.xticks(XtickS, fontsize=fs1-4); plt.yticks([]); plt.text(0.5,Te1,r'  A$_3$', fontsize=fs)
plt.legend([r'p53$_{A}$','p53$_{M}$'], loc='center right', prop={'size': Legfs},facecolor=None)
plt.text(0.525,1000,'Active state', fontsize=Legfs)
plt.text(1.325,1000,'Cancer state', fontsize=Legfs)

DATAX=Data[3]
Ind1=0; Ind2=0; Ind3=0
for INDEX, data in zip(range(len(DATAX[:,1])),DATAX[:,1]):
    if data<(518*0.95): 
        Ind1=INDEX        
        print INDEX, DATAX[:,0][INDEX]; break

for INDEX, dataA,dataM in zip(range(len(DATAX[:,1])),DATAX[:,1],DATAX[:,2]):
    if dataM>dataA: 
        Ind2=INDEX            
        print INDEX, DATAX[:,0][INDEX]; break

for INDEX, dataA,dataM in zip(range(len(DATAX[:,1])),DATAX[:,1],DATAX[:,2]):
    if dataM>max(DATAX[:,1]):
        Ind3=INDEX  
        print INDEX, DATAX[:,0][INDEX], 'Ind3'; break

plt.subplot(grid[107:194, 95:190])
plt.plot(Data[3][:,0],Data[3][:,1],'-g' ,Data[3][:,0],Data[3][:,2],'-r', linewidth=lw)
plt.fill_between(Data[3][:Ind1,0], [-50 for AA in range(len(Data[3][:Ind1,0]))], 2000, facecolor='y', alpha=ColStrength)
plt.fill_between(Data[3][Ind1-1:Ind2,0], [-50 for AA in range(len(Data[3][Ind1-1:Ind2,0]))], 2000, facecolor='c', alpha=ColStrength)
plt.fill_between(Data[3][Ind2-1:Ind3,0], [-50 for AA in range(len(Data[3][Ind2-1:Ind3,0]))], 2000, facecolor='m', alpha=ColStrength)
plt.fill_between(Data[3][Ind3-1:,0], [-50 for AA in range(len(Data[3][Ind3-1:,0]))], 2000, facecolor='grey', alpha=ColStrength)
plt.xlim(XliM); plt.ylim(YliM); plt.xticks(XtickS, fontsize=fs1-4); plt.yticks(YtickS, fontsize=fs1-4); plt.ylabel('Concentration', fontsize=fs1)
plt.text(Te0,Te1,r'  B$_1$', fontsize=fs); plt.legend([r'p53$_{A}$','p53$_{M}$'], loc='best', prop={'size': Legfs},facecolor=None)
plt.text(0.525,1000,'Active state', fontsize=Legfs)
plt.text(0.93,800,'Apoptotic state', fontsize=Legfs)
plt.text(1.8,1800,'Pre-malignant state', fontsize=Legfs,rotation=90)
plt.text(2.25,1500,'Cancer state', fontsize=Legfs,rotation=90)

DATAX=Data[4]
Ind1=0; Ind2=0; Ind3=0
for INDEX, data in zip(range(len(DATAX[:,1])),DATAX[:,1]):
    if data<(518*0.95): 
        Ind1=INDEX        
        print INDEX, DATAX[:,0][INDEX]; break

for INDEX, dataA,dataM in zip(range(len(DATAX[:,1])),DATAX[:,1],DATAX[:,2]):
    if dataM>dataA: 
        Ind2=INDEX            
        print INDEX, DATAX[:,0][INDEX]; break

for INDEX, dataA,dataM in zip(range(len(DATAX[:,1])),DATAX[:,1],DATAX[:,2]):
    if dataM>max(DATAX[:,1]):
        Ind3=INDEX  
        print INDEX, DATAX[:,0][INDEX], 'Ind3'; break

plt.subplot(grid[107:194, 200:295])
plt.plot(Data[4][:,0],Data[4][:,1],'-g' ,Data[4][:,0],Data[4][:,2],'-r', linewidth=lw)
plt.fill_between(Data[4][:Ind1,0], [-50 for AA in range(len(Data[4][:Ind1,0]))], 2000, facecolor='y', alpha=ColStrength)
plt.fill_between(Data[4][Ind1-1:Ind2,0], [-50 for AA in range(len(Data[4][Ind1-1:Ind2,0]))], 2000, facecolor='c', alpha=ColStrength)
plt.fill_between(Data[4][Ind2-1:Ind3,0], [-50 for AA in range(len(Data[4][Ind2-1:Ind3,0]))], 2000, facecolor='m', alpha=ColStrength)
plt.fill_between(Data[4][Ind3-1:,0], [-50 for AA in range(len(Data[4][Ind3-1:,0]))], 2000, facecolor='grey', alpha=ColStrength)
plt.xlim(XliM); plt.ylim(YliM); plt.xticks(XtickS, fontsize=fs1-4); plt.yticks([]); plt.text(Te0,Te1,r'  B$_2$', fontsize=fs)
plt.legend([r'p53$_{A}$','p53$_{M}$'], loc='best', prop={'size': Legfs},facecolor=None)
plt.text(0.525,1000,'Active state', fontsize=Legfs)
plt.text(1.25,1600,'Apoptotic state', fontsize=Legfs,rotation=90)
plt.text(1.45,1800,'Pre-malignant state', fontsize=Legfs,rotation=90)
plt.text(1.7,1000,'Cancer state', fontsize=Legfs,rotation=0)


DATAX=Data[5]
Ind2=0
for INDEX, dataA,dataM in zip(range(len(DATAX[:,1])),DATAX[:,1],DATAX[:,2]):
    if dataM>dataA: 
        Ind2=INDEX            
        print INDEX, DATAX[:,0][INDEX]; break

plt.subplot(grid[107:194, 305:400])
plt.plot(Data[5][:,0],Data[5][:,1],'-g' ,Data[5][:,0],Data[5][:,2],'-r', linewidth=lw)
plt.fill_between(Data[5][:Ind2,0], [-50 for AA in range(len(Data[5][:Ind2,0]))], 2000, facecolor='y', alpha=ColStrength)
plt.fill_between(Data[5][Ind2-1:,0], [-50 for AA in range(len(Data[5][Ind2-1:,0]))], 2000, facecolor='grey', alpha=ColStrength)
plt.xlim(XliM); plt.ylim(YliM); plt.xticks(XtickS, fontsize=fs1-4); plt.yticks([]); plt.text(0.5,Te1,r'  B$_3$', fontsize=fs)
plt.legend([r'p53$_{A}$','p53$_{M}$'], loc='center right', prop={'size': Legfs},facecolor=None)
plt.text(.525,1000,'Active state', fontsize=Legfs)
plt.text(1.25,1000,'Cancer state', fontsize=Legfs)

DATAX=Data[6]
Ind1=93; Ind2=0
for INDEX,dataA1,dataA2,dataM1,dataM2 in zip(range(len(DATAX[:,1])),DATAX[:,1],DATAX[:,2],DATAX[:,3],DATAX[:,4]):
    #print dataM, dataA
    if (dataM1+dataM2)/2.0>(dataA1+dataA2)/2.0: 
        Ind2=INDEX            
        print INDEX, DATAX[:,0][INDEX]; break

plt.subplot(grid[214:300, 95:190])
plt.plot(Data[6][:,0],(Data[6][:,1]+Data[6][:,2])/2.0,'-g' ,Data[6][:,0],(Data[6][:,3]+Data[6][:,4])/2.0,'-r', linewidth=lw)
plt.plot(Data[6][Ind1:Ind2,0],Data[6][Ind1:Ind2,5],'-b' ,Data[6][Ind1:Ind2,0],Data[6][Ind1:Ind2,6],'-k', linewidth=lw)
plt.fill_between(Data[6][:Ind2,0], [-50 for AA in range(len(Data[6][:Ind2,0]))], 2000, facecolor='y', alpha=ColStrength)
plt.fill_between(Data[6][Ind1:Ind2,0], Data[6][Ind1:Ind2,5], Data[6][Ind1:Ind2,6], facecolor='wheat', alpha=ColStrength)
plt.fill_between(Data[6][Ind2-1:,0], [-50 for AA in range(len(Data[6][Ind2-1:,0]))], 2000, facecolor='grey', alpha=ColStrength)
plt.ylim(YliM); plt.yticks(YtickS, fontsize=fs1-4); plt.xlabel(r'Magnitude of stress ($I$)', fontsize=fs); plt.ylabel('Concentration', fontsize=fs1)
plt.xticks([1.0,2.0,3.0,4.0,5.0,6.0,7.0,7.5], fontsize=fs1-4); plt.legend([r'p53$_{A}$','p53$_{M}$'], loc='best', prop={'size': Legfs},facecolor=None); plt.xlim(1.0,7.0);
plt.text(2.75,Te1,r' C$_1$', fontsize=fs)
plt.text(1.05,1000,'Active state', fontsize=Legfs)
plt.text(6.0,1500,'Cancer state', fontsize=Legfs,rotation=90)

DATAX=Data[7]
Ind1=98; Ind2=0
for INDEX,dataA1,dataA2,dataM1,dataM2 in zip(range(len(DATAX[:,1])),DATAX[:,1],DATAX[:,2],DATAX[:,3],DATAX[:,4]):
    #print dataM, dataA
    if (dataM1+dataM2)/2.0>(dataA1+dataA2)/2.0: 
        Ind2=INDEX            
        print INDEX, DATAX[:,0][INDEX]; break

plt.subplot(grid[214:300, 200:295])
plt.plot(Data[7][:,0],(Data[7][:,1]+Data[7][:,2])/2.0,'-g' ,Data[7][:,0],(Data[7][:,3]+Data[7][:,4])/2.0,'-r', linewidth=lw)
plt.plot(Data[7][Ind1:Ind2,0],Data[7][Ind1:Ind2,5],'-b' ,Data[7][Ind1:Ind2,0],Data[7][Ind1:Ind2,6],'-k', linewidth=lw)
plt.fill_between(Data[7][:Ind2,0], [-50 for AA in range(len(Data[7][:Ind2,0]))], 2000, facecolor='y', alpha=ColStrength)
plt.fill_between(Data[7][Ind1:Ind2,0], Data[7][Ind1:Ind2,5], Data[7][Ind1:Ind2,6], facecolor='wheat', alpha=ColStrength)
plt.fill_between(Data[7][Ind2-1:,0], [-50 for AA in range(len(Data[7][Ind2-1:,0]))], 2000, facecolor='grey', alpha=ColStrength)
plt.ylim(YliM); plt.yticks([]); plt.xlabel(r'Magnitude of stress ($I$)', fontsize=fs); plt.xticks([1.0,2.0,3.0,4.0,5.0,6.0,7.0,7.5], fontsize=fs1-4);
plt.xlim(1.0,7.0);
plt.text(2.75,Te1,r' C$_2$', fontsize=fs); plt.legend([r'p53$_{A}$','p53$_{M}$'], loc='best', prop={'size': Legfs},facecolor=None)
plt.text(1.05,1000,'Active state', fontsize=Legfs)
plt.text(4.0,1000,'Cancer state', fontsize=Legfs)

DATAX=Data[8]
Ind2=0
for INDEX,dataA1,dataA2,dataM1,dataM2 in zip(range(len(DATAX[:,1])),DATAX[:,1],DATAX[:,2],DATAX[:,3],DATAX[:,4]):
    #print dataM, dataA
    if (dataM1+dataM2)/2.0>(dataA1+dataA2)/2.0: 
        Ind2=INDEX            
        print INDEX, DATAX[:,0][INDEX]; break

plt.subplot(grid[214:300, 305:400])
plt.plot(Data[8][:,0],(Data[8][:,1]+Data[8][:,2])/2.0,'-g' ,Data[8][:,0],(Data[8][:,3]+Data[8][:,4])/2.0,'-r', linewidth=lw)
plt.fill_between(Data[8][:Ind2,0], [-50 for AA in range(len(Data[8][:Ind2,0]))], 2000, facecolor='y', alpha=0.25)
plt.fill_between(Data[8][Ind2-1:,0], [-50 for AA in range(len(Data[8][Ind2-1:,0]))], 2000, facecolor='grey', alpha=0.25)
plt.ylim(YliM); plt.yticks([]); plt.xlabel(r'Magnitude of stress ($I$)', fontsize=fs); plt.xticks([1.0,2.0,3.0,4.0,5.0,6.0,7.0,7.5], fontsize=fs1-4);
plt.xlim(1.0,7.0);
plt.text(1.0,Te1,r' C$_3$', fontsize=fs); plt.legend([r'p53$_{A}$','p53$_{M}$'], loc='center right', prop={'size': Legfs},facecolor=None)
plt.text(1.05,1000,'Active state', fontsize=Legfs)
plt.text(3.0,1000,'Cancer state', fontsize=Legfs)


# Stress ploting 
TIME = np.array(Time)/3600.0
plt.subplot(grid[12:42, 0:60])
plt.plot(TIME, [Stress(Ti)[0] for Ti in Time],'-b', linewidth=lw)
plt.xlim([0-max(TIME)*0.0125,max(TIME)*1.0125]); plt.ylim([-0.1, 1.1])
plt.yticks([]); plt.ylabel('Stress', fontsize=fs-2)
plt.text(0,1.25,'A', fontsize=fs); plt.xticks([0,40,80,120,160], fontsize=fs-4)
plt.title('Constant stress', fontsize=fs1-3)


plt.subplot(grid[118:148, 0:60])
plt.plot(TIME, [Stress(Ti)[1] for Ti in Time],'-b', linewidth=lw)
plt.xlim([0-max(TIME)*0.0125,max(TIME)*1.0125]); plt.ylim([-0.2, 2.2])
plt.yticks([]); plt.ylabel('Stress', fontsize=fs-2)
plt.text(0,2.50,'B', fontsize=fs); plt.xticks([0,40,80,120,160], fontsize=fs-4)
plt.title('Oscillatory stress', fontsize=fs1-3)


plt.subplot(grid[222:248, 0:60])
plt.plot(TIME, [Stress(Ti)[2] for Ti in Time],'-b', linewidth=lw)
plt.xlim([0-max(TIME)*0.0125,max(TIME)*1.0125]); plt.ylim([-0.1, 1.1])
plt.yticks([]); plt.xlabel('Time (hrs)', fontsize=fs-2); plt.ylabel('Stress', fontsize=fs-2)
plt.text(0,1.25,'C', fontsize=fs); plt.xticks([0,40,80,120,160], fontsize=fs-4)
plt.title('Decaying stress', fontsize=fs1-3)


plt.savefig('fig4.eps')
plt.savefig('fig4.pdf')
#plt.show()

end_time = datetime.now()
print 'MAIN PROGRAM IS COMPLETED||Duration||H:M:S||{}'.format(end_time - start_time), '\n'
###############################################################################################################
###############################################################################################################
