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

Data=[np.genfromtxt('./data/Ser{}.dat'.format(Ind)) for Ind in range(1,10)]

A = 1.0; T=3600*6; lemda=0.05/3600.0; Time=[time for index, time in zip(range(12*10**6), np.linspace(0, 6*10**5, 12*10**6)) if index%100==0]
def Stress(t): return [A, A + A*math.sin(2*math.pi*t/T), A*np.exp(-lemda*t)]

plt.figure(figsize=(16,7.5))
grid = plt.GridSpec(300, 400, left=0.023, right=0.990, top=0.985, bottom=0.080, hspace = 0.0, wspace=0.0)

lw=1.25; fs=18; Legfs=11; fs1=10; fs2=15
XtickS=[0,40,80,120,160]; YtickS=[0,600,1200,1800]
XliM=[-4.25,170]; YliM=[-50,2000]; Te0=0; Te1=1775

# Data ploting
plt.subplot(grid[0:92, 95:190])
plt.plot(Data[0][:,0],Data[0][:,1],'-g' ,Data[0][:,0],Data[0][:,2],'-r', linewidth=lw)
plt.xlim(XliM); plt.ylim(YliM); plt.xticks([]); plt.yticks(YtickS, fontsize=fs-5); plt.ylabel('Concentration', fontsize=fs)
plt.legend([r'p53$_{A}$','p53$_{M}$'], loc='best', prop={'size': Legfs},facecolor=None); plt.text(Te0,Te1,r'  A$_1$', fontsize=fs)
plt.text(25,1800,'Active state', fontsize=fs2)


plt.subplot(grid[0:92, 200:295])
plt.plot(Data[1][:,0],Data[1][:,1],'-g' ,Data[1][:,0],Data[1][:,2],'-r', linewidth=lw)
plt.xlim(XliM); plt.ylim(YliM); plt.xticks([]); plt.yticks([]); plt.text(Te0,Te1,r'  A$_2$', fontsize=fs)
plt.legend([r'p53$_{A}$','p53$_{M}$'], loc='best', prop={'size': Legfs},facecolor=None)
plt.text(25,1800,'Apoptotic state', fontsize=fs2)

plt.subplot(grid[0:92, 305:400])
plt.plot(Data[2][:,0],Data[2][:,1],'-g' ,Data[2][:,0],Data[2][:,2],'-r', linewidth=lw)
plt.xlim(XliM); plt.ylim(YliM); plt.xticks([]); plt.yticks([]); plt.text(Te0,Te1,r'  A$_3$', fontsize=fs)
plt.legend([r'p53$_{A}$','p53$_{M}$'], loc='center right', prop={'size': Legfs},facecolor=None)
plt.text(25,1800,'Cancer state', fontsize=fs2)

plt.subplot(grid[104:196, 95:190])
plt.plot(Data[3][:,0],Data[3][:,1],'-g' ,Data[3][:,0],Data[3][:,2],'-r', linewidth=lw)
plt.xlim(XliM); plt.ylim(YliM); plt.xticks([]); plt.yticks(YtickS, fontsize=fs-5); plt.ylabel('Concentration', fontsize=fs)
plt.text(Te0,Te1,r'  B$_1$', fontsize=fs); plt.legend([r'p53$_{A}$','p53$_{M}$'], loc='best', prop={'size': Legfs},facecolor=None)
plt.text(25,1800,'Active state', fontsize=fs2)

plt.subplot(grid[104:196, 200:295])
plt.plot(Data[4][:,0],Data[4][:,1],'-g' ,Data[4][:,0],Data[4][:,2],'-r', linewidth=lw)
plt.xlim(XliM); plt.ylim(YliM); plt.xticks([]); plt.yticks([]); plt.text(Te0,Te1,r'  B$_2$', fontsize=fs)
plt.legend([r'p53$_{A}$','p53$_{M}$'], loc='best', prop={'size': Legfs},facecolor=None)
plt.text(25,1800,'Apoptotic state', fontsize=fs2)

plt.subplot(grid[104:196, 305:400])
plt.plot(Data[5][:,0],Data[5][:,1],'-g' ,Data[5][:,0],Data[5][:,2],'-r', linewidth=lw)
plt.xlim(XliM); plt.ylim(YliM); plt.xticks([]); plt.yticks([]); plt.text(Te0,Te1,r'  B$_3$', fontsize=fs)
plt.legend([r'p53$_{A}$','p53$_{M}$'], loc='center right', prop={'size': Legfs},facecolor=None)
plt.text(25,1800,'Cancer state', fontsize=fs2)

plt.subplot(grid[208:300, 95:190])
plt.plot(Data[6][:,0],Data[6][:,1],'-g' ,Data[6][:,0],Data[6][:,2],'-r', linewidth=lw)
plt.xlim(XliM); plt.ylim(YliM); plt.yticks(YtickS, fontsize=fs-5); plt.xlabel('Time (hrs)', fontsize=fs); plt.ylabel('Concentration', fontsize=fs)
plt.xticks(XtickS, fontsize=fs-5); plt.legend([r'p53$_{A}$','p53$_{M}$'], loc='best', prop={'size': Legfs},facecolor=None);plt.text(Te0,Te1,r'  C$_1$', fontsize=fs)
plt.text(25,1800,'Active state', fontsize=fs2)

plt.subplot(grid[208:300, 200:295])
plt.plot(Data[7][:,0],Data[7][:,1],'-g' ,Data[7][:,0],Data[7][:,2],'-r', linewidth=lw)
plt.xlim(XliM); plt.ylim(YliM); plt.yticks([]); plt.xlabel('Time (hrs)', fontsize=fs); plt.xticks(XtickS, fontsize=fs-5);
plt.text(Te0,Te1,r'  C$_2$', fontsize=fs); plt.legend([r'p53$_{A}$','p53$_{M}$'], loc='best', prop={'size': Legfs},facecolor=None)
plt.text(25,1800,'Active state', fontsize=fs2)

plt.subplot(grid[208:300, 305:400])
plt.plot(Data[8][:,0],Data[8][:,1],'-g' ,Data[8][:,0],Data[8][:,2],'-r', linewidth=lw)
plt.xlim(XliM); plt.ylim(YliM); plt.yticks([]); plt.xlabel('Time (hrs)', fontsize=fs); plt.xticks(XtickS, fontsize=fs-5);
plt.text(Te0,Te1,r'  C$_3$', fontsize=fs); plt.legend([r'p53$_{A}$','p53$_{M}$'], loc='center right', prop={'size': Legfs},facecolor=None)
plt.text(25,1800,'Cancer state', fontsize=fs2)

# Stress ploting 
TIME = np.array(Time)/3600.0
plt.subplot(grid[10:42, 0:60])
plt.plot(TIME, [Stress(Ti)[0] for Ti in Time],'-b', linewidth=lw)
plt.xlim([0-max(TIME)*0.0125,max(TIME)*1.0125]); plt.ylim([-0.1, 1.1])
plt.yticks([]); plt.ylabel('Stress', fontsize=fs-2)
plt.text(0,1.25,'A', fontsize=fs); plt.xticks(XtickS, fontsize=fs-5)
plt.title('Constant stress', fontsize=fs-4)


plt.subplot(grid[112:142, 0:60])
plt.plot(TIME, [Stress(Ti)[1] for Ti in Time],'-b', linewidth=lw)
plt.xlim([0-max(TIME)*0.0125,max(TIME)*1.0125]); plt.ylim([-0.2, 2.2])
plt.yticks([]); plt.ylabel('Stress', fontsize=fs-2)
plt.text(0,2.50,'B', fontsize=fs); plt.xticks(XtickS, fontsize=fs-5)
plt.title('Oscillatory stress', fontsize=fs-4)


plt.subplot(grid[216:242, 0:60])
plt.plot(TIME, [Stress(Ti)[2] for Ti in Time],'-b', linewidth=lw)
plt.xlim([0-max(TIME)*0.0125,max(TIME)*1.0125]); plt.ylim([-0.1, 1.1])
plt.yticks([]); plt.xlabel('Time (hrs)', fontsize=fs-2); plt.ylabel('Stress', fontsize=fs-2)
plt.text(0,1.25,'C', fontsize=fs); plt.xticks(XtickS, fontsize=fs-5)
plt.title('Decaying stress', fontsize=fs-4)


plt.savefig('fig2.pdf')
plt.savefig('fig2.eps')

#plt.show()

end_time = datetime.now()
print 'MAIN PROGRAM IS COMPLETED||Duration||H:M:S||{}'.format(end_time - start_time), '\n'
###############################################################################################################
###############################################################################################################
