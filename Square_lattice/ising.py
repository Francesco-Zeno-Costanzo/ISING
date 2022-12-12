import time
import glob
import numpy as np
import random as rn
import imageio as io
from PIL import Image
import matplotlib.pyplot as plt

start_time=time.time()
J=1
H=0
L=100
R=np.zeros((L, L))

def init():
    for i in range(L):
        for j in range(L):
            if i<L/2:
                R[i][j]=1
            else:
                R[i][j]=-1

def Metropolis(B):
    for i in range(L):
        for j in range(L):
            c = rn.randint(0, L-1)
            b = rn.randint(0, L-1)
            dE=2*J*(R[c][b])*( R[c][(b+1)%L] + R[(c+1)%L][b] + R[c][(b-1)%L] + R[(c-1)%L][b]) + H*(R[c][b])
            if dE<0.0:
                R[c][b]=-R[c][b]
            else:
                if(rn.random())<np.exp(-dE*B):
                    R[c][b]=-R[c][b]



N=200
beta=0.4
init()
fig = plt.figure(0)
plt.xlim(0, L)
plt.ylim(0, L)
for i in range(L):
    for j in range(L):
        if R[i][j]>0:
            plt.errorbar(i, j, fmt='o', color='k')
        else:
            plt.errorbar(i, j, fmt='o', color='w')
plt.title("Raggiungimento ergodicità \n Campo magnetico esterno B=%d \n "r"$\beta$=%f" %( H, beta))

plt.savefig(r'C:\Users\franc\Desktop\codici python\gif/%d'%(0))
plt.close(fig)

for t in range(1, N+1):
    Metropolis(beta)
    fig = plt.figure(t)
    plt.xlim(0, L)
    plt.ylim(0, L)
    for i in range(L):
        for j in range(L):
            if R[i][j]>0:
                plt.errorbar(i, j, fmt='o', color='k')
            else:
                plt.errorbar(i, j, fmt='o', color='w')
    plt.title("Raggiungimento ergodicità \n Campo magnetico esterno B=%d \n " r"$\beta$=%f" %( H, beta))
    plt.savefig(r'C:\Users\franc\Desktop\codici python\gif/%d'%(t))
    plt.close(fig)

b=(time.time() - start_time)
print("--- %s secondi per i plot ---" %b)

frames=[]
imgs=sorted(glob.glob(r'C:\Users\franc\Desktop\codici python\gif/*.png'))
imgs.sort(key=len)
for i in imgs:
    new_frame=Image.open(i)
    frames.append(new_frame)


frames[0].save('erg.gif',format='GIF',append_images=frames[:],save_all=True,duration=100,loop=0)

l=(time.time() - start_time-b)
print("--- %s secondi per creare l'animazione ---" %l)
