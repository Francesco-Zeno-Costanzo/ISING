import numpy as np
import matplotlib.pyplot as plt

B, E, M, C, X= np.loadtxt(r'C:\Users\franc\Desktop\C\ising20.txt', unpack = True)
B1, E1, M1, C1, X1= np.loadtxt(r'C:\Users\franc\Desktop\C\ising30.txt', unpack = True)
B2, E2, M2, C2, X2= np.loadtxt(r'C:\Users\franc\Desktop\C\ising40.txt', unpack = True)
B3, E3, M3, C3, X3= np.loadtxt(r'C:\Users\franc\Desktop\C\ising50.txt', unpack = True)
#B4, E4, M4, C4, X4= np.loadtxt(r'C:\Users\franc\Desktop\C\\ising10.txt', unpack = True)
H=0
n=1
b=1/8
g=7/4
a=0
bc=0.4406868

plt.figure(1)
plt.suptitle("Simulazione del modello di Ising 2D tramite Metropolis \n Campo magnetico esterno B=%d " %( H))

plt.subplot(221)
plt.xlabel(r"$\beta$", fontsize=15)
plt.ylabel("Energia", fontsize=15)
plt.grid()
plt.errorbar(B,  E,   fmt='.', color='black', label='L=20')
plt.errorbar(B1, E1,  fmt='v', color='black', label='L=30')
plt.errorbar(B2, E2,  fmt='*', color='black', label='L=40')
plt.errorbar(B3, E3,  fmt='^', color='black', label='L=50')
#plt.errorbar(B4, E4,  fmt='<', color='black', label='L=10')
plt.legend(loc='best')

plt.subplot(222)
plt.xlabel(r"$\beta$", fontsize=15)
plt.ylabel("Magnetizzazione", fontsize=15)
plt.grid()
plt.errorbar(B,  M,  fmt='.', color='black', label='L=20')
plt.errorbar(B1, M1, fmt='v', color='black', label='L=30')
plt.errorbar(B2, M2, fmt='*', color='black', label='L=40')
plt.errorbar(B3, M3, fmt='^', color='black', label='L=50')
#plt.errorbar(B4, M4, fmt='<', color='black', label='L=10')
plt.legend(loc='best')

plt.subplot(223)
plt.xlabel(r"$\beta$", fontsize=15)
plt.ylabel("Calore specifico", fontsize=15)
plt.grid()
plt.errorbar(B,  C,  fmt='.', color='black', label='L=20')
plt.errorbar(B1, C1, fmt='v', color='black', label='L=30')
plt.errorbar(B2, C2, fmt='*', color='black', label='L=40')
plt.errorbar(B3, C3, fmt='^', color='black', label='L=50')
#plt.errorbar(B4, C4, fmt='<', color='black', label='L=10')
plt.legend(loc='best')

plt.subplot(224)
plt.xlabel(r"$\beta$", fontsize=15)
plt.ylabel("Suscettività", fontsize=15)
plt.grid()
plt.errorbar(B,  X,  fmt='.', color='black', label='L=20')
plt.errorbar(B1, X1, fmt='v', color='black', label='L=30')
plt.errorbar(B2, X2, fmt='*', color='black', label='L=40')
plt.errorbar(B3, X3, fmt='^', color='black', label='L=50')
#plt.errorbar(B4, X4, fmt='<', color='black', label='L=10')
plt.legend(loc='best')

'''
fig2 = plt.figure(2)
plt.suptitle("Simulazione del modello di Ising 2D tramite Metropolis \n Campo magnetico esterno B=%d \n" r"$\alpha = %d, b= %f, \nu = %d, \gamma = %f, \beta_c= %f$" %( H, a, b, n, g, bc))
frame1=fig2.add_axes((.55,.1,.4,.35))
#frame1=fig1.add_axes((trasla lateralmente, trasla verticamente, larghezza, altezza))
plt.xlabel(r"$(\beta-\beta_c)L^{1/ \nu}$", fontsize=15)
plt.ylabel(r"$|M|/L^{-b/ \nu}$", fontsize=15)
plt.grid()
plt.errorbar((B-bc)*20**(1/n),  M/(20**(-b/n)),  fmt='.', color='black', label='L=20')
plt.errorbar((B1-bc)*30**(1/n), M1/(30**(-b/n)), fmt='v', color='black', label='L=30')
plt.errorbar((B2-bc)*40**(1/n), M2/(40**(-b/n)), fmt='*', color='black', label='L=40')
plt.errorbar((B3-bc)*50**(1/n), M3/(50**(-b/n)), fmt='^', color='black', label='L=50')
#plt.errorbar((B4-bc)*10**(1/n), M4/(10**(-b/n)), fmt='<', color='black', label='L=10')
plt.legend(loc='best')

frame2=fig2.add_axes((.06,.1,.42,.35))
plt.xlabel(r"$(\beta-\beta_c)L^{1/ \nu}$", fontsize=15)
plt.ylabel(r"$C/L^{\alpha/ \nu}$", fontsize=15)
plt.grid()
plt.errorbar((B-bc)*20**(1/n),  C/(20**(a/n)),  fmt='.', color='black', label='L=20')
plt.errorbar((B1-bc)*30**(1/n), C1/(30**(a/n)), fmt='v', color='black', label='L=30')
plt.errorbar((B2-bc)*40**(1/n), C2/(40**(a/n)), fmt='*', color='black', label='L=40')
plt.errorbar((B3-bc)*50**(1/n), C3/(50**(a/n)), fmt='^', color='black', label='L=50')
#plt.errorbar((B4-bc)*10**(1/n), C4/(10**(a/n)), fmt='<', color='black', label='L=10')
plt.legend(loc='best')

frame3=fig2.add_axes((.3,.55,.422,.35))
plt.xlabel(r"$(\beta-\beta_c)L^{1/ \nu}$", fontsize=15)
plt.ylabel(r"$ \chi /L^{\gamma/ \nu}$", fontsize=15)
plt.grid()
plt.errorbar((B-bc)*20**(1/n),  X/(20**(g/n)),  fmt='.', color='black', label='L=20')
plt.errorbar((B1-bc)*30**(1/n), X1/(30**(g/n)), fmt='v', color='black', label='L=30')
plt.errorbar((B2-bc)*40**(1/n), X2/(40**(g/n)), fmt='*', color='black', label='L=40')
plt.errorbar((B3-bc)*50**(1/n), X3/(50**(g/n)), fmt='^', color='black', label='L=50')
#plt.errorbar((B4-bc)*10**(1/n), X4/(10**(g/n)), fmt='<', color='black', label='L=10')
plt.legend(loc='best')
'''

fig2 = plt.figure(2)
plt.suptitle("Simulazione del modello di Ising 2D tramite Metropolis \n Campo magnetico esterno B=%d \n" r"$\alpha = %d, b= %f, \nu = %d, \gamma = %f, \beta_c= %f$" %( H, a, b, n, g, bc))
plt.xlabel(r"$(\beta-\beta_c)L^{1/ \nu}$", fontsize=15)
plt.ylabel(r"$|M|/L^{-b/ \nu}$", fontsize=15)
plt.grid()
plt.errorbar((1/B3-1/bc)/(1/bc), X3, fmt='^', color='black', label='L=50')

plt.legend(loc='best')
plt.show()
