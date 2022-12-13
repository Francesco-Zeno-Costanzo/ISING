import numpy as np
import matplotlib.pyplot as plt

import grafici

#il numero finale identifica la grandezza del reticolo
E10, M10, C10, X10, cb10, dE10, dM10, dC10, dX10, dcb10 = np.loadtxt(r'datiplot/dati10.dat', unpack=True)
E15, M15, C15, X15, cb15, dE15, dM15, dC15, dX15, dcb15 = np.loadtxt(r'datiplot/dati15.dat', unpack=True)
"""
E20, M20, C20, X20, cb20, dE20, dM20, dC20, dX20, dcb20 = np.loadtxt(r'datiplot/20.dat', unpack=True)
E25, M25, C25, X25, cb25, dE25, dM25, dC25, dX25, dcb25 = np.loadtxt(r'datiplot/25.dat', unpack=True)
E30, M30, C30, X30, cb30, dE30, dM30, dC30, dX30, dcb30 = np.loadtxt(r'datiplot/30.dat', unpack=True)
E35, M35, C35, X35, cb35, dE35, dM35, dC35, dX35, dcb35 = np.loadtxt(r'datiplot/35.dat', unpack=True)
E40, M40, C40, X40, cb40, dE40, dM40, dC40, dX40, dcb40 = np.loadtxt(r'datiplot/40.dat', unpack=True)
E45, M45, C45, X45, cb45, dE45, dM45, dC45, dX45, dcb45 = np.loadtxt(r'datiplot/45.dat', unpack=True)
E50, M50, C50, X50, cb50, dE50, dM50, dC50, dX50, dcb50 = np.loadtxt(r'datiplot/50.dat', unpack=True)
"""

H = 0


##Alcuni plot a titolo espositivo

#calcolo array temperature
par = np.loadtxt(r'init.txt', unpack=True)
bmin, bmax, npassi = par[3:6]

npassi = int(npassi)

B = np.zeros(npassi)
for i in range(1, npassi+1):
    B[i-1] = bmin + (i-1)*(bmax - bmin)/(npassi-1)

Title = f'Simulazione del modello di Ising 2D tramite Metropolis \n Campo magnetico esterno B={H}'
xlabel = r'$\beta$ [u.a.]'   
grafici.plot(B, [E10, E15], 
                [dE10, dE15],
                [10, 15], 1, Title, xlabel, "Energia [u.a.]")
                
grafici.plot(B, [M10, M15],
                [dM10, dM15],
                [10, 15], 2, Title, xlabel, "Magetizzazione [u.a.]")
                
grafici.plot(B, [C10, C15], 
                [dC10, dC15],
                [10, 15], 3, Title, xlabel, "Calore specifico [u.a.]")
                 
grafici.plot(B, [X10, X15], 
                [dX10, dX15], 
                [10, 15], 4, Title, xlabel, "Suscettività [u.a.]") 
 
grafici.plotbinder(B, [cb10, cb15], 
                      [dcb10, dcb15],
                      [10, 15], 5, Title, xlabel, "Cumulante di binder") 

plt.show()
