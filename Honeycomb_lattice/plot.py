"""
In this code we use temperature instead of beta
"""
import numpy as np
import matplotlib.pyplot as plt

import grafici

path = r'dataplot/'  # path dei dati

J1  = 1
J2  = ["0", "0025","005", "0075", "01", "0125", "015", "0175", "02"]
J2p = [0, -0.025, -0.05, -0.075, -0.1, -0.125, -0.15, -0.175, -0.2]

# liste i cui elementi saranno le curve delle quantità termodinamiche
E  = [] ; dE  = []  # energia del sistema ed errore
M  = [] ; dM  = []  # magnetizzazione del sistema ed errore
C  = [] ; dC  = []  # calore specifico del sistema ed errore
X  = [] ; dX  = []  # suscettività del sistema ed errore
CB = [] ; dCB = []  # cumulante di binder

for i, j in enumerate(J2):
    # Leggo i dati al variare di L e li conservo
    Data = np.loadtxt(path+f'data_L_50_J2_{j}.dat', unpack=True)
    ene, mag, cal, chi, cb, dene, dmag, dcal, dchi, dcb = Data
    E.append(ene) ; dE.append(dene)
    M.append(mag) ; dM.append(dmag)
    C.append(cal) ; dC.append(dcal)
    X.append(chi) ; dX.append(dchi)
    CB.append(cb) ; dCB.append(dcb)


# calcolo array temperature
par = np.loadtxt(r'init.txt', max_rows=6, unpack=True)
Tmin, Tmax, npassi = par[3:6]

npassi = int(npassi)
T = np.zeros(npassi)

for i in range(1, npassi+1):
    T[i-1] = Tmin + (i-1)*(Tmax - Tmin)/(npassi-1)

C = [C[i]/T**2 for i in range(len(C))]
X = [X[i]/T**2 for i in range(len(X))]

Title = f'Ising 2D su reticolo esagonale'
xlabel = r'$\beta$ [u.a.]'
leg = r"$J_2$"

grafici.plot(T, E,  dE,  J2p, leg, 1, Title, xlabel, "Energia [u.a.]")
grafici.plot(T, M,  dM,  J2p, leg, 2, Title, xlabel, "Magetizzazione [u.a.]")
grafici.plot(T, C,  dC,  J2p, leg, 3, Title, xlabel, "Calore specifico [u.a.]")
grafici.plot(T, X,  dX,  J2p, leg, 4, Title, xlabel, "Suscettività [u.a.]")
grafici.plot(T, CB, dCB, J2p, leg, 5, Title, xlabel, "Cumulante [u.a.]")

#===================================================================

TMAX = []
for i in range(len(J2)):
    idx = np.where(C[i]==max(C[i]))[0][0]
    TMAX.append(T[idx])

TMAX = np.array(TMAX)

plt.figure(6)
plt.title("Phase diagram", fontsize=15)
plt.ylabel(r"$T_c$", fontsize=15)
plt.xlabel(r"$J_2$", fontsize=15)
plt.plot(J2p, TMAX, 'bo')#, marker='o')
plt.plot([0, -1/4], [2/np.log(2+np.sqrt(3)), 0], "rp")
plt.ylim(-0.05, 1.6)
plt.xlim(-0.3, 0.01)
plt.grid()
plt.show()
