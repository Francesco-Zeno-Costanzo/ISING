"""
Nel seguente codice l'analisi è fatta in temperatura e non beta
"""
import numpy as np
import matplotlib.pyplot as plt

import grafici

path = r'dataplot/'  # path dei dati

J1  = 1
J2  = ["0", "01"]
J2p = [0, -0.1]
# liste i cui elementi saranno le curve delle quantità termodinamiche
E = [] ; dE = []  # energia del sistema ed errore
M = [] ; dM = []  # magnetizzazione del sistema ed errore
C = [] ; dC = []  # calore specifico del sistema ed errore
X = [] ; dX = []  # suscettività del sistema ed errore

for i, j in enumerate(J2):
    # Leggo i dati al variare di L e li conservo
    Data = np.loadtxt(path+f'data_L_10_J2_{j}.dat', unpack=True)
    ene, mag, cal, chi, dene, dmag, dcal, dchi = Data
    E.append(ene) ; dE.append(dene)
    M.append(mag) ; dM.append(dmag)
    C.append(cal) ; dC.append(dcal)
    X.append(chi) ; dX.append(dchi)


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
xlabel = r'$T$ [u.a.]'

grafici.plot(T, E, dE, J2, 1, Title, xlabel, "Energia [u.a.]")
grafici.plot(T, M, dM, J2, 2, Title, xlabel, "Magetizzazione [u.a.]")
grafici.plot(T, C, dC, J2, 3, Title, xlabel, "Calore specifico [u.a.]")
grafici.plot(T, X, dX, J2, 4, Title, xlabel, "Suscettività [u.a.]")

#===================================================================

TMAX = []
for i in range(len(J2)):
    idx = np.where(C[i]==max(C[i]))[0][0]
    TMAX.append(T[idx])

plt.figure(5)
plt.plot(J2p, TMAX, 'bo')
plt.grid()
plt.show()
