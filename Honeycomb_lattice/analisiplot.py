import numpy as np
import matplotlib.pyplot as plt

import grafici

path = r'datiplot/'        #path dei dati

L = np.arange(10, 50+1, 5) #reticoli considerati

#liste i cui elementi saranno le curve delle quantità termodinamiche
E = [] ; dE = []           #energia del sistema ed errore
M = [] ; dM = []           #magnetizzazione del sistema ed errore
C = [] ; dC = []           #calore specifico del sistema ed errore
X = [] ; dX = []           #suscettività del sistema ed errore

CB = [] ; dCB = []         #cumulante di binder

for i, l in enumerate(L):
    #Leggo i dati al variare di L e li conservo
    Data = np.loadtxt(path+f'dati{l}.dat', unpack=True)
    ene, mag, cal, chi, cb, dene, dmag, dcal, dchi, dcb = Data
    E.append(ene) ; dE.append(dene)
    M.append(mag) ; dM.append(dmag)
    C.append(cal) ; dC.append(dcal)
    X.append(chi) ; dX.append(dchi)
    CB.append(cb) ; dCB.append(dcb)



#-------------------------------------------------------------------------
#Alcuni plot a titolo espositivo
#-------------------------------------------------------------------------


#calcolo array temperature
par = np.loadtxt(r'init.txt', max_rows=6, unpack=True)
bmin, bmax, npassi = par[3:6]

npassi = int(npassi)
B = np.zeros(npassi)

for i in range(1, npassi+1):
    B[i-1] = bmin + (i-1)*(bmax - bmin)/(npassi-1)

Title = f'Ising 2D su reticolo esagonale'
xlabel = r'$\beta$ [u.a.]' 

#grafici.plot(B, E, dE, L, 1, Title, xlabel, "Energia [u.a.]")

#grafici.plot(B, M, dM, L, 2, Title, xlabel, "Magetizzazione [u.a.]")
      
#grafici.plot(B, C, dC, L, 3, Title, xlabel, "Calore specifico [u.a.]")
                 
#grafici.plot(B, X, dX, L, 4, Title, xlabel, "Suscettività [u.a.]") 
 
#grafici.plotbinder(B, CB, dCB, L, 5, Title, xlabel, "Cumulante di binder") 

#valore di beta ricavato
bc = 0.661
dbc = 0.001
print(f"beta critico = {bc:.3f} +- {dbc:.3f} \n ")

def F(x, m, g, q):
    '''
    funzione modello per i fit successivi
    '''
    return m*x**g + q


#funzione per estrarre il massimo della suscettività e relativo errore
def MS(*arg):
    '''
    *arg permette alla funzione di prendere in input infiniti argomenti.
    La funzione è scritta affinchè la prima metà degli argomenti passati siano
    gli array dei valori centrali mentre la seconda parte quelli degli errori.
    Restitusce due array contenenti il massimo dei valori centrali e relativi errori
    '''
    Q = len(arg)//2 #lunghezza degli array
    l = np.zeros(Q)
    dl = np.zeros(Q)
    k = np.zeros(Q, dtype=np.int64)#verra usato come indice quindi deve essere intero
	
    for i in range(len(arg)):

        if i < Q:
            l[i] = np.max(arg[i]) #trovo il massimo
            k[i] = np.where(arg[i] == l[i])[0][0] #conservo il rispettivo indice
			
        else:
            #trovo il valore dell'errore del punto usando l'indice
            dl[i-Q] = arg[i][k[i-Q]]
			
    return l, dl



#-------------------------------------------------------------------------
#Stima di gamma/nu
#-------------------------------------------------------------------------

print('stima di gamma/nu \n')

TOT = X + dX
MX, dMX = MS(*TOT)

Title = 'Massimo della suscettività al variare di L'
ylabel = r'$\chi_{max}$ [a.u.]'
xlabel = r'L [a.u.]'
init = np.array([1, 1, 1])
 

pars, dpars = grafici.fit(F, L, MX, dMX, init, 6, Title, xlabel, ylabel)


#-------------------------------------------------------------------------
#Stima di gamma
#-------------------------------------------------------------------------

print('----------------------------------------------- \n')
print('stima di gamma \n')

# dato che prima abbiamo stimato gamma/nu con
# questa stima possiamo stimare sia gamma che nu

index1 = 57                     #seleziono solo un range dei dati
index2 = 87
x  =      B[index1:index2]-bc   #temperatura ridotta
y  =  X[-1][index1:index2]      #suscettività ultimo reticolo
dy = dX[-1][index1:index2]      #errore associato

#valori che mi aspetto per i parametri ottimali
init = np.array([0.00352, -7/4, -0.013]) #aiutano la convergenza del fit

Title = 'Suscettivita reticolo 50X50 dopo il punto critico'
ylabel = r'$\chi$ [a.u.]'
xlabel = r'$\beta - \beta_c$ [u.a.]'

pars1, dpars1 = grafici.fit(F, x, y, dy, init, 7, Title, xlabel, ylabel)


#-------------------------------------------------------------------------
#Stima di beta
#-------------------------------------------------------------------------

print('----------------------------------------------- \n')
print('stima di beta \n')

def G(x, m, g):
    '''
    funzione modello per il fit di beta senza intercetta
    '''
    return m*x**g


#Per stimare questo indice va calcolato nelle altre osservabili l'indice di beta pseudo_crit
indxcrit = 51 

#selezione della manetizzazione al punto critico per i vari reticoli
y  = np.array([ M[i][indxcrit] for i in range(len(L))])
dy = np.array([dM[i][indxcrit] for i in range(len(L))])

#valori che mi aspetto per i parametri ottimali
init = np.array([0.9, -1/8]) #aiutano la convergenza del fit

Title = 'Magnetizzazione al punto critico al variare di L'
ylabel = r'$M_{\beta_c}$ [a.u.]'
xlabel = r'L [a.u.]'

pars2, dpars2 = grafici.fit(G, L, y, dy, init, 8, Title, xlabel, ylabel)



#-------------------------------------------------------------------------
#Stima di alpha
#-------------------------------------------------------------------------

print('----------------------------------------------- \n')
print('stima di alpha \n')

def Fa(x, m, g, q):
    '''funzione modello per i fit successivi
    '''
    return m*x**g*(1 + q*np.log(x) + np.log(2))

TOT = C + dC
CM, dCM = MS(*TOT)


#valori che mi aspetto per i parametri ottimali
init = np.array([0.5, 0.1, 2]) #aiutano la convergenza del fit

Title = 'Massimo del calore specifico al variare di L'
ylabel = r'$C_{max}$ [a.u.]'
xlabel = r'$\beta - \beta_c$ [u.a.]'

pars3, dpars3 = grafici.fit(Fa, L, CM, dCM, init, 9, Title, xlabel, ylabel)
print('----------------------------------------------- \n')


#-------------------------------------------------------------------------
#Calcolo deigli indici critici ed associati errori
#-------------------------------------------------------------------------


g = -pars1[1]
dg = dpars[1]

n = -pars1[1]/pars[1]
dn = np.sqrt((dpars[1]/pars[1])**2 + (dpars[1]/pars1[1])**2)*n

b = abs(pars2[1])
db = dpars2[1]

a = abs(pars3[1])
da = dpars3[1]


print('\n')
print(f"gamma = {g:.5f} +- {dg:.5f}")
print(f"beta  = {b:.5f} +- {db:.5f}")
print(f"nu    = {n:.5f} +- {dn:.5f}")
print(f"alpha = {a:.5f} +- {da:.5f}")

#distanza dai valori teorici
print('dg/nu = %.5f' %(pars[1]-7/4))
print('db    = %.5f' %(pars2[1]-1/8))
print('dg    = %.5f' %(pars1[1]+7/4))
print('da    = %.5f' %(pars3[1]-0))




#-------------------------------------------------------------------------
#Grafici del finite size scaling
#-------------------------------------------------------------------------



Title = 'Finite size scaling della magnetizzazione'
xlabel = r'$(\beta-\beta_c)L^{1/ \nu}$'
ylabel = r'$|M|/L^{-b/ \nu}$'
    
grafici.FSS(B, 1/n, -b/n, bc, M, dM, L, 10, Title, xlabel, ylabel)

Title = 'Finite size scaling della suscettività'
xlabel = r'$(\beta-\beta_c)L^{1/ \nu}$'
ylabel = r'$ \chi /L^{\gamma/ \nu}$'

grafici.FSS(B, 1/n, g/n, bc, X, dX, L, 11, Title, xlabel, ylabel) 

Title = 'Finite size scaling del calore specifico'
xlabel = r'$(\beta-\beta_c)L^{1/ \nu}$'
ylabel = r'$ C /L^{\alpha/ \nu}$'

#Per tenere conto delle correzzioni logaritmiche
for i in range(len(L)):
    for j in range(len(B)):
        C[i][j]  /= np.log(2*L[i]**2)
        dC[i][j] /= np.log(2*L[i]**2)

grafici.FSS(B, 1/n, a/n, bc, C, dC, L, 12, Title, xlabel, ylabel) 



plt.show()
