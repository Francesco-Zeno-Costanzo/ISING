import numpy as np
import matplotlib.pyplot as plt

import grafici

#il numero finale identifica la grandezza del reticolo
E10, M10, C10, X10, cb10, dE10, dM10, dC10, dX10, dcb10 = np.loadtxt(r'datiplot/dati10.dat', unpack=True)
E20, M20, C20, X20, cb20, dE20, dM20, dC20, dX20, dcb20 = np.loadtxt(r'datiplot/dati20.dat', unpack=True)
E30, M30, C30, X30, cb30, dE30, dM30, dC30, dX30, dcb30 = np.loadtxt(r'datiplot/dati30.dat', unpack=True)
E40, M40, C40, X40, cb40, dE40, dM40, dC40, dX40, dcb40 = np.loadtxt(r'datiplot/dati40.dat', unpack=True)
E50, M50, C50, X50, cb50, dE50, dM50, dC50, dX50, dcb50 = np.loadtxt(r'datiplot/dati50.dat', unpack=True)

ret = [10, 20, 30, 40, 50]
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
grafici.plot(B, [E10, E20, E30, E40, E50], 
                [dE10, dE20, dE30, dE40, dE50],
                ret, 1, Title, xlabel, "Energia [u.a.]")
                
grafici.plot(B, [M10, M20, M30, M40, M50],
                [dM10, dM20, dM30, dM40, dM50],
                ret, 2, Title, xlabel, "Magetizzazione [u.a.]")
                
grafici.plot(B, [C10, C20, C30, C40, C50], 
                [dC10, dC20, dC30, dC40, dC50],
                ret, 3, Title, xlabel, "Calore specifico [u.a.]")
                 
grafici.plot(B, [X10, X20, X30, X40, X50], 
                [dX10, dX20, dX30, dX40, dX50], 
                ret, 4, Title, xlabel, "Suscettività [u.a.]") 
 
grafici.plotbinder(B, [cb10, cb20, cb30, cb40, cb50], 
                      [dcb10, dcb20, dcb30, dcb40, dcb50],
                      ret, 5, Title, xlabel, "Cumulante di binder") 

#valore di beta ricavato
bc = 0.665
dbc = 0.001
print(f"beta critico = {bc:.3f} +- {dbc:.3f}")

def F(x, m, g, q):
    '''funzione modello per i fit successivi
    '''
    return m*x**g + q

##stima di gamma/nu

print('\n')

print('stima di gamma/nu')

#funzione per estrarre il massimo della suscettività e relativo errore
def MS(*arg):
    '''
    *arg permette alla funzione di prendere in input infiniti argomenti.
    La funzione è scritta affinche la prima metà degli argomenti passati siano
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
	
MX, dMX = MS(X10, X20, X30, X40, X50, dX10, dX20, dX30, dX40, dX50)
N = np.arange(10, 51, 10)
Title = 'Massimo della suscettività al variare di L'
ylabel = r'$\chi_{max} [a.u.]$'
xlabel = r'L [a.u.]'
init = np.array([1, 1, 1])
  
pars, dpars = grafici.fit(F, N, MX, dMX, init, 6, Title, xlabel, ylabel)

##stima di gamma
print('stima di gamma')

#dato che prima abbiamo stimato gamma/nu con questa stima possiamo stimare sia gamma che nu


#seleziono solo un range dei dati
x = B[33:]-bc
y = X50[33:]
dy = dX50[33:]
#valori che mi aspetto per i parametri ottimali
init = np.array([0.00352, -7/4, -0.013]) #aiutano la convergenza del fit

Title = 'Suscettivita reticolo 40X40 dopo il punto critico'
ylabel = r'$chi [a.u.]$'
xlabel = r'$\beta - \beta_c$ [u.a.]$'

pars1, dpars1 = grafici.fit(F, x, y, dy, init, 7, Title, xlabel, ylabel)

print('\n')
print('stima di beta')

#seleziono solo un range dei dati
x = B[31:38] - bc
y = M50[31:38]
dy = dM50[31:38]

#valori che mi aspetto per i parametri ottimali
init = np.array([0.93, 1/8, 0.64]) #aiutano la convergenza del fit

Title = 'Magnetizzazione reticolo 50X50 dopo il punto critico'
ylabel = r'M [a.u.]'
xlabel = r'$\beta - \beta_c$ [u.a.]$'

pars2, dpars2 = grafici.fit(F, x, y, dy, init, 8, Title, xlabel, ylabel)

##stima di alpha
print('\n')
print('stima di alpha')

def Fa(x, m, g):
    '''funzione modello per i fit successivi
    '''
    return m*x**g 
    
CM, dCM = MS(C10, C20, C30, C40, C50, dC10, dC20, dC30, dC40, dC50)
N = np.arange(10, 51, 10)


#valori che mi aspetto per i parametri ottimali
init = np.array([0.1, 0]) #aiutano la convergenza del fit

Title = 'Massimo del calore specifico al variare di L'
ylabel = r'Calore specifico [a.u.]'
xlabel = r'$\beta - \beta_c$ [u.a.]$'

pars3, dpars3 = grafici.fit(Fa, N, CM, dCM, init, 9, Title, xlabel, ylabel)

#caclolo deigli indici critici ed associati errori


g = -pars1[1]
dg = dpars[1]

n = -pars1[1]/pars[1]
dn = np.sqrt((dpars[1]/pars[1])**2 + (dpars[1]/pars1[1])**2)*n

b = pars2[1]
db = dpars2[1]

a = pars3[1]
da = dpars3[1]


print('\n')
print(f"gamma = {g:.5f} +- {dg:.5f}")
print(f"beta  = {b:.5f} +- {db:.5f}")
print(f"nu    = {n:.5f} +- {dn:.5f}")
print(f"alpha = {a:.5f} +- {da:.5f}")

print('dg/nu = %.5f' %(pars[1]-7/4))
print('db = %.5f' %(pars2[1]-1/8))
print('dg = %.5f' %(pars1[1]+7/4))
print('da = %.5f' %(pars3[1]-0))

##grafici del finite size scaling

Title = 'Finite size scaling della magnetizzazione'
xlabel = r'$(\beta-\beta_c)L^{1/ \nu}$'
ylabel = r'$|M|/L^{-b/ \nu}$'
    
grafici.FSS(B, 1/n, -b/n, bc, 
            [M20, M30, M40, M50],
            [dM20, dM30, dM40, dM50], 
            [20, 30, 40, 50], 10, Title, xlabel, ylabel)

Title = 'Finite size scaling della suscettività'
xlabel = r'$(\beta-\beta_c)L^{1/ \nu}$'
ylabel = r'$ \chi /L^{\gamma/ \nu}$'

grafici.FSS(B, 1/n, g/n, bc, 
            [X20, X30, X40, X50],
            [dX20, dX30, dX40, dX50],
            [20, 30, 40, 50], 11, Title, xlabel, ylabel) 

plt.show()
