import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
       
def plot(x, args, d_args, L, str_leg, k, title, xlabel, ylabel):
    """
    Funzione che esegue dei plot di dati
    su uno stesso grafico
    
    Parameters
    ----------
    x : 1darray
        array sulle x, uguali per tutti
    args : list
        lista degli array sulle y
    dargs : lista
        lista degli errori sulle y
    str_leg : string
        string of legend
    L : lista
        lista che la bella i dati, da mettere nella legenda
    k : int
        numero del grafico
    title : string
        titolo del grafico
    xlabel : string
        label asse sulle x
    ylabel : string
        label asse sulle y
    """
    plt.figure(k)
    plt.title(title)
    plt.xlabel(xlabel, fontsize=15)
    plt.ylabel(ylabel, fontsize=15)
    plt.grid()
  
    colors = plt.cm.jet(np.linspace(0, 1, len(args)))

    for p, dp, i in zip(args, d_args, range(len(args))):

        plt.errorbar(x, p, dp, fmt='.', linestyle='--', color=colors[i], label=f'{str_leg}={L[i]}')
                   
    plt.legend(loc='best')


def plotbinder(B, args, d_args, L, k, title, xlabel, ylabel):
    """
    Funzione che esegue il plot per il cumulante
    di binder su uno stesso grafico con un sottografico
    per evidenziare meglio l'intersezione delle curve
    
    Parameters
    ----------
    B : 1darray
        array sulle x, uguali per tutti
    args : list
        lista degli array sulle y
    dargs : lista
        lista degli errory sulle y
    L : lista
        lista che la bella i dati, da mettere nella legenda
    k : int
        numero del grafico
    title : string
        titolo del grafico
    xlabel : string
        label asse sulle x
    ylabel : string
        label asse sulle y
    """
    fig = plt.figure(k)
    main_ax = fig.add_subplot() #creo le variabili per il garfico
    main_ax.set_title(title)
    main_ax.set_xlabel(xlabel, fontsize=15)
    main_ax.set_ylabel(ylabel, fontsize=15)
    main_ax.grid()
    
    right_inset_ax = fig.add_axes([.55, .5, .3, .3])
    right_inset_ax.grid()
    right_inset_ax.set_xlim(0.645, 0.675)
    right_inset_ax.set_ylim(0.99, 1.4)
    
    colors = plt.cm.jet(np.linspace(0, 1, len(args)))
    
    for p, dp, i in zip(args, d_args, range(len(args))):

        main_ax.errorbar(B, p, dp, fmt='.', color=colors[i], label=f'L={L[i]}')
        right_inset_ax.errorbar(B, p, dp, fmt='.', linestyle='--', color=colors[i], label=f'L={L[i]}')
        
    main_ax.legend(loc='best')
    right_inset_ax.legend(loc='best')
    

def FSS(B, q1, q2, bc, args, d_args, L, k, title, xlabel, ylabel):
    """
    Funzione che esegue dei plot di dati
    su uno stesso grafico per il finite size scaling
    
    Parameters
    ----------
    B : 1darray
        array sulle x, uguali per tutti
    q1 : float
        indice critico per scaling sulle x
    q2 : float
        indice critico per scaling sulle y
    bc : float
        temperatura critica del sistema
    args : list
        lista degli array sulle y
    dargs : lista
        lista degli errory sulle y
    L : lista
        lista che la bella i dati, da mettere nella legenda
    k : int
        numero del grafico
    title : string
        titolo del grafico
    xlabel : string
        label asse sulle x
    ylabel : string
        label asse sulle y
    """     
    plt.figure(k)
    plt.title(title)
    plt.xlabel(xlabel, fontsize=15)
    plt.ylabel(ylabel, fontsize=15)
    plt.grid()
    
    for p, dp, i in zip(args, d_args, range(len(args))):

        plt.errorbar((B-bc)*L[i]**(q1), p/L[i]**q2, dp/L[i]**q2, fmt='.', color='k', label=f'L={L[i]}')
        
    plt.legend(loc='best')
    
    
def fit(F, x, y, dy, init, k, Title, xlabel, ylabel, plot=True):
    """
    Funzione che esegue un fit
    
    Parameters
    ----------
    F : callable
        funzione di fit
    x : 1darray
        array sulle x
    y : 1darray
        array sulle y
    dy : 1darray
        listra degli errory sulle y
    init : 1darray
        array di init per aiutare la convergenza del fit
    k : int
        numero del grafico
    title : string
        titolo del grafico
    xlabel : string
        label asse sulle x
    ylabel : string
        label asse sulle y
    plot : bool
        se True il grafico viene mostrato
    
    Return
    ------
    pars : 1darray
        array parametri ottimali
    err : 1darray
        array errori sui parametri
    """
    pars, covm = curve_fit(F, x, y, init, sigma=dy)  
    err = np.sqrt(covm.diagonal())
    for p, dp, i in zip(pars, err, range(len(pars))):
        print(f"pars{i} = {p:.5f} +- {dp:.5f}")
    

    chisq = sum(((y - F(x, *pars))/dy)**2.)
    ndof = len(y) - len(pars)
    print('chi quadro = %.3f (%d dof)' % (chisq, ndof))


    c = np.zeros((len(pars), len(pars)))

    for i in range(0, len(pars)):
        for j in range(0, len(pars)):
            c[i][j]=(covm[i][j])/(np.sqrt(covm.diagonal()[i])*np.sqrt(covm.diagonal()[j]))
    print(c) #matrice di correlazione

    if plot :
        fig = plt.figure(k)
        frame1 = fig.add_axes((.1,.35,.8,.6))
        frame1.set_title(Title, fontsize=20)
        frame1.set_ylabel(ylabel, fontsize=15)
        frame1.grid()

        frame1.errorbar(x, y, dy, fmt='.', color='black', label='dati')
        t = np.linspace(np.min(x), np.max(x), 1000)
        frame1.plot(t, F(t, *pars), color='blue', label='best fit')   
        frame1.legend(loc='best')
        frame2 = fig.add_axes((.1,.1,.8,.2))

        ff = (y - F(x, *pars))/dy
        frame2.set_ylabel('Residui Normalizzati')
        frame2.set_xlabel(xlabel, fontsize=15)

        frame2.plot(t, 0*t, color='red', linestyle='--', alpha=0.5) #grafico la retta costantemente zero
        frame2.plot(x, ff, '.', color='black') #grafico i residui normalizzati
        frame2.grid()
    
    return pars, err
    
    
def plotfit(F, Pars, X, args, d_args, L, k, title, xlabel, ylabel):
    """
    Funzione che esegue dei plot di fit
    su uno stesso grafico
    
    Parameters
    ----------
    F : callable
        funzione di fit
    Pars : lista
        lista dei parametri ottimali    
    X : lista
        lista degli array sulle x
    args : list
        lista degli array sulle y
    dargs : lista
        lista degli errori sulle y
    L : lista
        lista che la bella i dati, da mettere nella legenda
    k : int
        numero del grafico
    title : string
        titolo del grafico
    xlabel : string
        label asse sulle x
    ylabel : string
        label asse sulle y
        array di init per aiutare la convergenza del fit
    """
    fig = plt.figure(k)
    
    frame1 = fig.add_axes((.1,.35,.8,.6))
    frame1.set_title(title, fontsize=20)
    frame1.set_ylabel(ylabel, fontsize=15)
    frame1.grid()
    
    frame2 = fig.add_axes((.1,.1,.8,.2))
    frame2.set_ylabel('Residui Normalizzati')
    frame2.set_xlabel(xlabel, fontsize=15)
    frame2.grid()
    
    colors = plt.cm.jet(np.linspace(0, 1, len(args)))
    
    for x, y, dy, pars, i in zip(X, args, d_args, Pars, range(len(args))):
      
        frame1.errorbar(x, y, dy, fmt='.', color=colors[i], label=f'dati per L={L[i]}')
        t = np.linspace(np.min(x), np.max(x), 1000)
        frame1.plot(t, F(t, *pars), color=colors[i])   
        frame1.legend(loc='best')
        

        ff = (y - F(x, *pars))/dy
        

        frame2.plot(t, 0*t, color='red', linestyle='--', alpha=0.5) #grafico la retta costantemente zero
        frame2.plot(x, ff, '.', color=colors[i]) #grafico i residui normalizzati
        

