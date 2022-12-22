# Simulazione

Prendiamo 100 mila misure, con indice di decorrelazione 10. Il range di temperature è poi 0.4-0.9, con 100 temperature.
| Reticolo | Elapsed time |                 
| -------- | ------------ |    
| L = 10   | eta: 11 min |
| L = 15   | eta: 23 min |
| L = 20   | eta: 44 min |
| L = 25   | eta: 1 h 5 min |
| L = 30   | eta: 1 h 30 min |
| L = 35   | eta: 2 h 20 min |
| L = 40   | eta: 3 h |
| L = 45   | eta: 3 h 45 min |
| L = 50   | eta: 4 h 37 min |

# Analisi

| Reticolo | Elapsed time |
| -------- | ------------ |
| L = 10   | eta: 26 s    |
| L = 15   | eta: 22 s    |
| L = 20   | eta: 21 s    |
| L = 25   | eta: 22 s    |
| L = 30   | eta: 22 s    |
| L = 35   | eta: 22 s    |
| L = 40   | eta: 21 s    |
| L = 45   | eta: 23 s    |
| L = 50   | eta: 22 s    |


# Analisi Python
I parametri che restituiscono gli indici critici nel programma sono i pars1

$\beta_c = 0.661(1)$

| Indici critici | Valore atteso | Valore calcolato |
| ----------- | ----------- | ----------- |
| $\beta$      |     $\frac18$ $\sim$ 0.125  | 0.139(5) |
| $\gamma$ | $\frac74$ $\sim$ 1.75 | 1.75(2) |
| $\nu$ | 1 | 1.00(2) |
| $\alpha$ | 0 | 0.01(2) |
