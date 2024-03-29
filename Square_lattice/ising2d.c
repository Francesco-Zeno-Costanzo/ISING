#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<math.h>
#include<time.h>

#define L       50
#define J       1.0
#define H       0.0
#define step    70
#define Bm      0.35
#define BM      0.5
#define decorr  100
#define N       10000

int Reticolo[L][L]; // Rappresentazione del reticolo bidimensionale

// Funzione per inizializzare il reticolo con spin casuali (+1 o -1)
void init (void) {
    for(int i=0;i<L;i++) {
        for(int j=0;j<L; j++) {
            Reticolo[i][j] = ((rand()%2)*2)-1; // Assegna casualmente +1 o -1
        }
    }
}

// Implementazione dell'algoritmo di Metropolis per l'evoluzione del reticolo in una data temperatura B
void Metropolis (double B) {
    int c=0;
    int b=0;
    double dE=0;
    for(int i=0; i<L; i++){
        for(int j=0; j<L; j++){
            c = rand()%L;
            b = rand()%L;
            dE=2*J*(Reticolo[c][b])*( Reticolo[c][(b+1)%L] + Reticolo[(c+1)%L][b]
            + Reticolo[c[(b+L-1)%L] + Reticolo[(c+L-1)%L][b]) + H*(Reticolo[c][b]);
            if(dE<0.0){
                Reticolo[c][b]=Reticolo[c][b]*(-1); // Accetta mossa
            }
            else{
                if((rand()/(RAND_MAX + 1.0))<exp(-dE*B)){
                    Reticolo[c][b]=Reticolo[c][b]*(-1); // Accetta mossa con probabilità di Boltzmann
                }
            }

        }
    }
}

// Calcola l'energia totale del sistema
double Energia (void) {
    double E = 0;
    for(int i=0; i<L; i++){
        for(int j=0; j<L; j++){
            E += - J*(Reticolo[i][j])*( Reticolo[i][(j+1)%L] + Reticolo[(i+1)%L][j]
                 + Reticolo[i][(j+L-1)%L] + Reticolo[(i+L-1)%L][j])*0.5 - H*(Reticolo[i][j]);
        }
    }
    return E;
}

// Calcola la magnetizzazione totale del sistema
double Magnetizzazione(void){
    double mag = 0;
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            mag += Reticolo[i][j] ;
        }
    }
    return mag;
}

// Funzione principale
int main (void){
    clock_t start = clock(); // Misura il tempo di esecuzione
    srand(time(NULL)); // Inizializza il generatore di numeri casuali con il tempo corrente

    // Dichiarazione e inizializzazione di variabili
    double B[step]={0}; // Temperatura
    double E[step]={0}; // Energia media
    double C[step]={0}; // Capacità termica
    double M[step]={0}; // Magnetizzazione media
    double X[step]={0}; // Suscettibilità magnetica
    double En = 0;      // Energia
    double Mag = 0;     // Magnetizzazione
    double En2 = 0;     // Energia quadratica media
    double Mag2 = 0;    // Magnetizzazione quadratica media
    double Ene = 0;     // Energia temporanea
    double mag = 0;     // Magnetizzazione temporanea

    double n1=1/(N*L*L*1.0); // Costante di normalizzazione
    double n2=n1/(N*1.0);    // Costante di normalizzazione

    init(); // Inizializza il reticolo

    // Ciclo su diverse temperature
    for(int j=0; j<step; j++){
        B[j]=Bm + j*((BM - Bm)/(step*1.0)); // Assegna valore alla temperatura
        En=Mag=En2=Mag2=Ene=mag=0;

        // Simulazione Monte Carlo
        for(int t=0; t<N; t++){
            for(int d=0; d<decorr; d++){
                Metropolis(B[j]); // Evoluzione del reticolo
            }
            Ene = Energia(); // Calcola l'energia
            mag = Magnetizzazione(); // Calcola la magnetizzazione

            En += Ene;
            if(mag>0){
                Mag += mag;
            }
            else{
                Mag += -mag;
            }
            En2  += Ene*Ene;
            Mag2 += mag*mag;
        }

        // Calcola le grandezze termodinamiche
        E[j] = n1*En;
        M[j] = n1*Mag;
        C[j] = (n1*En2  - n2*(En*En));
        X[j] = (n1*Mag2 - n2*(Mag*Mag));
    }

    // Salva i risultati su file
    FILE *fd;
    fd=fopen("ising50.txt", "w");
    if( fd==NULL ) {
        perror("Errore in apertura del file");
      }
    for(int q=0; q<step; q++){
        fprintf(fd, "%f \t %f \t %f \t %f \t %f \n", B[q], E[q], M[q], C[q], X[q]);
    }
      fclose(fd);

      clock_t end = clock();
      printf("Tempo di esecuzione =  %f secondi \n", ((double)(end - start)) / CLOCKS_PER_SEC); // Stampa il tempo di esecuzione
    return 0;
}
