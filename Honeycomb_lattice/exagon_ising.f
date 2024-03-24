      program honeycomb_ising
      
C===================================================================================
C     Il codice simula un modello di ising calssico su reticolo esagone con
C     accoppiamenti fino a terzi vicini; l'hamiltoniana del sistema è:
C
C     H = - J_0*sum(s_i*s_j) - J_1*sum(s_i*s_j)  - J_2*sum(s_i*s_j) + bext*sum(s_i)
C               primi vicini       secondi vicini      terzi vicini   campo esterno
C
C     Facciamo un esempio per chiarire, supponendo L=5 e scegliando
C     un sito a caso del reticolo (e.g. 17) si avrà la struttura:
C
C                   25          27
C             24          26          28
C
C             15          17          19
C                   16          18
C          
C                   7           9
C                         8
C
C     I primi   vicini di 17 sono: 16, 18 e 26
C     I secondi vicini di 17 sono: 7, 9, 25, 27, 15 e 19
C     I terzi   vicini di 17 sono: 8, 24 e 28
C
C     Per funzionare il codice neccessita di alcuni file:
C     1) Un file chiamto randomseed, contente almeno un numero, necassario per
C        l'utilizzo di ran2(). La subroutine ranfinish lo aggiornerà poi, quindi
C        creato una volta non ce ne si deve più preoccupare
C
C     2) Un file chiamto init.txt contente i parametri della simulazione.
C        dai commenti sottostanti si evince cosa siano queste varibili,
C        sottolineamo però che i_start se:
C        1) Vale 0 si ha partenza caldo (i.e. siti random)
C        2) Vale 1 si ha partenza a freddo (i.e. tutti i siti = 1)
C        ( scegliere se 0 o 1 può essere utile se si utilizza beta o T )
C
C===================================================================================

      include "parameter.f"

C===================================================================
C     File che permette di cambiare più facilmente i parametri della
C     simulazione facendolo una sola volta piuttosto che diverse
C     volte all'interno del codice
C===================================================================

      character :: cr  ! Variaboile per il caricamento
      cr = char(13)    ! serve per sovrascivere sulla shell

      call cpu_time(start)
      call ranstart

      open(1, file='init.txt',status='old')                 ! File coi parametri
      open(2, file='data/data_L_10_J2_01.dat',status='unknown')   ! File coi risultati

      read(1,*) misure      ! Numero di misure
      read(1,*) i_dec       ! Updating fra una misura e l'altra
      read(1,*) bext        ! Valore del campo esterno
      read(1,*) bmin        ! Temperatura inversa minima
      read(1,*) bmax        ! Temperatura inversa missima
      read(1,*) npassi      ! Numero di temperature
      read(1,*) J_1         ! Costante di accoppiamento primi   vicini
      read(1,*) J_2         ! Costante di accoppiamento secondi vicini
      read(1,*) J_3         ! Costante di accoppiamento terzi   vicini
      read(1,*) i_start     ! Flag per decidere inizializzazione del reticolo

      write(2, *) misure    ! Scrivo quantità che seriviranno nell'analisi
      write(2, *) nvol
      write(2, *) npassi

      call bordi()          ! Condizioni di bordo continuo
      call init(i_start)    ! Inizializzo il sistema

      ! La scrittura avviene su un unico file che poi verrà letto a blocchi.
      ! Ogni blocco corrisponde ad una temperatura diversa
      ! ci sono npassi blocchi, ciascuno lungo misure

      write(*,*) 'Taglia del reticolo: ',nlatt
      write(*,*) 'J_1, J_2, J_3: ', J_1, J_2, J_3

      do k = 1, npassi  ! Ciclo sulle temperature

          beta = bmin + (k-1)*(bmax-bmin)/float(npassi-1)

          write(*, '(A1)', advance='no')  cr  ! Scrivo cr necessario per pulire la shell

          do i = 1, misure                    ! Ciclo sulle misure a fissa temperatura

              do j = 1, i_dec
                  call metropolis(beta)       ! Decorrela il sistema
              enddo

              call magnetizzazione(mag)       ! Misurazione delle osservabili
              call energia(ene)

              write(2,*) abs(mag), ene        ! Salvataggio risultati

          enddo

          !Sezione di caricamento, scrivo percentuale
          write(*,'(f8.1, a2)', advance='NO') k/float(npassi)*100," %"
          flush(6) !pulisco la riga

      enddo

      call ranfinish

      call cpu_time(finish)
      
      elapsed = finish - start            ! Elapsed time in secodni
      i_h = elapsed/3600                  ! Ore trascorse
      i_m = modulo(int(elapsed/60), 60)   ! Minuti trascorsi
      i_s = modulo(int(elapsed), 60)      ! Secondi trascorsi
      
      write(*,*)
      print*, 'Tempo di esecuzione :'
      print '(I3," ore ", I3," min", I3, " sec ")', i_h, i_m, i_s

      end program honeycomb_ising

C============================================================================
C Subroutine per condizioni al bordo per reticolo esagonale
C============================================================================
      
      subroutine bordi()

      include "parameter.f"

      !*************************** Primi vicini ***************************
      if (J_1 /= 0) then

          do i = 1, nvol

              nn1(i) = i + 1                          ! Codizioni per ogni sito
              nn2(i) = i - 1

              if (modulo(i, 2) == 0) then             ! Condizioni che regolano
                  nn3(i) = i - 2*nlatt + 1            ! l'esagonalità del reticolo.
              else                                    ! Il reticolo è un array ma le
                  nn3(i) = i + 2*nlatt - 1            ! interazioni con questi vicini
              endif                                   ! lo mappano in un reticolo esagonale

              if ((nn3(i) < 0).or.(nn3(i)>nvol)) then ! Condizioni, periodico-elicoidali
                  nn3(i) = modulo(nn3(i), nvol)	      ! sopra e sotto identificati, mentre
              endif                                   ! la riga i-esima si attacca alla i+1-esima

          enddo
          nn1(nvol) = 1                               ! Problemi agli estremi
          nn2(1) = nvol
      endif

      ! *************************** Secondi vicini ***************************
      if (J_2 /= 0) then

          do i = 1, nvol

              ! Dentro il reticolo
              nnn1(i) = i - 2*nlatt
              nnn2(i) = i - 2*nlatt + 2
              nnn3(i) = i + 2*nlatt - 2
              nnn4(i) = i + 2*nlatt
              nnn5(i) = i + 2
              nnn6(i) = i - 2

              ! Bordi del reticolo
              if ((nnn1(i) < 0).or.(nnn1(i)>nvol)) then
                  nnn1(i) = modulo(nnn1(i), nvol)
              endif

              if ((nnn2(i) < 0).or.(nnn2(i)>nvol)) then   
                  nnn2(i) = modulo(nnn2(i), nvol)
              endif

              if ((nnn3(i) < 0).or.(nnn3(i)>nvol)) then   
                  nnn3(i) = modulo(nnn3(i), nvol)
              endif

              if ((nnn4(i) < 0).or.(nnn4(i)>nvol)) then   
                  nnn4(i) = modulo(nnn4(i), nvol)
              endif

              if ((nnn5(i) < 0).or.(nnn5(i)>nvol)) then   
                  nnn5(i) = modulo(nnn5(i), nvol)
              endif

              if ((nnn6(i) < 0).or.(nnn6(i)>nvol)) then   
                  nnn6(i) = modulo(nnn6(i), nvol)
              endif
          enddo
      endif

      !*************************** Terzi vicini ***************************
      if (J_3 /= 0) then

          do i = 1, nvol

              ! Dentro il reticolo
              if (modulo(i, 2) == 0) then
              nnnn1(i) = i + 2*nlatt - 1
              nnnn2(i) = i - 2*nlatt - 1
              nnnn3(i) = i - 2*nlatt + 3
              else
              nnnn1(i) = i - 2*nlatt + 1
              nnnn2(i) = i + 2*nlatt + 1
              nnnn3(i) = i + 2*nlatt - 3
              endif

              ! Bordo del reticolo
              if ((nnnn1(i) < 0).or.(nnnn1(i)>nvol)) then  
                  nnnn1(i) = modulo(nnnn1(i), nvol)
              endif

              if ((nnnn2(i) < 0).or.(nnnn2(i)>nvol)) then  
                  nnnn2(i) = modulo(nnnn2(i), nvol)
              endif

              if ((nnnn3(i) < 0).or.(nnnn3(i)>nvol)) then  
                  nnnn3(i) = modulo(nnnn3(i), nvol)
              endif
          enddo
      endif

      return
      end
 
 
C============================================================================
C Inizializzazione del reticolo
C============================================================================

      subroutine init(i_start)
      include "parameter.f"

      if (i_start == 0) then
          do i = 1, nvol           ! Iniazzializzazione del reticolo.
              x = ran2()           ! Idealmente corrisponde ad avere T infinita
              if(x<0.5) then       ! perchè il varole +-1 del sito è causale
                  campo(i) = 1
              else
                  campo(i) = -1
              endif
          enddo
      elseif (i_start == 1) then
          do i = 1, nvol           ! Iniazzializzazione del reticolo.
              campo(i) = 1         ! Tutti gli spin a 1 come se fosse T = 0
          enddo
      endif
      return
      end

C============================================================================
C Subroutine che implementa il metropolis
C============================================================================

      subroutine metropolis(beta)

      include "parameter.f"
      
      do i = 1, nvol                 ! Ciclo su tutti i siti

          F4 = 0                     ! Initialize "force"

          if (J_1 /= 0) then

              F1 = campo(nn1(i)) + campo(nn2(i)) + campo(nn3(i))
              F4 = F4 - J_1*F1

          endif
          if (J_2 /= 0) then

              F2 = campo(nnn1(i)) + campo(nnn2(i)) + campo(nnn3(i)) +
     &             campo(nnn4(i)) + campo(nnn5(i)) + campo(nnn6(i))
              F4 = F4 - J_2*F2

          endif
          if (J_3 /= 0) then

              F3 = campo(nnnn1(i)) + campo(nnnn2(i)) + campo(nnnn3(i))
              F4 = F4 - J_3*F3

          endif

          F = (F4 + bext)/beta       ! Aggiungo eventuale campo esterno

          p = exp(-2.0*campo(i)*F)   ! Probabilità di accettare la mossa

          x = ran2()                 ! Numero casuale per il test

          if(x < p) then             ! Test di accettanza
              campo(i) = -campo(i)   ! se F è negativo il test è
          endif                      ! passato di default
      enddo
      return
      end

C============================================================================
C Subroutine misura della magnetizzazione
C============================================================================

      subroutine magnetizzazione(mag)

      include "parameter.f"

      mag = 0                 ! Inizializzo la variabile

      if (J_2 >= 0) then
          do i = 1, nvol      ! Ciclo su tutto il reticolo

              mag = mag + campo(i)  ! e sommo ogni sito

          enddo
      else 
          do i = 1, nvol      ! Ciclo su tutto il reticolo

              mag = mag + ((-1)**i)*campo(i)  !e sommo ogni sito

          enddo
      endif
      return
      end


C============================================================================
C Subroutine misura dell'energia
C============================================================================

      subroutine energia(ene)

      include "parameter.f"

      ene = 0                      ! Inzializzo la variabile

      do i = 1, nvol               ! Ciclo su tutti i siti
          F4 = 0
          if (J_1 /= 0) then

              F1 = campo(nn1(i)) + campo(nn2(i)) + campo(nn3(i))
              F4 = F4 - J_1*F1

          endif
          if (J_2 /= 0) then

              F2 = campo(nnn1(i)) + campo(nnn2(i)) + campo(nnn3(i)) +
     &             campo(nnn4(i)) + campo(nnn5(i)) + campo(nnn6(i))
              F4 = F4 - J_2*F2

          endif
          if (J_3 /= 0) then

              F3 = campo(nnnn1(i)) + campo(nnnn2(i)) + campo(nnnn3(i))
              F4 = F4 - J_3*F3

          endif
 
          ene = ene - 0.5*F4*campo(i) ! 0.5 per non sovracontare
          ene = ene - bext*campo(i)   ! Eventuale campo esterno

      enddo
      return
      end

c============================================================================
c  RANDOM NUMBER GENERATOR: standard ran2 from numerical recipes
c============================================================================
      function ran2()
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      real ran2,am,eps,rnmx
      parameter(im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1,
     &          ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,
     &          ir2=3791,ntab=32,ndiv=1+imm1/ntab,eps=1.2e-7,
     &          rnmx=1.-eps)
      integer idum2,j,k,iv,iy
      common /dasav/ idum,idum2,iv(ntab),iy
c      save iv,iy,idum2
c      data idum2/123456789/, iv/NTAB*0/, iy/0/

      if(idum.le.0) then
         idum=max0(-idum,1)
         idum2=idum
         do j=ntab+8,1,-1
            k=idum/iq1
            idum=ia1*(idum-k*iq1)-k*ir1
            if(idum.lt.0) idum=idum+im1
            if(j.le.ntab) iv(j)=idum
         enddo
         iy=iv(1)
      endif
      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
      if(idum.lt.0) idum=idum+im1
      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
      if(idum2.lt.0) idum2=idum2+im2
      j=1+iy/ndiv
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1) iy=iy+imm1
      ran2=min(am*iy,rnmx)

      return
      end

c=============================================================================
      subroutine ranstart
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      common /dasav/ idum,idum2,iv(32),iy

      open(unit=23, file='randomseed', status='unknown')
      read(23,*) idum
      read(23,*,end=117) idum2
      do i=1,32
         read(23,*) iv(i)
      enddo
      read(23,*) iy
      close(23)
      goto 118                          !!takes account of the first start
 117  if(idum.ge.0) idum = -idum -1     !!
      close(23)
 118  continue                          !!

      return
      end

c=============================================================================
      subroutine ranfinish
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      common /dasav/ idum,idum2,iv(32),iy

      open(unit=23, file='randomseed', status='unknown')
      write(23,*) idum
      write(23,*) idum2
      do i=1,32
         write(23,*) iv(i)
      enddo
      write(23,*) iy
      close(23)

      return
      end
c=============================================================================

