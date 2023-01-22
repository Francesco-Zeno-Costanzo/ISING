      program honeycomb_ising
      
C===============================================================
C Il codice simula un modello di ising su reticolo esagone con
C accoppiamenti fino a terzi vicini; l'hamiltoniana è:
C
C H = - J_0*sum(s_i*s_j) + J_1*sum(s_i*s_j)  - J_2*sum(s_i*s_j)
C           primi vicini       secondi vicini      terzi vicini
C
C Le costanti di accoppiamento sia assumono essere positive.
C Facciamo un esempio per chiarire, supponendo L=5 e scegliando
C un sito a caso del reticolo (e.g. 17) si avrà la struttura:
C
C                25          27
C          24          26          28
C
C          15          17          19
C                16          18
C          
C                7           9
C                      8
C
C I primi   vicini di 17 sono: 16, 18 e 26
C I secondi vicini di 17 sono: 7, 9, 25, 27, 15 e 19
C I terzi   vicini di 17 sono: 8, 24 e 28
C===============================================================	
  
      include "parameter.f"
C===============================================================
C file che permette di cambiare più facilmente i parametri della
C simulazione facendolo una sola volta piuttosto che diverse
C volte all'interno del codice 
C===============================================================
      
      character :: cr  !variaboile per il caricamento
      cr = char(13)    !serve per sovrascivere sulla shell

      call cpu_time(start)
      call ranstart

      open(1, file='init.txt',status='old')		      !file coi parametri
      open(2, file='anti_dati/dati10.dat',status='unknown')	!file coi risultati
	
      read(1,*) misure           !numero di misure
      read(1,*) i_dec            !updating fra una misura e l'altra
      read(1,*) bext      	   !valore del campo esterno
      read(1,*) bmin 		   !temperatura inversa minima
      read(1,*) bmax		   !temperatura inversa missima
      read(1,*) npassi	         !numero di temperature
      read(1,*) J_0              !costante di accoppiamento primi   vicini
      read(1,*) J_1              !costante di accoppiamento secondi vicini
      read(1,*) J_2              !costante di accoppiamento terzi   vicini

      write(2, *) misure	   !scrivo quantità che seriviranno
      write(2, *) nvol	         !nell'analisi
      write(2, *) npassi
     
      call bordi()		   !condizioni di bordo continuo
      call init()		         !inizializzo la matrice
	
      !la scrittura avviene su un unico file che poi verrà letto a blocchi
      !ogni blocco corrisponde ad una temperatura diversa
      !ci sono npassi blocchi, ciascuno lungo misure
      
      write(*,*) 'Taglia del reticolo: ',nlatt
	
      do k = 1, npassi		 !ciclo sulle temperature
	
          beta = bmin + (k-1)*(bmax-bmin)/float(npassi-1)
          write(*, '(A1)', advance='no')  cr !scrivo cr necessario per pulire la shell
		
          do i = 1, misure 	               !ciclo sulle misure a fissa temperatura

              do j = 1, i_dec
                  call metropolis(beta)      !decorrela la matrice  
              enddo

              call magnetizzazione(mag)	   !misurazione delle osservabili
              call energia(ene)

              write(2,*) abs(mag), ene	         !salvataggio risultati

          enddo

          !Sezione di caricamento, scrivo percentuale
          write(*,'(f8.1, a2)', advance='NO') k/float(npassi)*100, " %"
          flush(6) !pulisco la riga

      enddo
	
      call ranfinish
	
      call cpu_time(finish)
      
      elapsed = finish - start            !elapsed time in secodni
      i_h = elapsed/3600                  !ore trascorse
      i_m = modulo(int(elapsed/60), 60)   !minuti trascorsi
      i_s = modulo(int(elapsed), 60)      !secondi trascorsi
      
      write(*,*)
      print*, 'Tempo di esecuzione :'
      print '(I3," ore ", I3," min", I3, " sec ")', i_h, i_m, i_s
	
      end program honeycomb_ising


C============================================================================
C Subroutine per condizioni al bordo per reticolo esagonale
C============================================================================

      
      subroutine bordi() 		
	   
      include "parameter.f"
      
      !primi vicini
      if (J_0 /= 0) then
      	
          do i = 1, nvol
      
              nn1(i) = i + 1		                  !codizioni per ogni sito
              nn2(i) = i - 1
      	  
              if (modulo(i, 2) == 0) then             !Condizioni che regolano
      	      nn3(i) = i - 2*nlatt + 1            !l'esagonalità del reticolo.
              else                                    !Il reticolo è un array ma le
      	      nn3(i) = i + 2*nlatt - 1            !interazioni con questi vicini
              endif                                   !lo mappano in un reticolo esagonale
      	  
              if ((nn3(i) < 0).or.(nn3(i)>nvol)) then !condizioni, periodico-elicoidali
                  nn3(i) = modulo(nn3(i), nvol)	      !sopra e sotto identificati, mentre
              endif                                   !la riga i-esima si attacca alla i+1 -esima
      	  
          enddo
      
          nn1(nvol) = 1			              !problemi agli estremi
          nn2(1) = nvol
      endif
      
      !secondi vicini
      if (J_1 /= 0) then
      	
          do i = 1, nvol
          
              nnn1(i) = i - 2*nlatt
              nnn2(i) = i - 2*nlatt + 2
              nnn3(i) = i + 2*nlatt - 2
              nnn4(i) = i + 2*nlatt
              nnn5(i) = i + 2
              nnn6(i) = i - 2
              
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
      
      !terzi vicini
      if (J_2 /= 0) then
      	
          do i = 1, nvol
          
              if (modulo(i, 2) == 0) then
      	      nnnn1(i) = i + 2*nlatt - 1
      	      nnnn2(i) = i - 2*nlatt - 1
      	      nnnn3(i) = i - 2*nlatt + 3
              else
      	      nnnn1(i) = i - 2*nlatt + 1
      	      nnnn2(i) = i + 2*nlatt + 1
      	      nnnn3(i) = i + 2*nlatt - 3
              endif
              
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
C Inizializzazione del reticolo a caldo
C============================================================================
      	
      subroutine init()
      	
      include "parameter.f"
      	
      do i = 1, nvol           !iniazzializzazione del reticolo
          x = ran2()           !idealmente corrispone ad avere T infinita
          if(x<0.5) then       !perchè il varole +-1 del sito è causale
              campo(i) = 1
          else
              campo(i) = -1
          endif
      enddo
    
      return
      end


C============================================================================
C Subroutine che implementa il metropolis
C============================================================================

      
      subroutine metropolis(beta)
	
      include "parameter.f"
      
      do j = 1, nvol				      !ciclo su tutti i siti
	    
	    F4 = 0
          i = int(nvol*ran2()+1)	            !scelta random di un sito
	    
	    if (J_0 /= 0) then
              
              i1 = nn1(i)			      !calcolo dei primi vicini
              i2 = nn2(i)
              i3 = nn3(i)
              F1 = campo(i1) + campo(i2) + campo(i3)
              F4 = F4 -J_0*F1
              
          endif
          if (J_1 /= 0) then
          
              l1 = nnn1(i)			      !calcolo dei secondi vicini
              l2 = nnn2(i)
              l3 = nnn3(i)
              l4 = nnn4(i)
              l5 = nnn5(i)
              l6 = nnn6(i)
              F2 = campo(l1) + campo(l2) + campo(l3) + campo(l4) +
     &             campo(l5) + campo(l6)
              F4 = F4 + J_1*F2
     
          endif
          if (J_2 /= 0) then
          
              k1 = nnnn1(i)			      !calcolo dei terzi vicini
              k2 = nnnn2(i)
              k3 = nnnn3(i)
              F3 = campo(k1) + campo(k2) + campo(k3)
              F4 = F4 - J_2*F3
              
          endif
          
          F = beta*(F4 + bext)		            !aggiungo eventuale campo esterno
     		
          ispin = campo(i)
     		
          p =  exp(-2.0*ispin*F)		!probabilità di accettare la mossa
     		
          x = ran2()			      !numero casuale per il test
     		
          if(x < p) then			!test di accettanza
              campo(i) = -ispin	      !se F è negativo il test è
     	    endif				      !passato di default
	
      enddo
      
      return
      end
	

C============================================================================
C Subroutine misura della magnetizzazione
C============================================================================


      subroutine magnetizzazione(mag)
	
      include "parameter.f"
      
      mag = 0			  !inizializzo la variabile
      
      if (J_1 <= 0) then 
          do i = 1, nvol		  !ciclo su tutto il reticolo
      
              mag = mag + campo(i)  !e sommo ogni sito

          enddo
      else 
          do i = 1, nvol		  !ciclo su tutto il reticolo
      
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
      
      ene = 0					      !inzializzo la variabile
      
      do i = 1, nvol				      !ciclo su tutti i siti

          F4 = 0
	    
	    if (J_0 /= 0) then
              
              i1 = nn1(i)			      !calcolo dei primi vicini
              i2 = nn2(i)
              i3 = nn3(i)
              F1 = campo(i1) + campo(i2) + campo(i3)
              F4 = F4 - J_0*F1
              
          endif
          if (J_1 /= 0) then
          
              l1 = nnn1(i)			      !calcolo dei secondi vicini
              l2 = nnn2(i)
              l3 = nnn3(i)
              l4 = nnn4(i)
              l5 = nnn5(i)
              l6 = nnn6(i)
              F2 = campo(l1) + campo(l2) + campo(l3) + campo(l4) +
     &             campo(l5) + campo(l6)
              F4 = F4 + J_1*F2
     
          endif
          if (J_2 /= 0) then
          
              k1 = nnnn1(i)			      !calcolo dei terzi vicini
              k2 = nnnn2(i)
              k3 = nnnn3(i)
              F3 = campo(k1) + campo(k2) + campo(k3)
              F4 = F4 - J_2*F3
              
          endif
 
          ene = ene - 0.5*F4*campo(i)	!0.5 per non sovracontare
          ene = ene - bext*campo(i)      	!eventuale campo esterno

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

