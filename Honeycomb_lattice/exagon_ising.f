      program honeycomb_ising
	      
      include "parameter.f"
C===============================================================
C file che permette di cambiare più facilmente i parametri della
C simulazione facendolo una sola volta piuttosto che diverse
C volte all'interno del codice 
C===============================================================

      call cpu_time(start)
      call ranstart

      open(1, file='init.txt',status='old')		!file coi parametri
      open(2, file='dati/dati15.dat',status='unknown')	!file coi risultati
	
      read(1,*) misure           !numero di misure
      read(1,*) i_dec            !updating fra una misura e l'altra
      read(1,*) bext      	 !valore del campo esterno
      read(1,*) bmin 		 !temperatura inversa minima
      read(1,*) bmax		 !temperatura inversa missima
      read(1,*) npassi	         !numero di temperature
      read(1,*) J_acc            !costante di accoppiamento
       
      write(2, *) misure	 !scrivo quantità che seriviranno
      write(2, *) nvol	         !nell'analisi
      write(2, *) npassi
     
      call bordi()		 !condizioni di bordo continuo
      call init()		 !inizializzo la matrice
	
      !la scrittura avviene su un unico file che poi verra letto a blocchi
      !ogni blocco corrisponde ad una temperatura diversa
      !ci sono npassi blocchi, ciascuno lungo misure
	
      do k = 1, npassi		 !ciclo sulle temperature
	
          beta = bmin + (k-1)*(bmax-bmin)/float(npassi-1)
		
	  do i = 1, misure 	!ciclo sulle misure a fissa temperatura

              do j = 1, i_dec
		  call metropolis(beta)   !decorrela la matrice
	      enddo
			
	      call magnetizzazione(mag)	!misurazione delle osservabili
	      call energia(ene)
	      
	      write(2,*) abs(mag), ene	!salvataggio risultati
	  
	  enddo
      enddo
	
      call ranfinish
	
      call cpu_time(finish)
      print '("tempo di esecuzione= ", f16.8," secondi.")', finish-start
	
      end program honeycomb_ising

C============================================================================

      subroutine bordi() 		
	   
      include "parameter.f"
      	
      do i = 0, nvol
      
          n1(i) = i + 1		!codizioni per ogni sito
      	  n2(i) = i - 1
      	  
      	  if (modulo(i, 2) == 0) then
      	      n3(i) = i + 2*nlatt - 1
      	  else
      	      n3(i) = i - 2*nlatt + 1
      	  endif
      	  
      	  if (n3(i) < 0) then
      	      n3(i) = nvol - n3(i)
      	  endif
      	  
      enddo
      
      n1(nvol - 1) = 0
      n2(0) = nvol - 1
      	
      return
      end
 
C============================================================================
      	
      subroutine init()
      	
      include "parameter.f"
      	
      do i = 0, nvol
          x=ran2()
      	  if(x<0.5) then
      	      campo(i) = 1
      	  else
      	      campo(i) = -1
      	  endif
      enddo
    
      return
      end

C============================================================================

      subroutine metropolis(beta)
	
      include "parameter.f"
      
      do j = 1, nvol				!ciclo su tutti i siti
		
          i = int(nvol*ran2())	                !scelta random di un sito
		
	  i1 = n1(i)			        !calcolo dei primi vicini
	  i2 = n2(i)
	  i3 = n3(i)
		
	  F = campo(i1) + campo(i2) + campo(i3) !sommo i primi vicini
     	  F = beta*(F + bext)		        !aggiungo eventuale campo esterno
     		
     	  ispin = campo(i)
     		
     	  p =  exp(-2.0*J_acc*ispin*F)		!probabilità di accettare la mossa
     		
     	  x = ran2()			        !numero casuale per il test
     		
     	  if(x < p) then			!test di accettanza
     	      campo(i) = -ispin	                !se F è negativo il test è
     	  endif				        !passato di default
	
      enddo
      return
      end
	

C============================================================================

      subroutine magnetizzazione(mag)
	
      include "parameter.f"
      
      mag = 0			!inzializzo la variabile
      
      do i = 0, nlatt		!ciclo su tutto il retico
      
          mag = mag + campo(i)  !e sommo ogni sito

      enddo
      
      return
      end

C============================================================================

      subroutine energia(ene)
	
      include "parameter.f"
      
      ene=0					!inzializzo la variabile
      do i = 1, nvol				!ciclo su tutti i siti
	
	  i1 = n1(i)			        !calcolo dei primi vicini
	  i2 = n2(i)
	  i3 = n3(i)
		
	  F = campo(i1) + campo(i2) + campo(i3) !sommo i primi vicini
 
     	  ene = ene - 0.5*J_acc*F*campo(i)	        !0.5 per non sovracontare
          ene = ene - bext*campo(i)      	!eventuale campo esterno

      enddo
	
      return
      end


C============================================================================


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

