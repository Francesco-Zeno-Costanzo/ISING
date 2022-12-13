      program analisi
	
      real*8, dimension(:), allocatable :: Ene, Mag
      real*8:: daver_e, daver_m, daver_c, daver_x, dcb
      real*8:: aver_e, aver_m, aver_c, aver_x, cb
      integer :: R
      
      common N, R, nvol
      call cpu_time(start)
      call ranstart
	
      !apro file da cui leggere i dati da analizzare e file su
      !cui scrivere i risultati con relativi errori
      open(unit=0, file="dati/dati15.dat", status="old", action="read")
      open(unit=2, file='datiplot/dati15.dat',status='unknown')
	
      R = 100				!numero di ricampionamenti
      
      read(0, *) N			!leggo i primi valori che sono 
      read(0, *) nvol			!necessari per l'analisi
      read(0, *) npassi
	
      allocate(Ene(N), Mag(N))

      do j=1, npassi			!leggo a blocchi il file
          do i = 1, N		        !ogni blocco una temperatura diversa
	      read(0, *) Mag(i), Ene(i)
	  enddo
	  
	  aver_e = sum(Ene)/float(N*nvol) 	        !Energia media
	  aver_m = sum(Mag)/float(N*nvol)		!magnetizzazione media
		
	  c = (sum(Ene**2)/float(N) - (sum(Ene)/float(N))**2)
	  x = (sum(Mag**2)/float(N) - (sum(Mag)/float(N))**2)
		
	  aver_c = c/float(nvol)			!calore specifico medio
	  aver_x = x/float(nvol)			!suscettività media
		
	  cb = (sum(Mag**4)/(sum(Mag**2)**2))*N

	  !calcolo errore con bootstrap, se l'ultimo parametro è 1
	  !viene calcolato anche l'errore sul cumulate di binder
		
   	  call errore(Ene, daver_e, daver_c, dcb, 0)
   	  call errore(Mag, daver_m, daver_x, dcb, 1)

    	  write(2,*) aver_e, aver_m, aver_c, aver_x, cb,	!salvo su file
     &	             daver_e, daver_m, daver_c, daver_x, dcb
		
      enddo
	
      call ranfinish
	
      call cpu_time(finish)
      print '("tempo di esecuzione= ", f8.4," secondi.")', finish-start

      end program analisi
	
C=============================================================================

      subroutine errore(x, dx, dx1, db, B)
      common N, R, nvol
      
      real*8, dimension(:), allocatable :: z, u, a, q
      real*8, dimension(N) :: x
      integer(16) :: nb, Dd, i, j, l
      real*8 :: media_x, media_dx, dx, dx1, media_b, db
      integer :: g, R, B
      
      dx = 0
      dx1 = 0
      db = 0
      media_x = 0
      media_dx = 0
      media_db = 0

      allocate(z(N), u(R), a(R))
      if(B==1) then
	  allocate(q(R))
      endif
	
      !il calcolo dell'errore avviene tramite il binned bootstrap poichè
      !a priori non si può sapere qual è il tempo di decorellazione migliore
      !da inserire nella simulazione, e per non sprecare tempo macchina,
      !bisogna tenere conto di una possibile correlazione frai i dati
	
      Dd = 1e4			!dimensione dei blocchi	
      nb = N/Dd			!numero dei blocchi
      do l = 1, R		!ciclo sui ricampionamenti
          do i = 1, nb
		
	      j = int(ran2()*N +1) !scelgo sito a caso
	      do g= 1, Dd		 	
	          z((i-1)*Dd+g) = x(mod(j+g-2,N)+1) !ricampiono a blocchi	
	      enddo
			
	  enddo
		
	  !calcolo della media e della varianza dei ricampionamenti
		
	  a(l) = sum(z)/float(N*nvol)
	  u(l) = (sum(z**2)/float(N) - (sum(z)/float(N))**2)/float(nvol)	
	  if(B==1) then
	      q(l) = (sum(z**4)/(sum(z**2)**2))*N
	  endif	
      enddo
	
      media_dx = sum(u)/float(R)		!calcolo la media degli estimatori
      media_x  = sum(a)/float(R)
	
      if(B==1) then
          media_b = sum(q)/float(R)
      endif
	
      do i=1, R
          dx1 = dx1 + (u(i) - media_dx)**2 !calcolo scarto quadratico
	  dx  = dx  + (a(i) - media_x )**2
	  if(B==1) then
	      db = db + (q(i) - media_b)**2
	  endif
      enddo
	
      dx  = sqrt(dx/float(R - 1))		!prendo l'errore sul campione
      dx1 = sqrt(dx1/float(R - 1))		!perchè  sono stai effettuati
      if(B==1) then				!R ricampionamenti
	  db = sqrt(db/float(R - 1))	        !quindi divido solo per R-1
      endif					!e non R(R-1)
	
      return
      end
	
C=============================================================================
      function ran2()
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      real*4 ran2,am,eps,rnmx
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
