      program analisi

      real*8, dimension(:), allocatable :: Ene, Mag
      real*8:: daver_e, daver_m, daver_c, daver_x, dcb
      real*8:: aver_e, aver_m, aver_c, aver_x, cb
      integer :: R
      integer(16) :: DB
      
      call cpu_time(start)
      call ranstart

      !apro file da cui leggere i dati da analizzare e file su
      !cui scrivere i risultati con relativi errori
      open(unit=0, file="data/data_L_10_J2_01.dat", status="old")
      open(unit=2, file='dataplot/data_L_10_J2_01.dat',status='unknown')

      R  = 100              ! Numero di ricampionamenti
      DB = 1000             ! Dimesione blocchi
      
      read(0, *) N          ! Leggo i primi valori che sono
      read(0, *) nvol       ! necessari per l'analisi
      read(0, *) npassi

      allocate(Ene(N), Mag(N))

      do j=1, npassi                   ! leggo a blocchi il file
          do i = 1, N                  ! ogni blocco una temperatura diversa
              read(0, *) Mag(i), Ene(i)
          enddo

          aver_e = sum(Ene)/float(N*nvol)    ! Energia media
          aver_m = sum(Mag)/float(N*nvol)    ! Magnetizzazione media

          c = (sum(Ene**2)/float(N) - (sum(Ene)/float(N))**2)
          x = (sum(Mag**2)/float(N) - (sum(Mag)/float(N))**2)

          aver_c = c/float(nvol)            ! Calore specifico medio
          aver_x = x/float(nvol)            ! Suscettività media

          cb = (sum(Mag**4)/(sum(Mag**2)**2))*N !cumulante di binder

          !calcolo errore con bootstrap, se l'ultimo parametro è 1
          !viene calcolato anche l'errore sul cumulate di binder

          call bootstrap(N, Ene, daver_e, 0, R, nvol, DB)
          call bootstrap(N, Mag, daver_m, 0, R, nvol, DB)
          call bootstrap(N, Ene, daver_c, 1, R, nvol, DB)
          call bootstrap(N, Mag, daver_x, 1, R, nvol, DB)

          write(2,*) aver_e, aver_m, aver_c, aver_x,  !salvo su file
     &             daver_e, daver_m, daver_c, daver_x

      enddo

      call ranfinish

      call cpu_time(finish)
      print '("tempo di esecuzione= ", f8.4," secondi.")', finish-start

      end program analisi

C=============================================================================
C Subroutine per il calcolo degli errori tramite binned bootstrap
C=============================================================================

      subroutine bootstrap(N, x, dx, ics, R, nvol, Db)
C=============================================================================
C     Subroutine for calculating errors using binned bootstrap
C
C     Parameters
C     N : int
C         size of data, length(x)
C     x : one dimensional array of  length(x) = N
C         data, markov chain
C     dx : float
C         variable to which the error will be written
C     ics : int
C         flag if 0 computere error on mean, if 1 compute error on variance
C     R : int
C         number of resampling
C     nvol : float
C         volume of lattice
C     Db : int
C         seize of the blocks
C=============================================================================
      real*8, dimension(:), allocatable :: z, a   ! axuliar array
      real*8, dimension(N) :: x                   ! initial array
      integer(16) :: nb, Db, i, j, l              ! parameter and indices
      real*8 :: media_x, dx                       ! final results
      integer :: g, R, ics                        ! others, pamaters

      allocate(z(N), a(R))

      ! The calculation of the error is done through the binned bootstrap since
      ! it is not possible to know a priori which is the best decorellation time
      ! to insert in the simulation, and in order not to waste machine time,
      ! a possible correlation between the data must be taken into account

      nb = N/Db         ! number of blocks

      do l = 1, R       ! loop of resampling

          do i = 1, nb
              j = int(ran2()*N +1) ! I choose site at random
              do g = 1, Db
                  z((i-1)*Db+g) = x(mod(j+g-2,N)+1) ! block's resampling
              enddo
          enddo

          ! calculation of the mean or the variance of the resamplings

          if (ics==0) then
              a(l) = sum(z)/float(N*nvol)
          endif

          if (ics==1) then
              a(l) = sum(z**2)/float(N) - (sum(z)/float(N))**2
              a(l) = a(l)/float(nvol)
          endif

      enddo

      media_x  = sum(a)/float(R)         ! mean
      dx = 0
      do i=1, R
          dx = dx + (a(i) - media_x )**2 ! standard deviation
      enddo

      dx  = sqrt(dx/float(R - 1))
      ! I take the error on the sample because we made R
      ! resamples so I only divide by R-1 and not R(R-1)

      return
      end
c============================================================================
c  RANDOM NUMBER GENERATOR: standard ran2 from numerical recipes
c============================================================================

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
