      program correlation
      
      real*8, dimension(:), allocatable :: Ene, Mag
      real*8, dimension(:), allocatable :: Corr_ene, Corr_mag
      
      common N, i_corr
      
      open(unit=0, file="daticorr/dati10.dat", status="old")
      open(unit=2, file='daticorrplot/dati10.dat',status='unknown')
      
      read(0, *) N			!leggo i primi valori che sono 
      read(0, *) nvol			!necessari per l'analisi
      read(0, *) npassi
      
      i_corr = N/2
	
      allocate(Ene(N), Mag(N))
      allocate(Corr_ene(i_corr), Corr_mag(i_corr))

      
      do i = 1, N
	    read(0, *) Mag(i), Ene(i)
	enddo
	
	call corr(Mag, Corr_mag, tau_mag)
	call corr(Ene, Corr_ene, tau_ene)
	
	do j = 1, i_corr
	    write(2,*) Corr_mag(j), Corr_ene(j)
	enddo
	
      print '("tempo di auto-corr_mag = ", f4.1," iterate.")', tau_mag
      print '("tempo di auto-corr_ene = ", f4.1," iterate.")', tau_ene
      
      end program correlation
      
C======================================================================
C SUBROUTINE CORRELAZIONE
C======================================================================

      subroutine corr(x, c, tau)
      
      common N, i_corr
      
      real*8, dimension(i_corr) :: c
      real*8, dimension(N) :: x
      
      media = sum(x)/float(N)
      
      do k = 1, i_corr
          
          yc = 0.0
          do i = 1, N-k
              yc = yc + (x(i)-media)*(x(i+k)-media)
          enddo
          
			
          yc = yc/float(N-k)
 	    c(k) = yc
      enddo
      c = c/(c(1))

      tau = sum(c)
      
      return
      end
