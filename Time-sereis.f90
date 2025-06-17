PROGRAM PredatorPrey
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  	double precision, allocatable :: x(:),a1(:),a(:)
    COMMON/C1/ar,ak,abeta,abeta0,amu,eps,ags,agi,d,T
    COMMON/C3/N

    OPEN(1, FILE = 'ts_0.8.dat')
    OPEN(2, FILE = 'hs_0.8.dat')
    


    NRUN = 50000000
	transcient=30000000
	NT = 100
	!iteration=30000
	
	allocate(x(nt),a1(nrun-transcient),a(nrun-transcient))	
    ! Set parameters
    ar = 2.0D0
    ak = 1000.0D0
    abeta = 6.9D0
    abeta0 = 0.8d0
    amu = 0.001D0
    eps = 0.2D0
    ags = 9.1D0
    agi = 1.0D0
    d = 0.1D0
    z1=0
	z2=0
	z3=0
	!n1=100
	sum_x = 0.0d0 
	sum_y = 0.0d0


    ! Initial conditions
    X(1) = 0.05011D0
    X(2) = 0.0204D0
    X(3) = 0.3040D0

    ! Simulation parameters
    T = 0.0D0		
    H = 0.01D0
    N = 3

   	jj = 0
	kk = 0
	sum_x = 0.0d0
	sum_y = 0.0d0

    ! Runge-Kutta integration
    DO I = 1, NRUN
        T = T + H
        CALL RAKU4(X, H)
            pi=4*atan(1.0d0)
        beta = abeta*(1.0d0 + abeta0*(sin(0.01d0*T)))
        IF(I .gt. transcient)WRITE(1,*) T, X(1), X(2), X(3)
        z3=z2
    	z2=z1 
	    z1=x(3)
        IF(I .gt. transcient) then
	        if (z2 .gt. z1 .and.  z2 .gt. z3)then 
        	jj=jj+1
        	a1(jj)=x(3)
        	sum_y=sum_y+x(3)
	        end if 
	
	        kk=kk+1
	        a(kk)=x(3)
	        sum_x=sum_x+x(3)
	    end if
	!end if 
	end do
	
	!sum_z=0.0
	!pss=sum_x/kk
!	std_dev = (x(2)-pss)**2  
!	sum_z=sum_z+std_dev
	
	!write(*,*) ps,pss
	!do i=1,kk
	
	                             
	sum_x = sum_x/float(kk)
	
!----------------------------------------------------
	sum_z = 0.0d0
	

	do i = 1,kk
	    std_dev = (a(i)-sum_x)**2
	    sum_z = sum_z + std_dev
	end do
	sum_z = sqrt(sum_z/float(kk))
	sum_y=sum_y/float(jj)     
	hs=sum_y+8.0*sum_z
	!write(*,*)d,hs
	write(2,*)hs


    STOP
END PROGRAM PredatorPrey

SUBROUTINE RAKU4(X, H)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    DIMENSION X(3), TEMP(3), RK1(3), RK2(3), RK3(3), RK4(3), PRIME(3)
    COMMON/C3/N

    ! First step
    CALL EQUATION(X, PRIME)
    DO I = 1, N
        RK1(I) = H * PRIME(I)
        TEMP(I) = X(I) + 0.5D0 * RK1(I)
    END DO

    ! Second step
    CALL EQUATION(TEMP, PRIME)
    DO I = 1, N
        RK2(I) = H * PRIME(I)
        TEMP(I) = X(I) + 0.5D0 * RK2(I)
    END DO

    ! Third step
    CALL EQUATION(TEMP, PRIME)
    DO I = 1, N
        RK3(I) = H * PRIME(I)
        TEMP(I) = X(I) + RK3(I)
    END DO

    ! Fourth step
    CALL EQUATION(TEMP, PRIME)
    DO I = 1, N
        RK4(I) = H * PRIME(I)
        X(I) = X(I) + (RK1(I) + 2.0D0 * (RK2(I) + RK3(I)) + RK4(I)) / 6.0D0
    END DO
END SUBROUTINE RAKU4

SUBROUTINE EQUATION(X, PRIME)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    DIMENSION X(3), PRIME(3)
    COMMON/C1/ar,ak,abeta,abeta0,amu,eps,ags,agi,d,T

    pi=4*atan(1.0d0)
    beta = abeta*(1.0d0 + abeta0*(sin(0.01d0*T)))

    PRIME(1) = ar * (X(1) + X(2)) * (1.0D0 - (X(1) + X(2)) / ak) - ags * X(1) * X(3) - beta * X(1) * X(2)
    PRIME(2) = beta * X(1) * X(2) - agi * X(2) * X(3) - amu * X(2)
    PRIME(3) = eps * (ags * X(1) + agi * X(2)) * X(3) - d * X(3)
END SUBROUTINE EQUATION

