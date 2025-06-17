C	*********************************************************************
C	PROGRAM FOR CALCULATING LYAPUNOV EXPONENT FOR LIENARD OSCILLATOR
C	*********************************************************************

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION X(100),CUM(100),ZNORM(100),GSC(100)

	 COMMON/C1/ar,ak,abeta,abet0,amu,ags,agi,d,T,eps,beta
	 COMMON/C3/N
	
      open(14,file='test_avg_ly.dat')


       
		
      ar = 2.0d0
      ak = 1000.0d0
      abeta = 6.90d0
      amu = 0.001d0
      ags = 9.10d0
      agi = 1.0d0
      eps = 0.2d0
      d = 0.1d0
      T = 0.0D0
      N=3
      NN=(N*N)+N

      NPAS = 2000
      NRUN = 500
      NT = 500
      NLYA = 1000*NRUN
      IT = NT*NPAS

	H = 0.01D0
 !      call random_number(x(1))
!      X(1) = sqrt(-2*log(u1))*cos(2*pi*u2)
      	X(1) = 0.051101d0
	X(2) = 0.0204d0
	X(3) = 0.3040d0

!       write (*,12)xi,x(1)
!12     format(2F9.3)
      DO I = 1,IT
        CALL RK1(X,H)
      END DO
C
C     INITIAL CONDITION FOR LINEAR SYSTEM

    	 X(4) = 0.051100d0
	X(5) = 0.0204d0
	X(6) = 0.3040d0
	X(7) = 0.051100d0
	X(8) = 0.0204d0
	X(9) = 0.3040d0
	X(10) = 0.051100d0
  
      DO I=1,N
      X((N+1)*I)=1.0D0
      CUM(I)=0.0D0
      END DO
c
      do abet0 =0.0d0,3.0d0,0.005d0
      DO I = 1,NLYA
        t=t+h
      CALL RK4(X,H)
c
C     NORMALIZE FIRST VECTOR

      ZNORM(1)=0.0D0
      DO J=1,N
      ZNORM(1)=ZNORM(1)+X(N*J+1)**2
      END DO
      ZNORM(1)=SQRT(ZNORM(1))
      DO J=1,N
      X(N*J+1)=X(N*J+1)/ZNORM(1)
      END DO
C
C     GENERATE THE NEW ORTHONORMAL SET OF VECTORS
      DO 40 J=2,N
C
C     GENERATE J-1 GSR COEFFICIENTS
      DO 10 K=1,J-1
      GSC(K)=0.0D0
      DO 10 L=1,N
      GSC(K)=GSC(K)+X(N*L+J)*X(N*L+K)
10    CONTINUE
C
C     CONSTRUCT A NEW VECTOR
      DO 20 K=1,N
      DO 20 L=1,J-1
      X(N*K+J)=X(N*K+J)-GSC(L)*X(N*K+L)
20    CONTINUE
C
C     CALCULATE THE VECTOR'S NORM
      ZNORM(J)=0.0D0
      DO 30 K=1,N
      ZNORM(J)=ZNORM(J)+X(N*K+J)**2
30    CONTINUE
      ZNORM(J)=SQRT(ZNORM(J))
C
C     NORMALIZE THE NEW VECTOR
      DO 40 K=1,N
      X(N*K+J)=X(N*K+J)/ZNORM(J)
40    CONTINUE
                                                                                
C     UPDATE RUNNING VECTOR MAGNITUDES
      DO K=1,N
      CUM(K)=CUM(K)+DLOG(ZNORM(K))/DLOG(2.0D0)
      END DO
      END DO

      DO K=1,N
      CUM(K)=CUM(K)/(FLOAT(NLYA)*H)
      END DO
C
       
       
15	format(1x,4f10.6)
      WRITE(14,*)abet0,CUM(1),CUM(2),CUM(3)
C      WRITE(*,*)aep_v,CUM(1),CUM(2),CUM(3)
      print 15,abet0,CUM(1),CUM(2),CUM(3)
      END DO
	STOP
      END
C	********************************************************************
      SUBROUTINE RK4(X,H)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(100),TEMP(100),AK1(100),AK2(100),AK3(100)
      DIMENSION AK4(100),PRIME(100)

      COMMON/C3/N

      NN=(N*N)+N

      DO I=1,NN
      TEMP(I)=X(I)
      END DO

      CALL DERIVE(TEMP,PRIME)
      DO I=1,NN
      AK1(I)=H*PRIME(I)
      TEMP(I)=X(I)+AK1(I)/2.0D0
      END DO

      CALL DERIVE(TEMP,PRIME)
      DO I=1,NN
      AK2(I)=H*PRIME(I)
      TEMP(I)=X(I)+AK2(I)/2.0D0
      END DO

      CALL DERIVE(TEMP,PRIME)
      DO I=1,NN
      AK3(I)=H*PRIME(I)
      TEMP(I)=X(I)+AK3(I)
      END DO

      CALL DERIVE(TEMP,PRIME)
      DO I=1,NN
      AK4(I)=H*PRIME(I)
      X(I)=X(I)+1/6.0D0*(AK1(I)+2.0D0*(AK2(I)+AK3(I))+AK4(I))
      END DO

      RETURN
      END
C
C 
      SUBROUTINE RK1(X,H)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(100),TEMP(100),AK1(100),AK2(100),AK3(100)
      DIMENSION AK4(100),PRIME(100)
      COMMON/C3/N

      DO I=1,N
      TEMP(I)=X(I)
      END DO

      CALL DER(TEMP,PRIME)
      DO I=1,N
      AK1(I)=H*PRIME(I)
      TEMP(I)=X(I)+AK1(I)/2.0D0
      END DO

      CALL DER(TEMP,PRIME)
      DO I=1,N
      AK1(I)=H*PRIME(I)
      TEMP(I)=X(I)+AK1(I)/2.0D0
      END DO

      CALL DER(TEMP,PRIME)
      DO I=1,N
      AK2(I)=H*PRIME(I)
      TEMP(I)=X(I)+AK2(I)/2.0D0
      END DO

      CALL DER(TEMP,PRIME)
      DO I=1,N
      AK3(I)=H*PRIME(I)
      TEMP(I)=X(I)+AK3(I)
      END DO

      CALL DER(TEMP,PRIME)
      DO I=1,N
      AK4(I)=H*PRIME(I)
                                                                                
      X(I)=X(I)+1/6.0D0*(AK1(I)+2.0D0*(AK2(I)+AK3(I))+AK4(I))
      END DO

      RETURN
      END

      SUBROUTINE DER(X,PRIME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(100),PRIME(100)
        
	
	COMMON/C1/ar,ak,abeta,abet0,amu,ags,agi,d,T,eps,beta
	COMMON/C3/N
	
      pi=4*atan(1.0d0)
       beta = abeta*(1.0d0 + abet0*(sin(0.01d0*T)))
	
	!write(*,*)amu,xigma,abeta,agma,at,aep
	!write(*,*)abeta0,T
		
	PRIME(1) = ar*(x(1)+x(2))*(1.0d0-(x(1)+x(2))/ak)
     & -ags*x(1)*x(3)-beta*x(1)*x(2)
	PRIME(2) = beta*x(1)*x(2) -agi*x(2)*x(3)-amu*x(2)
	PRIME(3) = eps*(ags*x(1)+agi*x(2))*x(3)-d*x(3)

      RETURN
      END
C
C
      SUBROUTINE DERIVE(X,PRIME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(100),PRIME(100)
        
	COMMON/C1/ar,ak,abeta,abet0,amu,ags,agi,d,T,eps,beta
	COMMON/C3/N
	
      pi=4*atan(1.0d0)
       beta = abeta*(1.0d0 + abet0*(sin(0.01d0*T)))
	
	!write(*,*)ar,ak,beta,abet0,amu,ags,agi,d,T,eps
	!write(*,*)abeta0,T
		
	PRIME(1) = ar*(x(1)+x(2))*(1.0d0-(x(1)+x(2))/ak)
     & -ags*x(1)*x(3)-beta*x(1)*x(2)
	PRIME(2) = beta*x(1)*x(2) -agi*x(2)*x(3)-amu*x(2)
	PRIME(3) = eps*(ags*x(1)+agi*x(2))*x(3)-d*x(3)
      

      DO I= 0,N-1
      pi=4*atan(1.0d0)
       beta = abeta*(1.0d0+abet0*(sin(0.01d0*T)))
	!write(*,*) abeta0
	
       PRIME(4+I) = (ar*(1.0d0-((x(1)+x(2))/ak))-ar*((x(1)+x(2))/ak)
     & -x(3)*ags-x(2)*beta)*x(4+I)+(ar*(1.0d0-((x(1)+x(2))/ak))
     & -(ar*((x(1)+x(2))/ak))-beta*x(1))*x(7+I)-x(1)*ags*x(10+I) 
      PRIME(7+I) = x(2)*beta*x(4+I) + (-x(3)*agi-amu+x(1)*beta)*x(7+I)
     & -x(2)*agi*x(10+I)
      PRIME(10+I)= x(3)*eps*ags*x(4+I)+eps*x(3)*agi*x(7+I)
     & +((x(1)*ags+x(2)*agi)*eps-d)*x(10+I)
    
	END DO
      RETURN
      END
                                                                      
