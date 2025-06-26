      PROGRAM PredatorPrey
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      ! Local parameters and loop indices
      INTEGER :: NRUN, transcient, NT, I, jj, kk
      INTEGER :: alloc_size1, alloc_size2
      
      ! Model variables and accumulators
      DOUBLE PRECISION :: T, H, pi, beta
      DOUBLE PRECISION :: z1, z2, z3
      DOUBLE PRECISION :: sum_x, sum_y, std_dev, sum_z, hs
      
      ! State and result arrays (allocatable)
      DOUBLE PRECISION, ALLOCATABLE :: x(:), a1(:), a(:)
      
      ! Common blocks for parameters
      COMMON /C1/ ar, ak, abeta, abeta0, amu, eps, ags, agi, d, T
      COMMON /C3/ N

      !--------------------------------------------------------------------
      ! Open output files: time series and “hs” statistic
      !--------------------------------------------------------------------
      OPEN(1, FILE = 'ts_0.8.dat')
      OPEN(2, FILE = 'hs_0.8.dat')

      !--------------------------------------------------------------------
      ! Simulation settings
      !--------------------------------------------------------------------
      NRUN      = 50000000             ! Total RK steps
      transcient= 30000000             ! Transient steps to discard
      NT        =       100            ! Number of points per block
      
      alloc_size1 = NRUN - transcient  ! Size for post-transient arrays
      alloc_size2 = NRUN - transcient

      ALLOCATE(x(    NT),              &
               a1( alloc_size1 ),     &
               a ( alloc_size2  ))

      !--------------------------------------------------------------------
      ! Model parameters (COMMON /C1/)
      !--------------------------------------------------------------------
      ar     =  2.0D0
      ak     =  1000.0D0
      abeta  =  6.9D0
      abeta0 =  0.8D0
      amu    =  0.001D0
      eps    =  0.2D0
      ags    =  9.1D0
      agi    =  1.0D0
      d      =  0.1D0

      !--------------------------------------------------------------------
      ! Initialize auxiliary variables
      !--------------------------------------------------------------------
      z1 = 0.0D0
      z2 = 0.0D0
      z3 = 0.0D0

      sum_x = 0.0D0
      sum_y = 0.0D0

      !--------------------------------------------------------------------
      ! Initial conditions for prey, infected prey, predator
      !--------------------------------------------------------------------
      x(1) = 0.05011D0
      x(2) = 0.0204D0
      x(3) = 0.3040D0

      !--------------------------------------------------------------------
      ! Time‐integration parameters
      !--------------------------------------------------------------------
      T = 0.0D0                ! Start time
      H = 0.01D0               ! Time step
      N = 3                    ! Dimension of system

      jj = 0                   ! Peak counter
      kk = 0                   ! Sample counter

      !--------------------------------------------------------------------
      ! Main RK4 integration loop
      !--------------------------------------------------------------------
      DO I = 1, NRUN
          T = T + H
          CALL RAKU4(x, H)         ! Advance state by one RK4 step

          pi   = 4.0D0 * ATAN(1.0D0)
          beta = abeta * (1.0D0 + abeta0 * SIN(0.01D0 * T))

          ! Record full state after transient
          IF (I .GT. transcient) THEN
              WRITE(1,*) T, x(1), x(2), x(3)
          END IF

          ! Track recent predator values for peak detection
          z3 = z2
          z2 = z1
          z1 = x(3)

          IF (I .GT. transcient) THEN
              ! Detect local maxima of x(3)
              IF (z2 .GT. z1 .AND. z2 .GT. z3) THEN
                  jj = jj + 1
                  a1(jj) = x(3)
                  sum_y = sum_y + x(3)
              END IF

              ! Accumulate all post-transient samples
              kk = kk + 1
              a(kk) = x(3)
              sum_x = sum_x + x(3)
          END IF
      END DO

      !--------------------------------------------------------------------
      ! Compute statistics on predator time series
      !--------------------------------------------------------------------
      sum_x = sum_x / DBLE(kk)          ! Mean of all samples

      sum_z = 0.0D0
      DO I = 1, kk
          std_dev = (a(I) - sum_x)**2
          sum_z = sum_z + std_dev
      END DO
      sum_z = SQRT(sum_z / DBLE(kk))    ! Standard deviation

      sum_y = sum_y / DBLE(jj)          ! Mean of peaks
      hs    = sum_y + 8.0D0 * sum_z     ! Combined statistic

      WRITE(2,*) hs                     ! Output hs to file

      STOP
      END PROGRAM PredatorPrey


      SUBROUTINE RAKU4(X, H)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER :: I
      DOUBLE PRECISION :: TEMP(3), RK1(3), RK2(3), RK3(3), RK4(3), PRIME(3)
      COMMON /C3/ N

      ! 4th‐order Runge–Kutta integration for 3-dimensional system
      CALL EQUATION(X, PRIME)
      DO I = 1, N
          RK1(I) = H * PRIME(I)
          TEMP(I) = X(I) + 0.5D0 * RK1(I)
      END DO

      CALL EQUATION(TEMP, PRIME)
      DO I = 1, N
          RK2(I) = H * PRIME(I)
          TEMP(I) = X(I) + 0.5D0 * RK2(I)
      END DO

      CALL EQUATION(TEMP, PRIME)
      DO I = 1, N
          RK3(I) = H * PRIME(I)
          TEMP(I) = X(I) +       RK3(I)
      END DO

      CALL EQUATION(TEMP, PRIME)
      DO I = 1, N
          RK4(I) = H * PRIME(I)
          X(I)  = X(I) + (RK1(I) + 2.0D0*(RK2(I) + RK3(I)) + RK4(I)) / 6.0D0
      END DO
      END SUBROUTINE RAKU4


      SUBROUTINE EQUATION(X, PRIME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER :: I
      DOUBLE PRECISION :: pi, beta
      DOUBLE PRECISION :: xloc(3), primeloc(3)
      COMMON /C1/ ar, ak, abeta, abeta0, amu, eps, ags, agi, d, T
      COMMON /C3/ N

      ! Copy into locals for clarity
      DO I = 1, N
          xloc(I) = X(I)
      END DO

      ! Compute time-varying transmission rate β
      pi   = 4.0D0 * ATAN(1.0D0)
      beta = abeta * (1.0D0 + abeta0 * SIN(0.01D0 * T))

      ! Equations of motion
      PRIME(1) = ar * (xloc(1) + xloc(2)) * (1.0D0 - (xloc(1) + xloc(2)) / ak)  &
                 - ags * xloc(1) * xloc(3)                                   &
                 - beta * xloc(1) * xloc(2)
      PRIME(2) = beta * xloc(1) * xloc(2)                                      &
                 - agi  * xloc(2) * xloc(3)                                   &
                 - amu  * xloc(2)
      PRIME(3) = eps * (ags * xloc(1) + agi * xloc(2)) * xloc(3)               &
                 - d    * xloc(3)
      END SUBROUTINE EQUATION
