      PROGRAM PredatorPrey
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER N, I, NRUN
      DOUBLE PRECISION X(3), H, T
      DOUBLE PRECISION z1, z2, z3

      ! Model parameters (COMMON block C1)
      COMMON /C1/ ar, ak, abeta, abeta0, amu, eps, ags, agi, d, T
      ! Dimension parameter (COMMON block C3)
      COMMON /C3/ N

      OPEN(1, FILE = 'bif.dat')      ! Output file for bifurcation data

      N = 3                          ! Number of equations

      ! Loop over forcing amplitude abeta0
      DO abeta0 = 0.0D0, 3.0D0, 0.005D0
          
          ! Set fixed parameter values
          ar    = 2.0D0
          ak    = 1000.0D0
          abeta = 6.9D0
          amu   = 0.001D0
          eps   = 0.2D0
          ags   = 9.1D0
          agi   = 1.0D0
          d     = 0.1D0

          ! Initialize time‐dependent forcing variables
          z1 = 0.0D0
          z2 = 0.0D0
          z3 = 0.0D0

          ! Initial conditions for prey, infected prey, predator
          X(1) = 0.051100D0
          X(2) = 0.0204D0
          X(3) = 0.3040D0

          ! Integration settings
          T    = 0.0D0            ! Start time
          H    = 0.01D0           ! Time step
          NRUN = 1600000          ! Total number of steps

          ! Time‐stepping via 4th‐order Runge–Kutta
          DO I = 1, NRUN
              T = T + H
              CALL RAKU4(X, H)

              ! Shift recent z‐values for local maximum detection
              z3 = z2
              z2 = z1
              z1 = X(3)

              ! After transient period, record peaks of X(3)
              IF (I .GT. 1500000) THEN
                  IF (z2 .GT. z1 .AND. z2 .GT. z3) THEN
                      WRITE(1,*) abeta0, X(3)
                      WRITE(*,*) abeta0, X(3)
                  END IF
              END IF
          END DO

      END DO

      STOP
      END PROGRAM PredatorPrey


      SUBROUTINE RAKU4(X, H)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER I, N
      DOUBLE PRECISION X(3), TEMP(3), RK1(3), RK2(3), RK3(3), RK4(3), PRIME(3)

      COMMON /C3/ N

      ! First sub‐step
      CALL EQUATION(X, PRIME)
      DO I = 1, N
          RK1(I) = H * PRIME(I)
          TEMP(I) = X(I) + 0.5D0 * RK1(I)
      END DO

      ! Second sub‐step
      CALL EQUATION(TEMP, PRIME)
      DO I = 1, N
          RK2(I) = H * PRIME(I)
          TEMP(I) = X(I) + 0.5D0 * RK2(I)
      END DO

      ! Third sub‐step
      CALL EQUATION(TEMP, PRIME)
      DO I = 1, N
          RK3(I) = H * PRIME(I)
          TEMP(I) = X(I) + RK3(I)
      END DO

      ! Fourth sub‐step and combine
      CALL EQUATION(TEMP, PRIME)
      DO I = 1, N
          RK4(I) = H * PRIME(I)
          X(I) = X(I) + (RK1(I) + 2.0D0*(RK2(I) + RK3(I)) + RK4(I)) / 6.0D0
      END DO

      END SUBROUTINE RAKU4


      SUBROUTINE EQUATION(X, PRIME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER I
      DOUBLE PRECISION X(3), PRIME(3)
      DOUBLE PRECISION pi, beta

      COMMON /C1/ ar, ak, abeta, abeta0, amu, eps, ags, agi, d, T

      pi = 4.0D0 * ATAN(1.0D0)

      ! Time‐varying transmission rate with sinusoidal forcing
      beta = abeta * (1.0D0 + abeta0 * SIN(0.01D0 * T))

      ! Prey (X1) growth, competition, predation, and infection loss
      PRIME(1) = ar*(X(1)+X(2))*(1.0D0 - (X(1)+X(2))/ak) 
     &           - ags*X(1)*X(3) 
     &           - beta*X(1)*X(2)

      ! Infected prey (X2) dynamics: infection gain, predation, and mortality
      PRIME(2) = beta*X(1)*X(2) 
     &           - agi*X(2)*X(3) 
     &           - amu*X(2)

      ! Predator (X3) growth from consumption and natural death
      PRIME(3) = eps*(ags*X(1) + agi*X(2))*X(3) 
     &           - d*X(3)

      END SUBROUTINE EQUATION
