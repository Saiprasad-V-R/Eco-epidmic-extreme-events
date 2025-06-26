      PROGRAM Eco-Epidemic
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER N, I, K, L, NPAS, NRUN, NT, NLYA, IT, NN
      DOUBLE PRECISION X(100), CUM(100), ZNORM(100), GSC(100)
      DOUBLE PRECISION H, T

      COMMON /C1/ ar, ak, abeta, abet0, amu, ags, agi, d, T, eps, beta
      COMMON /C3/ N

      ! Open output file for average Lyapunov exponents
      OPEN(14, FILE='test_avg_ly.dat')

      ! Model parameters
      ar    = 2.0D0
      ak    = 1000.0D0
      abeta = 6.90D0
      amu   = 0.001D0
      ags   = 9.10D0
      agi   = 1.0D0
      eps   = 0.2D0
      d     = 0.1D0
      T     = 0.0D0
      N     = 3
      NN    = N*N + N        ! Total dimension for state + variational

      ! Simulation settings
      NPAS = 2000
      NRUN = 500
      NT   = 500
      NLYA = 1000 * NRUN     ! Number of steps for Lyapunov calc
      IT   = NT * NPAS       ! Transient spin-up steps
      H    = 0.01D0          ! Time step

      ! Initial conditions for main system
      X(1) = 0.051101D0
      X(2) = 0.0204D0
      X(3) = 0.3040D0

      ! --- Remove transients via simple RK stepping ---
      DO I = 1, IT
          CALL RK1(X, H)
      END DO

      ! Initialize variational vectors and cumulants
      DO I = 1, N
          X(N+I) = X(I)       ! Copy state into first perturbation block
          CUM(I)   = 0.0D0    ! Reset cumulant sums
      END DO

      ! --- Loop over forcing amplitude abet0 ---
      DO abet0 = 0.0D0, 3.0D0, 0.005D0

          T = 0.0D0          ! Reset time for this run

          ! Reset perturbation vectors to identity basis
          DO I = 1, N
              DO K = 1, N
                  IF (K .EQ. I) THEN
                      X(N*K + I) = 1.0D0
                  ELSE
                      X(N*K + I) = 0.0D0
                  END IF
              END DO
              CUM(I) = 0.0D0
          END DO

          ! --- Main RK4 integration for Lyapunov exponents ---
          DO I = 1, NLYA
              T = T + H
              CALL RK4(X, H)

              ! Orthonormalize perturbation vectors (Gram–Schmidt)
              ! First vector
              ZNORM(1) = 0.0D0
              DO K = 1, N
                  ZNORM(1) = ZNORM(1) + X(N*K+1)**2
              END DO
              ZNORM(1) = SQRT(ZNORM(1))
              DO K = 1, N
                  X(N*K+1) = X(N*K+1) / ZNORM(1)
              END DO

              ! Remaining vectors
              DO J = 2, N
                  ! Compute projections onto earlier vectors
                  DO K = 1, J-1
                      GSC(K) = 0.0D0
                      DO L = 1, N
                          GSC(K) = GSC(K) + X(N*L+J)*X(N*L+K)
                      END DO
                  END DO

                  ! Subtract projections
                  DO K = 1, N
                      DO L = 1, J-1
                          X(N*K+J) = X(N*K+J) - GSC(L)*X(N*K+L)
                      END DO
                  END DO

                  ! Normalize J-th vector
                  ZNORM(J) = 0.0D0
                  DO K = 1, N
                      ZNORM(J) = ZNORM(J) + X(N*K+J)**2
                  END DO
                  ZNORM(J) = SQRT(ZNORM(J))
                  DO K = 1, N
                      X(N*K+J) = X(N*K+J) / ZNORM(J)
                  END DO
              END DO

              ! Accumulate log norms for exponents
              DO K = 1, N
                  CUM(K) = CUM(K) + LOG(ZNORM(K)) / LOG(2.0D0)
              END DO
          END DO

          ! Compute average Lyapunov exponents
          DO K = 1, N
              CUM(K) = CUM(K) / (NLYA * H)
          END DO

          ! Output for current abet0
          WRITE(14,*) abet0, CUM(1), CUM(2), CUM(3)
          PRINT *,  abet0, CUM(1), CUM(2), CUM(3)
      END DO

      STOP
      END PROGRAM LienardLyapunov


      SUBROUTINE RK4(X, H)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER I, NN
      DOUBLE PRECISION X(100), TEMP(100), AK1(100), AK2(100), AK3(100), AK4(100), PRIME(100)
      COMMON /C3/ N

      NN = N*N + N     ! Total system dimension

      ! Copy state to TEMP
      DO I = 1, NN
          TEMP(I) = X(I)
      END DO

      ! Fourth‐order Runge–Kutta steps
      CALL DERIVE(TEMP, PRIME)
      DO I = 1, NN
          AK1(I) = H*PRIME(I)
          TEMP(I) = X(I) + AK1(I)/2.0D0
      END DO

      CALL DERIVE(TEMP, PRIME)
      DO I = 1, NN
          AK2(I) = H*PRIME(I)
          TEMP(I) = X(I) + AK2(I)/2.0D0
      END DO

      CALL DERIVE(TEMP, PRIME)
      DO I = 1, NN
          AK3(I) = H*PRIME(I)
          TEMP(I) = X(I) + AK3(I)
      END DO

      CALL DERIVE(TEMP, PRIME)
      DO I = 1, NN
          AK4(I) = H*PRIME(I)
          X(I) = X(I) + (AK1(I) + 2.0D0*(AK2(I)+AK3(I)) + AK4(I))/6.0D0
      END DO
      END SUBROUTINE RK4


      SUBROUTINE RK1(X, H)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER I
      DOUBLE PRECISION X(100), TEMP(100), AK1(100), AK2(100), AK3(100), AK4(100), PRIME(100)
      COMMON /C3/ N

      ! Simple RK4 for base system (first N variables)
      DO I = 1, N
          TEMP(I) = X(I)
      END DO

      CALL DER(TEMP, PRIME)
      DO I = 1, N
          AK1(I) = H*PRIME(I)
          TEMP(I) = X(I) + AK1(I)/2.0D0
      END DO

      CALL DER(TEMP, PRIME)
      DO I = 1, N
          AK2(I) = H*PRIME(I)
          TEMP(I) = X(I) + AK2(I)/2.0D0
      END DO

      CALL DER(TEMP, PRIME)
      DO I = 1, N
          AK3(I) = H*PRIME(I)
          TEMP(I) = X(I) + AK3(I)
      END DO

      CALL DER(TEMP, PRIME)
      DO I = 1, N
          AK4(I) = H*PRIME(I)
          X(I) = X(I) + (AK1(I) + 2.0D0*(AK2(I)+AK3(I)) + AK4(I))/6.0D0
      END DO
      END SUBROUTINE RK1


      SUBROUTINE DER(X, PRIME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER I
      DOUBLE PRECISION X(100), PRIME(100), pi
      COMMON /C1/ ar, ak, abeta, abet0, amu, ags, agi, d, T, eps, beta
      COMMON /C3/ N

      ! Compute current transmission β
      pi   = 4.0D0 * ATAN(1.0D0)
      beta = abeta * (1.0D0 + abet0 * SIN(0.01D0 * T))

      ! Base dynamics (prey, infected prey, predator)
      PRIME(1) = ar*(X(1)+X(2))*(1.0D0-(X(1)+X(2))/ak) - ags*X(1)*X(3) - beta*X(1)*X(2)
      PRIME(2) = beta*X(1)*X(2) - agi*X(2)*X(3) - amu*X(2)
      PRIME(3) = eps*(ags*X(1)+agi*X(2))*X(3) - d*X(3)
      END SUBROUTINE DER


      SUBROUTINE DERIVE(X, PRIME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER I, NN
      DOUBLE PRECISION X(100), PRIME(100), pi
      COMMON /C1/ ar, ak, abeta, abet0, amu, ags, agi, d, T, eps, beta
      COMMON /C3/ N

      NN = N*N + N     ! Dimension including variational part

      ! Compute current transmission β
      pi   = 4.0D0 * ATAN(1.0D0)
      beta = abeta * (1.0D0 + abet0 * SIN(0.01D0 * T))

      ! Base system equations
      PRIME(1) = ar*(X(1)+X(2))*(1.0D0-(X(1)+X(2))/ak) - ags*X(1)*X(3) - beta*X(1)*X(2)
      PRIME(2) = beta*X(1)*X(2) - agi*X(2)*X(3) - amu*X(2)
      PRIME(3) = eps*(ags*X(1)+agi*X(2))*X(3) - d*X(3)

      ! Variational (Jacobian) equations
      DO I = 0, N-1
          PRIME(4+I) = (ar*(1.0D0-(X(1)+X(2))/ak) - ar*((X(1)+X(2))/ak)
     &                   - ags*X(3) - beta*X(2))*X(4+I)
     &               + (ar*(1.0D0-(X(1)+X(2))/ak) - ar*((X(1)+X(2))/ak)
     &                   - beta*X(1))*X(7+I)
     &               - ags*X(1)*X(10+I)

          PRIME(7+I) = beta*X(1)*X(4+I)
     &               + (-agi*X(3)-amu+beta*X(1))*X(7+I)
     &               - agi*X(2)*X(10+I)

          PRIME(10+I)= eps*ags*X(3)*X(4+I)
     &               + eps*agi*X(3)*X(7+I)
     &               + (eps*(ags*X(1)+agi*X(2)) - d)*X(10+I)
      END DO
      END SUBROUTINE DERIVE
