C#######################################################################
C           Implementation of Algorith described in the paper
C  "Simulating Coulomb and Log-Gases with Hybrid Monte Carlo Algorithms"
C                 by Djalil Chafaï and Gregoire Ferré
C            DOI: https://doi.org/10.1007/s10955-018-2195-6
C
C     Author: Pimenta, J. V. A.
C#######################################################################
C#######################################################################
      PROGRAM HMC

      IMPLICIT REAL*8 (A-H,O-Z)

C#######################################################################
C#######################################################################
C     Considerations:
C
C     We onsider gases that in a subspace S of R^n, with dimension m
C
C     The external potential is given by V : S -> R
C
C     The interaction potential is given by W : S -> (-inf, inf]
C
C     We take N particles in S, and for any N >= 2, we define P_n to be
C     P_N = 1/Z_N exp(-beta H_N(x_1, ..., x_N)) dx_1 ... dx_N
C     .: Z_N is the partition function, beta is the inverse temperature
C     and H_N is the Hamiltonian
C
C     The Hamiltonian is given by
C     H_N(x_1,...,x_N) = 1/N sum_{i=1}^N V(x_i) + 
C                        1/2N^2 sum_{i \neq j} W(x_i - x_j)
C
C     Note that the Hamiltonian is invariant by permutation and func of
C     \mu_N = 1/N sum_{i=1}^N \delta_{x_i} 
C
C     It is know that the empirical measure converges to an non random
C     measure, i.e. \mu_N -> arg inf \epsilon(\mu) where
C     \epsilon(\mu) = \int V(x)d\mu(x) 
C                   + 1/2 \int\int W(x-y)d\mu(x)d\mu(y)
C
C ----------------------------------------------------------------------
C
C     The Hybrid Monte Carlo algorithm is described as follows:
C
C     Define a space E = R^mn and let U_N : E -> R be smooth such that
C     exp{-\beta_N U_N} is a density w.r.t. the Lebesgue measure on E
C
C     Let (X_t, V_t) be a Markov chain on E x E with invariant measure
C     \pi_N(dx dv) = exp{-\beta_N U_N(x)} dx dp solution of
C
C     dX_t = \alpha_N \Nabla U_N(Y_t) dt
C     dY_t = - \alpha_N \Nabla H_N(X_t) dt 
C            - \gamma_N \alpha_N \Nabla U_N(X_t) dt
C            + \sqrt{2 ((\gamma_N \alpha_N) / \beta_N) dB_t}
C
C     where (B_t)_{t>0} is a Brownian motion on E
C     \gamma_N is a coeficient that represents the friction
C     \alpha_N is a coeficient that represents the step size
C
C#######################################################################
C#######################################################################
C     Parameters:
C
C     N: number of particles, scalar
C     m: dimension of the space, scalar
C     beta: inverse temperature, scalar
C     alpha: step size, scalar
C     gamma: friction, scalar
C     nsteps: number of steps, scalar
C     eta: coeficient exp{-\gamma_N \alpha_N \delta t}, scalar
C
C     xk: position at k, vector of size (N,m)
C     xtilde: candidate positions for k+1, vector of size (N,m)
C     xk1: position at k + 1, vector of size (N,m)
C     vk: velocities at k, vector of size (N,m)
C     vtilde: gaussian modified velocities, vector of size (N,m)
C     vtildehalf: candidate velocities for k+1/2, vector of size (N,m)
C     vtildek1: candidate velocities for k+1, vector of size (N,m)
C     vk1: velocitie at k + 1, vector of size (N,m)
C     F: force, vector of size (N,m)
C     GVe: gradient of the external potential, vector of size (m)
C     GW: gradient of the interaction potential, vector of size (m)

      PARAMETER (N = 8, m = 2, beta = 1, alpha = 0.1, gamma = 0.1)

      DIMENSION xk(N,m), xtilde(N,m),
     &      vk(N,m),vtildehalf(N,m),vtilde(N,m),vtildek1(N,m),
     &      F(N,m), GVe(m), GW(m)

      COMMON /X/ xk, xtilde
      COMMON /V/ vk, vtildehalf, vtilde, vtildek1
      COMMON /G/ F, GVe, GW

      nsteps = 1000   
      tstep = 0.1

      eta = EXP(-gamma * alpha * tstep)

C#######################################################################
C#######################################################################
C     Step 1: Initialization of (x0, v0)
C     We take x0 = 0 and v0 = 0
C     We also initialize xk = x0 and vk = v0

      CALL INIT()

C ---------------------------------------------------------------------
      DO 10 k = 1, nsteps

C ---------------------------------------------------------------------
C     Step 2: Update the velocities with
C     vtilde_k = \eta vk_k + \sqrt{(1 - \eta^2) / \beta_N} \G_k
C     where \G_k is a standard Gaussian vector
C     and \eta = \exp{-\gamma_N \alpha_N \delta t}

      CALL GaussianUPT_P(eta, beta)

C ---------------------------------------------------------------------
C     Step 3: Calculate 
C     vtildehalf_k = vtilde_k - \alpha_N \Nabla H_N(xk) timestep/2
C     xtilde_k = xk + \alpha_N vtildehalf_k timestep
C     vtilde = vtildehalf_k - \alpha_N \Nabla H_N(xtilde_k) timestep/2
C
C     Note that we update a half step of the velocities then we update
C     the positions and then we update the velocities again

      CALL UPDATE_halfV(tstep, alpha)
      CALL UPDATE_X(tstep, alpha)
      CALL UPDATE_V(tstep)

C ---------------------------------------------------------------------
C     Step 4: define the acceptance probability
C     prob=1 ^ exp{-\beta_N (H_N(xtilde_k)-H_N(xk)+vtilde_k^2/2-vk^2/2)}

      p = CALC_PROB(beta)
      p = MIN(1.0, p)

C ---------------------------------------------------------------------
C     Step 5: accept or reject the candidate with probability p

      IF (p > RAND(0)) THEN
            xk = xtilde
            vk = vtildek1
      ELSE
            xk = xk
            vk = -vtildek1
      END IF

C ---------------------------------------------------------------------

 10   END DO
C#######################################################################

      END PROGRAM HMC

C#######################################################################
C     Subroutines:

C     INIT: initialization of (x0, v0)
c           modifies xk, vk and F
      SUBROUTINE INIT()
            PARAMETER(N = 8, m = 2)
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION xk(N,m), xtilde(N,m),
     &      vk(N,m),vtildehalf(N,m),vtilde(N,m),vtildek1(N,m),
     &      F(N,m), GVe(m), GW(m)
            COMMON /X/ xk, xtilde
            COMMON /V/ vk, vtildehalf, vtilde, vtildek1
            COMMON /G/ F, GVe, GW

            DO i = 1, N
            DO j = 1, m

                  xk(i,j) = -1.0 + 2*RAND(0)
                  vk(i,j) = 0.0
                  F(i,j) = 0.0
                  
            END DO
            END DO

      END SUBROUTINE INIT

C     GaussianUPT_P: update the velocities with the gaussian variable
c           modifies vtilde
      SUBROUTINE GaussianUPT_P(eta, beta)
            PARAMETER(N = 8, m = 2)
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION 
     &      vk(N,m),vtildehalf(N,m),vtilde(N,m),vtildek1(N,m)
            COMMON /V/ vk, vtildehalf, vtilde, vtildek1

            DO i = 1, N
            DO j = 1, m
                  v_tilde = eta*vk(i,m)+SQRT((1 - eta**2)/beta)*Gauss()
                  vtilde(i,m) = v_tilde
            END DO
            END DO

      END SUBROUTINE GaussianUPT_P

C     UPDATE_halfV: update the velocities a half step
c           modifies vtildehalf
      SUBROUTINE UPDATE_halfV(tstep, alpha)
            PARAMETER(N = 8, m = 2)
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION xk(N,m), xtilde(N,m),
     &      vk(N,m),vtildehalf(N,m),vtilde(N,m),vtildek1(N,m),
     &      F(N,m), GVe(m), GW(m)
            COMMON /X/ xk, xtilde
            COMMON /V/ vk, vtildehalf, vtilde, vtildek1
            COMMON /G/ F, GVe, GW

            Call GRAD_H(.FALSE.)

            DO i = 1, N
            DO j = 1, m

                  vtildehalf(i,j) = vtilde(i,j) - alpha*F(i,j)*tstep/2

            END DO
            END DO

      END SUBROUTINE UPDATE_halfV

C     UPDATE_X: update the positions
c           modifies xtilde
      SUBROUTINE UPDATE_X(tstep, alpha)
            PARAMETER(N = 8, m = 2)
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION xk(N,m), xtilde(N,m),
     &      vk(N,m),vtildehalf(N,m),vtilde(N,m),vtildek1(N,m)
            COMMON /X/ xk, xtilde
            COMMON /V/ vk, vtildehalf, vtilde, vtildek1

            DO i = 1, N
            DO j = 1, m

                  xtilde(i,j) = xk(i,j) +
     &                          alpha * vtildehalf(i,j) * tstep

            END DO
            END DO

      END SUBROUTINE UPDATE_X

C     UPDATE_V: complete the update of the velocities
c           modifies vtildek1
      SUBROUTINE UPDATE_V(tstep)
            PARAMETER(N = 8, m = 2)
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION xk(N,m), xtilde(N,m),
     &      vk(N,m),vtildehalf(N,m),vtilde(N,m),vtildek1(N,m),
     &      F(N,m), GVe(m), GW(m)
            COMMON /X/ xk, xtilde
            COMMON /V/ vk, vtildehalf, vtilde, vtildek1
            COMMON /G/ F, GVe, GW

            CALL GRAD_H(.TRUE.)

            DO i = 1, N
            DO j = 1, m

                  vtildek1(i,j) = vtildehalf(i,j) - alpha*F(i,j)*tstep/2

            END DO
            END DO

      END SUBROUTINE UPDATE_V

C     GRAD_H: gradient of the Hamiltonian (force)
c           modifies F, GVe and GW
      SUBROUTINE GRAD_H(next)
            PARAMETER(N = 8, m = 2)
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION x(N,m), xk(N,m), xtilde(N,m)
            DIMENSION F(N,m), F_aux(N,N,m), GVe(m), GW(m)
            LOGICAL next
            COMMON /G/ F, GVe, GW
            COMMON /X/ xk, xtilde

            if (next) then
                  x = xtilde
            else
                  x = xk
            end if

            DO i = 1, N
                  DO j = 1, i-1
                        CALL GRAD_W(x(i,1:m), x(j,1:m))
                        F_aux(i,j,:) = - GW
                  END DO
            END DO

            DO i = 1, N

                  F(i,:) = 0.0
                  DO j = i+1, N
                        F(i,:) = F(i,:) - F_aux(i,j,:)
                  END DO

                  DO j = 1, i-1
                        F(i,:) = F(i,:) + F_aux(i,j,:)
                  END DO

                  F(i,:) = F(i,:) / N

                  CALL GRAD_Ve(x(i,1:m))
                  F(i,:) = F(i,:) - GVe

                  F(i,:) = F(i,:) / N

            END DO

      END SUBROUTINE GRAD_H

C     GRAD_Ve: gradient of the external potential
c           modifies GVe
      SUBROUTINE GRAD_Ve(x)
            PARAMETER(N = 8, m = 2)
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION x(m), GVe(m), GW(m), F(N,m)
            COMMON /G/ F, GVe, GW

c                 Gradient of V(x) = ||x||^2 / (2 * beta)
            GVe = x / beta

      END SUBROUTINE GRAD_Ve

C     GRAD_W: gradient of the interaction potential
c           modifies GW
      SUBROUTINE GRAD_W(x, y)
            PARAMETER(N = 8, m = 2)
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION x(m), y(m), v(m), GW(m), F(N,m), GVe(m)
            COMMON /G/ F, GVe, GW

c                 Gradient of W(x) = -log(||x||)
            v = x-y
            GW = -v / NORM2(v)**2

      END SUBROUTINE GRAD_W
C#######################################################################
C#######################################################################
C     Functions:
C     1-FUNCTION ARE CALLED IN MAIN
C     2-FUNCTION ARE CALLED IN SUBROUTINES
C     3-FUNCTION ARE CALLED IN FUNCTIONS

C     (2-FUNCTION) G: gaussian variable
c           return a standard gaussian variable, scalar
      FUNCTION Gauss()
            IMPLICIT REAL*8 (A-H,O-Z)
            mu = 0
            sigma = 1
            Gauss = mu + sigma * SQRT(-2*LOG(RAND(0)))*COS(2*PI*RAND(0))
            RETURN
      END FUNCTION Gauss

C     (1-FUNCTION) CALC_PROB: calculate the acceptance probability
c           return the acceptance probability, scalar
      FUNCTION CALC_PROB(beta)
            PARAMETER(N = 8, m = 2)
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION vk(N,m),vtildehalf(N,m),vtilde(N,m),vtildek1(N,m)
            COMMON /V/ vk, vtildehalf, vtilde, vtildek1

            CALC_PROB = 1
            DO i = 1, N
            DO j = 1, m

            CALC_PROB = CALC_PROB * EXP(-beta * (H(.TRUE.)-H(.FALSE.) +
     &                  vtilde(i,j)**2 / 2 - vk(i,j)**2 / 2))
      
            END DO
            END DO
       
            RETURN
            END FUNCTION CALC_PROB
      
C     (3-FUNCTION) H: Hamiltonian
c           return the Hamiltonian, scalar
      FUNCTION H(next)
            PARAMETER(N = 8, m = 2)
            IMPLICIT REAL*8 (A-H,O-Z)
            LOGICAL next
            DIMENSION x(N,m), xk(N,m), xtilde(N,m)
            COMMON /X/ xk, xtilde

            if (next) then
                  x = xtilde
            else
                  x = xk
            end if

            H = 0.0

            DO i = 1, N
                  H = H + Ve(x(i,:))
                  DO j = i+1, N
                        H = H + W(x(i,:), x(j,:)) / N
                  END DO
            END DO

            H = H / N

            RETURN
      END FUNCTION H

C     (3-FUNCTION) V: external potential
c           return the external potential, scalar
      FUNCTION Ve(x)
            PARAMETER(N = 8, m = 2)
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION x(m)

c                 V(x) = ||x||^2 / (2 * beta)
            Ve = DOT_PRODUCT(x,x) / (2 * beta)

            RETURN
      END FUNCTION Ve

C     (3-FUNCTION) W: interaction potential
c           return the interaction potential, scalar
      FUNCTION W(x, y)
            PARAMETER(N = 8, m = 2)
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION x(m), y(m)

c           W(x) = -log(||x||)
            W = -LOG(NORM2(x-y))

            RETURN
      END FUNCTION W

C#######################################################################
