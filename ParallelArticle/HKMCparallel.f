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

      USE OMP_LIB

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
C     sdn: coeficient \sqrt{(1 - \eta^2) / \beta_N}, scalar
C     pi: coeficient \pi, scalar
C     tstep: time step, scalar
C     niter: number of iterations to register, scalar
C     p: acceptance probability, scalar
C
C     xk: position at k, vector of size (N,m)
C     xtildek1: candidate positions for k+1, vector of size (N,m)
C     vk: velocities at k, vector of size (N,m)
C     vtilde: gaussian modified velocities, vector of size (N,m)
C     vtildek1: candidate velocities for k+1, vector of size (N,m)
C     F: force, vector of size (N,m)
C     GVe: gradient of the external potential, vector of size (m)
C     GW: gradient of the interaction potential, vector of size (m)

      PARAMETER (N = 50, m = 1, beta = 2.0, alpha = 0.1)

      DIMENSION xk(N,m), xtildek1(N,m),
     &      vk(N,m),vtilde(N,m),vtildek1(N,m),
     &      F(N,m), GVe(m), GW(m)

      COMMON /X/ xk, xtildek1
      COMMON /V/ vk, vtilde, vtildek1
      COMMON /G/ F, GVe, GW

      nsteps = 10000
      niter = 500
      tstep = 0.1
      gamma = 10.0
      t = 1.0
      a = 5.0

      eta = EXP(-gamma * alpha * tstep)
      sdn = SQRT((1 - eta**2) / beta) / N
      pi = ACOS(-1.d0)

      call srand(987991650)

C#######################################################################
C#######################################################################
C     Step 1: Initialization of (x0, v0)
C     We take x0 = 0 and v0 = 0
C     We also initialize xk = x0 and vk = v0

      CALL INIT()

C ---------------------------------------------------------------------

      DO 10 k = 1, INT(nsteps/2)

C ---------------------------------------------------------------------
C     Step 2: Update the velocities with
C     vtilde_k = \eta vk_k + \sqrt{(1 - \eta^2) / \beta_N} \G_k
C     where \G_k is a standard Gaussian vector
C     and \eta = \exp{-\gamma_N \alpha_N \delta t}
      
      CALL GaussianV(eta, sdn, pi)

C ---------------------------------------------------------------------
C     Step 3: Calculate 
C     vtildehalf_k = vtilde_k - \alpha_N \Nabla H_N(xk) timestep/2
C     xtildek1_k = xk + \alpha_N vtildehalf_k timestep
C     vtilde = vtildehalf_k - \alpha_N \Nabla H_N(xtildek1_k) timestep/2
C
C     Note that we update a half step of the velocities then we update
C     the positions and then we update the velocities again

      CALL UPDATE(tstep, alpha, beta, t, a)

C ---------------------------------------------------------------------
C     Step 4: define the acceptance probability
C     prob=1 ^ exp{-\beta_N (H_N(xtildek1_k)-H_N(xk)+vtilde_k^2/2-vk^2/2)}

      p = PROB(beta, t, a)

C ---------------------------------------------------------------------
C     Step 5: accept or reject the candidate with probability p

      IF (RAND(0) <= p) THEN
            xk = xtildek1
            vk = vtildek1
      ELSE
            vk = -vtildek1
            CALL GRAD_H(.FALSE., beta, t, a)
      END IF

C ---------------------------------------------------------------------

 10   END DO

C#######################################################################
C     We do a parallel run of the algorithm to save some data

      OPEN(1,FILE='dataX.txt',STATUS='UNKNOWN')
      OPEN(2,FILE='dataV.txt',STATUS='UNKNOWN')

      PRINT *, OMP_IN_PARALLEL(), OMP_GET_MAX_THREADS()

      !$OMP PARALLEL
      print *, 'hello from thread:', OMP_GET_THREAD_NUM()
      !$OMP END PARALLEL
 
      DO 20 k = INT(nsteps/2)+1, nsteps
      
      CALL GaussianV(eta, sdn, pi)

      CALL UPDATE(tstep, alpha, beta, t, a)

      p = PROB(beta, t, a)

      IF (RAND(0) <= p) THEN
            xk = xtildek1
            vk = vtildek1
      ELSE
            vk = -vtildek1
            CALL GRAD_H(.FALSE., beta, t, a)
      END IF

C ---------------------------------------------------------------------
c     Save some data every 1000 steps
      IF (MOD(k,niter) == 0) THEN
            WRITE(1,*) xk / SQRT(2*beta)
            WRITE(2,*) vk
            WRITE(*,*) k, XK(1,1), VK(1,1)
      END IF

C ---------------------------------------------------------------------

 20   END DO

      CLOSE(1)
      CLOSE(2)

C#######################################################################

      END PROGRAM HMC

C#######################################################################
C     Subroutines:

C     INIT: initialization of (x0, v0)
c           modifies xk, vk, F
      SUBROUTINE INIT()
            PARAMETER(N = 50, m = 1)
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION xk(N,m), xtildek1(N,m),
     &                vk(N,m),vtilde(N,m),vtildek1(N,m),
     &                F(N,m), GVe(m), GW(m)
            COMMON /X/ xk, xtildek1
            COMMON /V/ vk, vtilde, vtildek1
            COMMON /G/ F, GVe, GW

            DO i = 1, N
            DO j = 1, m

                  xk(i,j) = -1.0 + 2*RAND(0)
                  vk(i,j) = 0.0
                  F(i,j) = 0.0                  

            END DO
            END DO

      END SUBROUTINE INIT

C     GaussianV: update the velocities with the gaussian variable
c           modifies vtilde
      SUBROUTINE GaussianV(eta, sdn, pi)
            PARAMETER(N = 50, m = 1)
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION vk(N,m),vtilde(N,m),vtildek1(N,m)
            COMMON /V/ vk, vtilde, vtildek1

            DO i = 1, N
            DO j = 1, m
                  vtilde(i,j) = eta * vk(i,j) + sdn * Gauss(pi)       
            END DO
            END DO

      END SUBROUTINE GaussianV

C     UPDATE: update the positions and velocities
c           modifies xtildek1 and vtildek1
      SUBROUTINE UPDATE(tstep, alpha, beta, t, a)
            PARAMETER(N = 50, m = 1)
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION xk(N,m), xtildek1(N,m),
     &                vk(N,m),vtilde(N,m),vtildek1(N,m),
     &                F(N,m), GVe(m), GW(m)
            COMMON /X/ xk, xtildek1
            COMMON /V/ vk, vtilde, vtildek1
            COMMON /G/ F, GVe, GW

            DO i = 1, N
                  vtilde(i,:) = vtilde(i,:) + alpha*F(i,:)*tstep/2
                  xtildek1(i,:) = xk(i,:) + alpha*vtilde(i,:)*tstep
            END DO

            CALL GRAD_H(.TRUE., beta, t, a)

            DO i = 1, N
                  vtildek1(i,:) = vtilde(i,:) + alpha*F(i,:)*tstep/2
            END DO

      END SUBROUTINE UPDATE  

C     GRAD_H: gradient of the Hamiltonian (force)
c           modifies F, GVe and GW
      SUBROUTINE GRAD_H(next, beta, t, a)
            PARAMETER(N = 50, m = 1)
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION x(N,m), xk(N,m), xtildek1(N,m),
     &                F_aux(N,N,m), F(N,m), GVe(m), GW(m)
            LOGICAL next
            COMMON /G/ F, GVe, GW
            COMMON /X/ xk, xtildek1

            IF (next) THEN
                  x = xtildek1
            ELSE
                  x = xk
            END IF

            DO i = 1, N
                  DO j = 1, i-1
                        CALL GRAD_W(x(i,:), x(j,:))
                        F_aux(i,j,:) = - GW
                  END DO
            END DO

            DO i = 1, N

                  F(i,:) = 0.0

                  DO j = 1,i-1
                        F(i,:) = F(i,:) + F_aux(i,j,:)
                  END DO

                  DO j = i+1,N
                        F(i,:) = F(i,:) - F_aux(j,i,:)
                  END DO

                  F(i,:) = F(i,:) / N

                  CALL GRAD_Ve(x(i,:), beta, t, a)
                  F(i,:) = F(i,:) - GVe

                  F(i,:) = F(i,:) / N

            END DO

      END SUBROUTINE GRAD_H

C     GRAD_Ve: gradient of the external potential
c           modifies GVe
      SUBROUTINE GRAD_Ve(x, beta, t, a)
            PARAMETER(N = 50, m = 1)
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION x(m), F(N,m), GVe(m), GW(m)
            COMMON /G/ F, GVe, GW

            Gve = 0
c                 Gradient of V(x) = ||x||^2 / (2 * beta) Beta - Hermite
            ! GVe = x / beta
c                 Gradient of V(x) = 1/4 x^4 + 1/2 x^2 t
            ! DO i = 1, m
            !       GVe(i) = GVe(i) + x(i)**3 + x(i)*t
            ! END DO
c                 Gradient of V(x) = t/(2α) x^(2α)
            DO i = 1, m
                  GVe(i) = GVe(i) + t/(2*a) * x(i)**(2*a-1)
            END DO

      END SUBROUTINE GRAD_Ve

C     GRAD_W: gradient of the interaction potential
c           modifies GW
      SUBROUTINE GRAD_W(x, y)
            PARAMETER(N = 50, m = 1)
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION x(m), y(m), v(m),
     &                F(N,m), GW(m), GVe(m)
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

C     (2-FUNCTION) Gauss: gaussian variable
c           return a standard gaussian variable, scalar
      FUNCTION Gauss(pi)
            IMPLICIT REAL*8 (A-H,O-Z)
            Gauss = 1000
            DO WHILE (Gauss > 100)
                  Gauss = sqrt(-2.*LOG(RAND()))*COS(2.*pi*RAND())
            END DO
            RETURN
      END FUNCTION Gauss

C     (1-FUNCTION) PROB: calculate the acceptance probability
c           return the acceptance probability, scalar
      FUNCTION PROB(beta, t, a)
            PARAMETER(N = 50, m = 1)
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION vk(N,m),vtilde(N,m),vtildek1(N,m)
            COMMON /V/ vk, vtilde, vtildek1

            Hi = H(.FALSE.,beta, t, a)
            Hf = H(.TRUE.,beta, t, a)

            PROB = EXP(-beta * (Hf-Hi))

            DO i = 1, N
                  dtilde = DOT_PRODUCT(vtildek1(i,:),vtildek1(i,:))
                  dk = DOT_PRODUCT(vk(i,:),vk(i,:))
                  PROB = PROB * EXP(-beta * (dtilde - dk) / 2)
            END DO
       
            RETURN 
            END FUNCTION PROB
      
C     (3-FUNCTION) H: Hamiltonian
c           return the Hamiltonian, scalar
      FUNCTION H(next, beta, t, a)
            PARAMETER(N = 50, m = 1)
            IMPLICIT REAL*8 (A-H,O-Z)
            LOGICAL next
            DIMENSION x(N,m), xk(N,m), xtildek1(N,m)
            COMMON /X/ xk, xtildek1

            if (next) then
                  x = xtildek1
            else
                  x = xk
            end if

            H = 0.0

            DO i = 1, N
                  H = H + Ve(x(i,:), beta, t, a)
                  DO j = i+1, N
                        H = H + W(x(i,:), x(j,:)) / (2*N)
                  END DO
            END DO

            H = H / N

            RETURN
      END FUNCTION H

C     (3-FUNCTION) Ve: external potential
c           return the external potential, scalar
      FUNCTION Ve(x, beta, t, a)
            PARAMETER(N = 50, m = 1)
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION x(m)

c                 V(x) = ||x||^2 / (2*beta) Beta - Hermite
            ! Ve = DOT_PRODUCT(x,x) / (2*beta)
c                 V(x) = 1/4 x^4 + 1/2 x^2 t
            ! Ve = 0.25*DOT_PRODUCT(x,x)**4 + 0.5*DOT_PRODUCT(x,x)**2 *t
c                 V(x) = t/(2α) x^(2α)
            Ve = t/(2*a) * DOT_PRODUCT(x,x)**(2*a)

            RETURN
      END FUNCTION Ve

C     (3-FUNCTION) W: interaction potential
c           return the interaction potential, scalar
      FUNCTION W(x, y)
            PARAMETER(N = 50, m = 1)
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION x(m), y(m)

c                 W(x) = -log(||x||)
            W = -LOG(NORM2(x-y))

            RETURN
      END FUNCTION W

C#######################################################################