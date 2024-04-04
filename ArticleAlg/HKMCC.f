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
C     We consider gases that in a subspace S of R^n, with dimension m
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
C     \mu_N = 1/N sum_{i=1}^{N} \delta_{x_i} 
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
C     alpha: time scale, scalar
C     gamma: friction, scalar
C     nsteps: number of steps, scalar
C     eta: coeficient exp{-\gamma_N \alpha_N \delta t}, scalar
C     sdn: coeficient \sqrt{(1 - \eta^2) / \beta_N}, scalar
C     pi: coeficient \pi, scalar
C     tstep: time step, scalar
C     niter: number of iterations to register, scalar
C     p: acceptance probability, scalar
C     t: coeficient of the external potential, scalar
C     a: coeficient of the external potential, scalar
C     Hi: Hamiltonian at k, scalar
C     Hf: Hamiltonian at k+1, scalar
C     accSteps: accepted steps, scalar
C
C     xk: position at k, vector of size (N,m)
C     xtildek1: candidate positions for k+1, vector of size (N,m)
C     vk: velocities at k, vector of size (N,m)
C     vtilde: gaussian modified velocities, vector of size (N,m)
C     vtildek1: candidate velocities for k+1, vector of size (N,m)
C     GH: gradient of Hamiltonian, vector of size (N,m)
C     GVe: gradient of the external potential, vector of size (m)
C     GW: gradient of the interaction potential, vector of size (m)
C     GHold: gradient of Hamiltonian at k-1, vector of size (N,m)

      PARAMETER (N = 100, m = 2, beta = 2.0, alpha = 1.0)

C#######################################################################

      DIMENSION xk(N,m), xtildek1(N,m),
     &      vk(N,m),vtilde(N,m),vtildek1(N,m),
     &      GH(N,m), GHold(N, m), GVe(m), GW(m)

      COMMON /X/ xk, xtildek1
      COMMON /V/ vk, vtilde, vtildek1
      COMMON /G/ GH, GHold, GVe, GW
      COMMON /Hk/ Hi, Hf
      COMMON /r/ accSteps

C#######################################################################

      Hi = 0.0
      Hf = 0.0
      accSteps = 0.0

      nsteps = 2000000
      niter = 500
      tstep = 0.5
      gamma = 1.0

      t = 2.0
      a = 2.0

C#######################################################################

      eta = EXP(-gamma * alpha * tstep)
      sdn = SQRT((1 - eta**2) / beta) / N
      pi = ACOS(-1.d0)

      call srand(987991650)

C#######################################################################
C#######################################################################
C     Step 1: Initialization of (x0, v0)
C     We take x0 = U(-1,1) and v0 = 0
C     We also initialize xk = x0 and vk = v0

      CALL INIT()
      Hi = H(.FALSE.,beta, t, a)

C ---------------------------------------------------------------------

      OPEN(1,FILE='Complex/a2t2.dat',STATUS='UNKNOWN')
      OPEN(2,FILE='./H.dat',STATUS='UNKNOWN')

C ---------------------------------------------------------------------

      DO 10 k = 1, nsteps

c     to ensure ergocity we need to pick randomly Dt
      ttstep = tstep + 0.01*(RAND(0) - 0.5)

C ---------------------------------------------------------------------
C     Step 2: Update the velocities with
C     vtilde_k = \eta vk_k + \sqrt{(1 - \eta^2) / \beta_N} \G_k
C     where \G_k is a standard Gaussian vector
C     and \eta = \exp{-\gamma_N \alpha_N \delta t}
      
      CALL L2_ORNSUHLEN(eta, sdn, pi)

C ---------------------------------------------------------------------
C     Step 3: Calculate 
C     vtildehalf_k = vtilde_k - \alpha_N \Nabla H_N(xk) timestep/2
C     xtildek1_k = xk + \alpha_N vtildehalf_k timestep
C     vtilde = vtildehalf_k - \alpha_N \Nabla H_N(xtildek1_k) timestep/2
C
C     Note that we update a half step of the velocities then we update
C     the positions and then we update the velocities again

      CALL L1_VERLET(ttstep, alpha, beta, t, a)

C ---------------------------------------------------------------------
C     Step 4: define the acceptance probability
C     prob=1 ^ exp{-\beta_N (H_N(xtildek1_k)-H_N(xk)+vtilde_k^2/2-vk^2/2)}
C     accept or reject the candidate with probability p

      CALL METROPOLIS(beta, t, a) 

C ---------------------------------------------------------------------
c     Save some data every niter steps, after nsteps/2
      IF (MOD(k,niter) == 0) THEN
            WRITE(*,*) k, XK(N,:), Hi, accSteps/k
c            Hi, accSteps/k
            IF (k > nsteps/5) THEN
                  WRITE(1,*) xk
                  WRITE(2,*) Hi
            END IF
      END IF

C ---------------------------------------------------------------------

 10   END DO
      
      CLOSE(1)
      CLOSE(2)

C#######################################################################

      END PROGRAM HMC

C#######################################################################
C     Subroutines:
C     INIT: initialization of (x0, v0)
c           modifies xk, vk, GH
      SUBROUTINE INIT()
            PARAMETER(N = 100, m = 2)
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION xk(N,m), xtildek1(N,m),
     &                vk(N,m), vtilde(N,m), vtildek1(N,m),
     &                GH(N,m), GHold(N, m), GVe(m), GW(m)
            COMMON /X/ xk, xtildek1
            COMMON /V/ vk, vtilde, vtildek1
            COMMON /G/ GH, GHold, GVe, GW

            DO i = 1, N
            DO j = 1, m

                  xk(i,j) = -1.0 + 2*RAND(0)
                  vk(i,j) = 0.0
                  GH(i,j) = 0.0                  

            END DO
            END DO

      END SUBROUTINE INIT

C     L2_ORNSUHLEN: update the velocities with the gaussian variable
c           modifies vtilde
      SUBROUTINE L2_ORNSUHLEN(eta, sdn, pi)
            PARAMETER(N = 100, m = 2)
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION vk(N,m),vtilde(N,m),vtildek1(N,m)
            COMMON /V/ vk, vtilde, vtildek1

            DO i = 1, N
            DO j = 1, m
                  vtilde(i,j) = eta * vk(i,j) + sdn * Gauss(pi)       
            END DO
            END DO

      END SUBROUTINE L2_ORNSUHLEN

C     L1_VERLET: update the positions and velocities
c           modifies xtildek1 and vtildek1
      SUBROUTINE L1_VERLET(tstep, alpha, beta, t, a)
            PARAMETER(N = 100, m = 2)
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION xk(N,m), xtildek1(N,m),
     &                vk(N,m),vtilde(N,m),vtildek1(N,m),
     &                GH(N,m), GHold(N, m), GVe(m), GW(m)
            COMMON /X/ xk, xtildek1
            COMMON /V/ vk, vtilde, vtildek1
            COMMON /G/ GH, GHold, GVe, GW

            DO i = 1, N
                  vtilde(i,:) = vtilde(i,:) + alpha*GH(i,:)*tstep/2
                  xtildek1(i,:) = xk(i,:) + alpha*vtilde(i,:)*tstep
            END DO

            GHold = GH
            CALL GRAD_H(beta, t, a)
      
            DO i = 1, N
                  vtildek1(i,:) = vtilde(i,:) + alpha*GH(i,:)*tstep/2
            END DO

      END SUBROUTINE L1_VERLET  

C     GRAD_H: gradient of the Hamiltonian (force)
c           modifies GH, GVe and GW
      SUBROUTINE GRAD_H(beta, t, a)
            PARAMETER(N = 100, m = 2)
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION xk(N,m), xtildek1(N,m),
     &                GH_aux(N,N,m), GH(N,m), GHold(N, m), GVe(m), GW(m)
            COMMON /G/ GH, GHold, GVe, GW
            COMMON /X/ xk, xtildek1

            DO i = 1, N
                  DO j = 1, i-1
                        CALL GRAD_W(xtildek1(i,:), xtildek1(j,:))
                        GH_aux(i,j,:) = - GW(:)
                  END DO
            END DO

            DO i = 1, N

                  GH(i,:) = 0.0

                  DO j = 1,i-1
                        GH(i,:) = GH(i,:) + GH_aux(i,j,:)
                  END DO

                  DO j = i+1,N
                        GH(i,:) = GH(i,:) - GH_aux(j,i,:)
                  END DO

                  GH(i,:) = GH(i,:) / N 

                  CALL GRAD_Ve(xtildek1(i,:), beta, t, a)
                  GH(i,:) = GH(i,:) - GVe(:)

                  GH(i,:) = GH(i,:) / N

            END DO

      END SUBROUTINE GRAD_H

C     GRAD_Ve: gradient of the external potential
c           modifies GVe
      SUBROUTINE GRAD_Ve(x, beta, t, a)
            PARAMETER(N = 100, m = 2)
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION x(m), GH(N,m), GHold(N, m), GVe(m), GW(m)
            COMPLEX*16 Xc, Xcc, VGVe
            COMMON /G/ GH, GHold, GVe, GW

            GVe = 0
c           -----------------------------------------------------------
c                 Gradient of V(x) = ||x||^2 / (2 * beta) Beta - Hermite
c                 Quadratic potential
            ! GVe = x / beta
c           -----------------------------------------------------------
c                 Gradient of V(x) = 1/4 x^4 + 1/2 x^2 t
c                 Quartic potential
            ! DO i = 1, m
            !       GVe(i) = GVe(i) + (x(i)**3 + x(i)*t)/beta
            ! END DO
c           -----------------------------------------------------------
c                 Gradient of V(x) = t x^(2α)
c                 Monic potential
            ! DO i = 1, m
            !        GVe(i) = GVe(i) + (2*a*t * x(i)**(2*a-1)) / beta
            ! END DO
c           -----------------------------------------------------------
c                 Gradient of V(z) = ||z||^2a - Re(t z^a)
c                 Complex potential

            Xc = CMPLX(x(1), x(2), KIND=16)
            Xcc = CONJG(Xc)
            VGVe = (2*a*(Xc*Xcc)**(a-1)*Xc + t*a*Xcc**(a-1)) / beta
            GVe(1) = REAL(VGVe, KIND=8)
            GVe(2) = AIMAG(VGVe)

c           -----------------------------------------------------------

            ! CALL NUMERICAL_GRAD_Ve(x, beta, t, a)

      END SUBROUTINE GRAD_Ve

C     NUMERICAL_GRAD_Ve: numerical gradient of the external potential
c           modifies GVe
      SUBROUTINE NUMERICAL_GRAD_Ve(x, beta, t, a)
            PARAMETER(N = 100, m = 2)
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION x(m), GH(N,m), GHold(N, m), GVe(m), GW(m),
     +               x1(m), x2(m)
            COMMON /G/ GH, GHold, GVe, GW

            eps = 1e-6

            DO i = 1, m
                  x1 = x
                  x1(i) = x1(i) + eps
                  x2 = x
                  x2(i) = x2(i) - eps

            GVe(i) = GVe(i) + (Ve(x1,beta,t,a)-Ve(x2,beta,t,a)) /(2*eps)
            END DO

      END SUBROUTINE NUMERICAL_GRAD_Ve

C     GRAD_W: gradient of the interaction potential
c           modifies GW
      SUBROUTINE GRAD_W(x, y)
            PARAMETER(N = 100, m = 2)
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION x(m), y(m), v(m),
     &                GH(N,m), GHold(N,m), GW(m), GVe(m)
            COMMON /G/ GH, GHold, GVe, GW

c          -----------------------------------------------------------
c                 Gradient of W(x) = -log(||x||)
c                 Logarithmic interaction - 2D restricted to 1D
            v = x-y
            GW = -v / NORM2(v)**2
c          -----------------------------------------------------------

      END SUBROUTINE GRAD_W

C     METROPOLIS: define the acceptance probability
c           modifies xk, vk
      SUBROUTINE METROPOLIS(beta, t, a)
            PARAMETER(N = 100, m = 2)
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION xk(N,m), xtildek1(N,m),
     &                vk(N,m),vtilde(N,m),vtildek1(N,m)
            COMMON /Hk/ Hi, Hf
            COMMON /X/ xk, xtildek1
            COMMON /V/ vk, vtilde, vtildek1
            COMMON /r/ accSteps

            p = PROBLOG(beta, t, a)
            
            IF (LOG(RAND(0)) <= p) THEN
                  xk = xtildek1
                  vk = vtildek1
                  Hi = Hf
                  accSteps = accSteps + 1
            ELSE
                  vk = -vtildek1
                  GH = GHold
            END IF

      END SUBROUTINE METROPOLIS


C#######################################################################
C#######################################################################
C     Functions:
C     1-FUNCTION ARE CALLED IN MAIN
C     2-FUNCTION ARE CALLED IN SUBROUTINES
C     3-FUNCTION ARE CALLED IN FUNCTIONS

C     (2-FUNCTION) Gauss: gaussian variable
c         Im using the Box-Muller method, described first
c         in the paper "A note on the generation of random normal
c         deviates" by George Box and Mervin Muller.
c           return a standard gaussian variable, scalar
      FUNCTION Gauss(pi)
            IMPLICIT REAL*8 (A-H,O-Z)
            REAL*8 u1, u2

            u1 = 0
            DO WHILE (u1 == 0)
                  u1 = RAND()
            END DO
            u2 = RAND()

            Gauss = sqrt(-2.*LOG(u1))*COS(2.*pi*u2)

            RETURN
      END FUNCTION Gauss

C     (1-FUNCTION) PROBLOG: calculate the log of acceptance probability
c           return the log of acceptance probability, scalar
      FUNCTION PROBLOG(beta, t, a)
            PARAMETER(N = 100, m = 2)
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION vk(N,m),vtilde(N,m),vtildek1(N,m)
            COMMON /V/ vk, vtilde, vtildek1
            COMMON /Hk/ Hi, Hf

            Hf = H(.TRUE.,beta, t, a)

            PROBLOG = - beta * (Hf-Hi)

            DO i = 1, N
C                 Energia cinética
                  aux_Kf = DOT_PRODUCT(vtildek1(i,:),vtildek1(i,:))
                  aux_Ki = DOT_PRODUCT(vk(i,:),vk(i,:))
                  PROBLOG = PROBLOG - beta * (aux_Kf - aux_Ki)/2
            END DO
       
            RETURN 
            END FUNCTION PROBLOG
      
C     (3-FUNCTION) H: Hamiltonian
c           return the Hamiltonian, scalar
      FUNCTION H(next, beta, t, a)
            PARAMETER(N = 100, m = 2)
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
                        H = H + W(x(i,:), x(j,:)) / N
                  END DO
            END DO

            H = H / N

            RETURN
      END FUNCTION H

C     (3-FUNCTION) Ve: external potential
c           return the external potential, scalar
      FUNCTION Ve(x, beta, t, a)
            PARAMETER(N = 100, m = 2)
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION x(m)
            COMPLEX*16 Xc
            Ve = 0.0
            ! COMPLEX*16 Xc

c           -----------------------------------------------------------            
c                 V(x) = ||x||^2 / (2 * beta)
c                 Quadratic potential
            ! Ve = DOT_PRODUCT(x,x) / (2*beta)
c           -----------------------------------------------------------            
c                 V(x) = 1/4 x^4 + 1/2 x^2 t
c                 Quartic potential
            ! DO i = 1, m
            !       Ve = Ve + (x(i)**4 / 4 + x(i)**2 * t / 2)/beta
            ! END DO
c           -----------------------------------------------------------            
c                 V(x) = t x^(2α)
c                 Monic potential
            ! DO i = 1, m
            !       Ve = Ve + (t * x(i)**(2*a))/(beta) 
            ! END DO
c           ----------------------------------------------------------- 
c                 V(z) = ||z||^2a - Re(t z^a)
c                 Complex potential
            Xc = CMPLX(x(1), x(2), KIND=16)
            Ve = ABS(Xc)**(2*a) - t*REAL(Xc**a)
c           -----------------------------------------------------------

            RETURN
      END FUNCTION Ve

C     (3-FUNCTION) W: interaction potential
c           return the interaction potential, scalar
      FUNCTION W(x, y)
            PARAMETER(N = 100, m = 2)
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION x(m), y(m)

c           -----------------------------------------------------------
c                 W(x) = -log(||x||)
c                 Logarithmic interaction - 2D restricted to 1D
            W = -LOG(NORM2(x-y))
c           -----------------------------------------------------------

            RETURN
      END FUNCTION W

C#######################################################################