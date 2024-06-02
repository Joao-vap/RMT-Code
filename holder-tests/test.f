      PROGRAM test
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (N = 100, t0 = 0.1, t1 = 1.0, a = 2.0)
      DIMENSION x(2), GVe(2)
      COMMON /G/ GVe

      OPEN(1,FILE='../Complex/testr.dat',STATUS='UNKNOWN')
      OPEN(2,FILE='../Complex/testi.dat',STATUS='UNKNOWN')
      OPEN(3,FILE='../Complex/test.dat',STATUS='UNKNOWN')

C     loop grid and calculate gradient
      DO i = 1, N
            DO j = 1, N
                  x(1) = -10.0 + 20.0*(i-1)/(N-1)
                  x(2) = -10.0 + 20.0*(j-1)/(N-1)
                  CALL NUMERICAL_GRAD_Ve(x,t0,t1,a)
                  WRITE(1,*) GVe(1)
                  WRITE(2,*) GVe(2)
                  WRITE(3,*) Ve(x,a,t0, t1)
            END DO
      END DO

      END PROGRAM test


      SUBROUTINE NUMERICAL_GRAD_Ve(x, t0, t1, a)
            PARAMETER(N = 100)
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION x(2), GVe(2), x1(2), x2(2)
            COMMON /G/ GVe

            eps = 1e-6

            DO i = 1, 2
                  x1 = x
                  x2 = x
                  x1(i) = x1(i) + eps
                  x2(i) = x2(i) - eps

                  GVe(i) = (Ve(x1,a,t0,t1)-Ve(x2,a,t0,t1))/(2*eps)
            END DO

      END SUBROUTINE NUMERICAL_GRAD_Ve

      SUBROUTINE GRAD_Ve(x, t, a)
            PARAMETER(N = 100)
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION x(2), GVe(2)
            COMPLEX*16 Xc, Xcc, VGVe
            COMMON /G/ GVe

      ! GRadient pf Ve = |z|^(2a) - t*Re(z^a)
            Xc = CMPLX(x(1), x(2), KIND=16)
            Xcc = CONJG(Xc)
            VGVe = 2*a*(Xc*Xcc)**(a-1)*Xc + t*a*Xcc**(a-1)
            GVe(1) = REAL(VGVe, KIND=8)
            GVe(2) = AIMAG(VGVe)

      END SUBROUTINE GRAD_Ve

C     (3-FUNCTION) Ve: external potential
c           return the external potential, scalar
      FUNCTION Ve(x, a, t0, t1)
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION x(m)
            COMPLEX*16 Xc

            Xc = CMPLX(x(1), x(2), KIND=16)
            Ve=1/t0*(ABS(Xc)**(a)-2*REAL(Xc**(a+1)/(a+1)+t1*Xc, KIND=8))

            RETURN
      END FUNCTION Ve