* test program
      PROGRAM TDAMED
       IMPLICIT NONE
       INTEGER N, J
       PARAMETER (N = 10)
       DOUBLE PRECISION X(N), AVE, DEV
       DOUBLE PRECISION Y(N), YMED, DYDX(N), TAU, DMDX
       INTEGER IERR

       DATA X / 3., 8., 3., -1., 0., 5., 4., 8., 9., 2. /
       DATA DYDX / 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. /

       AVE = 0.D0
       DO J = 1, N
         Y(J) = X(J)
         AVE = AVE + Y(J)
       END DO
       AVE = AVE / DBLE(N)
       DEV = 0.D0
       DO J = 1, N
         DEV = DEV + (Y(J) - AVE)**2
       END DO
       DEV = SQRT(DEV / DBLE(N))

       PRINT*, 'x is: ', X
       TAU = DEV / 10.D0
       PRINT *, 'using tau: ', TAU
       CALL DAMED (Y, DYDX, N, TAU, YMED, DMDX, IERR)
       PRINT *, 'damed returns #', IERR, ' median: ', YMED
       PRINT*, 'derivative: ', DMDX

      END  ! of test program
