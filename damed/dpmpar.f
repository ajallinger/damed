*-----------------------------------------------------------------------
* wrapper around FORTRAN 95 double precision machine constants
*-----------------------------------------------------------------------
      FUNCTION DPMPAR (INPT)
       IMPLICIT NONE
       INTEGER INPT
       DOUBLE PRECISION DPMPAR, D
       SAVE D
       DATA D / 1.D0 /

       DPMPAR = 0.D0
       IF (INPT .EQ. 1) THEN
         DPMPAR = EPSILON(D)  ! machine precision
       ELSE IF (INPT .EQ. 2) THEN
         DPMPAR = TINY(D)     ! smallest positive number
       ELSE IF (INPT .EQ. 3) THEN
         DPMPAR = HUGE(D)     ! largest positive number
       END IF
       RETURN
      END  ! of DPMPAR
