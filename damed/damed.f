*-----------------------------------------------------------------------
* damed - differentiable approximation to the median
* by Andy Allinger, 2023, released to the public domain
*
*    Permission  to  use, copy, modify, and distribute this software and
*    its documentation  for  any  purpose  and  without  fee  is  hereby
*    granted,  without any conditions or restrictions.  This software is
*    provided "as is" without express or implied warranty.
*
* Intended for use in unconstrained optimization.
*
*___Name______Type______In/Out___Description____________________________
*   Y(N)      Double    In       Array of functions of x
*   DYDX(N)   Double    In       Partial derivatives dy[i]/dx
*   N         Integer   In       Length of Y
*   TAU       Double    In       Small number
*   YMED      Double    Out      Smooth approximate median of Y
*   DMDX      Double    Out      Derivative of the smooth median
*   IERR      Integer   Out      Error code
*                                 0 = no errors
*                                 1 = division by zero
*                                 2 = errors from root finding
*-----------------------------------------------------------------------
      SUBROUTINE DAMED (Y, DYDX, N, TAU, YMED, DMDX, IERR)
       IMPLICIT NONE
       INTEGER N, IERR
       DOUBLE PRECISION Y(N), DYDX(N), TAU, YMED, DMDX

       DOUBLE PRECISION ZERO, HALF, ONE, BIG
       PARAMETER (ZERO = 0.D0, HALF = 0.5D0, ONE = 1.D0, BIG = 1.D36)

       INTEGER J, IFLAG
       DOUBLE PRECISION ATOL, B, C
       DOUBLE PRECISION BEL, OOT, R, S, Z, ZMIN

       DOUBLE PRECISION COUNTF, DPMPAR
       LOGICAL DAFDIV
       EXTERNAL COUNTF

*-----------------------------------------------------------------------
*             Begin.
*-----------------------------------------------------------------------
       IERR = 0
       B = Y(1)                                  ! range of Y
       C = Y(1)
       DO J = 2, N
         IF (Y(J) < B) B = Y(J)
         IF (Y(J) > C) C = Y(J)
       END DO
       R = (C + B) * HALF
       S = (C - B) * HALF
       ATOL = DPMPAR(1) * MAX(R, S)

       CALL DISECT (COUNTF, B, C, ATOL, IFLAG, Y, N, TAU)
       IF (IFLAG .NE. 0) THEN
         IERR = 2
         RETURN
       END IF
       YMED = HALF * (B + C)

       OOT = ONE / TAU
       ZMIN = BIG
       DO J = 1, N
         Z = ((Y(J) - YMED) * OOT)**2
         IF (Z < ZMIN) ZMIN = Z
       END DO

       R = ZERO
       S = ZERO
       DO J = 1, N
         BEL = EXP(ZMIN - ((Y(J) - YMED) * OOT)**2 )
         R = R + BEL * DYDX(J)
         S = S + BEL
       END DO

       DMDX = ZERO
       IF (DAFDIV(R, S)) THEN
         DMDX = R / S
       ELSE
         IERR = 1
       END IF

       RETURN
      END  ! of damed


*-----------------------------------------------------------------------
*             The counting function, using approximate absolute value
*-----------------------------------------------------------------------
      FUNCTION COUNTF (YMED, Y, N, TAU)
       IMPLICIT NONE
       INTEGER N
       DOUBLE PRECISION YMED, Y(N), TAU, COUNTF

       DOUBLE PRECISION ZERO, ONE
       PARAMETER (ZERO = 0.D0, ONE = 1.D0)
       INTEGER J
       DOUBLE PRECISION C, S, T, OOT, XJ

       INTRINSIC ERF

       OOT = ONE / TAU
       C = ZERO                                 ! compensated addition
       S = ZERO
       DO J = 1, N
         XJ = ERF( (Y(J) - YMED) * OOT )
         T = S + XJ
         IF (ABS(S) .GE. ABS(XJ)) THEN
           C = C + ((S - T) + XJ)
         ELSE
           C = C + ((XJ - T) + S)
         END IF
         S = T
       END DO
       COUNTF = S + C

       RETURN
      END  ! of countf
