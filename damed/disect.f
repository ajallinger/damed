*-----------------------------------------------------------------------
* disect - solve an equation by bisection, allowing an interval f = 0.
* by Andy Allinger, 2018-2023, released to the public domain.
*
*    Permission  to  use, copy, modify, and distribute this software and
*    its documentation  for  any  purpose  and  without  fee  is  hereby
*    granted,  without any conditions or restrictions.  This software is
*    provided "as is" without express or implied warranty.
*
*___Name______Type_______In/Out___Description___________________________
*   FUNC      Function   In       The function to be solved
*                                   It has the calling sequence:
*                                   F = FUNC (X, DARG, N, TAU)
*   B         Double     Both     One end of the interval to search
*   C         Double     Both     The other end of the interval
*                                   On return, B and C have been changed
*                                   to bounds on the solution.
*   ATOL      Double     In       Absolute error tolerance for solution
*   IERR      Integer    Out      Error code
*                                   0 = no errors
*                                   1 = no zero crossing
*                                   2 = too many iterations
*   DARG(N)   Double     Thru     Array passed to user function
*   N         Integer    Thru     Length of DARG
*   TAU       Double     Thru     Argument passed to user function
*-----------------------------------------------------------------------
      SUBROUTINE DISECT (FUNC, B, C, ATOL, IERR, DARG, N, TAU)
       IMPLICIT NONE
       INTEGER IERR, N
       DOUBLE PRECISION FUNC, B, C, ATOL, DARG(N), TAU
       EXTERNAL FUNC

*           Constants
       INTEGER MAXIT
       DOUBLE PRECISION ZERO, HALF
       PARAMETER ( MAXIT = 9000,
     $             ZERO = 0.0D0,
     $             HALF = 0.5D0 )

*           Local variables
       INTEGER KOUNT
       DOUBLE PRECISION F, FB, FC,               ! f values
     $                  XL, XG, ROOTL, ROOTG     ! x values

*-----------------------------------------------------------------------
*        Begin.
*-----------------------------------------------------------------------
       IERR = 0
       IF (B > C) THEN  ! ensure that B < C
         F = B
         B = C
         C = F
       END IF

*          Check that a root is present
       IF (C - B < ATOL) RETURN
       FB = FUNC (B, DARG, N, TAU)
       FC = FUNC (C, DARG, N, TAU)
       KOUNT = 2

*         No change of sign
       IF (.NOT. ((FB < ZERO .AND. FC > ZERO)
     $      .OR.  (FB > ZERO .AND. FC < ZERO)) ) THEN
         IERR = 1
         RETURN
       END IF

*-----------------------------------------------------------------------
*        Find ROOTL, the least value in [B,C] such that f(ROOTL) = 0
*-----------------------------------------------------------------------
       XL = B
       XG = C

  10   ROOTL = HALF * (XL + XG)
       IF (XG - XL .LE. ATOL) GO TO 20
       F = FUNC (ROOTL, DARG, N, TAU)

       KOUNT = KOUNT + 1
       IF (KOUNT > MAXIT) THEN
         IERR = 2
         RETURN
       END IF

       IF (FB < ZERO) THEN  ! search for f(x) < 0
         IF (F < ZERO) THEN
           XL = ROOTL
         ELSE
           XG = ROOTL
         END IF
       ELSE                 ! search for f(x) > 0
         IF (F > ZERO) THEN
           XL = ROOTL
         ELSE
           XG = ROOTL
         END IF
       END IF
       GO TO 10

*-----------------------------------------------------------------------
*        Find ROOTG, the greatest value in [B,C] such that f(ROOTG) = 0
*-----------------------------------------------------------------------
  20   XL = B
       XG = C

  30   ROOTG = HALF * (XL + XG)
       IF (XG - XL .LE. ATOL) GO TO 40
       F = FUNC (ROOTG, DARG, N, TAU)

       KOUNT = KOUNT + 1
       IF (KOUNT > MAXIT) THEN
         IERR = 2
         RETURN
       END IF

       IF (FC < ZERO) THEN  ! search for f(x) < 0
         IF (F < ZERO) THEN
           XG = ROOTG
         ELSE
           XL = ROOTG
         END IF
       ELSE                ! search for f(x) > 0
         IF (F > ZERO) THEN
           XG = ROOTG
         ELSE
           XL = ROOTG
         END IF
       END IF
       GO TO 30

*             Final answer
  40   B = ROOTL
       C = ROOTG
       RETURN
      END  ! of disect
