      SUBROUTINE SLARTG( F, G, CS, SN, R )
*
*  -- LAPACK auxiliary routine (version 1.0b) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      REAL               CS, F, G, R, SN
*     ..
*
*  Purpose
*  =======
*
*  SLARTG generate a plane rotation so that
*
*     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
*     [ -SN  CS  ]     [ G ]     [ 0 ]
*
*  This is a faster version of the BLAS1 routine SROTG, except for
*  the following differences:
*     F and G are unchanged on return.
*     If G=0, then CS=1 and SN=0.
*     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
*        floating point operations (saves work in SBDSQR when
*        there are zeros on the diagonal).
*
*  Arguments
*  =========
*
*  F       (input) REAL
*          The first component of vector to be rotated.
*
*  G       (input) REAL
*          The second component of vector to be rotated.
*
*  CS      (output) REAL
*          The cosine of the rotation.
*
*  SN      (output) REAL
*          The sine of the rotation.
*
*  R       (output) REAL
*          The nonzero component of the rotated vector.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )
      REAL               ONE
      PARAMETER          ( ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      REAL               T, TT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
*     ..
*     .. Executable Statements ..
*
      IF( G.EQ.ZERO ) THEN
         CS = ONE
         SN = ZERO
         R = F
      ELSE IF( F.EQ.ZERO ) THEN
         CS = ZERO
         SN = ONE
         R = G
      ELSE
         IF( ABS( F ).GT.ABS( G ) ) THEN
            T = G / F
            TT = SQRT( ONE+T*T )
            CS = ONE / TT
            SN = T*CS
            R = F*TT
         ELSE
            T = F / G
            TT = SQRT( ONE+T*T )
            SN = ONE / TT
            CS = T*SN
            R = G*TT
         END IF
      END IF
      RETURN
*
*     End of SLARTG
*
      END
