      REAL*4 FUNCTION CHISQR(N, O, C, E)
*
*  Computes reduced Chi-squared
*
      REAL*4 O(1), C(1), E(1)
      CHISQR = 0.
      NSUM = 0
      DO I = 1, N
        IF(E(I).GT.0.) THEN
          ADD = (O(I) - C(I)) / E(I)
          CHISQR = CHISQR + ADD*ADD
          NSUM = NSUM + 1
        END IF
      END DO
      IF(NSUM.GT.0) CHISQR = CHISQR/NSUM
      RETURN
      END
