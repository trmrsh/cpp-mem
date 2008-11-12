*
      SUBROUTINE MOVER( SOURCE, DEST, N)
      REAL SOURCE(1), DEST(1)
      DO I=1,N
        DEST(I) = SOURCE(I)
      END DO
      RETURN
      END
