      PROGRAM TESTF
      
      INTEGER NDATA
      PARAMETER (NDATA = 80)
      REAL MODEL(NDATA), DATA(NDATA), X(NDATA)
      REAL E(NDATA), IMAGE(NDATA), GAUSS2
      INTEGER SEED/-8269/
      REAL  SIGMA/0.01/
      REAL CAIM, RMAX, C, TEST, CNEW, S, RNEW, ACC, SNEW, SUMF
      INTEGER MODE/10/
      INTEGER LEVEL/20/
      COMMON /MECOMP/ NJ,MJ,NK,MK,KA(40),KB(40),KC(40),KD(40),
     * IOUT,IIN,L0,L1,M0,M10,M11,M20,M21,M3
      COMMON /MECOMS/ ST(100000)
      DATA IOUT/6/

      CALL MEMCORE(100000,NDATA,1,1,NDATA,IER)
      IF(IER.NE.0) STOP 'FAILED TO FORMAT CORE'

      DO I=1,NDATA
         X(I) = I+1
         MODEL(I) = 1.
      END DO
      
      MODEL(40) = 3.
      MODEL(41) = 3.
      MODEL(39) = 2.
      MODEL(42) = 2.
      MODEL(38) = 1.5
      MODEL(43) = 1.5
  
      CALL OP(MODEL,DATA,NDATA)

      DO I = 1, NDATA
         DATA(I)  = GAUSS2(DATA(I),SIGMA,SEED)
         E(I)     = SIGMA
         IMAGE(I) = 0.5
      END DO

      DO I = 1, NDATA
         ST( KB(1)+I-1) = IMAGE(I)
         ST(KB(21)+I-1) = DATA(I)
         ST(KB(22)+I-1) = 2./(E(I)*E(I)*NDATA)
      END DO

      CAIM = 1.
      RMAX = 0.2
      
      DO I = 1, 10
         CALL MEMPRM( MODE, LEVEL, CAIM, RMAX, 1., ACC, C, TEST,
     &        CNEW, S, RNEW, SNEW, SUMF)
      END DO

      DO I = 1, NDATA
         IMAGE(I) = ST(KB(1)+I)
      END DO
      END

      SUBROUTINE OPUS(J, K)
      COMMON /MECOMP/ NJ,MJ,NK,MK,KA(40),KB(40),KC(40),KD(40),
     * IOUT,IIN,L0,L1,M0,M10,M11,M20,M21,M3
      COMMON /MECOMS/ ST(100000)
      INTEGER NDATA
      PARAMETER (NDATA = 80)

      WRITE(*,*) '     OPUS ',J,' ---> ',K
      CALL OP(ST(KB(J)),ST(KB(K)), NDATA)
      RETURN
      END


      SUBROUTINE OP(MODEL, DATA, NPT) 
      INTEGER NSP, NPT
      REAL DATA(NPT), MODEL(NPT)
      PARAMETER (NSP=9)
      REAL WGT(NSP)
      DATA WGT/1.,2.,4.,5.,6.,5.,4.,2.,1./
      INTEGER K
      REAL SUM, WNORM
      
      WNORM = 0
      DO I=1,NSP
         WNORM = WNORM + WGT(I)
      END DO
      DO I=1,NPT
         SUM = 0.
         K = I - NSP/2
         DO J = 1, NSP
            IF(K.LT.1) THEN
               SUM = SUM + WGT(J)*MODEL(1)
            ELSE IF(K .GT. NPT) THEN
               SUM = SUM + WGT(J)*MODEL(NPT)
            ELSE
               SUM = SUM + WGT(J)*MODEL(K);
            END IF
            K = K + 1
         END DO
         DATA(I) = SUM/WNORM
      END DO
      RETURN
      END

      SUBROUTINE TROPUS(K, J)
      COMMON /MECOMP/ NJ,MJ,NK,MK,KA(40),KB(40),KC(40),KD(40),
     * IOUT,IIN,L0,L1,M0,M10,M11,M20,M21,M3
      COMMON /MECOMS/ ST(100000)
      INTEGER NDATA
      PARAMETER (NDATA = 80)

      WRITE(*,*) '   TROPUS ',J,' <--- ',K
      CALL TR(ST(KB(K)),ST(KB(J)), NDATA)
      RETURN
      END


      SUBROUTINE TR(DATA, MODEL, NPT) 
      INTEGER NSP, NPT
      REAL DATA(NPT), MODEL(NPT)
      PARAMETER (NSP=9)
      REAL WGT(NSP)
      DATA WGT/1.,2.,4.,5.,6.,5.,4.,2.,1./
      INTEGER K
      REAL SUM, WNORM
      
      WNORM = 0
      DO I=1,NSP
         WNORM = WNORM + WGT(I)
      END DO
      DO I=1,NPT
         MODEL(I) = 0.
      END DO
      DO I=1,NPT
         K = I - NSP/2
         DO J = 1, NSP
            IF(K.LT.1) THEN
               MODEL(I) = MODEL(I) + WGT(J)*DATA(1)
            ELSE IF(K .GT. NPT) THEN
               MODEL(I) = MODEL(I) + WGT(J)*DATA(NPT)
            ELSE
               MODEL(I) = MODEL(I) + WGT(J)*DATA(K)
            END IF
            K = K + 1
         END DO
         MODEL(I) = MODEL(I)/WNORM
      END DO
      RETURN
      END


