      program test
      double precision d(6,6),v(6),vec(6,6)
      d(1,1) = 1.;
      d(1,2) = 67788.2
      d(2,1) = d(1,2)
      d(1,3) = -0.024412
      d(3,1) = d(1,3)
      d(2,2) = 1.
      d(2,3) = 6628.52
      d(3,2) = d(2,3)
      d(3,3) = 1.
      call meml33(3,d,v,vec)
      print *, (v(i),i=1,3)
      end
      
      
      SUBROUTINE MEML33(ND,AMAT,VAL,VEC)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     
*     Purpose:
*     Diagonalisation of symmetric matrix
*     k      k         k
*     Solves      AMAT(i,j) * VEC(j)  = VAL  * VEC(i)
*     for eigenvalues VAL and eigenvectors VEC, for k=1,2,...,ND
*     
*     Parameters:
*     ARGUMENT TYPE  I/O  DIMENSION   DESCRIPTION
*     ND      I     I      -         Number of dimensions
*     AMAT    R*8   I     ND,ND      Matrix (from 6*6 array)
*     VAL     R*8     O    ND        Eigenvalues
*     VEC     R*8     O   ND,ND      Eigenvectors
*     
*     Globals:
*     VARIABLE  COMMON BLOCK  I/O  DESCRIPTION
*     -
*     
*     External calls:
*     -
*     
*     History:
*     John Skilling    8 Nov 1985     Initial release
*     
*     
*     Notes:
*     (1) Eigenvalues VAL are returned in increasing order
*     (2) Eigenvectors VEC are normalised to V.AMAT.V=1
*     (3) Input matrix AMAT is preserved
*     (4) Only the upper triangle j>=i of AMAT(i,j) is read
*     
*     (5) Algorithm is to repeatedly square AMAT until the largest
*     eigenvalue dominates, then to subtract off that eigenvector,
*     and repeat until all the eigenvectors have been found.
*     This ND**4 algorithm would be too slow for large matrices,
*     but is robust for small.
*     
*-----------------------------------------------------------------------
      IMPLICIT CHARACTER*1 (A-H,O-Z)
      REAL*8 AMAT(6,6),VAL(6),VEC(6,6)
      REAL*8 Z(6,6),P(6,6),Q(6,6),X(6),Y(6),T,A,B,C,EPS
*     EPS is related to machine accuracy
      PARAMETER (EPS=3.0E-16)
*     Copy AMAT to a positive definite matrix Z (by adding C to eigenvals)
      A=0.0
      B=0.0
      T=0.0
      DO 2 I=1,ND
         A=A+AMAT(I,I)
         DO 1  J=1,ND
            C=AMAT(I,J)
            IF(ABS(C).LT.1.0E-18) C=0.0
            B=B+C**2
 1       CONTINUE
         T=T+1.0
    2 CONTINUE
      C=MAX(B-A*A/T,0.0D0)
      C=A/T-SQRT(C)-EPS*SQRT(B)
      print *,c
      DO 4 I=1,ND
         DO 3 J=I,ND
            Z(I,J)=AMAT(I,J)
            Z(J,I)=AMAT(I,J)
 3       CONTINUE
         Z(I,I)=Z(I,I)-C
    4 CONTINUE
* LMAX = maximum number of inner loop iterates
      T=-LOG(EPS)/(EPS*LOG(2.0))
      T=LOG(T)/LOG(2.0)
      LMAX=IFIX(SNGL(T+2.0))
*
      N=ND
* Outer loop over N for successively smaller eigenvalues
    5 CONTINUE
        T=0.0
        DO 7 I=1,ND
          DO 6 J=1,ND
            P(I,J)=Z(I,J)
    6     CONTINUE
          T=T+P(I,I)
    7   CONTINUE
        
        print *,"t = ",t
*
        L=0
* Inner loop over L for squaring P and setting T=trace
    8   CONTINUE
          L=L+1
          T=1.0/T
          DO 10 I=1,ND
            DO  9 J=1,ND
              Q(I,J)=P(I,J)*T
              IF(ABS(Q(I,J)).LT.1.0E-18) Q(I,J)=0.0
    9       CONTINUE
   10     CONTINUE
          T=0.0
          DO 13 I=1,ND
            DO 12 J=I,ND
              A=0.0
              DO 11 K=1,ND
                A=A+Q(I,K)*Q(K,J)
   11         CONTINUE
              IF(ABS(A).LT.1.0E-18) A=0.0
              P(I,J)=A
              P(J,I)=A
   12       CONTINUE
            T=T+P(I,I)
   13     CONTINUE
          print *,"end t = ",t
        IF( T.LT.1.0-EPS*10.0 .AND. L.LE.LMAX )  GOTO 8
* End inner loop when P is dyadic
*
* K = largest column = estimate of current largest eigenvector
        K=1
        A=0.0
        DO 14 I=1,ND
          IF(P(I,I).GT.A) THEN
            A=P(I,I)
            K=I
          ENDIF
   14   CONTINUE
*
* P(P(largest column)) = better estimate X of eigenvector
        DO 16 I=1,ND
          A=0.0
          DO 15  J=1,ND
            A=A+P(I,J)*P(J,K)
   15     CONTINUE
          IF(ABS(A).LT.1.0E-18) A=0.0
          Y(I)=A
   16   CONTINUE
        T=0.0
        DO 18 I=1,ND
          A=0.0
          DO 17 J=1,ND
            A=A+P(I,J)*Y(J)
   17     CONTINUE
          IF(ABS(A).LT.1.0E-18) A=0.0
          T=T+A*A
          X(I)=A
   18   CONTINUE
* Repeat..  Orthogonalise X to previous eigenvectors
        K=ND
        IF(K.EQ.N) GOTO 22
   19     A=0.0
          DO 20 I=1,ND
            A=A+X(I)*VEC(I,K)
   20     CONTINUE
          DO 21 I=1,ND
            B=X(I)-A
            IF(ABS(B).LT.1.0E-18) B=0.0
            X(I)=B
   21     CONTINUE
          K=K-1
          IF(K.GT.N) GOTO 19
* ..while
* Normalise eigenvector X
   22   A=0.0
        DO 23 I=1,ND
          A=A+X(I)*X(I)
   23   CONTINUE
        A=1.0/SQRT(A)
* Copy eigenvector X into output array VEC
        DO 24 I=1,ND
          X(I)=A*X(I)
          VEC(I,N)=X(I)
   24   CONTINUE
* Set eigenvalue VAL directly from the eigenvector
        A=0.0
        DO 26 I=1,ND
          DO 25 J=1,ND
            B=X(I)*X(J)
            IF(ABS(B).LT.1.0E-18) B=0.0
            A=A+Z(I,J)*B
   25     CONTINUE
   26   CONTINUE
        VAL(N)=A+C
* Finish ?   (if full set of eigenvectors has been found)
        N=N-1
        IF(N.LE.0) RETURN
* Otherwise, remove dyadic  X.Xtranspose  from matrix Z
        DO 28 I=1,ND
          DO 27 J=I,ND
            B=X(I)*X(J)
            IF(ABS(B).LT.1.0E-18) B=0.0
            Z(I,J)=Z(I,J)-A*B
            IF(ABS(Z(I,J)).LT.1.0E-18) Z(I,J)=0.0
            Z(J,I)=Z(I,J)
            print *, "z ", z(i,j)
   27     CONTINUE
   28   CONTINUE
      GOTO 5
* End outer loop
*
      END



