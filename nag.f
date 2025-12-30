      SUBROUTINE F01LBF(N,M1,M2,A,IA,AL,IL,IN,IV,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     MARK 14 REVISED. IER-727 (DEC 1989).
C
C     PDK, DNAC, NPL, TEDDINGTON. JUNE 1976.
C     NPL DNAC LIBRARY SUBROUTINE BANDLU.
C     REVISED BY N.A.G. CENTRAL OFFICE, 1979.
C
C     LU DECOMPOSITION OF A BAND MATRIX A WITH M1 SUB AND M2 SUPER-
C     DIAGONALS. THE MATRIX A IS PRESENTED AS AN IW*N ARRAY
C     WITH EACH ROW TOP JUSTIFIED, WHERE IW=MIN(N,M1+M2+1).
C     L AND U ARE FOUND USING GAUSSIAN
C     ELIMINATION WITH PARTIAL PIVOTING. U IS OVERWRITTEN ON A, THE
C     DIAGONAL ELEMENTS BEING STORED AS THEIR RECIPROCALS AND L IS
C     STORED AS A SEPARATE N*M1 MATRIX AL. THE ARRAY AL MUST
C     HAVE DIMENSIONS AT LEAST N*1 IN THE CALLING PROGRAM
C     IF M1 IS ZERO. DETAILS OF THE PIVOTING ARE
C     STORED IN THE VECTOR IN WHICH IS SUCH THAT IN(I)=K IF ROWS I
C     AND K WERE INTERCHENGED AT THE ITH MAJOR STEP.
C     IA AND IL ARE THE ROW DIMENSIONS OF A AND AL IN THE
C     CALLING PROGRAM.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01LBF')
C     .. Scalar Arguments ..
      INTEGER           IA, IFAIL, IL, IV, M1, M2, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,*), AL(IL,*)
      INTEGER           IN(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  EPS, ONE, X, Y, ZERO
      INTEGER           I, IK, IR, ISAVE, IW, IWW, J, JR, K, M
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN
C     .. Data statements ..
      DATA              ZERO/0.0D0/, ONE/1.0D0/
C     .. Executable Statements ..
      ISAVE = IFAIL
      IFAIL = 1
      IV = 0
      EPS = X02AJF()
      IW = MIN(M1+M2+1,N)
      IF (N.LT.1 .OR. M1.LT.0 .OR. M1.GT.N-1 .OR. M2.LT.0 .OR. M2.GT.N-
     *    1 .OR. IA.LT.IW .OR. IL.LT.1 .OR. IL.LT.M1) GO TO 380
      IFAIL = 2
      M = M2 + 1
      K = IW - M2 - 1
      IF (K.LE.0) GO TO 60
      DO 40 I = 1, K
         M = M + 1
         DO 20 J = M, IW
            A(J,I) = ZERO
   20    CONTINUE
   40 CONTINUE
   60 M = N - IW + M1 + 2
      J = IW + 1
      IF (M.GT.N) GO TO 120
      DO 100 I = M, N
         J = J - 1
         DO 80 K = J, IW
            A(K,I) = ZERO
   80    CONTINUE
  100 CONTINUE
C
C     ZEROS INSERTED.
C
  120 DO 180 I = 1, N
         X = ZERO
         DO 140 J = 1, IW
            X = X + ABS(A(J,I))
  140    CONTINUE
         IF (X.GT.ZERO) GO TO 160
         IR = I
         GO TO 360
  160    AL(1,I) = ONE/X
  180 CONTINUE
C
C     ROW NORMS OF A CALCULATED AND THEIR RECIPROCALS
C     STORED IN FIRST COLUMN OF AL.
C
      IFAIL = 3
      DO 340 IR = 1, N
         X = ZERO
         M = MIN(N,IR+M1)
         IWW = MIN(IW,N-IR+1)
         DO 200 I = IR, M
            Y = ABS(A(1,I))*AL(1,I)
            IF (Y.LE.X) GO TO 200
            X = Y
            J = I
  200    CONTINUE
         IF (X.LT.EPS) GO TO 360
         IN(IR) = J
C
C        (IR)TH PIVOT ELEMENT SELECTED.
C
         IF (J.EQ.IR) GO TO 240
         DO 220 I = 1, IWW
            X = A(I,IR)
            A(I,IR) = A(I,J)
            A(I,J) = X
  220    CONTINUE
         AL(1,J) = AL(1,IR)
C
C        ROWS IR AND J INTERCHANGED.
C
  240    JR = IR + 1
         Y = ONE/A(1,IR)
         IF (JR.GT.M) GO TO 320
         DO 300 I = JR, M
            X = A(1,I)*Y
            IF (IWW.LT.2) GO TO 280
            DO 260 J = 2, IWW
               A(J-1,I) = A(J,I) - X*A(J,IR)
  260       CONTINUE
  280       IK = I - IR
            AL(IK,IR) = X
            A(IWW,I) = ZERO
  300    CONTINUE
  320    A(1,IR) = Y
  340 CONTINUE
C
C     ELIMINATION COMPLETED.
C
      IFAIL = 0
      RETURN
  360 A(1,IR) = ZERO
      IV = IR
  380 IFAIL = P01ABF(ISAVE,IFAIL,SRNAME,0,P01REC)
      RETURN
      END
      SUBROUTINE F01NAF(N,ML,MU,A,NRA,TOL,IN,SCALE,IFAIL)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12B REVISED. IER-543 (FEB 1987).
C
C     F01NAF FACTORIZES THE N BY N COMPLEX BAND MATRIX C, WITH ML
C     SUB-DIAGONALS AND MU SUPER-DIAGONALS, AS
C
C     C = P*L*U,
C
C     WHERE P IS A PERMUTATION MATRIX, L IS A UNIT LOWER TRIANGULAR
C     MATRIX WITH AT MOST ML NON-ZERO SUB-DIAGONAL ELEMENTS PER
C     COLUMN AND U IS AN UPPER TRIANGULAR BAND MATRIX WITH
C     ( ML + MU ) SUPER-DIAGONALS.
C
C     C MUST BE SUPPLIED, COLUMN BY COLUMN, IN THE FIRST ( ML + MU + 1 )
C     ROWS OF THE ( 2*ML + MU + 1 ) BY N ARRAY A.
C
C     FOR A DESCRIPTION OF THE PARAMETERS AND USE OF THIS ROUTINE SEE
C     THE NAG LIBRARY MANUAL.
C
C     -- WRITTEN ON 6-DECEMBER-1983.  S.J.HAMMARLING.
C
C     NAG FORTRAN 66 GENERAL PURPOSE ROUTINE.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01NAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TOL
      INTEGER           IFAIL, ML, MU, N, NRA
C     .. Array Arguments ..
      COMPLEX*16        A(NRA,N)
      DOUBLE PRECISION  SCALE(N)
      INTEGER           IN(N)
C     .. Arrays in Common ..
      DOUBLE PRECISION  WMACH(15)
C     .. Local Scalars ..
      COMPLEX*16        CTEMP, CZERO, PIV
      DOUBLE PRECISION  APIV, COLNM, ONE, TEMP, TL, ZERO
      INTEGER           I, I1, I2, IA, II, J, J1, J2, K, KPIV, LA, NM1,
     *                  WIDTH
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          F06FBF, F06HBF, ZSCAL, ZAXPY, X02ZAZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, DIMAG, DBLE
C     .. Common blocks ..
      COMMON            /AX02ZA/WMACH
C     .. Save statements ..
      SAVE              /AX02ZA/
C     .. Data statements ..
      DATA              ZERO/0.0D+0/, ONE/1.0D+0/
      DATA              CZERO/(0.0D+0,0.0D+0)/
C     .. Executable Statements ..
C
C     CHECK INPUT PARAMETERS AND INITIALIZE.
C
      WIDTH = ML + 1 + MU
      IF (MAX(ML,MU).LT.N .AND. MIN(ML,MU).GE.0 .AND. NRA.GE.(ML+WIDTH))
     *    GO TO 20
      IFAIL = P01ABF(IFAIL,1,SRNAME,0,P01REC)
      RETURN
   20 CONTINUE
C
      CALL X02ZAZ
C
      TL = MAX(TOL,WMACH(3))
      NM1 = N - 1
C
C     SHIFT A DOWN ML ROWS TO ALLOW FOR THE ADDITIONAL ML
C     SUPER-DIAGONALS THAT WILL BE CREATED BY THE PARTIAL PIVOTING.
C
      IF (ML.EQ.0) GO TO 80
      DO 60 J = 1, N
         I1 = MAX(MU+1-J,0) + 1
         I2 = MIN(MU+1-J+N,WIDTH)
         IA = I2 + ML
         I = I2
         DO 40 II = I1, I2
            A(IA,J) = A(I,J)
            I = I - 1
            IA = IA - 1
   40    CONTINUE
   60 CONTINUE
   80 CONTINUE
C
C     COMPUTE THE ROW SCALE FACTORS OF A AND INSERT ZEROS INTO THE ML
C     ADDITIONAL SUB-DIAGONALS.
C
      CALL F06FBF(N,ZERO,SCALE,1)
C
      DO 120 J = 1, N
         I1 = MAX(WIDTH-J,0) + 1
C
         IF (I1.LE.ML) CALL F06HBF(ML-I1+1,CZERO,A(I1,J),1)
C
         I1 = MAX(WIDTH-J,ML) + 1
         I2 = MIN(N-J,ML) + WIDTH
         IA = I1 + J - WIDTH
         DO 100 I = I1, I2
            SCALE(IA) = SCALE(IA) + (ABS(DBLE(A(I,J)))+ABS(DIMAG(A(I,J))
     *                  ))
            IA = IA + 1
  100    CONTINUE
  120 CONTINUE
      DO 140 I = 1, N
         IF (SCALE(I).GT.ZERO) SCALE(I) = ONE/SCALE(I)
  140 CONTINUE
C
C     START THE ELIMINATION. ( N - 1 ) MAJOR STEPS.
C
      IF (ML.EQ.0) GO TO 320
      J2 = 1
      DO 300 K = 1, NM1
         LA = MIN(N-K,ML)
C
C           SELECT THE PIVOT ROW.
C
         KPIV = K
         APIV = SCALE(K)*(ABS(DBLE(A(WIDTH,K)))+ABS(DIMAG(A(WIDTH,K))))
         I1 = WIDTH + 1
         I2 = WIDTH + LA
         IA = K
         DO 180 I = I1, I2
            IA = IA + 1
            TEMP = SCALE(IA)*(ABS(DBLE(A(I,K)))+ABS(DIMAG(A(I,K))))
            IF (TEMP.LE.APIV) GO TO 160
            APIV = TEMP
            KPIV = IA
  160       CONTINUE
  180    CONTINUE
         IN(K) = KPIV
C
C           IF NECESSARY INTERCHANGE THE SCALE FACTORS AND THE PIVOT
C           ELEMENT.
C
         IA = WIDTH + KPIV - K
         PIV = A(IA,K)
         IF (KPIV.EQ.K) GO TO 200
         TEMP = SCALE(K)
         SCALE(K) = SCALE(KPIV)
         SCALE(KPIV) = TEMP
         A(IA,K) = A(WIDTH,K)
         A(WIDTH,K) = PIV
  200    CONTINUE
C
C           FORM THE MULTIPLIERS.
C
         IF (DBLE(PIV).EQ.ZERO .AND. DIMAG(PIV).EQ.ZERO) GO TO 280
C
         CALL ZSCAL(LA,ONE/PIV,A(WIDTH+1,K),1)
C
C              ELIMINATE.
C
         J1 = K + 1
         J2 = MIN(N,MAX(KPIV+MU,J2))
C
C              J2 = K ONLY IF KPIV = K AND MU = 0.
C
         IF (J2.EQ.K) GO TO 260
         DO 240 J = J1, J2
            IA = WIDTH + KPIV - J
            I1 = WIDTH + K - J
            CTEMP = A(IA,J)
            IF (KPIV.EQ.K) GO TO 220
            A(IA,J) = A(I1,J)
            A(I1,J) = CTEMP
  220       CONTINUE
C
            CALL ZAXPY(LA,-CTEMP,A(WIDTH+1,K),1,A(I1+1,J),1)
C
  240    CONTINUE
  260    CONTINUE
  280    CONTINUE
  300 CONTINUE
      GO TO 360
  320 IF (N.EQ.1) GO TO 360
      DO 340 K = 1, NM1
         IN(K) = K
  340 CONTINUE
  360 CONTINUE
C
C     TEST FOR NEAR SINGULARITY.
C
      IN(N) = 0
      DO 440 J = 1, N
         COLNM = ZERO
         I1 = MAX(WIDTH-J,0) + 1
         IF (I1.EQ.WIDTH) GO TO 400
         I2 = WIDTH - 1
         IA = J + I1 - WIDTH
         DO 380 I = I1, I2
            COLNM = COLNM + SCALE(IA)*(ABS(DBLE(A(I,J)))
     *              +ABS(DIMAG(A(I,J))))
            IA = IA + 1
  380    CONTINUE
  400    CONTINUE
         TEMP = SCALE(J)*(ABS(DBLE(A(WIDTH,J)))+ABS(DIMAG(A(WIDTH,J))))
         IF (TEMP.GT.TL*(COLNM+TEMP)) GO TO 420
         IN(N) = J
         GO TO 460
  420    CONTINUE
  440 CONTINUE
  460 CONTINUE
C
      IFAIL = 0
      RETURN
C
C     END OF F01NAF.
C
      END
      SUBROUTINE F03AFF(N,EPS,A,IA,D1,ID,P,IFAIL)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 3 REVISED.
C     MARK 4.5 REVISED
C     MARK 11 REVISED. VECTORISATION (JAN 1984).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12 REVISED (LEVEL 2 BLAS) (MAR 1986)
C     MARK 12 REVISED. EXTENDED BLAS (JUNE 1986)
C
C     UNSYMDET
C     THE UNSYMMETRIC MATRIX, A, IS STORED IN THE N*N ARRAY A(I,J),
C     I=1,N, J=1,N. THE DECOMPOSITION A=LU, WHERE L IS A
C     LOWER TRIANGULAR MATRIX AND U IS A UNIT UPPER TRIANGULAR
C     MATRIX, IS PERFORMED AND OVERWRITTEN ON A, OMITTING THE UNIT
C     DIAGONAL OF U. A RECORD OF ANY INTERCHANGES MADE TO THE ROWS
C     OF A IS KEPT IN P(I), I=1,N, SUCH THAT THE I-TH ROW AND
C     THE P(I)-TH ROW WERE INTERCHANGED AT THE I-TH STEP. THE
C     DETERMINANT, D1 * 2.0**ID, OF A IS ALSO COMPUTED. THE
C     SUBROUTINE
C     WILL FAIL IF A, MODIFIED BY THE ROUNDING ERRORS, IS SINGULAR
C     OR ALMOST SINGULAR. SETS IFAIL = 0 IF SUCCESSFUL ELSE IFAIL =
C     1.
C     1ST DECEMBER 1971
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F03AFF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  D1, EPS
      INTEGER           IA, ID, IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), P(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  X, Y
      INTEGER           I, J, K, L
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          DGEMV, DTRSV
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
      DO 20 I = 1, N
         P(I) = 0.0D0
   20 CONTINUE
      DO 60 J = 1, N
         DO 40 I = 1, N
            P(I) = P(I) + A(I,J)**2
   40    CONTINUE
   60 CONTINUE
      DO 80 I = 1, N
         IF (P(I).LE.0.0D0) GO TO 240
         P(I) = 1.0D0/SQRT(P(I))
   80 CONTINUE
      D1 = 1.0D0
      ID = 0
      DO 220 K = 1, N
         L = K
         X = 0.0D0
         DO 100 I = K, N
            Y = ABS(A(I,K)*P(I))
            IF (Y.LE.X) GO TO 100
            X = Y
            L = I
  100    CONTINUE
         IF (L.EQ.K) GO TO 140
         D1 = -D1
         DO 120 J = 1, N
            Y = A(K,J)
            A(K,J) = A(L,J)
            A(L,J) = Y
  120    CONTINUE
         P(L) = P(K)
  140    P(K) = L
         D1 = D1*A(K,K)
         IF (X.LT.8.0D0*EPS) GO TO 240
  160    IF (ABS(D1).LT.1.0D0) GO TO 180
         D1 = D1*0.0625D0
         ID = ID + 4
         GO TO 160
  180    IF (ABS(D1).GE.0.0625D0) GO TO 200
         D1 = D1*16.0D0
         ID = ID - 4
         GO TO 180
  200    IF (K.LT.N) THEN
            CALL DTRSV('L','N','N',K,A,IA,A(1,K+1),1)
            CALL DGEMV('N',N-K,K,-1.0D0,A(K+1,1),IA,A(1,K+1),1,1.0D0,
     *                 A(K+1,K+1),1)
         END IF
  220 CONTINUE
      IFAIL = 0
      RETURN
  240 IFAIL = P01ABF(IFAIL,1,SRNAME,0,P01REC)
      RETURN
      END
      SUBROUTINE F03AHF(N,A,IA,DETR,DETI,IDETE,RINT,IFAIL)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 3 REVISED.
C     MARK 4.5 REVISED
C     MARK 7 REVISED IER-143 (DEC 1978)
C     MARK 8 REVISED. IER-236 (APR 1980).
C     MARK 11 REVISED. VECTORISATION (JAN 1984).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12 REVISED. EXTENDED BLAS (JUNE 1986)
C
C     COMPDET
C     THE COMPLEX UNSYMMETRIC MATRIX, A, IS STORED IN THE ARRAY
C     A(N,N). THE DECOMPOSITION A = LU, WHERE L IS A LOWER
C     TRIANGULAR MATRIX AND U IS A UNIT UPPER TRIANGULAR MATRIX,
C     IS PERFORMED AND OVERWRITTEN ON A,
C     OMITTING THE UNIT DIAGONAL OF U. A RECORD OF ANY
C     INTERCHANGES MADE TO THE ROWS OF A IS KEPT IN RINT(I), I=1,N,
C     SUCH THAT THE I-TH ROW AND THE RINT(I)-TH ROW WERE
C     INTERCHANGED AT THE I-TH STEP. THE DETERMINANT, (DETR + I *
C     DETI) * 2.0**IDETE, OF A IS ALSO COMPUTED. THE SUBROUTINE
C     WILL
C     FAIL IF A, MODIFIED BY THE ROUNDING ERRORS, IS SINGULAR.
C     1ST DECEMBER 1971
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F03AHF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DETI, DETR
      INTEGER           IA, IDETE, IFAIL, N
C     .. Array Arguments ..
      COMPLEX*16        A(IA,N)
      DOUBLE PRECISION  RINT(N)
C     .. Local Scalars ..
      COMPLEX*16        CZ
      DOUBLE PRECISION  W, X, Y, Z
      INTEGER           I, J, K, L
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          ZGEMV, ZTRSV
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DIMAG, DBLE
C     .. Executable Statements ..
      DO 20 I = 1, N
         RINT(I) = 0.0D0
   20 CONTINUE
      DO 60 J = 1, N
         DO 40 I = 1, N
            RINT(I) = RINT(I) + DBLE(A(I,J))**2 + DIMAG(A(I,J))**2
   40    CONTINUE
   60 CONTINUE
      DO 80 I = 1, N
         IF (RINT(I).LE.0.0D0) GO TO 240
   80 CONTINUE
      DETR = 1.0D0
      DETI = 0.0D0
      IDETE = 0
      DO 220 K = 1, N
         L = K
         Z = 0.0D0
         DO 100 I = K, N
            X = DBLE(A(I,K))
            Y = DIMAG(A(I,K))
            X = (X*X+Y*Y)/RINT(I)
            IF (X.LE.Z) GO TO 100
            Z = X
            L = I
  100    CONTINUE
         IF (L.EQ.K) GO TO 140
         DETR = -DETR
         DETI = -DETI
         DO 120 J = 1, N
            CZ = A(K,J)
            A(K,J) = A(L,J)
            A(L,J) = CZ
  120    CONTINUE
         RINT(L) = RINT(K)
  140    RINT(K) = L
         X = DBLE(A(K,K))
         Y = DIMAG(A(K,K))
         Z = X*X + Y*Y
         W = X*DETR - Y*DETI
         DETI = X*DETI + Y*DETR
         DETR = W
         IF (ABS(DETR).LE.ABS(DETI)) W = DETI
         IF (W.EQ.0.0D0) GO TO 240
  160    IF (ABS(W).LT.1.0D0) GO TO 180
         W = W*0.0625D0
         DETR = DETR*0.0625D0
         DETI = DETI*0.0625D0
         IDETE = IDETE + 4
         GO TO 160
  180    IF (ABS(W).GE.0.0625D0) GO TO 200
         W = W*16.0D0
         DETR = DETR*16.0D0
         DETI = DETI*16.0D0
         IDETE = IDETE - 4
         GO TO 180
  200    IF (K.LT.N) THEN
            CALL ZTRSV('L','N','N',K,A,IA,A(1,K+1),1)
            CALL ZGEMV('N',N-K,K,(-1.0D0,0.0D0),A(K+1,1),IA,A(1,K+1),1,
     *                 (1.0D0,0.0D0),A(K+1,K+1),1)
         END IF
  220 CONTINUE
      IFAIL = 0
      RETURN
  240 IDETE = 0
      IFAIL = P01ABF(IFAIL,1,SRNAME,0,P01REC)
      RETURN
      END
      SUBROUTINE F04ADF(A,IA,B,IB,N,M,C,IC,WKSPCE,IFAIL)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 3 REVISED.
C     MARK 4.5 REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     APPROXIMATE SOLUTION OF A SET OF COMPLEX LINEAR
C     EQUATIONS WITH MULTIPLE RIGHT HAND SIDES.
C     1ST AUGUST 1971
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F04ADF')
C     .. Scalar Arguments ..
      INTEGER           IA, IB, IC, IFAIL, M, N
C     .. Array Arguments ..
      COMPLEX*16        A(IA,N), B(IB,M), C(IC,M)
      DOUBLE PRECISION  WKSPCE(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  DETI, DETR
      INTEGER           I, ISAVE, J
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          F03AHF, F04AKF
C     .. Executable Statements ..
      ISAVE = IFAIL
      IFAIL = 1
      DO 40 I = 1, N
         DO 20 J = 1, M
            C(I,J) = B(I,J)
   20    CONTINUE
   40 CONTINUE
      CALL F03AHF(N,A,IA,DETR,DETI,I,WKSPCE,IFAIL)
      IF (IFAIL.EQ.0) GO TO 60
      IFAIL = P01ABF(ISAVE,IFAIL,SRNAME,0,P01REC)
      RETURN
   60 CALL F04AKF(N,M,A,IA,WKSPCE,C,IC)
      RETURN
      END
      SUBROUTINE F04AEF(A,IA,B,IB,N,M,C,IC,WKSPCE,AA,IAA,BB,IBB,IFAIL)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 3 REVISED.
C     MARK 4.5 REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     ACCURATE SOLUTION OF A SET OF REAL LINEAR EQUATIONS
C     WITH MULTIPLE RIGHT HAND SIDES.
C     1ST AUGUST 1971
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F04AEF')
C     .. Scalar Arguments ..
      INTEGER           IA, IAA, IB, IBB, IC, IFAIL, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), AA(IAA,N), B(IB,M), BB(IBB,M), C(IC,M),
     *                  WKSPCE(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  D1, EPS, XXXX
      INTEGER           I, ISAVE, J
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F03AFF, F04AHF
C     .. Executable Statements ..
      ISAVE = IFAIL
      IFAIL = 1
      EPS = X02AJF()
C     COPY A TO WORKING ARRAY AA
      DO 40 I = 1, N
         DO 20 J = 1, N
            AA(I,J) = A(I,J)
   20    CONTINUE
   40 CONTINUE
      CALL F03AFF(N,EPS,AA,IAA,D1,I,WKSPCE,IFAIL)
      IF (IFAIL.EQ.0) GO TO 60
      IFAIL = P01ABF(ISAVE,IFAIL,SRNAME,0,P01REC)
      RETURN
   60 IFAIL = 1
      CALL F04AHF(N,M,A,IA,AA,IAA,WKSPCE,B,IB,EPS,C,IC,BB,IBB,I,IFAIL)
      IF (IFAIL.EQ.0) GO TO 80
      IFAIL = P01ABF(ISAVE,2,SRNAME,0,P01REC)
   80 RETURN
      END
      SUBROUTINE F04AHF(N,IR,A,IA,AA,IAA,P,B,IB,EPS,X,IX,BB,IBB,L,IFAIL)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 3 REVISED.
C     MARK 4 REVISED.
C     MARK 4.5 REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     UNSYMACCSOLVE
C     SOLVES AX=B WHERE A IS AN N*N UNSYMMETRIC MATRIX AND B IS AN
C     N*IR MATRIX OF RIGHT HAND SIDES, USING THE SUBROUTINE F04AJF.
C     THE SUBROUTINE MUST BY PRECEDED BY F03AFF IN WHICH L AND U
C     ARE PRODUCED IN AA(I,J) AND THE INTERCHANGES IN P(I). THE
C     RESIDUALS BB=B-AX ARE CALCULATED AND AD=BB IS SOLVED, OVER-
C     WRITING D ON BB. THE REFINEMENT IS REPEATED, AS LONG AS THE
C     MAXIMUM CORRECTION AT ANY STAGE IS LESS THAN HALF THAT AT THE
C     PREVIOUS STAGE, UNTIL THE MAXIMUM CORRECTION IS LESS THAN 2
C     EPS TIMES THE MAXIMUM X. SETS IFAIL = 1 IF THE SOLUTION FAILS
C     TO IMPROVE, ELSE IFAIL = 0. L IS THE NUMBER OF ITERATIONS.
C     ADDITIONAL PRECISION INNERPRODUCTS ARE ABSOLUTELY NECESSARY.
C     1ST DECEMBER 1971
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F04AHF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPS
      INTEGER           IA, IAA, IB, IBB, IFAIL, IR, IX, L, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), AA(IAA,N), B(IB,IR), BB(IBB,IR), P(N),
     *                  X(IX,IR)
C     .. Local Scalars ..
      DOUBLE PRECISION  BBMAX, D0, D1, D11, D2, XMAX
      INTEGER           I, ID2, IFAIL1, ISAVE, J
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          F04AJF, X03AAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
      ISAVE = IFAIL
      IFAIL1 = 0
      DO 40 J = 1, IR
         DO 20 I = 1, N
            X(I,J) = 0.0D0
            BB(I,J) = B(I,J)
   20    CONTINUE
   40 CONTINUE
      L = 0
      D0 = 0.0D0
   60 CALL F04AJF(N,IR,AA,IAA,P,BB,IBB)
      L = L + 1
      ID2 = 0
      D1 = 0.0D0
      DO 100 J = 1, IR
         DO 80 I = 1, N
            X(I,J) = X(I,J) + BB(I,J)
   80    CONTINUE
  100 CONTINUE
      DO 140 J = 1, IR
         XMAX = 0.0D0
         BBMAX = 0.0D0
         DO 120 I = 1, N
            IF (ABS(X(I,J)).GT.XMAX) XMAX = ABS(X(I,J))
            IF (ABS(BB(I,J)).GT.BBMAX) BBMAX = ABS(BB(I,J))
            CALL X03AAF(A(I,1),N*IA-I+1,X(1,J),(IR-J+1)
     *                  *IX,N,IA,1,-B(I,J),0.0D0,D11,D2,.TRUE.,IFAIL1)
            BB(I,J) = -D11
  120    CONTINUE
         IF (BBMAX.GT.D1*XMAX) D1 = BBMAX/XMAX
         IF (BBMAX.GT.2.0D0*EPS*XMAX) ID2 = 1
  140 CONTINUE
      IF ((D1.GT.0.67D0*D0) .AND. (L.NE.1)) GO TO 160
      D0 = D1
      IF (ID2.EQ.1) GO TO 60
      IFAIL = 0
      RETURN
  160 IFAIL = P01ABF(ISAVE,1,SRNAME,0,P01REC)
      RETURN
      END
      SUBROUTINE F04AJF(N,IR,A,IA,P,B,IB)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 4 REVISED.
C     MARK 4.5 REVISED
C     MARK 11 REVISED. VECTORISATION (JAN 1984).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12 REVISED. EXTENDED BLAS (JUNE 1986)
C
C     UNSYMSOL
C     SOLVES AX=B, WHERE A IS AN UNSYMMETRIC MATRIX AND B IS AN
C     N*IR
C     MATRIX OF IR RIGHT-HAND SIDES. THE SUBROUTINE F04AJF MUST BY
C     PRECEDED BY F03AFF IN WHICH L AND U ARE PRODUCED IN A(I,J),
C     FROM A, AND THE RECORD OF THE INTERCHANGES IS PRODUCED IN
C     P(I). AX=B IS SOLVED IN THREE STEPS, INTERCHANGE THE
C     ELEMENTS OF B, LY=B AND UX=Y. THE MATRICES Y AND THEN X ARE
C     OVERWRITTEN ON B.
C     1ST AUGUST 1971
C
C     .. Scalar Arguments ..
      INTEGER           IA, IB, IR, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), B(IB,IR), P(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  X
      INTEGER           I, I1, K
C     .. External Subroutines ..
      EXTERNAL          DTRSV
C     .. Executable Statements ..
C     INTERCHANGING ELEMENTS OF B
      DO 40 I = 1, N
         I1 = P(I) + 0.5D0
         IF (I1.EQ.I) GO TO 40
         DO 20 K = 1, IR
            X = B(I,K)
            B(I,K) = B(I1,K)
            B(I1,K) = X
   20    CONTINUE
   40 CONTINUE
      DO 60 K = 1, IR
C        SOLUTION OF LY= B
         CALL DTRSV('L','N','N',N,A,IA,B(1,K),1)
C        SOLUTION OF UX= Y
         CALL DTRSV('U','N','U',N,A,IA,B(1,K),1)
   60 CONTINUE
      RETURN
      END
      SUBROUTINE F04AKF(N,IR,A,IA,RINT,B,IB)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 4.5 REVISED
C     MARK 11 REVISED. VECTORISATION (JAN 1984).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12 REVISED. EXTENDED BLAS (JUNE 1986)
C
C     COMPSOL
C     SOLVES AX=B, WHERE A IS A COMPLEX UNSYMMETRIC MATRIX AND B
C     IS A COMPLEX MATRIX OF IR RIGHT-HAND SIDES.
C     THE SUBROUTINE F04AKF MUST BE PROCEDED BY
C     F03AHF IN WHICH L AND U ARE PRODUCED IN THE ARRAY A(N,N),
C     FROM A, AND THE RECORD OF THE INTERCHANGES IS
C     PRODUCED IN RINT(N). AX=B IS SOLVED IN THREE STEPS,
C     INTERCHANGE ELEMENTS OF B, LY=B AND UX=Y. THE MATRICES Y AND
C     THEN X ARE OVERWRITTEN ON B.
C     1ST AUGUST 1971
C
C     .. Scalar Arguments ..
      INTEGER           IA, IB, IR, N
C     .. Array Arguments ..
      COMPLEX*16        A(IA,N), B(IB,IR)
      DOUBLE PRECISION  RINT(N)
C     .. Local Scalars ..
      COMPLEX*16        CX
      INTEGER           I, I1, K
C     .. External Subroutines ..
      EXTERNAL          ZTRSV
C     .. Executable Statements ..
C     INTERCHANGE ELEMENTS OF B
      DO 40 I = 1, N
         I1 = RINT(I) + 0.5D0
         IF (I1.EQ.I) GO TO 40
         DO 20 K = 1, IR
            CX = B(I,K)
            B(I,K) = B(I1,K)
            B(I1,K) = CX
   20    CONTINUE
   40 CONTINUE
      DO 60 K = 1, IR
C        SOLUTION OF LY= B
         CALL ZTRSV('L','N','N',N,A,IA,B(1,K),1)
C        SOLUTION OF UX= Y
         CALL ZTRSV('U','N','U',N,A,IA,B(1,K),1)
   60 CONTINUE
      RETURN
      END
      SUBROUTINE F04ATF(A,IA,B,N,C,AA,IAA,WKS1,WKS2,IFAIL)
C     MARK 4 RELEASE NAG COPYRIGHT 1974.
C     MARK 4.5 REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     ACCURATE SOLUTION OF A SET OF REAL LINEAR EQUATIONS
C     WITH ONE RIGHT SIDE.
C     1ST. APRIL 1973
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F04ATF')
C     .. Scalar Arguments ..
      INTEGER           IA, IAA, IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), AA(IAA,N), B(N), C(N), WKS1(N), WKS2(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  D1, XXXX
      INTEGER           I, IT, J
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F03AFF, F04AHF
C     .. Executable Statements ..
      DO 40 I = 1, N
         DO 20 J = 1, N
            AA(I,J) = A(I,J)
   20    CONTINUE
   40 CONTINUE
      IT = 1
      CALL F03AFF(N,X02AJF(),AA,IAA,D1,I,WKS1,IT)
      IF (IT.EQ.0) GO TO 60
      IFAIL = P01ABF(IFAIL,1,SRNAME,0,P01REC)
      RETURN
   60 IT = 1
      CALL F04AHF(N,1,A,IA,AA,IAA,WKS1,B,N,X02AJF(),C,N,WKS2,N,I,IT)
      IF (IT.EQ.0) GO TO 80
      IFAIL = P01ABF(IFAIL,2,SRNAME,0,P01REC)
      RETURN
   80 IFAIL = 0
      RETURN
      END
      SUBROUTINE F04LDF(N,M1,M2,IR,A,IA,AL,IL,IN,B,IB,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14 REVISED. IER-742 (DEC 1989).
C
C     PDK, DNAC, NPL, TEDDINGTON. JUNE 1976.
C     NPL DNAC LIBRARY SUBROUTINE BANSOL.
C     REVISED BY N.A.G. CENTRAL OFFICE, 1979.
C
C     SOLUTION OF THE BAND EQUATIONS AX=B WHERE A,AL AND IN
C     ARE AS SUPPLIED BY SUBROUTINE F01LBF AND B IS AN N*IR MATRIX.
C     ALTHOUGH B MUST BE A TWO-DIMENSIONAL MATRIX IT IS PERMISSIBLE
C     TO PUT IR=1. IA,IL AND IB ARE THE ROW DIMENSIONS OF A, AL AND
C     B IN THE CALLING PROGRAM. X IS OVERWRITTEN ON B.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F04LDF')
C     .. Scalar Arguments ..
      INTEGER           IA, IB, IFAIL, IL, IR, M1, M2, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,*), AL(IL,*), B(IB,IR)
      INTEGER           IN(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  X, Y
      INTEGER           I, II, IK, IW, J, JJ, K, KK, M
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
      IW = MIN(M1+M2+1,N)
      IF (N.LT.1 .OR. M1.LT.0 .OR. M1.GT.N-1 .OR. M2.LT.0 .OR. M2.GT.N-
     *    1 .OR. IA.LT.IW .OR. IL.LT.1 .OR. IL.LT.M1 .OR. IB.LT.N)
     *    GO TO 160
      M = M1
      DO 60 K = 1, N
         IK = K + 1
         M = MIN(N,M+1)
         J = IN(K)
         DO 40 JJ = 1, IR
            X = B(J,JJ)
            B(J,JJ) = B(K,JJ)
            B(K,JJ) = X
            IF (IK.GT.M) GO TO 40
            DO 20 I = IK, M
               II = I - K
               B(I,JJ) = B(I,JJ) - X*AL(II,K)
   20       CONTINUE
   40    CONTINUE
   60 CONTINUE
C
C     FORWARD SUBSTITUTION COMPLETED.
C
      DO 140 K = 1, N
         M = MIN(IW,K)
         I = N + 1 - K
         II = I - 1
         Y = A(1,I)
         DO 120 JJ = 1, IR
            X = B(I,JJ)
            IF (M.EQ.1) GO TO 100
            DO 80 J = 2, M
               KK = J + II
               X = X - A(J,I)*B(KK,JJ)
   80       CONTINUE
  100       B(I,JJ) = X*Y
  120    CONTINUE
  140 CONTINUE
C
C     BACKWARD SUBSTITUTION COMPLETED.
C
      IFAIL = 0
      RETURN
  160 IFAIL = P01ABF(IFAIL,1,SRNAME,0,P01REC)
      RETURN
      END
      SUBROUTINE F04NAF(JOB,N,ML,MU,A,NRA,IN,B,TOL,IFAIL)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     F04NAF SOLVES ONE OF THE SYSTEMS OF EQUATIONS
C
C     C*X = B,   ( C**H )*X = B,   U*X = B,
C
C     WHERE C IS AN N BY N COMPLEX BAND MATRIX, WITH ML SUB-DIAGONALS
C     AND MU SUPER-DIAGONALS, THAT HAS BEEN FACTORIZED AS
C
C     C = P*L*U,
C
C     BY ROUTINE F01NAF, WITH THE FACTORS RETURNED IN THE
C     ( 2*ML + MU + 1 ) BY N ARRAY A.
C
C     FOR A DESCRIPTION OF THE PARAMETERS AND USE OF THIS ROUTINE SEE
C     THE NAG LIBRARY MANUAL.
C
C     -- WRITTEN ON 8-DECEMBER-1983.  S.J.HAMMMARLING.
C
C     NAG FORTRAN 66 GENERAL PURPOSE ROUTINE.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F04NAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TOL
      INTEGER           IFAIL, JOB, ML, MU, N, NRA
C     .. Array Arguments ..
      COMPLEX*16        A(NRA,N), B(N)
      INTEGER           IN(N)
C     .. Arrays in Common ..
      DOUBLE PRECISION  WMACH(15)
C     .. Local Scalars ..
      COMPLEX*16        AWJ, CPERT, CTEMP
      DOUBLE PRECISION  PERTI, PERTR, TWO, ZERO
      INTEGER           I, I1, IA, IERR, J, JJ, K, KK, KPIV, NM1, WIDTH
      LOGICAL           FAIL
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      COMPLEX*16        ZDOTC, F06CLF
      INTEGER           P01ABF
      EXTERNAL          ZDOTC, F06CLF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          ZAXPY, X02ZAZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, DCMPLX, DCONJG, DIMAG, DBLE
C     .. Common blocks ..
      COMMON            /AX02ZA/WMACH
C     .. Save statements ..
      SAVE              /AX02ZA/
C     .. Data statements ..
      DATA              TWO/2.0D+0/, ZERO/0.0D+0/
C     .. Executable Statements ..
C
C     CHECK INPUT PARAMETERS AND INITIALIZE.
C
      WIDTH = ML + 1 + MU
      IF (MAX(ML,MU).LT.N .AND. MIN(ML,MU).GE.0 .AND. NRA.GE.(ML+WIDTH)
     *     .AND. ABS(JOB).LT.4 .AND. JOB.NE.0) GO TO 20
      IFAIL = P01ABF(IFAIL,1,SRNAME,0,P01REC)
      RETURN
   20 CONTINUE
C
      CALL X02ZAZ
C
      NM1 = N - 1
      IF (JOB.GT.0) GO TO 100
C
C        COMPUTE TOL.
C
      IF (TOL.GT.ZERO) GO TO 80
      DO 60 J = 1, N
         I1 = 1 + MAX(WIDTH-J,0)
         DO 40 I = I1, WIDTH
            TOL = MAX(TOL,ABS(DBLE(A(I,J)))+ABS(DIMAG(A(I,J))))
   40    CONTINUE
   60 CONTINUE
      TOL = TOL*WMACH(3)
      IF (TOL.EQ.ZERO) TOL = WMACH(3)
   80 CONTINUE
  100 CONTINUE
C
      IF (ABS(JOB).NE.1 .OR. ML.EQ.0) GO TO 160
      DO 140 K = 1, NM1
         KPIV = IN(K)
         CTEMP = B(KPIV)
         IF (KPIV.EQ.K) GO TO 120
         B(KPIV) = B(K)
         B(K) = CTEMP
  120    CONTINUE
C
         CALL ZAXPY(MIN(N-K,ML),-CTEMP,A(WIDTH+1,K),1,B(K+1),1)
C
  140 CONTINUE
  160 CONTINUE
C
      IF (JOB.NE.1 .AND. JOB.NE.3) GO TO 220
      J = N
      DO 200 JJ = 1, N
C
         B(J) = F06CLF(B(J),A(WIDTH,J),FAIL)
         IF (FAIL) GO TO 540
C
         I1 = 1 + MAX(J-WIDTH,0)
         IF (I1.GE.J) GO TO 180
         IA = WIDTH + I1 - J
C
         CALL ZAXPY(J-I1,-B(J),A(IA,J),1,B(I1),1)
C
  180    CONTINUE
         J = J - 1
  200 CONTINUE
      IFAIL = 0
      RETURN
  220 CONTINUE
C
      IF (JOB.NE.(-1) .AND. JOB.NE.(-3)) GO TO 320
      J = N
      DO 300 JJ = 1, N
         AWJ = A(WIDTH,J)
         PERTR = TOL
         IF (DBLE(AWJ).LT.ZERO) PERTR = -PERTR
         PERTI = TOL
         IF (DIMAG(AWJ).LT.ZERO) PERTI = -PERTI
         CPERT = DCMPLX(PERTR,PERTI)
         CTEMP = B(J)
C
         B(J) = F06CLF(CTEMP,AWJ,FAIL)
C
C        +          WHILE( FAIL )LOOP
  240    IF ( .NOT. FAIL) GO TO 260
         AWJ = AWJ + CPERT
         CPERT = TWO*CPERT
         B(J) = F06CLF(CTEMP,AWJ,FAIL)
         GO TO 240
C        +          END WHILE
  260    CONTINUE
         I1 = 1 + MAX(J-WIDTH,0)
         IF (I1.GE.J) GO TO 280
         IA = WIDTH + I1 - J
C
         CALL ZAXPY(J-I1,-B(J),A(IA,J),1,B(I1),1)
C
  280    CONTINUE
         J = J - 1
  300 CONTINUE
      IFAIL = 0
      RETURN
  320 CONTINUE
C
      IF (JOB.EQ.(-2)) GO TO 380
      DO 360 J = 1, N
         CTEMP = B(J)
         I1 = 1 + MAX(J-WIDTH,0)
         IF (I1.GE.J) GO TO 340
         IA = WIDTH + I1 - J
C
         CTEMP = CTEMP - ZDOTC(J-I1,A(IA,J),1,B(I1),1)
C
  340    CONTINUE
C
         B(J) = F06CLF(CTEMP,DCONJG(A(WIDTH,J)),FAIL)
         IF (FAIL) GO TO 540
C
  360 CONTINUE
      GO TO 480
  380 CONTINUE
C
      DO 460 J = 1, N
         AWJ = DCONJG(A(WIDTH,J))
         PERTR = TOL
         IF (DBLE(AWJ).LT.ZERO) PERTR = -PERTR
         PERTI = TOL
         IF (DIMAG(AWJ).LT.ZERO) PERTI = -PERTI
         CPERT = DCMPLX(PERTR,PERTI)
         CTEMP = B(J)
         I1 = 1 + MAX(J-WIDTH,0)
         IF (I1.GE.J) GO TO 400
         IA = WIDTH + I1 - J
C
         CTEMP = CTEMP - ZDOTC(J-I1,A(IA,J),1,B(I1),1)
C
  400    CONTINUE
C
         B(J) = F06CLF(CTEMP,AWJ,FAIL)
C
C        +          WHILE( FAIL )LOOP
  420    IF ( .NOT. FAIL) GO TO 440
         AWJ = AWJ + CPERT
         CPERT = TWO*CPERT
         B(J) = F06CLF(CTEMP,AWJ,FAIL)
         GO TO 420
C        +          END WHILE
  440    CONTINUE
  460 CONTINUE
  480 CONTINUE
C
      IF (ML.EQ.0) GO TO 520
      K = NM1
      DO 500 KK = 1, NM1
C
         CTEMP = B(K) - ZDOTC(MIN(N-K,ML),A(WIDTH+1,K),1,B(K+1),1)
C
         KPIV = IN(K)
         IF (KPIV.NE.K) B(K) = B(KPIV)
         B(KPIV) = CTEMP
         K = K - 1
  500 CONTINUE
  520 CONTINUE
      IFAIL = 0
      RETURN
C
  540 IERR = J + 1
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,P01REC)
      RETURN
C
C     END OF F04NAF.
C
      END









      SUBROUTINE F06AAZ ( SRNAME, INFO )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      INTEGER            INFO
      CHARACTER*13       SRNAME
C     ..
C
C  Purpose
C  =======
C
C  F06AAZ  is an error handler for the Level 2 BLAS routines.
C
C  It is called by the Level 2 BLAS routines if an input parameter is
C  invalid.
C
C  Parameters
C  ==========
C
C  SRNAME - CHARACTER*13.
C           On entry, SRNAME specifies the name of the routine which
C           called F06AAZ.
C
C  INFO   - INTEGER.
C           On entry, INFO specifies the position of the invalid
C           parameter in the parameter-list of the calling routine.
C
C
C  Auxiliary routine for Level 2 Blas.
C
C  Written on 20-July-1986.
C
C     .. Local Scalars ..
      INTEGER            IFAIL
      CHARACTER*80       REC (1)
C     .. External Functions ..
      INTEGER            P01ABF
      EXTERNAL           P01ABF
C     ..
C     .. Executable Statements ..
      WRITE (REC (1),99999) SRNAME, INFO
      IFAIL = 0
      IFAIL = P01ABF (IFAIL, -1, SRNAME(1:6), 1, REC)
C
      RETURN
C
99999 FORMAT ( ' ** On entry to ', A13, ' parameter number ', I2,
     $         ' had an illegal value' )
C
C     End of F06AAZ.
C
      END
      COMPLEX*16       FUNCTION F06CLF( A, B, FAIL )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-605 (MAR 1988).
C     .. Scalar Arguments ..
      COMPLEX*16                        A, B
      LOGICAL                           FAIL
C     ..
C
C  F06CLF returns the value div given by
C
C     div = ( a/b      if a/b does not overflow,
C           (
C           ( 0.0      if a .eq. 0.0,
C           (
C           ( cflmax   if a .ne. 0.0 and a/b would overflow,
C
C  where
C
C     cflmax = ( flmax*sign( re( a/b ) ), flmax*sign( im( a/b ) ) )
C
C  and flmax is a large value, via the function name. In addition if a/b
C  would  overflow then  fail  is returned as  true, otherwise  fail  is
C  returned as false.
C
C  Note that when a and b are both zero, fail is returned as .true., but
C  div  is returned as  0.0. In all other cases of overflow  div is such
C  that abs( re( div ) ) and abs( im( div ) ) are flmax.
C
C  For  real  x and y,  if  y = 0,  sign( x/y )  is taken as  sign( x ).
C
C  Nag Fortran 77 O( 1 ) basic linear algebra routine.
C
C  -- Written on 27-April-1983.
C     Sven Hammarling, Nag Central Office.
C  -- Amended on 4-December-1987.
C     Sven Hammarling and Mick Pont, Nag Central Office.
C        To avoid extremely unlikely division by zero.
C
C
C     .. Parameters ..
      DOUBLE PRECISION         ONE
      PARAMETER              ( ONE  = 1.0D+0 )
      COMPLEX*16               ZERO
      PARAMETER              ( ZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Local Scalars ..
      COMPLEX*16               VALUE
      DOUBLE PRECISION         AI, AR, BI, BIG, BR, DIV, FLMAX, FLMIN,
     $                         NUMI, NUMR, TEMP
      LOGICAL                  FIRST
C     .. External Functions ..
      DOUBLE PRECISION         X02AMF
      EXTERNAL                 X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC                ABS, DBLE, DCMPLX, DIMAG, MAX, SIGN
C     .. Save statement ..
      SAVE                     BIG, FIRST, FLMAX
C     .. Data statements ..
      DATA                     FIRST/ .TRUE. /
C     ..
C     .. Executable Statements ..
      IF( A.EQ.ZERO )THEN
         VALUE = ZERO
         IF( B.EQ.ZERO )THEN
            FAIL = .TRUE.
         ELSE
            FAIL = .FALSE.
         END IF
      ELSE
C
         IF( FIRST )THEN
            FIRST = .FALSE.
            FLMIN =  X02AMF( )
            FLMAX =  1/FLMIN
            BIG   =  FLMAX/2
         END IF
C
         AR    =  DBLE ( A )
         AI    =  DIMAG( A )
         BR    =  DBLE ( B )
         BI    =  DIMAG( B )
         TEMP  =  MAX( ABS( AR ), ABS( AI ), ABS( BR ), ABS( BI ) )
         IF( TEMP.GE.BIG )THEN
            AR = AR/2
            AI = AI/2
            BR = BR/2
            BI = BI/2
         END IF
         IF( DCMPLX( BR, BI ).EQ.ZERO ) THEN
            VALUE =  DCMPLX( SIGN( FLMAX, DBLE ( A ) ),
     $                       SIGN( FLMAX, DIMAG( A ) )  )
            FAIL  = .TRUE.
         ELSE
            IF( ABS( BR ).GE.ABS( BI ) )THEN
               TEMP = BI/BR
               DIV  = BR     + TEMP*BI
               NUMR = AR     + TEMP*AI
               NUMI = AI     - TEMP*AR
            ELSE
               TEMP = BR/BI
               DIV  = BI      + TEMP*BR
               NUMR = AI      + TEMP*AR
               NUMI = TEMP*AI - AR
            END IF
            IF( ABS( DIV ).GE.ONE )THEN
               VALUE =  DCMPLX( NUMR/DIV, NUMI/DIV )
               FAIL  = .FALSE.
            ELSE
               TEMP  =  ABS( DIV )*FLMAX
               IF( ( ABS( NUMR ).LE.TEMP ).AND.
     $             ( ABS( NUMI ).LE.TEMP )      )THEN
                  VALUE =  DCMPLX( NUMR/DIV, NUMI/DIV )
                  FAIL  = .FALSE.
               ELSE
                  IF( DIV.GE.DBLE( ZERO ) )THEN
                     VALUE = DCMPLX( SIGN( FLMAX,  NUMR ),
     $                               SIGN( FLMAX,  NUMI )  )
                  ELSE
                     VALUE = DCMPLX( SIGN( FLMAX, -NUMR ),
     $                               SIGN( FLMAX, -NUMI )  )
                  END IF
                  FAIL = .TRUE.
               END IF
            END IF
         END IF
      END IF
C
      F06CLF = VALUE
      RETURN
C
C     End of F06CLF. ( CDIV )
C
      END
      SUBROUTINE F06FBF( N, CONST, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION   CONST
      INTEGER            INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
C     ..
C
C  F06FBF performs the operation
C
C     x = const*e,   e' = ( 1  1 ... 1 ).
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 22-September-1983.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      INTEGER            IX
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( CONST.NE.ZERO )THEN
            DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = CONST
   10       CONTINUE
         ELSE
            DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = ZERO
   20       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06FBF. ( SLOAD )
C
      END
      COMPLEX*16       FUNCTION F06GBF( N, X, INCX, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      COMPLEX*16                ZDOTC
      ENTRY                     ZDOTC ( N, X, INCX, Y, INCY )
C     .. Scalar Arguments ..
      INTEGER                           INCX, INCY, N
C     .. Array Arguments ..
      COMPLEX*16                        X( * ), Y( * )
C     ..
C
C  F06GBF returns the value
C
C     F06GBF = conjg( x' )*y
C
C
C  Nag Fortran 77 version of the Blas routine ZDOTC.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 21-September-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      COMPLEX*16               ZERO
      PARAMETER              ( ZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Local Scalars ..
      COMPLEX*16               SUM
      INTEGER                  I, IX, IY
C     .. Intrinsic Functions ..
      INTRINSIC                DCONJG
C     ..
C     .. Executable Statements ..
      SUM = ZERO
      IF( N.GT.0 )THEN
         IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
            DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
               SUM = SUM + DCONJG( X( IX ) )*Y( IX )
   10       CONTINUE
         ELSE
            IF( INCY.GE.0 )THEN
               IY = 1
            ELSE
               IY = 1 - ( N - 1 )*INCY
            END IF
            IF( INCX.GT.0 )THEN
               DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
                  SUM = SUM + DCONJG( X( IX ) )*Y( IY )
                  IY  = IY  + INCY
   20          CONTINUE
            ELSE
               IX = 1 - ( N - 1 )*INCX
               DO 30, I = 1, N
                  SUM = SUM + DCONJG( X( IX ) )*Y( IY )
                  IX  = IX  + INCX
                  IY  = IY  + INCY
   30          CONTINUE
            END IF
         END IF
      END IF
C
      F06GBF = SUM
      RETURN
C
C     End of F06GBF. ( ZDOTC )
C
      END
      SUBROUTINE F06GCF( N, ALPHA, X, INCX, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      ENTRY      ZAXPY ( N, ALPHA, X, INCX, Y, INCY )
C     .. Scalar Arguments ..
      COMPLEX*16         ALPHA
      INTEGER            INCX, INCY, N
C     .. Array Arguments ..
      COMPLEX*16         X( * ), Y( * )
C     ..
C
C  F06GCF performs the operation
C
C     y := alpha*x + y
C
C
C  Nag Fortran 77 version of the Blas routine ZAXPY.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- written on 28-April-1983.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Local Scalars ..
      INTEGER            I, IX, IY
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ALPHA.NE.ZERO )THEN
            IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
               DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
                  Y( IX ) = ALPHA*X( IX ) + Y( IX )
   10          CONTINUE
            ELSE
               IF( INCY.GE.0 )THEN
                  IY = 1
               ELSE
                  IY = 1 - ( N - 1 )*INCY
               END IF
               IF( INCX.GT.0 )THEN
                  DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
                     Y( IY ) = ALPHA*X( IX ) + Y( IY )
                     IY      = IY            + INCY
   20             CONTINUE
               ELSE
                  IX = 1 - ( N - 1 )*INCX
                  DO 30, I = 1, N
                     Y( IY ) = ALPHA*X( IX ) + Y( IY )
                     IX      = IX            + INCX
                     IY      = IY            + INCY
   30             CONTINUE
               END IF
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06GCF. ( ZAXPY )
C
      END
      SUBROUTINE F06GDF( N, ALPHA, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      ENTRY      ZSCAL ( N, ALPHA, X, INCX )
C     .. Scalar Arguments ..
      COMPLEX*16         ALPHA
      INTEGER            INCX, N
C     .. Array Arguments ..
      COMPLEX*16         X( * )
C     ..
C
C  F06GDF performs the operation
C
C     x := alpha*x
C
C
C  Nag Fortran 77 version of the Blas routine ZSCAL.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 26-November-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      COMPLEX*16         CZERO
      PARAMETER        ( CZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Local Scalars ..
      INTEGER            IX
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ALPHA.EQ.CZERO )THEN
            DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = CZERO
   10       CONTINUE
         ELSE
            DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = ALPHA*X( IX )
   20       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06GDF. ( ZSCAL )
C
      END
      SUBROUTINE F06HBF( N, CONST, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      COMPLEX*16         CONST
      INTEGER            INCX, N
C     .. Array Arguments ..
      COMPLEX*16         X( * )
C     ..
C
C  F06HBF performs the operation
C
C     x = const*e,   e' = ( 1  1 ... 1 ).
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 22-September-1983.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Local Scalars ..
      INTEGER            IX
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( CONST.NE.ZERO )THEN
            DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = CONST
   10       CONTINUE
         ELSE
            DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = ZERO
   20       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06HBF. ( CLOAD )
C
      END
      SUBROUTINE F06PAF( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     AXP4 VERSION FOR VECTOR MACHINES
C     .. Entry Points ..
      ENTRY      DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
C     ..
C
C  Purpose
C  =======
C
C  DGEMV  performs one of the matrix-vector operations
C
C     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
C
C  where alpha and beta are scalars, x and y are vectors and A is an
C  m by n matrix.
C
C  Parameters
C  ==========
C
C  TRANS  - CHARACTER*1.
C           On entry, TRANS specifies the operation to be performed as
C           follows:
C
C              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
C
C              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
C
C              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C           On entry, M specifies the number of rows of the matrix A.
C           M must be at least zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the number of columns of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - DOUBLE PRECISION.
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry, the leading m by n part of the array A must
C           contain the matrix of coefficients.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, m ).
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of DIMENSION at least
C           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
C           and at least
C           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
C           Before entry, the incremented array X must contain the
C           vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  BETA   - DOUBLE PRECISION.
C           On entry, BETA specifies the scalar beta. When BETA is
C           supplied as zero then Y need not be set on input.
C           Unchanged on exit.
C
C  Y      - DOUBLE PRECISION array of DIMENSION at least
C           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
C           and at least
C           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
C           Before entry with BETA non-zero, the incremented array Y
C           must contain the vector y. On exit, Y is overwritten by the
C           updated vector y.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C  -- DO-loops unrolled on 20-November-1986.
C     Peter Mayes, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP, TEMP1, TEMP2, TEMP3, TEMP4
      INTEGER            I, INFO, IY, J, JX, KX, KY, LENX, LENY, M4, N4
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.(TRANS.EQ.'N' .OR. TRANS.EQ.'n').AND.
     $         .NOT.(TRANS.EQ.'T' .OR. TRANS.EQ.'t').AND.
     $         .NOT.(TRANS.EQ.'C' .OR. TRANS.EQ.'c')      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06PAF/DGEMV ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
C
C     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
C     up the start points in  X  and  Y.
C
      IF( (TRANS.EQ.'N' .OR. TRANS.EQ.'n') )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
C
C     Start the operations. In this version the inner loops are all
C     equivalent to AXPY operations.
C
C     First form  y := beta*y.
C
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      JX = KX
      IF( (TRANS.EQ.'N' .OR. TRANS.EQ.'n') )THEN
C
C        Form  y := alpha*A*x + y.
C
         IF( INCY.EQ.1 )THEN
C**** U n r o l l   t o   d e p t h   4 ********************************
            N4 = 4*( N/4 )
            DO 60, J = 1, N4, 4
               TEMP1 = ALPHA*X( JX )
               TEMP2 = ALPHA*X( JX + INCX )
               TEMP3 = ALPHA*X( JX + 2*INCX )
               TEMP4 = ALPHA*X( JX + 3*INCX )
               IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.ZERO.OR.
     $             TEMP4.NE.ZERO )THEN
                  DO 50, I = 1, M
                     Y( I ) = ( ( ( ( Y( I ) + TEMP1*A( I, J ) )
     $                        + TEMP2*A( I, J + 1 ) )
     $                        + TEMP3*A( I, J + 2 ) )
     $                        + TEMP4*A( I, J + 3 ) )
   50             CONTINUE
               END IF
               JX = JX + 4*INCX
   60       CONTINUE
C**** Clean-up loop ****************************************************
            DO 80, J = N4 + 1, N, 1
               TEMP = ALPHA*X( JX )
               IF( TEMP.NE.ZERO )THEN
                  DO 70, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         ELSE
C**** U n r o l l   t o   d e p t h   4 ********************************
            N4 = 4*( N/4 )
            DO 100, J = 1, N4, 4
               TEMP1 = ALPHA*X( JX )
               TEMP2 = ALPHA*X( JX + INCX )
               TEMP3 = ALPHA*X( JX + 2*INCX )
               TEMP4 = ALPHA*X( JX + 3*INCX )
               IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.ZERO.OR.
     $             TEMP4.NE.ZERO )THEN
                  IY = KY
                  DO 90, I = 1, M
                     Y( IY ) = ( ( ( ( Y( IY ) + TEMP1*A( I, J ) )
     $                         + TEMP2*A( I, J + 1 ) )
     $                         + TEMP3*A( I, J + 2 ) )
     $                         + TEMP4*A( I, J + 3 ) )
                     IY = IY + INCY
   90             CONTINUE
               END IF
               JX = JX + 4*INCX
  100       CONTINUE
C**** Clean-up loop ****************************************************
            DO 120, J = N4 + 1, N, 1
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY = KY
                  DO 110, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY = IY + INCY
  110             CONTINUE
               END IF
               JX = JX + INCX
  120       CONTINUE
         END IF
      ELSE
C
C        Form  y := alpha*A'*x + y.
C
         IF( INCY.EQ.1 )THEN
C**** U n r o l l   t o   d e p t h   4 ********************************
            M4 = 4*( M/4 )
            DO 140, J = 1, M4, 4
               TEMP1 = ALPHA*X( JX )
               TEMP2 = ALPHA*X( JX + INCX )
               TEMP3 = ALPHA*X( JX + 2*INCX )
               TEMP4 = ALPHA*X( JX + 3*INCX )
               IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.ZERO.OR.
     $             TEMP4.NE.ZERO )THEN
                  DO 130, I = 1, N
                     Y( I ) = ( ( ( ( Y( I ) + TEMP1*A( J, I ) )
     $                        + TEMP2*A( J + 1, I ) )
     $                        + TEMP3*A( J + 2, I ) )
     $                        + TEMP4*A( J + 3, I ) )
  130             CONTINUE
               END IF
               JX = JX + 4*INCX
  140       CONTINUE
C**** Clean-up loop ****************************************************
            DO 160, J = M4 + 1, M, 1
               TEMP = ALPHA*X( JX )
               IF( TEMP.NE.ZERO )THEN
                  DO 150, I = 1, N
                     Y( I ) = Y( I ) + TEMP*A( J, I )
  150             CONTINUE
               END IF
               JX = JX + INCX
  160       CONTINUE
         ELSE
C**** U n r o l l   t o   d e p t h   4 ********************************
            M4 = 4*( M/4 )
            DO 180, J = 1, M4, 4
               TEMP1 = ALPHA*X( JX )
               TEMP2 = ALPHA*X( JX + INCX )
               TEMP3 = ALPHA*X( JX + 2*INCX )
               TEMP4 = ALPHA*X( JX + 3*INCX )
               IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.ZERO.OR.
     $             TEMP4.NE.ZERO )THEN
                  IY = KY
                  DO 170, I = 1, N
                     Y( IY ) = ( ( ( ( Y( IY ) + TEMP1*A( J, I ) )
     $                         + TEMP2*A( J + 1, I ) )
     $                         + TEMP3*A( J + 2, I ) )
     $                         + TEMP4*A( J + 3, I ) )
                     IY = IY + INCY
  170             CONTINUE
               END IF
               JX = JX + 4*INCX
  180       CONTINUE
C**** Clean-up loop ****************************************************
            DO 200, J = M4 + 1, M, 1
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY = KY
                  DO 190, I = 1, N
                     Y( IY ) = Y( IY ) + TEMP*A( J, I )
                     IY = IY + INCY
  190             CONTINUE
               END IF
               JX = JX + INCX
  200       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06PAF (DGEMV ).
C
      END
      SUBROUTINE F06PJF( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     AXP4 VERSION FOR VECTOR MACHINES
C     .. Entry Points ..
      ENTRY      DTRSV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
C     .. Scalar Arguments ..
      INTEGER            INCX, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * )
C     ..
C
C  Purpose
C  =======
C
C  DTRSV  solves one of the systems of equations
C
C     A*x = b,   or   A'*x = b,
C
C  where b and x are n element vectors and A is an n by n unit, or
C  non-unit, upper or lower triangular matrix.
C
C  No test for singularity or near-singularity is included in this
C  routine. Such tests must be performed before calling this routine.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the matrix is an upper or
C           lower triangular matrix as follows:
C
C              UPLO = 'U' or 'u'   A is an upper triangular matrix.
C
C              UPLO = 'L' or 'l'   A is a lower triangular matrix.
C
C           Unchanged on exit.
C
C  TRANS  - CHARACTER*1.
C           On entry, TRANS specifies the equations to be solved as
C           follows:
C
C              TRANS = 'N' or 'n'   A*x = b.
C
C              TRANS = 'T' or 't'   A'*x = b.
C
C              TRANS = 'C' or 'c'   A'*x = b.
C
C           Unchanged on exit.
C
C  DIAG   - CHARACTER*1.
C           On entry, DIAG specifies whether or not A is unit
C           triangular as follows:
C
C              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
C
C              DIAG = 'N' or 'n'   A is not assumed to be unit
C                                  triangular.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry with  UPLO = 'U' or 'u', the leading n by n
C           upper triangular part of the array A must contain the upper
C           triangular matrix and the strictly lower triangular part of
C           A is not referenced.
C           Before entry with UPLO = 'L' or 'l', the leading n by n
C           lower triangular part of the array A must contain the lower
C           triangular matrix and the strictly upper triangular part of
C           A is not referenced.
C           Note that when  DIAG = 'U' or 'u', the diagonal elements of
C           A are not referenced either, but are assumed to be unity.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, n ).
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the n
C           element right-hand side vector b. On exit, X is overwritten
C           with the solution vector x.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C  -- DO-loops unrolled on 20-November-1986.
C     Peter Mayes, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP1, TEMP2, TEMP3, TEMP4
      INTEGER            I, INFO, IX, J, JX, KX, N4
      LOGICAL            NOUNIT
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          MAX, MOD
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.(UPLO .EQ.'U' .OR. UPLO .EQ.'u').AND.
     $         .NOT.(UPLO .EQ.'L' .OR. UPLO .EQ.'l')      )THEN
         INFO = 1
      ELSE IF( .NOT.(TRANS.EQ.'N' .OR. TRANS.EQ.'n').AND.
     $         .NOT.(TRANS.EQ.'T' .OR. TRANS.EQ.'t').AND.
     $         .NOT.(TRANS.EQ.'C' .OR. TRANS.EQ.'c')      )THEN
         INFO = 2
      ELSE IF( .NOT.(DIAG .EQ.'U' .OR. DIAG .EQ.'u').AND.
     $         .NOT.(DIAG .EQ.'N' .OR. DIAG .EQ.'n')      )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06PJF/DTRSV ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 )
     $   RETURN
C
      NOUNIT = (DIAG.EQ.'N' .OR. DIAG.EQ.'n')
C
C     Set up the start point in X if the increment is not unity. This
C     will be  ( N - 1 )*INCX  too small for descending loops.
C
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
C
C     Start the operations. In this version the inner loops are all
C     equivalent to AXPY operations.
C
      IF( (TRANS.EQ.'N' .OR. TRANS.EQ.'n') )THEN
C
C        Form  x := inv( A )*x.
C
         IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
            IF( INCX.EQ.1 )THEN
C**** U n r o l l   t o   D e p t h  4 ********** L O O P  A ***********
               N4 = MOD( N, 4 ) + 1
               DO 20, J = N, N4, -4
C**** 4 by 4 triangle at bottom of block *******************************
                  IF( NOUNIT )
     $               X( J ) = X( J )/A( J, J )
                  X( J - 1 ) = X( J - 1 ) - X( J )*A( J - 1, J )
                  IF( NOUNIT )
     $               X( J - 1 ) = X( J - 1 )/A( J - 1, J - 1 )
                  X( J - 2 ) = X( J - 2 ) - X( J )*A( J - 2, J ) -
     $                         X( J - 1 )*A( J - 2, J - 1 )
                  IF( NOUNIT )
     $               X( J - 2 ) = X( J - 2 )/A( J - 2, J - 2 )
                  X( J - 3 ) = X( J - 3 ) - X( J )*A( J - 3, J ) -
     $                         X( J - 1 )*A( J - 3, J - 1 ) - X( J - 2 )
     $                         *A( J - 3, J - 2 )
                  IF( NOUNIT )
     $               X( J - 3 ) = X( J - 3 )/A( J - 3, J - 3 )
C**** Unrolled main loop ***********************************************
                  TEMP1 = X( J )
                  TEMP2 = X( J - 1 )
                  TEMP3 = X( J - 2 )
                  TEMP4 = X( J - 3 )
                  IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.
     $                ZERO.OR.TEMP4.NE.ZERO )THEN
                     DO 10, I = J - 4, 1, -1
                        X( I ) = ( ( ( ( X( I ) - TEMP1*A( I, J ) )
     $                           - TEMP2*A( I, J - 1 ) )
     $                           - TEMP3*A( I, J - 2 ) )
     $                           - TEMP4*A( I, J - 3 ) )
   10                CONTINUE
                  END IF
   20          CONTINUE
C**** Left-overs on top left corner ************************************
               IF( N4.GE.2 )THEN
                  IF( NOUNIT )
     $               X( N4 - 1 ) = X( N4 - 1 )/A( N4 - 1, N4 - 1 )
               END IF
               IF( N4.GE.3 )THEN
                  X( N4 - 2 ) = X( N4 - 2 ) - X( N4 - 1 )
     $                          *A( N4 - 2, N4 - 1 )
                  IF( NOUNIT )
     $               X( N4 - 2 ) = X( N4 - 2 )/A( N4 - 2, N4 - 2 )
               END IF
               IF( N4.GE.4 )THEN
                  X( N4 - 3 ) = X( N4 - 3 ) - X( N4 - 1 )
     $                          *A( N4 - 3, N4 - 1 ) - X( N4 - 2 )
     $                          *A( N4 - 3, N4 - 2 )
                  IF( NOUNIT )
     $               X( N4 - 3 ) = X( N4 - 3 )/A( N4 - 3, N4 - 3 )
               END IF
            ELSE
               JX = KX + ( N - 1 )*INCX
C**** U n r o l l   t o   D e p t h  4 ********** L O O P  B ***********
               N4 = MOD( N, 4 ) + 1
               DO 40, J = N, N4, -4
C**** 4 by 4 triangle at bottom of block *******************************
                  IF( NOUNIT )
     $               X( JX ) = X( JX )/A( J, J )
                  X( JX - INCX ) = X( JX - INCX ) - X( JX )
     $                             *A( J - 1, J )
                  IF( NOUNIT )
     $               X( JX - INCX ) = X( JX - INCX )/A( J - 1, J - 1 )
                  X( JX - 2*INCX ) = X( JX - 2*INCX ) - X( JX )
     $                               *A( J - 2, J ) - X( JX - INCX )
     $                               *A( J - 2, J - 1 )
                  IF( NOUNIT )
     $               X( JX - 2*INCX ) = X( JX - 2*INCX )
     $                                  /A( J - 2, J - 2 )
                  X( JX - 3*INCX ) = X( JX - 3*INCX ) - X( JX )
     $                               *A( J - 3, J ) - X( JX - INCX )
     $                               *A( J - 3, J - 1 ) -
     $                               X( JX - 2*INCX )*A( J - 3, J - 2 )
                  IF( NOUNIT )
     $               X( JX - 3*INCX ) = X( JX - 3*INCX )
     $                                  /A( J - 3, J - 3 )
C**** Unrolled main loop ***********************************************
                  TEMP1 = X( JX )
                  TEMP2 = X( JX - INCX )
                  TEMP3 = X( JX - 2*INCX )
                  TEMP4 = X( JX - 3*INCX )
                  IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.
     $                ZERO.OR.TEMP4.NE.ZERO )THEN
                     IX = JX - 3*INCX
                     DO 30, I = J - 4, 1, -1
                        IX = IX - INCX
                        X( IX ) = ( ( ( ( X( IX ) - TEMP1*A( I, J ) )
     $                            - TEMP2*A( I, J - 1 ) )
     $                            - TEMP3*A( I, J - 2 ) )
     $                            - TEMP4*A( I, J - 3 ) )
   30                CONTINUE
                  END IF
                  JX = JX - 4*INCX
   40          CONTINUE
C**** Left-overs on top left corner ************************************
               IF( N4.GE.2 )THEN
                  IF( NOUNIT )
     $               X( JX ) = X( JX )/A( N4 - 1, N4 - 1 )
               END IF
               IF( N4.GE.3 )THEN
                  X( JX - INCX ) = X( JX - INCX ) - X( JX )
     $                             *A( N4 - 2, N4 - 1 )
                  IF( NOUNIT )
     $               X( JX - INCX ) = X( JX - INCX )/A( N4 - 2, N4 - 2 )
               END IF
               IF( N4.GE.4 )THEN
                  X( JX - 2*INCX ) = X( JX - 2*INCX ) - X( JX )
     $                               *A( N4 - 3, N4 - 1 ) -
     $                               X( JX - INCX )*A( N4 - 3, N4 - 2 )
                  IF( NOUNIT )
     $               X( JX - 2*INCX ) = X( JX - 2*INCX )
     $                                  /A( N4 - 3, N4 - 3 )
               END IF
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
C**** U n r o l l   t o   D e p t h  4 ********** L O O P  C ***********
               N4 = 4*( N/4 )
               DO 60, J = 1, N4, 4
C**** 4 by 4 triangle at bottom of block *******************************
                  IF( NOUNIT )
     $               X( J ) = X( J )/A( J, J )
                  X( J + 1 ) = X( J + 1 ) - X( J )*A( J + 1, J )
                  IF( NOUNIT )
     $               X( J + 1 ) = X( J + 1 )/A( J + 1, J + 1 )
                  X( J + 2 ) = X( J + 2 ) - X( J )*A( J + 2, J ) -
     $                         X( J + 1 )*A( J + 2, J + 1 )
                  IF( NOUNIT )
     $               X( J + 2 ) = X( J + 2 )/A( J + 2, J + 2 )
                  X( J + 3 ) = X( J + 3 ) - X( J )*A( J + 3, J ) -
     $                         X( J + 1 )*A( J + 3, J + 1 ) - X( J + 2 )
     $                         *A( J + 3, J + 2 )
                  IF( NOUNIT )
     $               X( J + 3 ) = X( J + 3 )/A( J + 3, J + 3 )
C**** Unrolled main loop ***********************************************
                  TEMP1 = X( J )
                  TEMP2 = X( J + 1 )
                  TEMP3 = X( J + 2 )
                  TEMP4 = X( J + 3 )
                  IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.
     $                ZERO.OR.TEMP4.NE.ZERO )THEN
                     DO 50, I = J + 4, N
                        X( I ) = ( ( ( ( X( I ) - TEMP1*A( I, J ) )
     $                           - TEMP2*A( I, J + 1 ) )
     $                           - TEMP3*A( I, J + 2 ) )
     $                           - TEMP4*A( I, J + 3 ) )
   50                CONTINUE
                  END IF
   60          CONTINUE
C**** Left-overs on top left corner ************************************
               IF( N4 + 1.LE.N )THEN
                  IF( NOUNIT )
     $               X( N4 + 1 ) = X( N4 + 1 )/A( N4 + 1, N4 + 1 )
               END IF
               IF( N4 + 2.LE.N )THEN
                  X( N4 + 2 ) = X( N4 + 2 ) - X( N4 + 1 )
     $                          *A( N4 + 2, N4 + 1 )
                  IF( NOUNIT )
     $               X( N4 + 2 ) = X( N4 + 2 )/A( N4 + 2, N4 + 2 )
               END IF
               IF( N4 + 3.LE.N )THEN
                  X( N4 + 3 ) = X( N4 + 3 ) - X( N4 + 1 )
     $                          *A( N4 + 3, N4 + 1 ) - X( N4 + 2 )
     $                          *A( N4 + 3, N4 + 2 )
                  IF( NOUNIT )
     $               X( N4 + 3 ) = X( N4 + 3 )/A( N4 + 3, N4 + 3 )
               END IF
            ELSE
               JX = KX
C**** U n r o l l   t o   D e p t h  4 ********** L O O P  D ***********
               N4 = 4*( N/4 )
               DO 80, J = 1, N4, 4
C**** 4 by 4 triangle at bottom of block *******************************
                  IF( NOUNIT )
     $               X( JX ) = X( JX )/A( J, J )
                  X( JX + INCX ) = X( JX + INCX ) - X( JX )
     $                             *A( J + 1, J )
                  IF( NOUNIT )
     $               X( JX + INCX ) = X( JX + INCX )/A( J + 1, J + 1 )
                  X( JX + 2*INCX ) = X( JX + 2*INCX ) - X( JX )
     $                               *A( J + 2, J ) - X( JX + INCX )
     $                               *A( J + 2, J + 1 )
                  IF( NOUNIT )
     $               X( JX + 2*INCX ) = X( JX + 2*INCX )
     $                                  /A( J + 2, J + 2 )
                  X( JX + 3*INCX ) = X( JX + 3*INCX ) - X( JX )
     $                               *A( J + 3, J ) - X( JX + INCX )
     $                               *A( J + 3, J + 1 ) -
     $                               X( JX + 2*INCX )*A( J + 3, J + 2 )
                  IF( NOUNIT )
     $               X( JX + 3*INCX ) = X( JX + 3*INCX )
     $                                  /A( J + 3, J + 3 )
C**** Unrolled main loop ***********************************************
                  TEMP1 = X( JX )
                  TEMP2 = X( JX + INCX )
                  TEMP3 = X( JX + 2*INCX )
                  TEMP4 = X( JX + 3*INCX )
                  IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.
     $                ZERO.OR.TEMP4.NE.ZERO )THEN
                     IX = JX + 3*INCX
                     DO 70, I = J + 4, N
                        IX = IX + INCX
                        X( IX ) = ( ( ( ( X( IX ) - TEMP1*A( I, J ) )
     $                            - TEMP2*A( I, J + 1 ) )
     $                            - TEMP3*A( I, J + 2 ) )
     $                            - TEMP4*A( I, J + 3 ) )
   70                CONTINUE
                  END IF
                  JX = JX + 4*INCX
   80          CONTINUE
C**** Left-overs on top left corner ************************************
               IF( N4 + 1.LE.N )THEN
                  IF( NOUNIT )
     $               X( JX ) = X( JX )/A( N4 + 1, N4 + 1 )
               END IF
               IF( N4 + 2.LE.N )THEN
                  X( JX + INCX ) = X( JX + INCX ) - X( JX )
     $                             *A( N4 + 2, N4 + 1 )
                  IF( NOUNIT )
     $               X( JX + INCX ) = X( JX + INCX )/A( N4 + 2, N4 + 2 )
               END IF
               IF( N4 + 3.LE.N )THEN
                  X( JX + 2*INCX ) = X( JX + 2*INCX ) - X( JX )
     $                               *A( N4 + 3, N4 + 1 ) -
     $                               X( JX + INCX )*A( N4 + 3, N4 + 2 )
                  IF( NOUNIT )
     $               X( JX + 2*INCX ) = X( JX + 2*INCX )
     $                                  /A( N4 + 3, N4 + 3 )
               END IF
            END IF
         END IF
      ELSE
C
C        Form  x := inv( A' )*x.
C
         IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
            IF( INCX.EQ.1 )THEN
C**** U n r o l l   t o   D e p t h  4 ********** L O O P  E ***********
               N4 = 4*( N/4 )
               DO 100, J = 1, N4, 4
C**** 4 by 4 triangle at bottom of block *******************************
                  IF( NOUNIT )
     $               X( J ) = X( J )/A( J, J )
                  X( J + 1 ) = X( J + 1 ) - X( J )*A( J, J + 1 )
                  IF( NOUNIT )
     $               X( J + 1 ) = X( J + 1 )/A( J + 1, J + 1 )
                  X( J + 2 ) = X( J + 2 ) - X( J )*A( J, J + 2 ) -
     $                         X( J + 1 )*A( J + 1, J + 2 )
                  IF( NOUNIT )
     $               X( J + 2 ) = X( J + 2 )/A( J + 2, J + 2 )
                  X( J + 3 ) = X( J + 3 ) - X( J )*A( J, J + 3 ) -
     $                         X( J + 1 )*A( J + 1, J + 3 ) - X( J + 2 )
     $                         *A( J + 2, J + 3 )
                  IF( NOUNIT )
     $               X( J + 3 ) = X( J + 3 )/A( J + 3, J + 3 )
C**** Unrolled main loop ***********************************************
                  TEMP1 = X( J )
                  TEMP2 = X( J + 1 )
                  TEMP3 = X( J + 2 )
                  TEMP4 = X( J + 3 )
                  IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.
     $                ZERO.OR.TEMP4.NE.ZERO )THEN
                     DO 90, I = J + 4, N
                        X( I ) = ( ( ( ( X( I ) - TEMP1*A( J, I ) )
     $                           - TEMP2*A( J + 1, I ) )
     $                           - TEMP3*A( J + 2, I ) )
     $                           - TEMP4*A( J + 3, I ) )
   90                CONTINUE
                  END IF
  100          CONTINUE
C**** Left-overs on top left corner ************************************
               IF( N4 + 1.LE.N )THEN
                  IF( NOUNIT )
     $               X( N4 + 1 ) = X( N4 + 1 )/A( N4 + 1, N4 + 1 )
               END IF
               IF( N4 + 2.LE.N )THEN
                  X( N4 + 2 ) = X( N4 + 2 ) - X( N4 + 1 )
     $                          *A( N4 + 1, N4 + 2 )
                  IF( NOUNIT )
     $               X( N4 + 2 ) = X( N4 + 2 )/A( N4 + 2, N4 + 2 )
               END IF
               IF( N4 + 3.LE.N )THEN
                  X( N4 + 3 ) = X( N4 + 3 ) - X( N4 + 1 )
     $                          *A( N4 + 1, N4 + 3 ) - X( N4 + 2 )
     $                          *A( N4 + 2, N4 + 3 )
                  IF( NOUNIT )
     $               X( N4 + 3 ) = X( N4 + 3 )/A( N4 + 3, N4 + 3 )
               END IF
            ELSE
               JX = KX
C**** U n r o l l   t o   D e p t h  4 ********** L O O P  D ***********
               N4 = 4*( N/4 )
               DO 120, J = 1, N4, 4
C**** 4 by 4 triangle at bottom of block *******************************
                  IF( NOUNIT )
     $               X( JX ) = X( JX )/A( J, J )
                  X( JX + INCX ) = X( JX + INCX ) - X( JX )
     $                             *A( J, J + 1 )
                  IF( NOUNIT )
     $               X( JX + INCX ) = X( JX + INCX )/A( J + 1, J + 1 )
                  X( JX + 2*INCX ) = X( JX + 2*INCX ) - X( JX )
     $                               *A( J, J + 2 ) - X( JX + INCX )
     $                               *A( J + 1, J + 2 )
                  IF( NOUNIT )
     $               X( JX + 2*INCX ) = X( JX + 2*INCX )
     $                                  /A( J + 2, J + 2 )
                  X( JX + 3*INCX ) = X( JX + 3*INCX ) - X( JX )
     $                               *A( J, J + 3 ) - X( JX + INCX )
     $                               *A( J + 1, J + 3 ) -
     $                               X( JX + 2*INCX )*A( J + 2, J + 3 )
                  IF( NOUNIT )
     $               X( JX + 3*INCX ) = X( JX + 3*INCX )
     $                                  /A( J + 3, J + 3 )
C**** Unrolled main loop ***********************************************
                  TEMP1 = X( JX )
                  TEMP2 = X( JX + INCX )
                  TEMP3 = X( JX + 2*INCX )
                  TEMP4 = X( JX + 3*INCX )
                  IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.
     $                ZERO.OR.TEMP4.NE.ZERO )THEN
                     IX = JX + 3*INCX
                     DO 110, I = J + 4, N
                        IX = IX + INCX
                        X( IX ) = ( ( ( ( X( IX ) - TEMP1*A( J, I ) )
     $                            - TEMP2*A( J + 1, I ) )
     $                            - TEMP3*A( J + 2, I ) )
     $                            - TEMP4*A( J + 3, I ) )
  110                CONTINUE
                  END IF
                  JX = JX + 4*INCX
  120          CONTINUE
C**** Left-overs on top left corner ************************************
               IF( N4 + 1.LE.N )THEN
                  IF( NOUNIT )
     $               X( JX ) = X( JX )/A( N4 + 1, N4 + 1 )
               END IF
               IF( N4 + 2.LE.N )THEN
                  X( JX + INCX ) = X( JX + INCX ) - X( JX )
     $                             *A( N4 + 1, N4 + 2 )
                  IF( NOUNIT )
     $               X( JX + INCX ) = X( JX + INCX )/A( N4 + 2, N4 + 2 )
               END IF
               IF( N4 + 3.LE.N )THEN
                  X( JX + 2*INCX ) = X( JX + 2*INCX ) - X( JX )
     $                               *A( N4 + 1, N4 + 3 ) -
     $                               X( JX + INCX )*A( N4 + 2, N4 + 3 )
                  IF( NOUNIT )
     $               X( JX + 2*INCX ) = X( JX + 2*INCX )
     $                                  /A( N4 + 3, N4 + 3 )
               END IF
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
C**** U n r o l l   t o   D e p t h  4 ********** L O O P  G ***********
               N4 = MOD( N, 4 ) + 1
               DO 140, J = N, N4, -4
C**** 4 by 4 triangle at bottom of block *******************************
                  IF( NOUNIT )
     $               X( J ) = X( J )/A( J, J )
                  X( J - 1 ) = X( J - 1 ) - X( J )*A( J, J - 1 )
                  IF( NOUNIT )
     $               X( J - 1 ) = X( J - 1 )/A( J - 1, J - 1 )
                  X( J - 2 ) = X( J - 2 ) - X( J )*A( J, J - 2 ) -
     $                         X( J - 1 )*A( J - 1, J - 2 )
                  IF( NOUNIT )
     $               X( J - 2 ) = X( J - 2 )/A( J - 2, J - 2 )
                  X( J - 3 ) = X( J - 3 ) - X( J )*A( J, J - 3 ) -
     $                         X( J - 1 )*A( J - 1, J - 3 ) - X( J - 2 )
     $                         *A( J - 2, J - 3 )
                  IF( NOUNIT )
     $               X( J - 3 ) = X( J - 3 )/A( J - 3, J - 3 )
C**** Unrolled main loop ***********************************************
                  TEMP1 = X( J )
                  TEMP2 = X( J - 1 )
                  TEMP3 = X( J - 2 )
                  TEMP4 = X( J - 3 )
                  IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.
     $                ZERO.OR.TEMP4.NE.ZERO )THEN
                     DO 130, I = J - 4, 1, -1
                        X( I ) = ( ( ( ( X( I ) - TEMP1*A( J, I ) )
     $                           - TEMP2*A( J - 1, I ) )
     $                           - TEMP3*A( J - 2, I ) )
     $                           - TEMP4*A( J - 3, I ) )
  130                CONTINUE
                  END IF
  140          CONTINUE
C**** Left-overs on top left corner ************************************
               IF( N4.GE.2 )THEN
                  IF( NOUNIT )
     $               X( N4 - 1 ) = X( N4 - 1 )/A( N4 - 1, N4 - 1 )
               END IF
               IF( N4.GE.3 )THEN
                  X( N4 - 2 ) = X( N4 - 2 ) - X( N4 - 1 )
     $                          *A( N4 - 1, N4 - 2 )
                  IF( NOUNIT )
     $               X( N4 - 2 ) = X( N4 - 2 )/A( N4 - 2, N4 - 2 )
               END IF
               IF( N4.GE.4 )THEN
                  X( N4 - 3 ) = X( N4 - 3 ) - X( N4 - 1 )
     $                          *A( N4 - 1, N4 - 3 ) - X( N4 - 2 )
     $                          *A( N4 - 2, N4 - 3 )
                  IF( NOUNIT )
     $               X( N4 - 3 ) = X( N4 - 3 )/A( N4 - 3, N4 - 3 )
               END IF
            ELSE
               JX = KX + ( N - 1 )*INCX
C**** U n r o l l   t o   D e p t h  4 ********** L O O P  B ***********
               N4 = MOD( N, 4 ) + 1
               DO 160, J = N, N4, -4
C**** 4 by 4 triangle at bottom of block *******************************
                  IF( NOUNIT )
     $               X( JX ) = X( JX )/A( J, J )
                  X( JX - INCX ) = X( JX - INCX ) - X( JX )
     $                             *A( J, J - 1 )
                  IF( NOUNIT )
     $               X( JX - INCX ) = X( JX - INCX )/A( J - 1, J - 1 )
                  X( JX - 2*INCX ) = X( JX - 2*INCX ) - X( JX )
     $                               *A( J, J - 2 ) - X( JX - INCX )
     $                               *A( J - 1, J - 2 )
                  IF( NOUNIT )
     $               X( JX - 2*INCX ) = X( JX - 2*INCX )
     $                                  /A( J - 2, J - 2 )
                  X( JX - 3*INCX ) = X( JX - 3*INCX ) - X( JX )
     $                               *A( J, J - 3 ) - X( JX - INCX )
     $                               *A( J - 1, J - 3 ) -
     $                               X( JX - 2*INCX )*A( J - 2, J - 3 )
                  IF( NOUNIT )
     $               X( JX - 3*INCX ) = X( JX - 3*INCX )
     $                                  /A( J - 3, J - 3 )
C**** Unrolled main loop ***********************************************
                  TEMP1 = X( JX )
                  TEMP2 = X( JX - INCX )
                  TEMP3 = X( JX - 2*INCX )
                  TEMP4 = X( JX - 3*INCX )
                  IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.
     $                ZERO.OR.TEMP4.NE.ZERO )THEN
                     IX = JX - 3*INCX
                     DO 150, I = J - 4, 1, -1
                        IX = IX - INCX
                        X( IX ) = ( ( ( ( X( IX ) - TEMP1*A( J, I ) )
     $                            - TEMP2*A( J - 1, I ) )
     $                            - TEMP3*A( J - 2, I ) )
     $                            - TEMP4*A( J - 3, I ) )
  150                CONTINUE
                  END IF
                  JX = JX - 4*INCX
  160          CONTINUE
C**** Left-overs on top left corner ************************************
               IF( N4.GE.2 )THEN
                  IF( NOUNIT )
     $               X( JX ) = X( JX )/A( N4 - 1, N4 - 1 )
               END IF
               IF( N4.GE.3 )THEN
                  X( JX - INCX ) = X( JX - INCX ) - X( JX )
     $                             *A( N4 - 1, N4 - 2 )
                  IF( NOUNIT )
     $               X( JX - INCX ) = X( JX - INCX )/A( N4 - 2, N4 - 2 )
               END IF
               IF( N4.GE.4 )THEN
                  X( JX - 2*INCX ) = X( JX - 2*INCX ) - X( JX )
     $                               *A( N4 - 1, N4 - 3 ) -
     $                               X( JX - INCX )*A( N4 - 2, N4 - 3 )
                  IF( NOUNIT )
     $               X( JX - 2*INCX ) = X( JX - 2*INCX )
     $                                  /A( N4 - 3, N4 - 3 )
               END IF
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06PJF (DTRSV ).
C
      END
      SUBROUTINE F06SAF( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     DOT VERSION FOR VECTOR MACHINES
C     .. Entry Points ..
      ENTRY      ZGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
C     .. Scalar Arguments ..
      COMPLEX*16         ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
C     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), X( * ), Y( * )
C     ..
C
C  Purpose
C  =======
C
C  ZGEMV  performs one of the matrix-vector operations
C
C     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   or
C
C     y := alpha*conjg( A' )*x + beta*y,
C
C  where alpha and beta are scalars, x and y are vectors and A is an
C  m by n matrix.
C
C  Parameters
C  ==========
C
C  TRANS  - CHARACTER*1.
C           On entry, TRANS specifies the operation to be performed as
C           follows:
C
C              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
C
C              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
C
C              TRANS = 'C' or 'c'   y := alpha*conjg( A' )*x + beta*y.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C           On entry, M specifies the number of rows of the matrix A.
C           M must be at least zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the number of columns of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - COMPLEX*16      .
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  A      - COMPLEX*16       array of DIMENSION ( LDA, n ).
C           Before entry, the leading m by n part of the array A must
C           contain the matrix of coefficients.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, m ).
C           Unchanged on exit.
C
C  X      - COMPLEX*16       array of DIMENSION at least
C           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
C           and at least
C           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
C           Before entry, the incremented array X must contain the
C           vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  BETA   - COMPLEX*16      .
C           On entry, BETA specifies the scalar beta. When BETA is
C           supplied as zero then Y need not be set on input.
C           Unchanged on exit.
C
C  Y      - COMPLEX*16       array of DIMENSION at least
C           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
C           and at least
C           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
C           Before entry with BETA non-zero, the incremented array Y
C           must contain the vector y. On exit, Y is overwritten by the
C           updated vector y.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER        ( ONE  = ( 1.0D+0, 0.0D+0 ) )
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Local Scalars ..
      COMPLEX*16         TEMP
      INTEGER            I, INFO, IX, IY, J, JY, KX, KY, LENX, LENY
      LOGICAL            NOCONJ
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.(TRANS.EQ.'N' .OR. TRANS.EQ.'n').AND.
     $         .NOT.(TRANS.EQ.'T' .OR. TRANS.EQ.'t').AND.
     $         .NOT.(TRANS.EQ.'C' .OR. TRANS.EQ.'c')      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06SAF/ZGEMV ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
C
      NOCONJ = (TRANS.EQ.'T' .OR. TRANS.EQ.'t')
C
C     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
C     up the start points in  X  and  Y.
C
      IF( (TRANS.EQ.'N' .OR. TRANS.EQ.'n') )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
C
C     Start the operations. In this version the inner loops are all
C     equivalent to dot-product operations.
C
C     First form  y := beta*y.
C
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      JY = KY
      IF( (TRANS.EQ.'N' .OR. TRANS.EQ.'n') )THEN
C
C        Form  y := alpha*A*x + y.
C
         IF( INCX.EQ.1 )THEN
            DO 60, J = 1, M
               TEMP = ZERO
               DO 50, I = 1, N
                  TEMP = TEMP + A( J, I )*X( I )
   50          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
   60       CONTINUE
         ELSE
            DO 80, J = 1, M
               TEMP = ZERO
               IX   = KX
               DO 70, I = 1, N
                  TEMP = TEMP + A( J, I )*X( IX )
                  IX   = IX   + INCX
   70          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
   80       CONTINUE
         END IF
      ELSE
C
C        Form  y := alpha*A'*x + y  or  y := alpha*conjg( A' )*x + y.
C
         IF( INCX.EQ.1 )THEN
            DO 110, J = 1, N
               TEMP = ZERO
               IF( NOCONJ )THEN
                  DO 90, I = 1, M
                     TEMP = TEMP + A( I, J )*X( I )
   90             CONTINUE
               ELSE
                  DO 100, I = 1, M
                     TEMP = TEMP + DCONJG( A( I, J ) )*X( I )
  100             CONTINUE
               END IF
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  110       CONTINUE
         ELSE
            DO 140, J = 1, N
               TEMP = ZERO
               IX   = KX
               IF( NOCONJ )THEN
                  DO 120, I = 1, M
                     TEMP = TEMP + A( I, J )*X( IX )
                     IX   = IX   + INCX
  120             CONTINUE
               ELSE
                  DO 130, I = 1, M
                     TEMP = TEMP + DCONJG( A( I, J ) )*X( IX )
                     IX   = IX   + INCX
  130             CONTINUE
               END IF
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  140       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06SAF (ZGEMV ).
C
      END
      SUBROUTINE F06SJF( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     DOT VERSION FOR VECTOR MACHINES
C     .. Entry Points ..
      ENTRY      ZTRSV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
C     .. Scalar Arguments ..
      INTEGER            INCX, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
C     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), X( * )
C     ..
C
C  Purpose
C  =======
C
C  ZTRSV  solves one of the systems of equations
C
C     A*x = b,   or   A'*x = b,   or   conjg( A' )*x = b,
C
C  where b and x are n element vectors and A is an n by n unit, or
C  non-unit, upper or lower triangular matrix.
C
C  No test for singularity or near-singularity is included in this
C  routine. Such tests must be performed before calling this routine.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the matrix is an upper or
C           lower triangular matrix as follows:
C
C              UPLO = 'U' or 'u'   A is an upper triangular matrix.
C
C              UPLO = 'L' or 'l'   A is a lower triangular matrix.
C
C           Unchanged on exit.
C
C  TRANS  - CHARACTER*1.
C           On entry, TRANS specifies the equations to be solved as
C           follows:
C
C              TRANS = 'N' or 'n'   A*x = b.
C
C              TRANS = 'T' or 't'   A'*x = b.
C
C              TRANS = 'C' or 'c'   conjg( A' )*x = b.
C
C           Unchanged on exit.
C
C  DIAG   - CHARACTER*1.
C           On entry, DIAG specifies whether or not A is unit
C           triangular as follows:
C
C              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
C
C              DIAG = 'N' or 'n'   A is not assumed to be unit
C                                  triangular.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  A      - COMPLEX*16       array of DIMENSION ( LDA, n ).
C           Before entry with  UPLO = 'U' or 'u', the leading n by n
C           upper triangular part of the array A must contain the upper
C           triangular matrix and the strictly lower triangular part of
C           A is not referenced.
C           Before entry with UPLO = 'L' or 'l', the leading n by n
C           lower triangular part of the array A must contain the lower
C           triangular matrix and the strictly upper triangular part of
C           A is not referenced.
C           Note that when  DIAG = 'U' or 'u', the diagonal elements of
C           A are not referenced either, but are assumed to be unity.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, n ).
C           Unchanged on exit.
C
C  X      - COMPLEX*16       array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the n
C           element right-hand side vector b. On exit, X is overwritten
C           with the solution vector x.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Local Scalars ..
      COMPLEX*16         TEMP
      INTEGER            I, INFO, IX, J, JX, KX
      LOGICAL            NOCONJ, NOUNIT
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.(UPLO .EQ.'U' .OR. UPLO .EQ.'u').AND.
     $         .NOT.(UPLO .EQ.'L' .OR. UPLO .EQ.'l')      )THEN
         INFO = 1
      ELSE IF( .NOT.(TRANS.EQ.'N' .OR. TRANS.EQ.'n').AND.
     $         .NOT.(TRANS.EQ.'T' .OR. TRANS.EQ.'t').AND.
     $         .NOT.(TRANS.EQ.'C' .OR. TRANS.EQ.'c')      )THEN
         INFO = 2
      ELSE IF( .NOT.(DIAG .EQ.'U' .OR. DIAG .EQ.'u').AND.
     $         .NOT.(DIAG .EQ.'N' .OR. DIAG .EQ.'n')      )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06SJF/ZTRSV ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 )
     $   RETURN
C
      NOCONJ = (TRANS.EQ.'T' .OR. TRANS.EQ.'t')
      NOUNIT = (DIAG .EQ.'N' .OR. DIAG .EQ.'n')
C
C     Set up the start point in X if the increment is not unity. This
C     will be  ( N - 1 )*INCX  too small for descending loops.
C
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
C
C     Start the operations. In this version the inner loops are all
C     equivalent to dot-product operations.
C
      IF( (TRANS.EQ.'N' .OR. TRANS.EQ.'n') )THEN
C
C        Form  x := inv( A )*x.
C
         IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
            IF( INCX.EQ.1 )THEN
               DO 20, J = N, 1, -1
                  TEMP = X( J )
                  DO 10, I = N, J + 1, -1
                     TEMP = TEMP - A( J, I )*X( I )
   10             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( J ) = TEMP
   20          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 40, J = N, 1, -1
                  TEMP = X( JX )
                  IX   = KX
                  DO 30, I = N, J + 1, -1
                     TEMP = TEMP - A( J, I )*X( IX )
                     IX   = IX   - INCX
   30             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( JX ) = TEMP
                  JX      = JX   - INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 60, J = 1, N
                  TEMP = X( J )
                  DO 50, I = 1, J - 1
                     TEMP = TEMP - A( J, I )*X( I )
   50             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( J ) = TEMP
   60          CONTINUE
            ELSE
               JX = KX
               DO 80, J = 1, N
                  TEMP = X( JX )
                  IX   = KX
                  DO 70, I = 1, J - 1
                     TEMP = TEMP - A( J, I )*X( IX )
                     IX   = IX   + INCX
   70             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( JX ) = TEMP
                  JX      = JX   + INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
C
C        Form  x := inv( A' )*x  or  x := inv( conjg( A' ) )*x.
C
         IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
            IF( INCX.EQ.1 )THEN
               DO 110, J = 1, N
                  TEMP = X( J )
                  IF( NOCONJ )THEN
                     DO 90, I = 1, J - 1
                        TEMP = TEMP - A( I, J )*X( I )
   90                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/A( J, J )
                  ELSE
                     DO 100, I = 1, J - 1
                        TEMP = TEMP - DCONJG( A( I, J ) )*X( I )
  100                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/DCONJG( A( J, J ) )
                  END IF
                  X( J ) = TEMP
  110          CONTINUE
            ELSE
               JX = KX
               DO 140, J = 1, N
                  IX   = KX
                  TEMP = X( JX )
                  IF( NOCONJ )THEN
                     DO 120, I = 1, J - 1
                        TEMP = TEMP - A( I, J )*X( IX )
                        IX   = IX   + INCX
  120                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/A( J, J )
                  ELSE
                     DO 130, I = 1, J - 1
                        TEMP = TEMP - DCONJG( A( I, J ) )*X( IX )
                        IX   = IX   + INCX
  130                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/DCONJG( A( J, J ) )
                  END IF
                  X( JX ) = TEMP
                  JX      = JX   + INCX
  140          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 170, J = N, 1, -1
                  TEMP = X( J )
                  IF( NOCONJ )THEN
                     DO 150, I = N, J + 1, -1
                        TEMP = TEMP - A( I, J )*X( I )
  150                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/A( J, J )
                  ELSE
                     DO 160, I = N, J + 1, -1
                        TEMP = TEMP - DCONJG( A( I, J ) )*X( I )
  160                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/DCONJG( A( J, J ) )
                  END IF
                  X( J ) = TEMP
  170          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 200, J = N, 1, -1
                  IX   = KX
                  TEMP = X( JX )
                  IF( NOCONJ )THEN
                     DO 180, I = N, J + 1, -1
                        TEMP = TEMP - A( I, J )*X( IX )
                        IX   = IX   - INCX
  180                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/A( J, J )
                  ELSE
                     DO 190, I = N, J + 1, -1
                        TEMP = TEMP - DCONJG( A( I, J ) )*X( IX )
                        IX   = IX   - INCX
  190                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/DCONJG( A( J, J ) )
                  END IF
                  X( JX ) = TEMP
                  JX      = JX   - INCX
  200          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06SJF (ZTRSV ).
C
      END
      INTEGER FUNCTION P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-621 (APR 1988).
C     MARK 13B REVISED. IER-668 (AUG 1988).
C
C     P01ABF is the error-handling routine for the NAG Library.
C
C     P01ABF either returns the value of IERROR through the routine
C     name (soft failure), or terminates execution of the program
C     (hard failure). Diagnostic messages may be output.
C
C     If IERROR = 0 (successful exit from the calling routine),
C     the value 0 is returned through the routine name, and no
C     message is output
C
C     If IERROR is non-zero (abnormal exit from the calling routine),
C     the action taken depends on the value of IFAIL.
C
C     IFAIL =  1: soft failure, silent exit (i.e. no messages are 
C                 output)
C     IFAIL = -1: soft failure, noisy exit (i.e. messages are output)
C     IFAIL =-13: soft failure, noisy exit but standard messages from
C                 P01ABF are suppressed
C     IFAIL =  0: hard failure, noisy exit
C
C     For compatibility with certain routines included before Mark 12
C     P01ABF also allows an alternative specification of IFAIL in which
C     it is regarded as a decimal integer with least significant digits
C     cba. Then
C
C     a = 0: hard failure  a = 1: soft failure
C     b = 0: silent exit   b = 1: noisy exit
C
C     except that hard failure now always implies a noisy exit.
C
C     S.Hammarling, M.P.Hooper and J.J.du Croz, NAG Central Office.
C
C     .. Scalar Arguments ..
      INTEGER                 IERROR, IFAIL, NREC
      CHARACTER*(*)           SRNAME
C     .. Array Arguments ..
      CHARACTER*(*)           REC(*)
C     .. Local Scalars ..
      INTEGER                 I, NERR
      CHARACTER*72            MESS
C     .. External Subroutines ..
      EXTERNAL                P01ABZ, X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC               ABS, MOD
C     .. Executable Statements ..
      IF (IERROR.NE.0) THEN
C        Abnormal exit from calling routine
         IF (IFAIL.EQ.-1 .OR. IFAIL.EQ.0 .OR. IFAIL.EQ.-13 .OR.
     *       (IFAIL.GT.0 .AND. MOD(IFAIL/10,10).NE.0)) THEN
C           Noisy exit
            CALL X04AAF(0,NERR)
            DO 20 I = 1, NREC
               CALL X04BAF(NERR,REC(I))
   20       CONTINUE
            IF (IFAIL.NE.-13) THEN
               WRITE (MESS,FMT=99999) SRNAME, IERROR
               CALL X04BAF(NERR,MESS)
               IF (ABS(MOD(IFAIL,10)).NE.1) THEN
C                 Hard failure
                  CALL X04BAF(NERR,
     *                     ' ** NAG hard failure - execution terminated'
     *                        )
                  CALL P01ABZ
               ELSE
C                 Soft failure
                  CALL X04BAF(NERR,
     *                        ' ** NAG soft failure - control returned')
               END IF
            END IF
         END IF
      END IF
      P01ABF = IERROR
      RETURN
C
99999 FORMAT (' ** ABNORMAL EXIT from NAG Library routine ',A,': IFAIL',
     *  ' =',I6)
      END
      SUBROUTINE P01ABZ
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C
C     Terminates execution when a hard failure occurs.
C
C     ******************** IMPLEMENTATION NOTE ********************
C     The following STOP statement may be replaced by a call to an
C     implementation-dependent routine to display a message and/or
C     to abort the program.
C     *************************************************************
C     .. Executable Statements ..
      STOP
      END
      DOUBLE PRECISION FUNCTION X02AJF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS  (1/2)*B**(1-P)  IF ROUNDS IS .TRUE.
C     RETURNS  B**(1-P)  OTHERWISE
C
      DOUBLE PRECISION X02CON
C     DATA X02CON /Z'3CA0000000000001' /
      DATA X02CON /1.1102230246252D-12/
C     .. Executable Statements ..
      X02AJF = X02CON
      RETURN
      END
      DOUBLE PRECISION FUNCTION X02AMF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS THE 'SAFE RANGE' PARAMETER
C     I.E. THE SMALLEST POSITIVE MODEL NUMBER Z SUCH THAT
C     FOR ANY X WHICH SATISFIES X.GE.Z AND X.LE.1/Z
C     THE FOLLOWING CAN BE COMPUTED WITHOUT OVERFLOW, UNDERFLOW OR OTHER
C     ERROR
C
C        -X
C        1.0/X
C        SQRT(X)
C        LOG(X)
C        EXP(LOG(X))
C        Y**(LOG(X)/LOG(Y)) FOR ANY Y
C
      DOUBLE PRECISION X02CON
C     DATA X02CON /Z'0010000000000000' /
      DATA X02CON /2.2250738585072D-308/
C     .. Executable Statements ..
      X02AMF = X02CON
      RETURN
      END
      INTEGER FUNCTION X02BHF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS THE MODEL PARAMETER, B.
C
C     .. Executable Statements ..
      X02BHF =     2
      RETURN
      END
      INTEGER FUNCTION X02BJF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS THE MODEL PARAMETER, p.
C
C     .. Executable Statements ..
      X02BJF =    53
      RETURN
      END
      LOGICAL FUNCTION X02DAF(X)
C     MARK 8 RELEASE. NAG COPYRIGHT 1980.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     RETURNS .FALSE. IF THE SYSTEM SETS UNDERFLOWING QUANTITIES
C     TO ZERO, WITHOUT ANY ERROR INDICATION OR UNDESIRABLE WARNING
C     OR SYSTEM OVERHEAD.
C     RETURNS .TRUE. OTHERWISE, IN WHICH CASE CERTAIN LIBRARY
C     ROUTINES WILL TAKE SPECIAL PRECAUTIONS TO AVOID UNDERFLOW
C     (USUALLY AT SOME COST IN EFFICIENCY).
C
C     X IS A DUMMY ARGUMENT
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION        X
C     .. Executable Statements ..
      X02DAF =  .FALSE.
      RETURN
      END
      SUBROUTINE X02ZAZ
C     MARK 12 RE-ISSUE. NAG COPYRIGHT 1986.
C
C***********************************************************************
C
C     NAG version of the Stanford routine MCHPAR.
C     Sven Hammarling, NAG Central Office.
C
C     X02ZAZ sets machine parameters as follows:
C
C     WMACH(  1 ) = nbase  = base of floating-point arithmetic.
C     WMACH(  2 ) = ndigit = no. of base ( nbase ) digits in the
C                            mantissa
C     WMACH(  3 ) = eps    = relative machine accuracy. (X02AJF.)
C     WMACH(  4 ) = rteps  = sqrt( eps ).
C     WMACH(  5 ) = rmin   = small positive floating-point number whose
C                             reciprocal does not overflow.
C     WMACH(  6 ) = rtrmin = sqrt( rmin ).
C     WMACH(  7 ) = rmax   = 1/rmin
C     WMACH(  8 ) = rtrmax = sqrt( rmax ).
C     WMACH(  9 ) = undflw = 0 if underflow is not fatal, +ve otherwise.
C     WMACH( 10 ) = nin    = input  stream unit number. ( 5.)
C     WMACH( 11 ) = nout   = output stream unit number. ( X04ABF.)
C     WMACH( 12 )
C     WMACH( 13 )   Not currently used.
C     WMACH( 14 )
C     WMACH( 15 )
C
C     Note that constants that represent integers may hold a number just
C     less than the integer, so that the integer should be recovered by
C     adding, say, 0.25. e.g.
C
C     IBASE = WMACH( 1 ) + 0.25
C
C***********************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER       (ZERO=0.0D+0)
C     .. Arrays in Common ..
      DOUBLE PRECISION WMACH(15)
C     .. Local Scalars ..
      DOUBLE PRECISION EPS, RMAX, RMIN, UNDFLW
      INTEGER         NBASE, NDIGIT, NOUT
      LOGICAL         FIRST
C     .. External Functions ..
      DOUBLE PRECISION X02AJF, X02AMF
      INTEGER         X02BHF, X02BJF
      LOGICAL         X02DAF
      EXTERNAL        X02AJF, X02AMF, X02BHF, X02BJF, X02DAF
C     .. External Subroutines ..
      EXTERNAL        X04ABF
C     .. Intrinsic Functions ..
      INTRINSIC       SQRT
C     .. Common blocks ..
      COMMON          /AX02ZA/WMACH
C     .. Save statement ..
      SAVE            /AX02ZA/, FIRST
C     .. Data statements ..
      DATA            FIRST/.TRUE./
C     .. Executable Statements ..
C
      IF (FIRST) THEN
         FIRST = .FALSE.
C
         IF (X02DAF(ZERO)) THEN
            UNDFLW = 1
         ELSE
            UNDFLW = 0
         END IF
         NBASE = X02BHF()
         NDIGIT = X02BJF()
         EPS = X02AJF()
         RMIN = X02AMF()
         RMAX = 1/RMIN
C
         WMACH(1) = NBASE
         WMACH(2) = NDIGIT
         WMACH(3) = EPS
         WMACH(4) = SQRT(EPS)
         WMACH(5) = RMIN
         WMACH(6) = SQRT(RMIN)
         WMACH(7) = RMAX
         WMACH(8) = SQRT(RMAX)
         WMACH(9) = UNDFLW
      END IF
      CALL X04ABF(0,NOUT)
      WMACH(10) = 5
      WMACH(11) = NOUT
      RETURN
C
C     End of  X02ZAZ. (MCHPAR)
C
      END
      SUBROUTINE X03AAF(A,ISIZEA,B,ISIZEB,N,ISTEPA,ISTEPB,C1,C2,D1,D2,
     *                  SW,IFAIL)
C     MARK 10 RE-ISSUE. NAG COPYRIGHT 1982
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12 REVISED. IER-524 (AUG 1986).
C     DOUBLE PRECISION BASE VERSION
C
C     CALCULATES THE VALUE OF A SCALAR PRODUCT USING BASIC OR
C     ADDITIONAL PRECISION AND ADDS IT TO A BASIC OR ADDITIONAL
C     PRECISION INITIAL VALUE.
C
C     FOR THIS DOUBLE PRECISION VERSION, ALL ADDITIONAL (I.E.
C     QUADRUPLE) PRECISION COMPUTATION IS PERFORMED BY THE AUXILIARY
C     ROUTINE X03AAY.  SEE THE COMMENTS AT THE HEAD OF THAT ROUTINE
C     CONCERNING IMPLEMENTATION.
C
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='X03AAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  C1, C2, D1, D2
      INTEGER           IFAIL, ISIZEA, ISIZEB, ISTEPA, ISTEPB, N
      LOGICAL           SW
C     .. Array Arguments ..
      DOUBLE PRECISION  A(ISIZEA), B(ISIZEB)
C     .. Local Scalars ..
      DOUBLE PRECISION  SUM
      INTEGER           I, IERR, IS, IT
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          X03AAY
C     .. Executable Statements ..
      IERR = 0
      IF (ISTEPA.LE.0 .OR. ISTEPB.LE.0) IERR = 1
      IF (ISIZEA.LT.(N-1)*ISTEPA+1 .OR. ISIZEB.LT.(N-1)*ISTEPB+1)
     *    IERR = 2
      IF (IERR.EQ.0) GO TO 20
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,P01REC)
      RETURN
   20 IFAIL = 0
      IF (SW) GO TO 80
C
C     BASIC PRECISION CALCULATION
C
      SUM = C1
      IF (N.LT.1) GO TO 60
      IS = 1
      IT = 1
      DO 40 I = 1, N
         SUM = SUM + A(IS)*B(IT)
         IS = IS + ISTEPA
         IT = IT + ISTEPB
   40 CONTINUE
   60 D1 = SUM
      D2 = 0.0D0
      RETURN
C
C     ADDITIONAL PRECISION COMPUTATION
C
   80 CALL X03AAY(A,ISIZEA,B,ISIZEB,N,ISTEPA,ISTEPB,C1,C2,D1,D2)
      RETURN
      END
      SUBROUTINE X03AAY(A,ISIZEA,B,ISIZEB,N,ISTEPA,ISTEPB,C1,C2,D1,D2)
C     MARK 10 RELEASE. NAG COPYRIGHT 1982
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     DOUBLE PRECISION BASE VERSION
C
C     PERFORMS QUADRUPLE PRECISION COMPUTATION FOR X03AAF
C
C     **************************************************************
C     THIS FORTRAN CODE WILL NOT WORK ON ALL MACHINES. IT IS
C     PRESUMED TO WORK IF THE MACHINE SATISFIES ONE OF THE FOLLOWING
C     ASSUMPTIONS.
C
C     A.) THERE ARE AN EVEN NUMBER OF B-ARY DIGITS IN THE MANTISSA
C     -   OF A DOUBLE PRECISION NUMBER (WHERE B IS THE BASE FOR THE
C     -   REPRESENTATION OF FLOATING-POINT NUMBERS), AND THE
C     -   COMPUTED RESULT OF A DOUBLE PRECISION ADDITION,
C     -   SUBTRACTION OR MULTIPLICATION IS EITHER CORRECTLY ROUNDED
C     -   OR CORRECTLY CHOPPED.
C
C     B.) FLOATING-POINT NUMBERS ARE REPRESENTED TO THE BASE 2 (WITH
C     -   ANY NUMBER OF BITS IN THE MANTISSA OF A DOUBLE PRECISION
C     -   NUMBER), AND THE COMPUTED RESULT OF A DOUBLE PRECISION
C     -   ADDITION, SUBTRACTION OR MULTIPLICATION IS CORRECTLY
C     -   ROUNDED.
C
C     REFERENCES-
C
C     T.J. DEKKER  A FLOATING-POINT TECHNIQUE FOR EXTENDING THE
C     AVAILABLE PRECISION. NUMER. MATH. 18, 224-242 (1971)
C
C     S. LINNAINMAA  SOFTWARE FOR DOUBLED-PRECISION FLOATING-POINT
C     COMPUTATIONS. ACM TRANS. MATH. SOFTWARE 7, 272-283 (1981)
C
C     IF THE ABOVE ASSUMPTIONS ARE NOT SATISFIED, THIS ROUTINE MUST
C     BE IMPLEMENTED IN ASSEMBLY LANGUAGE. IN ANY CASE ASSEMBLY
C     LANGUAGE MAY BE PREFERABLE FOR GREATER EFFICIENCY.  CONSULT
C     NAG CENTRAL OFFICE.
C
C     THE ROUTINE MUST SIMULATE THE FOLLOWING QUADRUPLE PRECISION
C     CODING IN PSEUDO-FORTRAN, WHERE
C     - QEXTD CONVERTS FROM DOUBLE TO QUADRUPLE PRECISION
C     - DBLEQ CONVERTS FROM QUADRUPLE TO DOUBLE PRECISION.
C
C     QUADRUPLE PRECISION SUM
C     SUM = QEXTD(C1) + QEXTD(C2)
C     IF (N.LT.1) GO TO 80
C     IS = 1
C     IT = 1
C     DO 60 I = 1, N
C        SUM = SUM + QEXTD(A(IS))*QEXTD(B(IT))
C        IS = IS + ISTEPA
C        IT = IT + ISTEPB
C  60 CONTINUE
C  80 D1 = SUM -- ROUNDED TO DOUBLE PRECISION
C     D2 = DBLEQ(SUM-QEXTD(D1))
C
C     **************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  C1, C2, D1, D2
      INTEGER           ISIZEA, ISIZEB, ISTEPA, ISTEPB, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(ISIZEA), B(ISIZEB)
C     .. Local Scalars ..
      DOUBLE PRECISION  AA, BB, CONS, HA, HB, R, S, SUM, SUMM, TA, TB,
     *                  Z, ZZ
      INTEGER           I, IS, IT
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Data statements ..
C     ************* IMPLEMENTATION-DEPENDENT CONSTANT **************
C     CONS MUST BE SET TO  B**(T - INT(T/2)) + 1 , WHERE T IS THE
C     NUMBER OF B-ARY DIGITS IN THE MANTISSA OF A DOUBLE PRECISION
C     NUMBER.
C     FOR B = 16 AND T = 14 (E.G. IBM 370) OR
C     FOR B = 2 AND T = 56 (E.G. DEC VAX-11)
C     DATA CONS /268435457.0D0/
C     DATA CONS /Z'41A0000002000000'/ 
      DATA CONS /134217729.0D0/
C     **************************************************************
C     .. Executable Statements ..
      SUM = C1 + C2
      SUMM = (C1-SUM) + C2
      IF (N.LT.1) GO TO 80
      IS = 1
      IT = 1
      DO 60 I = 1, N
         AA = A(IS)
         BB = B(IT)
         Z = AA*CONS
         HA = (AA-Z) + Z
         TA = AA - HA
         Z = BB*CONS
         HB = (BB-Z) + Z
         TB = BB - HB
         Z = AA*BB
         ZZ = (((HA*HB-Z)+HA*TB)+TA*HB) + TA*TB
         R = Z + SUM
         IF (ABS(Z).GT.ABS(SUM)) GO TO 20
         S = (((SUM-R)+Z)+ZZ) + SUMM
         GO TO 40
   20    S = (((Z-R)+SUM)+SUMM) + ZZ
   40    SUM = R + S
         SUMM = (R-SUM) + S
         IS = IS + ISTEPA
         IT = IT + ISTEPB
   60 CONTINUE
C  80 D1 = SUM + (SUMM+SUMM)
C     *************** IMPLEMENTATION DEPENDENT CODE ****************
C     THE PREVIOUS STATEMENT ASSUMES THAT THE RESULT OF A DOUBLE
C     PRECISION ADDITION IS TRUNCATED.  IF IT IS ROUNDED, THEN
C     THE STATEMENT MUST BE CHANGED TO
   80 D1 = SUM + SUMM
C     **************************************************************
      D2 = (SUM-D1) + SUMM
      RETURN
      END
      SUBROUTINE X04AAF(I,NERR)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978
C     MARK 7C REVISED IER-190 (MAY 1979)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     IF I = 0, SETS NERR TO CURRENT ERROR MESSAGE UNIT NUMBER
C     (STORED IN NERR1).
C     IF I = 1, CHANGES CURRENT ERROR MESSAGE UNIT NUMBER TO
C     VALUE SPECIFIED BY NERR.
C
C     .. Scalar Arguments ..
      INTEGER           I, NERR
C     .. Local Scalars ..
      INTEGER           NERR1
C     .. Save statement ..
      SAVE              NERR1
C     .. Data statements ..
      DATA              NERR1/0/
C     .. Executable Statements ..
      IF (I.EQ.0) NERR = NERR1
      IF (I.EQ.1) NERR1 = NERR
      RETURN
      END
      SUBROUTINE X04ABF(I,NADV)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978
C     MARK 7C REVISED IER-190 (MAY 1979)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C      IF I = 0, SETS NADV TO CURRENT ADVISORY MESSAGE UNIT NUMBER
C     (STORED IN NADV1).
C     IF I = 1, CHANGES CURRENT ADVISORY MESSAGE UNIT NUMBER TO
C     VALUE SPECIFIED BY NADV.
C
C     .. Scalar Arguments ..
      INTEGER           I, NADV
C     .. Local Scalars ..
      INTEGER           NADV1
C     .. Save statement ..
      SAVE              NADV1
C     .. Data statements ..
      DATA              NADV1/6/
C     .. Executable Statements ..
      IF (I.EQ.0) NADV = NADV1
      IF (I.EQ.1) NADV1 = NADV
      RETURN
      END
      SUBROUTINE X04BAF(NOUT,REC)
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C
C     X04BAF writes the contents of REC to the unit defined by NOUT.
C
C     Trailing blanks are not output, except that if REC is entirely
C     blank, a single blank character is output.
C     If NOUT.lt.0, i.e. if NOUT is not a valid Fortran unit identifier,
C     then no output occurs.
C
C     .. Scalar Arguments ..
      INTEGER           NOUT
      CHARACTER*(*)     REC
C     .. Local Scalars ..
      INTEGER           I
C     .. Intrinsic Functions ..
      INTRINSIC         LEN
C     .. Executable Statements ..
      IF (NOUT.GE.0) THEN
C        Remove trailing blanks
         DO 20 I = LEN(REC), 2, -1
            IF (REC(I:I).NE.' ') GO TO 40
   20    CONTINUE
C        Write record to external file
   40    WRITE (NOUT,FMT=99999) REC(1:I)
      END IF
      RETURN
C
99999 FORMAT (A)
      END


