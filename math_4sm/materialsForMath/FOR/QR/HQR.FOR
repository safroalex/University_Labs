      SUBROUTINE HQR(NM,N,LOW,IGH,H,WR,WI,IERR)
C
      INTEGER I,J,K,L,M,N,EN,LL,MM,NA,IGH,ITS,LOW,MP2,ENM2,IERR
      REAL H(NM,N),WR(N),WI(N)
      REAL P,Q,R,S,T,W,X,Y,ZZ,MACHEP
      REAL SQRT,ABS,SIGN
      INTEGER MIN0
      LOGICAL NOTLAS
C
       MACHEP=1.0
   10 MACHEP=MACHEP/2.
      TOL1=1.+MACHEP
      IF(TOL1.GT.1.)GO TO 10
C     STORE ROOTS ISOLATED BY BALANC
      IERR=0
      DO 50 I=1,N
      IF(I.GE.LOW.AND.I.LE.IGH)GO TO 50
      WR(I)=H(I,I)
      WI(I)=0.0
   50 CONTINUE
C
      EN=IGH
      T=0.0
C     SEARCH FOR NEXT EIGENVALUES
   60 IF(EN.LT.LOW)GO TO 1001
      ITS=0
      NA=EN-1
      ENM2=NA-1
C     LOOK FOR SINGLE SMALL SUB-DIAGONAL FOR L=EN STEP-1 UNTIL LOW DO
   70 DO 80 LL=LOW,EN
      L=EN+LOW-LL
      IF(L.EQ.LOW)GO TO 100
      IF(ABS(H(L,L-1)).LE.MACHEP*(ABS(H(L-1,L-1))
     *   +ABS(H(L,L))))GO TO 100
   80 CONTINUE
C     FORM SHIFT
  100 X=H(EN,EN)
      IF(L.EQ.EN)GO TO 270
      Y=H(NA,NA)
      W=H(EN,NA)*H(NA,EN)
      IF(L.EQ.NA)GO TO 280
      IF(ITS.EQ.30)GO TO 1000
      IF(ITS.NE.10.AND.ITS.NE.20)GO TO 130
C     FORM EXEPTIONAL SHIFT
       T=T+X
C
       DO 120 I=LOW,EN
  120 H(I,I)=H(I,I)-X
C
      S=ABS(H(EN,NA))+ABS(H(NA,ENM2))
      X=0.75*S
      Y=X
      W=-0.4375*S*S
  130 ITS=ITS+1
C     LOOK FOR TWO CONSECUTIVE SMALL SUB-DIAGONAL ELEMENTS
C               FOR M=EN-2 STEP-1 UNTIL L DO...
      DO 140 MM=L,ENM2
      M=ENM2+L-MM
      ZZ=H(M,M)
      R=X-ZZ
      S=Y-ZZ
      P=(R*S-W)/H(M+1,M)+H(M,M+1)
      Q=H(M+1,M+1)-ZZ-R-S
      R=H(M+2,M+1)
      S=ABS(P)+ABS(Q)+ABS(R)
      P=P/S
      Q=Q/S
      R=R/S
      IF(M.EQ.L)GO TO 150
      IF(ABS(H(M,M-1))*(ABS(Q)+ABS(R)).LE.MACHEP*ABS(P)
     **(ABS(H(M-1,M-1))+ABS(ZZ)+ABS(H(M+1,M+1))))GO TO 150
  140 CONTINUE
C
  150 MP2=M+2
C
      DO 160 I=MP2,EN
      H(I,I-2)=0.0
      IF(I.EQ.MP2) GO TO 160
      H(I,I-3)=0.0
  160 CONTINUE
C     DOUBLE QR STEP INVOLVING ROWS TO EN AND COLUMNS M TO EN
      DO 260 K=M,NA
      NOTLAS=K.NE.NA
      IF(K.EQ.M) GO TO 170
      P=H(K,K-1)
      Q=H(K+1,K-1)
      R=0.0
      IF(NOTLAS) R=H(K+2,K-1)
      X=ABS(P)+ABS(Q)+ABS(R)
      IF(X.EQ.0.0) GO TO 260
      P=P/X
      Q=Q/X
      R=R/X
  170 S=SIGN(SQRT(P*P+Q*Q+R*R),P)
      IF(K.EQ.M) GO TO 180
      H(K,K-1)=-S*X
      GO TO 190
  180 IF(L.NE.M) H(K,K-1)=-H(K,K-1)
  190 P=P+S
      X=P/S
      Y=Q/S
      ZZ=R/S
      Q=Q/P
      R=R/P
C     ROW MODIFICATION
      DO 210 J=K,EN
      P=H(K,J)+Q*H(K+1,J)
      IF(.NOT.NOTLAS) GO TO 200
      P=P+R*H(K+2,J)
      H(K+2,J)=H(K+2,J)-P*ZZ
  200 H(K+1,J)=H(K+1,J)-P*Y
      H(K,J)=H(K,J)-P*X
  210 CONTINUE
C     COLUMN MODIFICATION
      J=MIN0(EN,K+3)
C
      DO 230 I=L,J
      P=X*H(I,K)+Y*H(I,K+1)
      IF(.NOT.NOTLAS) GO TO 220
      P=P+ZZ*H(I,K+2)
      H(I,K+2)=H(I,K+2)-P*R
  220 H(I,K+1)=H(I,K+1)-P*Q
      H(I,K)=H(I,K)-P
  230 CONTINUE
C
  260 CONTINUE
C
      GO TO 70
C     ONE ROOT FOUND
  270 WR(EN)=X+T
      WI(EN)=0.0
      EN=NA
      GO TO 60
C
  280 P=(Y-X)/2.0
      Q=P*P+W
      ZZ=SQRT(ABS(Q))
      X=X+T
      IF(Q.LT.0.0) GO TO 320
C     REAL PAIR
      ZZ=P+SIGN(ZZ,P)
      WR(NA)=X+ZZ
      WR(EN)=WR(NA)
      IF(ZZ.NE.0.0) WR(EN)=X-W/ZZ
      WI(NA)=0.0
      WI(EN)=0.0
      GO TO 330
C     COMPLEX PAIR
  320 WR(NA)=X+P
      WR(EN)=X+P
      WI(NA)=ZZ
      WI(EN)=-ZZ
  330 EN=ENM2
      GO TO 60
C     SET ERROR_NO CONVERGENCE TO EN EIGENVALUE AFTER 30 ITERATIONS
 1000 IERR=EN
 1001 RETURN
      END
