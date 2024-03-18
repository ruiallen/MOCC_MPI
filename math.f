
      subroutine gaussv (a,b,s, n,nop, n1,n4)
      implicit real*8 (a-h,o-z)
c*****solves n by n matrix equation  by gaussian elimination:
c...      a(p,i,j)x(p,j,k)=b(p,i,k)
c...   with vectorisation on passive index p<=nol
c...      s(nop) is scratch space
c*****answer is returned in b, a is destroyed.
       integer  p,n,n1,n4,nop
       double precision a(n4,n1,n1),b(n4,n1,n1), s(1)
      if(n.eq.1) go to 280
 2030  format("  gauss n nop n1 n4 ",5i5)
      ns = n - 1
      do 200 i=1,ns
      ig = i + 1
       do 1 p=1,nop
    1  s(p)=1.0d0/a(p,i,i)
      do 120 j=ig,n
       do 120 p=1,nop
  120  a(p,i,j)=a(p,i,j)*s(p)
      do 130 j=1,n
       do 130 p=1,nop
  130  b(p,i,j)=b(p,i,j)*s(p)
      do 150 k=1,n
      if(k-i) 151,150,151
  151  do 2 p=1,nop
    2  s(p)=a(p,k,i)
      do 140 l=ig,n
       do 140 p=1,nop
  140  a(p,k,l)=a(p,k,l) - a(p,i,l)*s(p)
      do 145 l=1,n
       do 145 p=1,nop
  145  b(p,k,l)=b(p,k,l) - b(p,i,l)*s(p)
150   continue
200   continue
       do 3 p=1,nop
    3  s(p)=1.0d0/a(p,n,n)
      do 230 j=1,n
       do 230 p=1,nop
  230  b(p,n,j)=b(p,n,j)*s(p)
      do 250 k=1,ns
       do 4 p=1,nop
    4  s(p)=a(p,k,n)
      do 250 l=1,n
       do 250 p=1,nop
  250  b(p,k,l)=b(p,k,l) - b(p,n,l)*s(p)
      return
  280  do 5 p=1,nop
    5  b(p,1,1)=b(p,1,1)/a(p,1,1)
      return
      end
c********************************************************************ZLB
      SUBROUTINE COULFG(XX,ETA1,XLMIN,XLMAX, FC,GC,FCP,GCP,
     X                  MODE1,KFN,IFAIL)
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C                                                                      C
C  Revised Coulomb wavefunction program using Steed's method           C
C                                                                      C
C  A. R. Barnett           Manchester  March   1981                    C
C                                                                      C
C  original program 'RCWFN'      in    CPC  8 (1974) 377-395           C
C                 + 'RCWFF'      in    CPC 11 (1976) 141-142           C
C  full description of algorithm in    CPC 21 (1981) 297-314           C
C  this version written up       in    CPC 27 (1982) 147-zzz           C
C                                                                      C
C  COULFG returns F,G,F',G', for real xx.GT.0,real eta1 (including 0), C
C   and real lambda(xlmin) .GT. -1 for integer-spaced lambda values    C
C   thus giving positive-energy solutions to the Coulomb Schrodinger   C
C   equation,to the Klein-Gordon equation and to suitable forms of     C
C   the Dirac equation ,also spherical & cylindrical Bessel equations  C
C                                                                      C
C  for a range of lambda values (xlmax - xlmin) must be an integer,    C
C  starting array element is m1 = max0(idint(xlmin+accur),0) + 1       C
C      see text for modifications for integer l-values                 C
C                                                                      C
C  if 'MODE' = 1  get F,G,F',G'   for integer-spaced lambda values     C
C            = 2      F,G      unused arrays must be dimensioned in    C
C            = 3      F               call to at least length (1)      C
C  if 'KFN'  = 0 real        Coulomb functions are returned            C
C            = 1 spherical   Bessel      "      "     "                C
C            = 2 cylindrical Bessel      "      "     "                C
C  the use of 'MODE' and 'KFN' is independent                          C
C                                                                      C
C  Precision:  results to within 2-3 decimals of 'machine accuracy'    C
C   in oscillating region x .GE. eta1 + sqrt(eta1**2 + xlm(xlm+1))     C
C   COULFG is coded for REAL*8 on IBM or equivalent  accur = 10**-16   C
C   use AUTODBL + extended precision on HX compiler  accur = 10**-33   C
C   for mantissas of 56 & 112 bits. For single precision CDC (48 bits) C
C   reassign SQRT=SQRT etc.  See text for complex arithmetic version  C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION    FC(1),GC(1),FCP(1),GCP(1)
      LOGICAL      ETANE0,XLTURN
C***  COULFG has calls to: SQRT,ABS,MOD,INT,SIGN,DFLOAT,MIN1
      DATA ZERO,ONE,TWO,TEN2,ABORT /0.0D0, 1.0D0, 2.0D0, 1.0D2, 2.0D4/
      DATA HALF,TM30 / 0.5D0, 1.0D-30 /
      DATA RT2DPI /0.79788 45608 02865 35587 98921 19868 76373 D0/
C *** This constant is  SQRT(TWO/PI):  use Q0 for IBM REAL*16: D0 for
C ***  REAL*8 & CDC double p:  E0 for CDC single p; and truncate value.
C
                        ACCUR = 1.0D-33
C ***            change accur to suit machine and precision required
      MODE  = 1
      IF(MODE1 .EQ. 2 .OR. MODE1 .EQ. 3 ) MODE = MODE1
      IFAIL = 0
      IEXP  = 1
      NPQ   = 0
      ETA   = ETA1
      GJWKB = ZERO
      PACCQ = ONE
      IF(KFN .NE. 0) ETA = ZERO
                 ETANE0  = ETA .NE. ZERO
      ACC   = ACCUR
      ACC4  = ACC*TEN2*TEN2
      ACCH  = dSQRT(ACC)
C ***    test range of XX, exit if.LE.SQRT(ACCUR) or if negative
C
      IF(XX .LE. ACCH)                          GO TO 120
      X     = XX
      XLM   = XLMIN
      IF(KFN .EQ. 2)  XLM = XLM - HALF
      IF(XLM .LE. -ONE .OR. XLMAX .LT. XLMIN)   GO TO 130
      E2MM1 = ETA*ETA + XLM*XLM + XLM
      XLTURN= X*(X - TWO*ETA) .LT. XLM*XLM + XLM
      DELL  = XLMAX - XLMIN + ACC
      IF(dABS(dMOD(DELL,ONE)) .GT. ACC) WRITE(6,1060)XLMAX,XLMIN,DELL
      LXTRA = idINT(DELL)
      XLL   = XLM + DFLOAT(LXTRA)
C ***       LXTRA is number of additional lambda values to be computed
C ***       XLL  is max lambda value, or 0.5 smaller for J,Y Bessels
C ***         determine starting array element (M1) from XLMIN
      M1  = MAX0(idINT(XLMIN + ACC),0) + 1
      L1  = M1 + LXTRA
C
C ***    evaluate cf1  =  f   =  fprime(xl,eta,x)/f(xl,eta,x)
C
      XI  = ONE/X
      FCL = ONE
      PK  = XLL + ONE
      PX  = PK  + ABORT
   10 EK  = ETA / PK
      F   = (EK + PK*XI)*FCL + (FCL - ONE)*XI
      PK1 =  PK + ONE
C ***   test ensures b1 .ne. zero for negative eta; fixup is exact.
             IF(dABS(ETA*X + PK*PK1) .GT. ACCH)  GO TO 20
             FCL  = (ONE + EK*EK)/(ONE + (ETA/PK1)**2)
             PK   =  TWO + PK
      GO TO 10
   20 D   =  ONE/((PK + PK1)*(XI + EK/PK1))
      DF  = -FCL*(ONE + EK*EK)*D
            IF(FCL .NE. ONE )  FCL = -ONE
            IF(D   .LT. ZERO)  FCL = -FCL
      F   =  F  + DF
C
C ***   begin cf1 loop on pk = k = lambda + 1
C
      P     = ONE
   30 PK    = PK1
        PK1 = PK1 + ONE
        EK  = ETA / PK
        TK  = (PK + PK1)*(XI + EK/PK1)
        D   =  TK - D*(ONE + EK*EK)
              IF(dABS(D) .GT. ACCH)             GO TO 40
              WRITE (6,1000) D,DF,ACCH,PK,EK,ETA,X
              P = P  +   ONE
              IF( P .GT. TWO )                  GO TO 140
   40 D     = ONE/D
              IF (D .LT. ZERO) FCL = -FCL
        DF  = DF*(D*TK - ONE)
        F   = F  + DF
              IF(PK .GT. PX)                    GO TO 140
      IF(dABS(DF) .GE. dABS(F)*ACC)             GO TO 30
                  NFP = PK - XLL - 1
      IF(LXTRA .EQ. 0)                          GO TO 60
C
C *** downward recurrence to lambda = XLM. Array GC,if present,stores RL
C
      FCL = FCL*TM30
      FPL = FCL*F
      IF(MODE .EQ. 1) FCP(L1) = FPL
                      FC (L1) = FCL
      XL  = XLL
      RL  = ONE
      EL  = ZERO
      DO 50  LP = 1,LXTRA
         IF(ETANE0) EL = ETA/XL
         IF(ETANE0) RL = dSQRT(ONE + EL*EL)
         SL    =  EL  + XL*XI
         L     =  L1  - LP
         FCL1  = (FCL *SL + FPL)/RL
         FPL   =  FCL1*SL - FCL *RL
         FCL   =  FCL1
         FC(L) =  FCL
         IF(MODE .EQ. 1) FCP(L)  = FPL
         IF(MODE .NE. 3 .AND. ETANE0) GC(L+1) = RL
   50 XL = XL - ONE
      IF(FCL .EQ. ZERO) FCL = ACC
      F  = FPL/FCL
C ***    now we have reached lambda = XLMIN = XLM
C ***    evaluate CF2 = P + I.Q  again using Steed's algorithm
C ***    see text for compact complex code for sp CDC or non-ANSI IBM
C
   60 IF( XLTURN ) CALL JWKB(X,ETA,dMAX1(XLM,ZERO),FJWKB,GJWKB,IEXP)
      IF( IEXP .GT. 1 .OR. GJWKB .GT. ONE/(ACCH*TEN2))  GO TO 80
          XLTURN = .FALSE.
      TA =  TWO*ABORT
      PK =  ZERO
      WI =  ETA + ETA
      P  =  ZERO
      Q  =  ONE - ETA*XI
      AR = -E2MM1
      AI =  ETA
      BR =  TWO*(X - ETA)
      BI =  TWO
      DR =  BR/(BR*BR + BI*BI)
      DI = -BI/(BR*BR + BI*BI)
      DP = -XI*(AR*DI + AI*DR)
      DQ =  XI*(AR*DR - AI*DI)
   70 P     = P  + DP
         Q  = Q  + DQ
         PK = PK + TWO
         AR = AR + PK
         AI = AI + WI
         BI = BI + TWO
         D  = AR*DR - AI*DI + BR
         DI = AI*DR + AR*DI + BI
         C  = ONE/(D*D + DI*DI)
         DR =  C*D
         DI = -C*DI
         A  = BR*DR - BI*DI - ONE
         B  = BI*DR + BR*DI
         C  = DP*A  - DQ*B
         DQ = DP*B  + DQ*A
         DP = C
         IF(PK .GT. TA)                         GO TO 150
      IF(dABS(DP)+dABS(DQ).GE.(dABS(P)+dABS(Q))*ACC)   GO TO 70
                      NPQ   = PK/TWO
                      PACCQ = HALF*ACC/dMIN1(dABS(Q),ONE)
                      IF(dABS(P) .GT. dABS(Q)) PACCQ = PACCQ*dABS(P)
C
C *** solve for FCM = F at lambda = XLM,then find norm factor W=W/FCM
C
      GAM = (F - P)/Q
            IF(Q .LE. ACC4*dABS(P))             GO TO 160
      W   = ONE/dSQRT((F - P)*GAM + Q)
            GO TO 90
C *** arrive here if G(XLM) .GT. 10**6 or IEXP .GT. 70 & XLTURN = .TRUE.
   80 W   = FJWKB
      GAM = GJWKB*W
      P   = F
      Q   = ONE
C
C *** normalise for spherical or cylindrical Bessel functions
C
   90                     ALPHA = ZERO
          IF(KFN  .EQ. 1) ALPHA = XI
          IF(KFN  .EQ. 2) ALPHA = XI*HALF
                          BETA  = ONE
          IF(KFN  .EQ. 1) BETA  = XI
          IF(KFN  .EQ. 2) BETA  = dSQRT(XI)*RT2DPI
      FCM  = dSIGN(W,FCL)*BETA
           FC(M1)  = FCM
                      IF(MODE .EQ. 3)           GO TO 100
           IF(.NOT. XLTURN)   GCL =  FCM*GAM
           IF(      XLTURN)   GCL =  GJWKB*BETA
           IF( KFN .NE. 0 )   GCL = -GCL
           GC(M1)  = GCL
           GPL =  GCL*(P - Q/GAM) - ALPHA*GCL
                      IF(MODE .EQ. 2)           GO TO 100
           GCP(M1) = GPL
           FCP(M1) = FCM*(F - ALPHA)
  100 IF(LXTRA .EQ. 0 ) RETURN
C *** upward recurrence from GC(M1),GCP(M1)  stored value is RL
C *** renormalise FC,FCP at each lambda and correct regular derivative
C ***    XL   = XLM here  and RL = ONE , EL = ZERO for Bessels
         W    = BETA*W/dABS(FCL)
         MAXL = L1 - 1
      DO 110 L = M1,MAXL
                      IF(MODE .EQ. 3)           GO TO 110
                      XL = XL + ONE
         IF(ETANE0)   EL = ETA/XL
         IF(ETANE0)   RL = GC(L+1)
                      SL = EL + XL*XI
         GCL1     = ((SL - ALPHA)*GCL - GPL)/RL
         GPL      =   RL*GCL -  (SL + ALPHA)*GCL1
         GCL      = GCL1
         GC(L+1)  = GCL1
                      IF(MODE .EQ. 2)           GO TO 110
         GCP(L+1) = GPL
         FCP(L+1) = W*(FCP(L+1) - ALPHA*FC(L+1))
  110 FC(L+1)     = W* FC(L+1)
      RETURN
 1000 FORMAT(/,' CF1 ACCURACY LOSS: D,DF,ACCH,K,ETA/K,ETA,X = ',1P,
     X       7D9.2,/)
C
C ***    error messages
C
  120 IFAIL = -1
      WRITE(6,1010) XX,ACCH
 1010 FORMAT(' FOR XX = ',1P,D12.3,' TRY SMALL-X  SOLUTIONS',
     X' OR X NEGATIVE',/ ,' SQUARE ROOT ACCURACY PARAMETER =  ',D12.3,/)
      RETURN
  130 IFAIL = -2
      WRITE (6,1020) XLMAX,XLMIN,XLM
 1020 FORMAT(/,' PROBLEM WITH INPUT ORDER VALUES:XLMAX,XLMIN,XLM = ',
     X1P,3D15.6,/)
      RETURN
  140 IFAIL =  1
      WRITE (6,1030) ABORT,F ,DF,PK,PX,ACC
 1030 FORMAT(' CF1 HAS FAILED TO CONVERGE AFTER ',F10.0,' ITERATIONS',/
     X,' F,DF,PK,PX,ACCUR =  ',1P,5D12.3,//)
      RETURN
  150 IFAIL =  2
      WRITE (6,1040) ABORT,P,Q,DP,DQ,ACC
 1040 FORMAT(' CF2 HAS FAILED TO CONVERGE AFTER ',F7.0,' ITERATIONS',/
     X,' P,Q,DP,DQ,ACCUR =  ',1P,4D17.7,D12.3,//)
      RETURN
  160 IFAIL =  3
      WRITE (6,1050) P,Q,ACC,DELL,LXTRA,M1
 1050 FORMAT(' FINAL Q.LE.ABS(P)*ACC*10**4 , P,Q,ACC = ',1P,3D12.3,4X,
     X' DELL,LXTRA,M1 = ',D12.3,2I5,/)
      RETURN
 1060 FORMAT(' XLMAX - XLMIN = DELL NOT AN INTEGER ',1P,3D20.10,/)
      END
c********************************************************************ZLB
      SUBROUTINE JWKB(XX,ETA1,XL,FJWKB,GJWKB,IEXP)
      implicit real*8 (a-h,o-z)
      double precision    XX,ETA1,XL,FJWKB,GJWKB,DZERO,zero
     . ,half,one,six,ten,rl35,aloge,x,eta,gh2,xll1,hll,hl,sl,
czlb     .  rl2,gh,phi,phi10,iexp
     .  rl2,gh,phi,phi10
C *** Computes JWKB approximations to Coulomb functions    for DL.GE. 0
C *** as modified by Biedenharn et al. Phys Rev 97 (1955) 542-554
C *** calls MAX1,SQRT,ALOG,EXO,ATAN2,FLOAT,INT        Barnett Feb 1981
      DATA   ZERO,HALF,ONE,SIX,TEN/ 0.0d0, 0.5d0, 1.0d0, 6.0d0, 10.0d0 /
      DATA  DZERO, RL35, ALOGE  /0.0D0, 35.0d0, 0.43429 45 d0 /
      X     = XX
      ETA   = ETA1
      GH2   = X*(ETA + ETA - X)
      XLL1  = dMAX1(XL*XL + XL,DZERO)
      IF(GH2 + XLL1 .LE. ZERO) RETURN
       HLL  = XLL1 + SIX/RL35
       HL   = dSQRT(HLL)
       SL   = ETA/HL + HL/X
       RL2  = ONE + ETA*ETA/HLL
       GH   = dSQRT(GH2 + HLL)/X
       PHI  = X*GH - HALF*( HL*dLOG((GH + SL)**2/RL2) - dLOG(GH) )
          IF(ETA .NE. ZERO) PHI = PHI - ETA*dATAN2(X*GH,X - ETA)
      PHI10 = -PHI*ALOGE
      IEXP  =  idINT(PHI10)
      IF(IEXP .GT. 70) GJWKB = TEN**(PHI10 - dble(IEXP))
      IF(IEXP .LE. 70) GJWKB = dEXP(-PHI)
      IF(IEXP .LE. 70) IEXP  = 0
      FJWKB = HALF/(GH*GJWKB)
      RETURN
      END
c********************************************************************ZLB
c P.C. Stancil 10-8-92
c Determines parameter array for spline
c From Numerical Recipes p.88 (1986)
c
      subroutine spline(x,y,n,yp1,ypn,y2)
      implicit real*8 (a-h,o-z)
c
      integer i,n,k,nmax
czlb      parameter (nmax=100)
      parameter (nmax=5000)
      double precision yp1,ypn,sig,p,qn,un
      double precision x(n),y(n),y2(n),u(nmax)
c
      if (yp1 .gt. 0.99d30) then
c      if (yp1 .gt. 0.99d60) then
         y2(1)=0.0d0
         u(1)=0.0d0
      else
         y2(1)=-0.5d0
         u(1)=(3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
c
      do 10 i=2,n-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+2.d0
         y2(i)=(sig-1.0d0)/p
         u(i)=(6.0d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     .         /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
 10   continue
c
c      if (ypn .gt. 0.99d30) then
      if (ypn .gt. 0.99d60) then
         qn=0.0d0
         un=0.0d0
      else
         qn=0.5d0
         un=(3.0d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
c
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0d0)
c
      do 20 k=n-1,1,-1
         y2(k)=y2(k)*y2(k+1)+u(k)
 20   continue
c
c      do i=1,n
c      write(11,999) x(i),y(i),y2(i)
c      end do
999   format(3g12.4)
      return
      end
c********************************************************************ZLB
c P.C. Stancil 10-8-92
c Evaluates a function using spline array from spline.f
c From Numerical Recipes p.89 (1986)
c
      subroutine splint(xa,ya,y2a,n,x,y)
      implicit real*8 (a-h,o-z)
c
      integer n,klo,khi,k
      double precision x,y,a,b,h
      double precision xa(n),ya(n),y2a(n)
c
      klo=1
      khi=n
c
 10   if ((khi-klo) .gt. 1) then
         k=(khi+klo)/2
         if (xa(k) .gt. x) then
            khi=k
         else
            klo=k
         endif
         goto 10
      endif
c
      h=xa(khi)-xa(klo)
      if (h .eq. 0.0d0) pause 'Bad XA input'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     .  ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*h**2/6.0d0
c
      return
      end
c********************************************************************ZLB
