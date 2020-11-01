      SUBROUTINE  DC3D0(ALPHA,X,Y,Z,DEPTH,DIP,POT1,POT2,POT3,POT4,      00010000
     *               UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET) 00020002
      IMPLICIT REAL*8 (A-H,O-Z)                                         00030000
      REAL*4   ALPHA,X,Y,Z,DEPTH,DIP,POT1,POT2,POT3,POT4,               00040000
     *         UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ             00050000
C                                                                       00060000
C********************************************************************   00070000
C*****                                                          *****   00080000
C*****    DISPLACEMENT AND STRAIN AT DEPTH                      *****   00090000
C*****    DUE TO BURIED POINT SOURCE IN A SEMIINFINITE MEDIUM   *****   00100000
C*****                         CODED BY  Y.OKADA ... SEP.1991   *****   00110002
C*****                         REVISED   Y.OKADA ... NOV.1991   *****   00120002
C*****                                                          *****   00130000
C********************************************************************   00140000
C                                                                       00150000
C***** INPUT                                                            00160000
C*****   ALPHA : MEDIUM CONSTANT  (LAMBDA+MYU)/(LAMBDA+2*MYU)           00170000
C*****   X,Y,Z : COORDINATE OF OBSERVING POINT                          00180000
C*****   DEPTH : SOURCE DEPTH                                           00190000
C*****   DIP   : DIP-ANGLE (DEGREE)                                     00200000
C*****   POT1-POT4 : STRIKE-, DIP-, TENSILE- AND INFLATE-POTENCY        00210000
C*****       POTENCY=(  MOMENT OF DOUBLE-COUPLE  )/MYU     FOR POT1,2   00220000
C*****       POTENCY=(INTENSITY OF ISOTROPIC PART)/LAMBDA  FOR POT3     00230000
C*****       POTENCY=(INTENSITY OF LINEAR DIPOLE )/MYU     FOR POT4     00240000
C                                                                       00250000
C***** OUTPUT                                                           00260000
C*****   UX, UY, UZ  : DISPLACEMENT ( UNIT=(UNIT OF POTENCY) /          00270000
C*****               :                     (UNIT OF X,Y,Z,DEPTH)**2  )  00280000
C*****   UXX,UYX,UZX : X-DERIVATIVE ( UNIT= UNIT OF POTENCY) /          00290000
C*****   UXY,UYY,UZY : Y-DERIVATIVE        (UNIT OF X,Y,Z,DEPTH)**3  )  00300000
C*****   UXZ,UYZ,UZZ : Z-DERIVATIVE                                     00310000
C*****   IRET        : RETURN CODE  ( =0....NORMAL,   =1....SINGULAR )  00320002
C                                                                       00330000
      COMMON /C1/DUMMY(8),R,dumm(15)                                    00340000
      DIMENSION  U(12),DUA(12),DUB(12),DUC(12)                          00350000
      DATA  F0/0.D0/                                                    00360000
C-----                                                                  00370000
      IF(Z.GT.0.) WRITE(*,'(''0** POSITIVE Z WAS GIVEN IN SUB-DC3D0'')')00380000
      DO 111 I=1,12                                                     00390000
        U(I)=F0                                                         00400000
        DUA(I)=F0                                                       00410000
        DUB(I)=F0                                                       00420000
        DUC(I)=F0                                                       00430000
  111 CONTINUE                                                          00440000
      AALPHA=ALPHA                                                      00450000
      DDIP=DIP                                                          00460000
      CALL DCCON0(AALPHA,DDIP)                                          00470000
C======================================                                 00480000
C=====  REAL-SOURCE CONTRIBUTION  =====                                 00490000
C======================================                                 00500000
      XX=X                                                              00510000
      YY=Y                                                              00520000
      ZZ=Z                                                              00530000
      DD=DEPTH+Z                                                        00540000
      CALL DCCON1(XX,YY,DD)                                             00550000
      IF(R.EQ.F0) GO TO 99                                              00560000
      PP1=POT1                                                          00570000
      PP2=POT2                                                          00580000
      PP3=POT3                                                          00590000
      PP4=POT4                                                          00600000
      CALL UA0(XX,YY,DD,PP1,PP2,PP3,PP4,DUA)                            00610000
C-----                                                                  00620000
      DO 222 I=1,12                                                     00630000
        IF(I.LT.10) U(I)=U(I)-DUA(I)                                    00640000
        IF(I.GE.10) U(I)=U(I)+DUA(I)                                    00650000
  222 CONTINUE                                                          00660000
C=======================================                                00670000
C=====  IMAGE-SOURCE CONTRIBUTION  =====                                00680000
C=======================================                                00690000
      DD=DEPTH-Z                                                        00700000
      CALL DCCON1(XX,YY,DD)                                             00710000
      CALL UA0(XX,YY,DD,PP1,PP2,PP3,PP4,DUA)                            00720000
      CALL UB0(XX,YY,DD,ZZ,PP1,PP2,PP3,PP4,DUB)                         00730000
      CALL UC0(XX,YY,DD,ZZ,PP1,PP2,PP3,PP4,DUC)                         00740000
C-----                                                                  00750000
      DO 333 I=1,12                                                     00760000
        DU=DUA(I)+DUB(I)+ZZ*DUC(I)                                      00770000
        IF(I.GE.10) DU=DU+DUC(I-9)                                      00780000
        U(I)=U(I)+DU                                                    00790000
  333 CONTINUE                                                          00800000
C=====                                                                  00810000
      UX=U(1)                                                           00820000
      UY=U(2)                                                           00830000
      UZ=U(3)                                                           00840000
      UXX=U(4)                                                          00850000
      UYX=U(5)                                                          00860000
      UZX=U(6)                                                          00870000
      UXY=U(7)                                                          00880000
      UYY=U(8)                                                          00890000
      UZY=U(9)                                                          00900000
      UXZ=U(10)                                                         00910000
      UYZ=U(11)                                                         00920000
      UZZ=U(12)                                                         00930000
      IRET=0                                                            00940002
      RETURN                                                            00950000
C=======================================                                00960000
C=====  IN CASE OF SINGULAR (R=0)  =====                                00970000
C=======================================                                00980000
   99 UX=F0                                                             00990000
      UY=F0                                                             01000000
      UZ=F0                                                             01010000
      UXX=F0                                                            01020000
      UYX=F0                                                            01030000
      UZX=F0                                                            01040000
      UXY=F0                                                            01050000
      UYY=F0                                                            01060000
      UZY=F0                                                            01070000
      UXZ=F0                                                            01080000
      UYZ=F0                                                            01090000
      UZZ=F0                                                            01100000
      IRET=1                                                            01110002
      RETURN                                                            01120000
      END                                                               01130000
      SUBROUTINE  UA0(X,Y,D,POT1,POT2,POT3,POT4,U)                      01140000
      IMPLICIT REAL*8 (A-H,O-Z)                                         01150000
      DIMENSION U(12),DU(12)                                            01160000
C                                                                       01170000
C********************************************************************   01180000
C*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-A)             *****   01190000
C*****    DUE TO BURIED POINT SOURCE IN A SEMIINFINITE MEDIUM   *****   01200000
C********************************************************************   01210000
C                                                                       01220000
C***** INPUT                                                            01230000
C*****   X,Y,D : STATION COORDINATES IN FAULT SYSTEM                    01240000
C*****   POT1-POT4 : STRIKE-, DIP-, TENSILE- AND INFLATE-POTENCY        01250000
C***** OUTPUT                                                           01260000
C*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     01270000
C                                                                       01280000
      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D  01290000
      COMMON /C1/P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,A3,A5,B3,C3,     01300000
     *           UY,VY,WY,UZ,VZ,WZ                                      01310000
      DATA F0,F1,F3/0.D0,1.D0,3.D0/                                     01320000
      DATA PI2/6.283185307179586D0/                                     01330000
C-----                                                                  01340000
      DO 111  I=1,12                                                    01350000
  111 U(I)=F0                                                           01360000
C======================================                                 01370000
C=====  STRIKE-SLIP CONTRIBUTION  =====                                 01380000
C======================================                                 01390000
      IF(POT1.NE.F0) THEN                                               01400000
        DU( 1)= ALP1*Q/R3    +ALP2*X2*QR                                01410000
        DU( 2)= ALP1*X/R3*SD +ALP2*XY*QR                                01420000
        DU( 3)=-ALP1*X/R3*CD +ALP2*X*D*QR                               01430000
        DU( 4)= X*QR*(-ALP1 +ALP2*(F1+A5) )                             01440000
        DU( 5)= ALP1*A3/R3*SD +ALP2*Y*QR*A5                             01450000
        DU( 6)=-ALP1*A3/R3*CD +ALP2*D*QR*A5                             01460000
        DU( 7)= ALP1*(SD/R3-Y*QR) +ALP2*F3*X2/R5*UY                     01470000
        DU( 8)= F3*X/R5*(-ALP1*Y*SD +ALP2*(Y*UY+Q) )                    01480000
        DU( 9)= F3*X/R5*( ALP1*Y*CD +ALP2*D*UY )                        01490000
        DU(10)= ALP1*(CD/R3+D*QR) +ALP2*F3*X2/R5*UZ                     01500000
        DU(11)= F3*X/R5*( ALP1*D*SD +ALP2*Y*UZ )                        01510000
        DU(12)= F3*X/R5*(-ALP1*D*CD +ALP2*(D*UZ-Q) )                    01520000
        DO 222 I=1,12                                                   01530000
  222   U(I)=U(I)+POT1/PI2*DU(I)                                        01540000
      ENDIF                                                             01550000
C===================================                                    01560000
C=====  DIP-SLIP CONTRIBUTION  =====                                    01570000
C===================================                                    01580000
      IF(POT2.NE.F0) THEN                                               01590000
        DU( 1)=            ALP2*X*P*QR                                  01600000
        DU( 2)= ALP1*S/R3 +ALP2*Y*P*QR                                  01610000
        DU( 3)=-ALP1*T/R3 +ALP2*D*P*QR                                  01620000
        DU( 4)=                 ALP2*P*QR*A5                            01630000
        DU( 5)=-ALP1*F3*X*S/R5 -ALP2*Y*P*QRX                            01640000
        DU( 6)= ALP1*F3*X*T/R5 -ALP2*D*P*QRX                            01650000
        DU( 7)=                          ALP2*F3*X/R5*VY                01660000
        DU( 8)= ALP1*(S2D/R3-F3*Y*S/R5) +ALP2*(F3*Y/R5*VY+P*QR)         01670000
        DU( 9)=-ALP1*(C2D/R3-F3*Y*T/R5) +ALP2*F3*D/R5*VY                01680000
        DU(10)=                          ALP2*F3*X/R5*VZ                01690000
        DU(11)= ALP1*(C2D/R3+F3*D*S/R5) +ALP2*F3*Y/R5*VZ                01700000
        DU(12)= ALP1*(S2D/R3-F3*D*T/R5) +ALP2*(F3*D/R5*VZ-P*QR)         01710000
        DO 333 I=1,12                                                   01720000
  333   U(I)=U(I)+POT2/PI2*DU(I)                                        01730000
      ENDIF                                                             01740000
C========================================                               01750000
C=====  TENSILE-FAULT CONTRIBUTION  =====                               01760000
C========================================                               01770000
      IF(POT3.NE.F0) THEN                                               01780000
        DU( 1)= ALP1*X/R3 -ALP2*X*Q*QR                                  01790000
        DU( 2)= ALP1*T/R3 -ALP2*Y*Q*QR                                  01800000
        DU( 3)= ALP1*S/R3 -ALP2*D*Q*QR                                  01810000
        DU( 4)= ALP1*A3/R3     -ALP2*Q*QR*A5                            01820000
        DU( 5)=-ALP1*F3*X*T/R5 +ALP2*Y*Q*QRX                            01830000
        DU( 6)=-ALP1*F3*X*S/R5 +ALP2*D*Q*QRX                            01840000
        DU( 7)=-ALP1*F3*XY/R5           -ALP2*X*QR*WY                   01850000
        DU( 8)= ALP1*(C2D/R3-F3*Y*T/R5) -ALP2*(Y*WY+Q)*QR               01860000
        DU( 9)= ALP1*(S2D/R3-F3*Y*S/R5) -ALP2*D*QR*WY                   01870000
        DU(10)= ALP1*F3*X*D/R5          -ALP2*X*QR*WZ                   01880000
        DU(11)=-ALP1*(S2D/R3-F3*D*T/R5) -ALP2*Y*QR*WZ                   01890000
        DU(12)= ALP1*(C2D/R3+F3*D*S/R5) -ALP2*(D*WZ-Q)*QR               01900000
        DO 444 I=1,12                                                   01910000
  444   U(I)=U(I)+POT3/PI2*DU(I)                                        01920000
      ENDIF                                                             01930000
C=========================================                              01940000
C=====  INFLATE SOURCE CONTRIBUTION  =====                              01950000
C=========================================                              01960000
      IF(POT4.NE.F0) THEN                                               01970000
        DU( 1)=-ALP1*X/R3                                               01980000
        DU( 2)=-ALP1*Y/R3                                               01990000
        DU( 3)=-ALP1*D/R3                                               02000000
        DU( 4)=-ALP1*A3/R3                                              02010000
        DU( 5)= ALP1*F3*XY/R5                                           02020000
        DU( 6)= ALP1*F3*X*D/R5                                          02030000
        DU( 7)= DU(5)                                                   02040000
        DU( 8)=-ALP1*B3/R3                                              02050000
        DU( 9)= ALP1*F3*Y*D/R5                                          02060000
        DU(10)=-DU(6)                                                   02070000
        DU(11)=-DU(9)                                                   02080000
        DU(12)= ALP1*C3/R3                                              02090000
        DO 555 I=1,12                                                   02100000
  555   U(I)=U(I)+POT4/PI2*DU(I)                                        02110000
      ENDIF                                                             02120000
      RETURN                                                            02130000
      END                                                               02140000
      SUBROUTINE  UB0(X,Y,D,Z,POT1,POT2,POT3,POT4,U)                    02150000
      IMPLICIT REAL*8 (A-H,O-Z)                                         02160000
      DIMENSION U(12),DU(12)                                            02170000
C                                                                       02180000
C********************************************************************   02190000
C*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-B)             *****   02200000
C*****    DUE TO BURIED POINT SOURCE IN A SEMIINFINITE MEDIUM   *****   02210000
C********************************************************************   02220000
C                                                                       02230000
C***** INPUT                                                            02240000
C*****   X,Y,D,Z : STATION COORDINATES IN FAULT SYSTEM                  02250000
C*****   POT1-POT4 : STRIKE-, DIP-, TENSILE- AND INFLATE-POTENCY        02260000
C***** OUTPUT                                                           02270000
C*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     02280000
C                                                                       02290000
      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D  02300000
      COMMON /C1/P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,A3,A5,B3,C3,     02310000
     *           UY,VY,WY,UZ,VZ,WZ                                      02320000
      DATA F0,F1,F2,F3,F4,F5,F8,F9                                      02330000
     *        /0.D0,1.D0,2.D0,3.D0,4.D0,5.D0,8.D0,9.D0/                 02340000
      DATA PI2/6.283185307179586D0/                                     02350000
C-----                                                                  02360000
      C=D+Z                                                             02370000
      RD=R+D                                                            02380000
      D12=F1/(R*RD*RD)                                                  02390000
      D32=D12*(F2*R+D)/R2                                               02400000
      D33=D12*(F3*R+D)/(R2*RD)                                          02410000
      D53=D12*(F8*R2+F9*R*D+F3*D2)/(R2*R2*RD)                           02420000
      D54=D12*(F5*R2+F4*R*D+D2)/R3*D12                                  02430000
C-----                                                                  02440000
      FI1= Y*(D12-X2*D33)                                               02450000
      FI2= X*(D12-Y2*D33)                                               02460000
      FI3= X/R3-FI2                                                     02470000
      FI4=-XY*D32                                                       02480000
      FI5= F1/(R*RD)-X2*D32                                             02490000
      FJ1=-F3*XY*(D33-X2*D54)                                           02500000
      FJ2= F1/R3-F3*D12+F3*X2*Y2*D54                                    02510000
      FJ3= A3/R3-FJ2                                                    02520000
      FJ4=-F3*XY/R5-FJ1                                                 02530000
      FK1=-Y*(D32-X2*D53)                                               02540000
      FK2=-X*(D32-Y2*D53)                                               02550000
      FK3=-F3*X*D/R5-FK2                                                02560000
C-----                                                                  02570000
      DO 111  I=1,12                                                    02580000
  111 U(I)=F0                                                           02590000
C======================================                                 02600000
C=====  STRIKE-SLIP CONTRIBUTION  =====                                 02610000
C======================================                                 02620000
      IF(POT1.NE.F0) THEN                                               02630000
        DU( 1)=-X2*QR  -ALP3*FI1*SD                                     02640000
        DU( 2)=-XY*QR  -ALP3*FI2*SD                                     02650000
        DU( 3)=-C*X*QR -ALP3*FI4*SD                                     02660000
        DU( 4)=-X*QR*(F1+A5) -ALP3*FJ1*SD                               02670000
        DU( 5)=-Y*QR*A5      -ALP3*FJ2*SD                               02680000
        DU( 6)=-C*QR*A5      -ALP3*FK1*SD                               02690000
        DU( 7)=-F3*X2/R5*UY      -ALP3*FJ2*SD                           02700000
        DU( 8)=-F3*XY/R5*UY-X*QR -ALP3*FJ4*SD                           02710000
        DU( 9)=-F3*C*X/R5*UY     -ALP3*FK2*SD                           02720000
        DU(10)=-F3*X2/R5*UZ  +ALP3*FK1*SD                               02730000
        DU(11)=-F3*XY/R5*UZ  +ALP3*FK2*SD                               02740000
        DU(12)= F3*X/R5*(-C*UZ +ALP3*Y*SD)                              02750000
        DO 222 I=1,12                                                   02760000
  222   U(I)=U(I)+POT1/PI2*DU(I)                                        02770000
      ENDIF                                                             02780000
C===================================                                    02790000
C=====  DIP-SLIP CONTRIBUTION  =====                                    02800000
C===================================                                    02810000
      IF(POT2.NE.F0) THEN                                               02820000
        DU( 1)=-X*P*QR +ALP3*FI3*SDCD                                   02830000
        DU( 2)=-Y*P*QR +ALP3*FI1*SDCD                                   02840000
        DU( 3)=-C*P*QR +ALP3*FI5*SDCD                                   02850000
        DU( 4)=-P*QR*A5 +ALP3*FJ3*SDCD                                  02860000
        DU( 5)= Y*P*QRX +ALP3*FJ1*SDCD                                  02870000
        DU( 6)= C*P*QRX +ALP3*FK3*SDCD                                  02880000
        DU( 7)=-F3*X/R5*VY      +ALP3*FJ1*SDCD                          02890000
        DU( 8)=-F3*Y/R5*VY-P*QR +ALP3*FJ2*SDCD                          02900000
        DU( 9)=-F3*C/R5*VY      +ALP3*FK1*SDCD                          02910000
        DU(10)=-F3*X/R5*VZ -ALP3*FK3*SDCD                               02920000
        DU(11)=-F3*Y/R5*VZ -ALP3*FK1*SDCD                               02930000
        DU(12)=-F3*C/R5*VZ +ALP3*A3/R3*SDCD                             02940000
        DO 333 I=1,12                                                   02950000
  333   U(I)=U(I)+POT2/PI2*DU(I)                                        02960000
      ENDIF                                                             02970000
C========================================                               02980000
C=====  TENSILE-FAULT CONTRIBUTION  =====                               02990000
C========================================                               03000000
      IF(POT3.NE.F0) THEN                                               03010000
        DU( 1)= X*Q*QR -ALP3*FI3*SDSD                                   03020000
        DU( 2)= Y*Q*QR -ALP3*FI1*SDSD                                   03030000
        DU( 3)= C*Q*QR -ALP3*FI5*SDSD                                   03040000
        DU( 4)= Q*QR*A5 -ALP3*FJ3*SDSD                                  03050000
        DU( 5)=-Y*Q*QRX -ALP3*FJ1*SDSD                                  03060000
        DU( 6)=-C*Q*QRX -ALP3*FK3*SDSD                                  03070000
        DU( 7)= X*QR*WY     -ALP3*FJ1*SDSD                              03080000
        DU( 8)= QR*(Y*WY+Q) -ALP3*FJ2*SDSD                              03090000
        DU( 9)= C*QR*WY     -ALP3*FK1*SDSD                              03100000
        DU(10)= X*QR*WZ +ALP3*FK3*SDSD                                  03110000
        DU(11)= Y*QR*WZ +ALP3*FK1*SDSD                                  03120000
        DU(12)= C*QR*WZ -ALP3*A3/R3*SDSD                                03130000
        DO 444 I=1,12                                                   03140000
  444   U(I)=U(I)+POT3/PI2*DU(I)                                        03150000
      ENDIF                                                             03160000
C=========================================                              03170000
C=====  INFLATE SOURCE CONTRIBUTION  =====                              03180000
C=========================================                              03190000
      IF(POT4.NE.F0) THEN                                               03200000
        DU( 1)= ALP3*X/R3                                               03210000
        DU( 2)= ALP3*Y/R3                                               03220000
        DU( 3)= ALP3*D/R3                                               03230000
        DU( 4)= ALP3*A3/R3                                              03240000
        DU( 5)=-ALP3*F3*XY/R5                                           03250000
        DU( 6)=-ALP3*F3*X*D/R5                                          03260000
        DU( 7)= DU(5)                                                   03270000
        DU( 8)= ALP3*B3/R3                                              03280000
        DU( 9)=-ALP3*F3*Y*D/R5                                          03290000
        DU(10)=-DU(6)                                                   03300000
        DU(11)=-DU(9)                                                   03310000
        DU(12)=-ALP3*C3/R3                                              03320000
        DO 555 I=1,12                                                   03330000
  555   U(I)=U(I)+POT4/PI2*DU(I)                                        03340000
      ENDIF                                                             03350000
      RETURN                                                            03360000
      END                                                               03370000
      SUBROUTINE  UC0(X,Y,D,Z,POT1,POT2,POT3,POT4,U)                    03380000
      IMPLICIT REAL*8 (A-H,O-Z)                                         03390000
      DIMENSION U(12),DU(12)                                            03400000
C                                                                       03410000
C********************************************************************   03420000
C*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-B)             *****   03430000
C*****    DUE TO BURIED POINT SOURCE IN A SEMIINFINITE MEDIUM   *****   03440000
C********************************************************************   03450000
C                                                                       03460000
C***** INPUT                                                            03470000
C*****   X,Y,D,Z : STATION COORDINATES IN FAULT SYSTEM                  03480000
C*****   POT1-POT4 : STRIKE-, DIP-, TENSILE- AND INFLATE-POTENCY        03490000
C***** OUTPUT                                                           03500000
C*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     03510000
C                                                                       03520000
      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D  03530000
      COMMON /C1/P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,A3,A5,B3,C3,um(6)03540000
      DATA F0,F1,F2,F3,F5,F7,F10,F15                                    03550000
     *        /0.D0,1.D0,2.D0,3.D0,5.D0,7.D0,10.D0,15.D0/               03560000
      DATA PI2/6.283185307179586D0/                                     03570000
C-----                                                                  03580000
      C=D+Z                                                             03590000
      Q2=Q*Q                                                            03600000
      R7=R5*R2                                                          03610000
      A7=F1-F7*X2/R2                                                    03620000
      B5=F1-F5*Y2/R2                                                    03630000
      B7=F1-F7*Y2/R2                                                    03640000
      C5=F1-F5*D2/R2                                                    03650000
      C7=F1-F7*D2/R2                                                    03660000
      D7=F2-F7*Q2/R2                                                    03670000
      QR5=F5*Q/R2                                                       03680000
      QR7=F7*Q/R2                                                       03690000
      DR5=F5*D/R2                                                       03700000
C-----                                                                  03710000
      DO 111  I=1,12                                                    03720000
  111 U(I)=F0                                                           03730000
C======================================                                 03740000
C=====  STRIKE-SLIP CONTRIBUTION  =====                                 03750000
C======================================                                 03760000
      IF(POT1.NE.F0) THEN                                               03770000
        DU( 1)=-ALP4*A3/R3*CD  +ALP5*C*QR*A5                            03780000
        DU( 2)= F3*X/R5*( ALP4*Y*CD +ALP5*C*(SD-Y*QR5) )                03790000
        DU( 3)= F3*X/R5*(-ALP4*Y*SD +ALP5*C*(CD+D*QR5) )                03800000
        DU( 4)= ALP4*F3*X/R5*(F2+A5)*CD   -ALP5*C*QRX*(F2+A7)           03810000
        DU( 5)= F3/R5*( ALP4*Y*A5*CD +ALP5*C*(A5*SD-Y*QR5*A7) )         03820000
        DU( 6)= F3/R5*(-ALP4*Y*A5*SD +ALP5*C*(A5*CD+D*QR5*A7) )         03830000
        DU( 7)= DU(5)                                                   03840000
        DU( 8)= F3*X/R5*( ALP4*B5*CD -ALP5*F5*C/R2*(F2*Y*SD+Q*B7) )     03850000
        DU( 9)= F3*X/R5*(-ALP4*B5*SD +ALP5*F5*C/R2*(D*B7*SD-Y*C7*CD) )  03860000
        DU(10)= F3/R5*   (-ALP4*D*A5*CD +ALP5*C*(A5*CD+D*QR5*A7) )      03870000
        DU(11)= F15*X/R7*( ALP4*Y*D*CD  +ALP5*C*(D*B7*SD-Y*C7*CD) )     03880000
        DU(12)= F15*X/R7*(-ALP4*Y*D*SD  +ALP5*C*(F2*D*CD-Q*C7) )        03890000
        DO 222 I=1,12                                                   03900000
  222   U(I)=U(I)+POT1/PI2*DU(I)                                        03910000
      ENDIF                                                             03920000
C===================================                                    03930000
C=====  DIP-SLIP CONTRIBUTION  =====                                    03940000
C===================================                                    03950000
      IF(POT2.NE.F0) THEN                                               03960000
        DU( 1)= ALP4*F3*X*T/R5          -ALP5*C*P*QRX                   03970000
        DU( 2)=-ALP4/R3*(C2D-F3*Y*T/R2) +ALP5*F3*C/R5*(S-Y*P*QR5)       03980000
        DU( 3)=-ALP4*A3/R3*SDCD         +ALP5*F3*C/R5*(T+D*P*QR5)       03990000
        DU( 4)= ALP4*F3*T/R5*A5              -ALP5*F5*C*P*QR/R2*A7      04000000
        DU( 5)= F3*X/R5*(ALP4*(C2D-F5*Y*T/R2)-ALP5*F5*C/R2*(S-Y*P*QR7)) 04010000
        DU( 6)= F3*X/R5*(ALP4*(F2+A5)*SDCD   -ALP5*F5*C/R2*(T+D*P*QR7)) 04020000
        DU( 7)= DU(5)                                                   04030000
        DU( 8)= F3/R5*(ALP4*(F2*Y*C2D+T*B5)                             04040000
     *                               +ALP5*C*(S2D-F10*Y*S/R2-P*QR5*B7)) 04050000
        DU( 9)= F3/R5*(ALP4*Y*A5*SDCD-ALP5*C*((F3+A5)*C2D+Y*P*DR5*QR7)) 04060000
        DU(10)= F3*X/R5*(-ALP4*(S2D-T*DR5) -ALP5*F5*C/R2*(T+D*P*QR7))   04070000
        DU(11)= F3/R5*(-ALP4*(D*B5*C2D+Y*C5*S2D)                        04080000
     *                                -ALP5*C*((F3+A5)*C2D+Y*P*DR5*QR7))04090000
        DU(12)= F3/R5*(-ALP4*D*A5*SDCD-ALP5*C*(S2D-F10*D*T/R2+P*QR5*C7))04100000
        DO 333 I=1,12                                                   04110000
  333   U(I)=U(I)+POT2/PI2*DU(I)                                        04120000
      ENDIF                                                             04130000
C========================================                               04140000
C=====  TENSILE-FAULT CONTRIBUTION  =====                               04150000
C========================================                               04160000
      IF(POT3.NE.F0) THEN                                               04170000
        DU( 1)= F3*X/R5*(-ALP4*S +ALP5*(C*Q*QR5-Z))                     04180000
        DU( 2)= ALP4/R3*(S2D-F3*Y*S/R2)+ALP5*F3/R5*(C*(T-Y+Y*Q*QR5)-Y*Z)04190000
        DU( 3)=-ALP4/R3*(F1-A3*SDSD)   -ALP5*F3/R5*(C*(S-D+D*Q*QR5)-D*Z)04200000
        DU( 4)=-ALP4*F3*S/R5*A5 +ALP5*(C*QR*QR5*A7-F3*Z/R5*A5)          04210000
        DU( 5)= F3*X/R5*(-ALP4*(S2D-F5*Y*S/R2)                          04220000
     *                               -ALP5*F5/R2*(C*(T-Y+Y*Q*QR7)-Y*Z)) 04230000
        DU( 6)= F3*X/R5*( ALP4*(F1-(F2+A5)*SDSD)                        04240000
     *                               +ALP5*F5/R2*(C*(S-D+D*Q*QR7)-D*Z)) 04250000
        DU( 7)= DU(5)                                                   04260000
        DU( 8)= F3/R5*(-ALP4*(F2*Y*S2D+S*B5)                            04270000
     *                -ALP5*(C*(F2*SDSD+F10*Y*(T-Y)/R2-Q*QR5*B7)+Z*B5)) 04280000
        DU( 9)= F3/R5*( ALP4*Y*(F1-A5*SDSD)                             04290000
     *                +ALP5*(C*(F3+A5)*S2D-Y*DR5*(C*D7+Z)))             04300000
        DU(10)= F3*X/R5*(-ALP4*(C2D+S*DR5)                              04310000
     *               +ALP5*(F5*C/R2*(S-D+D*Q*QR7)-F1-Z*DR5))            04320000
        DU(11)= F3/R5*( ALP4*(D*B5*S2D-Y*C5*C2D)                        04330000
     *               +ALP5*(C*((F3+A5)*S2D-Y*DR5*D7)-Y*(F1+Z*DR5)))     04340000
        DU(12)= F3/R5*(-ALP4*D*(F1-A5*SDSD)                             04350000
     *               -ALP5*(C*(C2D+F10*D*(S-D)/R2-Q*QR5*C7)+Z*(F1+C5))) 04360000
        DO 444 I=1,12                                                   04370000
  444   U(I)=U(I)+POT3/PI2*DU(I)                                        04380000
      ENDIF                                                             04390000
C=========================================                              04400000
C=====  INFLATE SOURCE CONTRIBUTION  =====                              04410000
C=========================================                              04420000
      IF(POT4.NE.F0) THEN                                               04430000
        DU( 1)= ALP4*F3*X*D/R5                                          04440000
        DU( 2)= ALP4*F3*Y*D/R5                                          04450000
        DU( 3)= ALP4*C3/R3                                              04460000
        DU( 4)= ALP4*F3*D/R5*A5                                         04470000
        DU( 5)=-ALP4*F15*XY*D/R7                                        04480000
        DU( 6)=-ALP4*F3*X/R5*C5                                         04490000
        DU( 7)= DU(5)                                                   04500000
        DU( 8)= ALP4*F3*D/R5*B5                                         04510000
        DU( 9)=-ALP4*F3*Y/R5*C5                                         04520000
        DU(10)= DU(6)                                                   04530000
        DU(11)= DU(9)                                                   04540000
        DU(12)= ALP4*F3*D/R5*(F2+C5)                                    04550000
        DO 555 I=1,12                                                   04560000
  555   U(I)=U(I)+POT4/PI2*DU(I)                                        04570000
      ENDIF                                                             04580000
      RETURN                                                            04590000
      END                                                               04600000
