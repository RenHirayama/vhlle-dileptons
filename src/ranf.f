c*********************************************************************72
c*********************************************************************72
c*** RANDOM NUMBER GENERATOR AND OTHER NUMERICAL TOOLS ***************72
c*********************************************************************72
c*********************************************************************72
c*********************************************************************72
c*********************************************************************72

      REAL*8 FUNCTION ranff(idum)

c*********************************************************************72
c
c  RANFF returns a unit pseudorandom R8.
c
c  Discussion:
c
c    This routine implements the recursion
c
c      seed = 16807 * seed mod ( 2**31 - 1 )
c      ranff = seed / ( 2**31 - 1 )
c
c    The integer arithmetic never requires more than 32 bits,
c    including a sign bit.
c
c    If the initial seed is 12345, then the first three computations are
c
c      Input     Output      R8_UNIFORM_01
c      SEED      SEED
c
c         12345   207482415  0.096616
c     207482415  1790989824  0.833995
c    1790989824  2035175616  0.947702
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 August 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Second Edition,
c    Springer, 1987,
c    ISBN: 0387964673,
c    LC: QA76.9.C65.B73.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, December 1986, pages 362-376.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley, 1998,
c    ISBN: 0471134031,
c    LC: T57.62.H37.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, Number 2, 1969, pages 136-143.
c
c  Parameters:
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R8_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.

c      implicit none

c      integer i4_huge
c      parameter ( i4_huge = 2147483647 )
c      integer k
c      integer seedin

c      write(0,*)'seedin',seedin ! Debug only
      
c      k = seedin / 127773

c      seedin = 16807 * ( seedin - k * 127773 ) - k * 2836

c      if ( seedin .lt. 0 ) then
c        seedin = seedin + i4_huge
c      end if

c      write(0,*)'new seedin',seedin ! Debug only
       

c..  Although SEED can be represented exactly as a 32 bit integer,
c..  it generally cannot be represented exactly as a 32 bit real number!

c      ranff = dble ( seedin ) * 4.656612875D-10
c      write(0,*)'***ranff***',ranff ! Debug only

c      return
c      end

      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
 11     continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ranff=min(AM*iy,RNMX)
      return
      END


c-----|----------------------------------------------------------------|X
 
       real*8 FUNCTION bessi1(x)
       implicit none
       double precision x

       REAL*8 ax
       DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,
     * q8,q9,y

       SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
       DATA p1,p2,p3,p4,p5,p6,p7/0.5d0,0.87890594d0,0.51498869d0,
     * 0.15084934d0,0.2658733d-1,0.301532d-2,0.32411d-3/
       DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,-0.3988024d-1,
     * -0.362018d-2,0.163801d-2,-0.1031555d-1,0.2282967d-1,
     * -0.2895312d-1,0.1787654d-1,-0.420059d-2/
       if (abs(x).lt.3.75d0) then
        y=(x/3.75d0)**2
        bessi1=x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
       else
        ax=abs(x)
        y=3.75d0/ax
        bessi1=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+
     * y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
       if(x.lt.0.0d0)bessi1=-bessi1
       endif
       return
       END

c-----|----------------------------------------------------------------|X

       real*8 function bessk1(x)
       implicit none
       real*8 x
c     uses bessi1
c     returns the modified Bessel function Ki(x) for positive real x
       real*8 bessi1
       double precision p1,p2,p3,p4,p5,p6,p7,q1,
     &  q2,q3,q4,q5,q6,q7,y
       save p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
       data p1,p2,p3,p4,p5,p6,p7/1.0d0,0.15443144d0,-0.67278579d0,
     &  -0.18156897d0,-0.1919402d-1,-0.110404d-2,-0.4686d-4/
       data q1,q2,q3,q4,q5,q6,q7/1.25331414d0,0.23498618d0,
     &  -0.3655620d-1,
     &  0.1504268d-1,-0.780353d-2,0.325614d-2,-0.68245d-3/
       if (x.le.2.0d0) then
          y=x*x/4.0d0
          bessk1=(log(x/2.0d0)*bessi1(x))+(1.0d0/x)*(p1+y*(p2+
     &      y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
       else
          y=2.0d0/x
          bessk1=(exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+
     &      y*(q6+y*q7))))))
       endif
       return
       end

c-----|----------------------------------------------------------------|X

        DOUBLE PRECISION FUNCTION BESSK(N,X)
C
C       =============================================================
C       Purpose: This program computes modified Bessel functions 
C                In(x) and Kn(x), and their derivatives using
C                subroutine IKNA
C       Input:   x --- Argument of In(x) and Kn(x) ( x � 0 )
C                n --- Order of In(x) and Kn(x)
C                      ( n = 0,1,���, n � 250 )
C       Output:  BI(n) --- In(x)
C                DI(n) --- In'(x)
C                BK(n) --- Kn(x)
C                DK(n) --- Kn'(x)
C       Example: Nmax = 5,    x = 10.0
C
C     n      In(x)          In'(x)         Kn(x)         Kn'(x)
C    ---------------------------------------------------------------
C     0   .2815717D+04   .2670988D+04   .1778006D-04  -.1864877D-04
C     1   .2670988D+04   .2548618D+04   .1864877D-04  -.1964494D-04
C     2   .2281519D+04   .2214685D+04   .2150982D-04  -.2295074D-04
C     3   .1758381D+04   .1754005D+04   .2725270D-04  -.2968563D-04
C     4   .1226491D+04   .1267785D+04   .3786144D-04  -.4239728D-04
C     5   .7771883D+03   .8378964D+03   .5754185D-04  -.6663236D-04
C       =============================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BI(0:250),DI(0:250),BK(0:250),DK(0:250)
c        WRITE(*,*)'  Please enter n, x '
c        READ(*,*)N,X
c        WRITE(*,25)N,X
c        WRITE(*,*)
c        IF (N.LE.10) THEN
           NS=1
c        ELSE
c           WRITE(*,*)'  Please enter order step Ns'
c           READ(*,*)NS
c        ENDIF
        CALL IKNA(N,X,NM,BI,DI,BK,DK)
c        WRITE(*,*)'  n      In(x)          In''(x) ',
c     &            '        Kn(x)         Kn''(x) '
c        WRITE(*,*)' -------------------------------',
c     &            '--------------------------------'
c        DO 10 K=0,NM,NS
c           WRITE(*,20)K,BI(K),DI(K),BK(K),DK(K)
c10      CONTINUE

        BESSK=BK(N)

20      FORMAT(1X,I3,4D15.7)
25      FORMAT(3X,6HNmax =,I3,',    ',3Hx =,F5.1)
        END


        SUBROUTINE IKNA(N,X,NM,BI,DI,BK,DK)
C
C       ========================================================
C       Purpose: Compute modified Bessel functions In(x) and
C                Kn(x), and their derivatives
C       Input:   x --- Argument of In(x) and Kn(x) ( x � 0 )
C                n --- Order of In(x) and Kn(x)
C       Output:  BI(n) --- In(x)
C                DI(n) --- In'(x)
C                BK(n) --- Kn(x)
C                DK(n) --- Kn'(x)
C                NM --- Highest order computed
C       Routines called:
C            (1) IK01A for computing I0(x),I1(x),K0(x) & K1(x)
C            (2) MSTA1 and MSTA2 for computing the starting 
C                point for backward recurrence
C       ========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BI(0:N),DI(0:N),BK(0:N),DK(0:N)
        NM=N
        IF (X.LE.1.0D-100) THEN
           DO 10 K=0,N
              BI(K)=0.0D0
              DI(K)=0.0D0
              BK(K)=1.0D+300
10            DK(K)=-1.0D+300
           BI(0)=1.0D0
           DI(1)=0.5D0
           RETURN
        ENDIF
        CALL IK01A(X,BI0,DI0,BI1,DI1,BK0,DK0,BK1,DK1)
        BI(0)=BI0
        BI(1)=BI1
        BK(0)=BK0
        BK(1)=BK1
        DI(0)=DI0
        DI(1)=DI1
        DK(0)=DK0
        DK(1)=DK1
        IF (N.LE.1) RETURN
        IF (X.GT.40.0.AND.N.LT.INT(0.25*X)) THEN
           H0=BI0
           H1=BI1
           DO 15 K=2,N
           H=-2.0D0*(K-1.0D0)/X*H1+H0
           BI(K)=H
           H0=H1
15         H1=H
        ELSE
           M=MSTA1(X,200)
           IF (M.LT.N) THEN
              NM=M
           ELSE
              M=MSTA2(X,N,15)
           ENDIF
           F0=0.0D0
           F1=1.0D-100
           F=1.0D-100
           DO 20 K=M,0,-1
              F=2.0D0*(K+1.0D0)*F1/X+F0
              IF (K.LE.NM) BI(K)=F
              F0=F1
20            F1=F
           S0=BI0/F
           DO 25 K=0,NM
25            BI(K)=S0*BI(K)
        ENDIF
        G0=BK0
        G1=BK1
        DO 30 K=2,NM
           G=2.0D0*(K-1.0D0)/X*G1+G0
           BK(K)=G
           G0=G1
30         G1=G
        DO 40 K=2,NM
           DI(K)=BI(K-1)-K/X*BI(K)
40         DK(K)=-BK(K-1)-K/X*BK(K)
        RETURN
        END


        SUBROUTINE IK01A(X,BI0,DI0,BI1,DI1,BK0,DK0,BK1,DK1)
C
C       =========================================================
C       Purpose: Compute modified Bessel functions I0(x), I1(1),
C                K0(x) and K1(x), and their derivatives
C       Input :  x   --- Argument ( x � 0 )
C       Output:  BI0 --- I0(x)
C                DI0 --- I0'(x)
C                BI1 --- I1(x)
C                DI1 --- I1'(x)
C                BK0 --- K0(x)
C                DK0 --- K0'(x)
C                BK1 --- K1(x)
C                DK1 --- K1'(x)
C       =========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(12),B(12),A1(8)
        PI=3.141592653589793D0
        EL=0.5772156649015329D0
        X2=X*X
        IF (X.EQ.0.0D0) THEN
           BI0=1.0D0
           BI1=0.0D0
           BK0=1.0D+300
           BK1=1.0D+300
           DI0=0.0D0
           DI1=0.5D0
           DK0=-1.0D+300
           DK1=-1.0D+300
           RETURN
        ELSE IF (X.LE.18.0) THEN
           BI0=1.0D0
           R=1.0D0
           DO 15 K=1,50
              R=0.25D0*R*X2/(K*K)
              BI0=BI0+R
              IF (DABS(R/BI0).LT.1.0D-15) GO TO 20
15         CONTINUE
20         BI1=1.0D0
           R=1.0D0
           DO 25 K=1,50
              R=0.25D0*R*X2/(K*(K+1))
              BI1=BI1+R
              IF (DABS(R/BI1).LT.1.0D-15) GO TO 30
25         CONTINUE
30         BI1=0.5D0*X*BI1
        ELSE
           DATA A/0.125D0,7.03125D-2,
     &            7.32421875D-2,1.1215209960938D-1,
     &            2.2710800170898D-1,5.7250142097473D-1,
     &            1.7277275025845D0,6.0740420012735D0,
     &            2.4380529699556D01,1.1001714026925D02,
     &            5.5133589612202D02,3.0380905109224D03/
           DATA B/-0.375D0,-1.171875D-1,
     &            -1.025390625D-1,-1.4419555664063D-1,
     &            -2.7757644653320D-1,-6.7659258842468D-1,
     &            -1.9935317337513D0,-6.8839142681099D0,
     &            -2.7248827311269D01,-1.2159789187654D02,
     &            -6.0384407670507D02,-3.3022722944809D03/
           K0=12
           IF (X.GE.35.0) K0=9
           IF (X.GE.50.0) K0=7
           CA=DEXP(X)/DSQRT(2.0D0*PI*X)
           BI0=1.0D0
           XR=1.0D0/X
           DO 35 K=1,K0
35            BI0=BI0+A(K)*XR**K
           BI0=CA*BI0
           BI1=1.0D0
           DO 40 K=1,K0
40            BI1=BI1+B(K)*XR**K
           BI1=CA*BI1
        ENDIF
        IF (X.LE.9.0D0) THEN
           CT=-(DLOG(X/2.0D0)+EL)
           BK0=0.0D0
           W0=0.0D0
           R=1.0D0
           WW=0.0D0
           DO 65 K=1,50
              W0=W0+1.0D0/K
              R=0.25D0*R/(K*K)*X2
              BK0=BK0+R*(W0+CT)
              IF (DABS((BK0-WW)/BK0).LT.1.0D-15) GO TO 70
65            WW=BK0
70         BK0=BK0+CT
        ELSE
           DATA A1/0.125D0,0.2109375D0,
     &             1.0986328125D0,1.1775970458984D01,
     &             2.1461706161499D02,5.9511522710323D03,
     &             2.3347645606175D05,1.2312234987631D07/
           CB=0.5D0/X
           XR2=1.0D0/X2
           BK0=1.0D0
           DO 75 K=1,8
75            BK0=BK0+A1(K)*XR2**K
           BK0=CB*BK0/BI0
        ENDIF
        BK1=(1.0D0/X-BI1*BK0)/BI0
        DI0=BI1
        DI1=BI0-BI1/X
        DK0=-BK1
        DK1=-BK0-BK1/X
        RETURN
        END


        INTEGER FUNCTION MSTA1(X,MP)
C
C       ===================================================
C       Purpose: Determine the starting point for backward  
C                recurrence such that the magnitude of    
C                Jn(x) at that point is about 10^(-MP)
C       Input :  x     --- Argument of Jn(x)
C                MP    --- Value of magnitude
C       Output:  MSTA1 --- Starting point   
C       ===================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        A0=DABS(X)
        N0=INT(1.1*A0)+1
        F0=ENVJ(N0,A0)-MP
        N1=N0+5
        F1=ENVJ(N1,A0)-MP
        DO 10 IT=1,20             
           NN=N1-(N1-N0)/(1.0D0-F0/F1)                  
           F=ENVJ(NN,A0)-MP
           IF(ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
 10        F1=F
 20     MSTA1=NN
        RETURN
        END


        INTEGER FUNCTION MSTA2(X,N,MP)
C
C       ===================================================
C       Purpose: Determine the starting point for backward
C                recurrence such that all Jn(x) has MP
C                significant digits
C       Input :  x  --- Argument of Jn(x)
C                n  --- Order of Jn(x)
C                MP --- Significant digit
C       Output:  MSTA2 --- Starting point
C       ===================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        A0=DABS(X)
        HMP=0.5D0*MP
        EJN=ENVJ(N,A0)
        IF (EJN.LE.HMP) THEN
           OBJ=MP
           N0=INT(1.1*A0)
        ELSE
           OBJ=HMP+EJN
           N0=N
        ENDIF
        F0=ENVJ(N0,A0)-OBJ
        N1=N0+5
        F1=ENVJ(N1,A0)-OBJ
        DO 10 IT=1,20
           NN=N1-(N1-N0)/(1.0D0-F0/F1)
           F=ENVJ(NN,A0)-OBJ
           IF (ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
10         F1=F
20      MSTA2=NN+10
        RETURN
        END

        REAL*8 FUNCTION ENVJ(N,X)
        DOUBLE PRECISION X
        ENVJ=0.5D0*DLOG10(6.28D0*N)-N*DLOG10(1.36D0*X/N)
        RETURN
        END

