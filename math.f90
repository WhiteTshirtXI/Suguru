      module mod_math
      use mod_common
      use mod_param
      public 
      contains
    





!-----------------------------------------------------------------------
      SUBROUTINE BLANK(A,N)
      CHARACTER*1 A(1)
      CHARACTER*1 BLNK
      SAVE        BLNK
      DATA        BLNK /' '/

      DO 10 I=1,N
         A(I)=BLNK
   10 CONTINUE
      RETURN
      END SUBROUTINE BLANK
!-----------------------------------------------------------------------
      SUBROUTINE VSQ (A,N)
      DIMENSION  A(1)


      DO 100 I = 1, N
 100     A(I) = A(I)**2
      RETURN
      END SUBROUTINE VSQ
!-----------------------------------------------------------------------
      SUBROUTINE VSQRT(A,N)
      DIMENSION  A(1)


      DO 100 I = 1, N
 100     A(I) = SQRT(A(I))
      RETURN
      END SUBROUTINE VSQRT
!-----------------------------------------------------------------------
      subroutine invers2(a,b,n)
      REAL A(1),B(1)


      DO 100 I=1,N
         A(I)=1./B(I)
 100  CONTINUE
      return
      END subroutine invers2
!-----------------------------------------------------------------------
      subroutine invcol1(a,n)
      REAL A(1)


      DO 100 I=1,N
         A(I)=1./A(I)
 100  CONTINUE
      return
      END subroutine invcol1
!-----------------------------------------------------------------------
      subroutine invcol2(a,b,n)

      REAL A(1),B(1)

      DO 100 I=1,N
         A(I)=A(I)/B(I)
 100  CONTINUE
      return
      END subroutine invcol2
!-----------------------------------------------------------------------
      subroutine invcol3(a,b,c,n)
      REAL A(1),B(1),C(1)


      DO 100 I=1,N
         A(I)=B(I)/C(I)
 100  CONTINUE
      return
      END subroutine invcol3
!-----------------------------------------------------------------------
      subroutine col4(a,b,c,d,n)
      REAL A(1),B(1),C(1),D(1)


      DO 100 I=1,N
         A(I)=B(I)*C(I)*D(I)
  100 CONTINUE
      return
      END subroutine col4
!-----------------------------------------------------------------------
      subroutine Xaddcol3(a,b,c,n)
      REAL A(1),B(1),C(1)


      DO 100 I=1,N
         A(I)=A(I)+B(I)*C(I)
  100 CONTINUE
      return
      END subroutine Xaddcol3
!-----------------------------------------------------------------------
      subroutine addcol4(a,b,c,d,n)
      REAL A(1),B(1),C(1),D(1)


      DO 100 I=1,N
         A(I)=A(I)+B(I)*C(I)*D(I)
  100 CONTINUE
      return
      END subroutine addcol4
!-----------------------------------------------------------------------
      subroutine ascol5 (a,b,c,d,e,n)
      REAL A(1),B(1),C(1),D(1),E(1)


      DO 100 I=1,N
         A(I) = B(I)*C(I)-D(I)*E(I)
 100  CONTINUE
      return
      END subroutine ascol5
!-----------------------------------------------------------------------
      subroutine sub2(a,b,n)
      REAL A(1),B(1)


      DO 100 I=1,N
         A(I)=A(I)-B(I)
 100  CONTINUE
      return
      END subroutine sub2
!-----------------------------------------------------------------------
      subroutine sub3(a,b,c,n)
      REAL A(1),B(1),C(1)


      DO 100 I=1,N
         A(I)=B(I)-C(I)
 100  CONTINUE
      return
      END subroutine sub3
!-----------------------------------------------------------------------
      subroutine subcol3(a,b,c,n)
      REAL A(1),B(1),C(1)


      DO 100 I=1,N
         A(I)=A(I)-B(I)*C(I)
  100 CONTINUE
      return
      END subroutine subcol3
!-----------------------------------------------------------------------
      subroutine subcol4(a,b,c,d,n)
      REAL A(1),B(1),C(1),D(1)

      DO 100 I=1,N
         A(I)=A(I)-B(I)*C(I)*D(I)
  100 CONTINUE
      return
      END subroutine subcol4
!-----------------------------------------------------------------------
      subroutine rzero(a,n)
      DIMENSION  A(1)
      DO 100 I = 1, N
 100     A(I ) = 0.0
      return
      END subroutine rzero
!-----------------------------------------------------------------------
      subroutine izero(a,n)
      INTEGER A(1)

      DO 100 I = 1, N
 100     A(I ) = 0
      return
      END subroutine izero
!-----------------------------------------------------------------------
      subroutine ione(a,n)
      INTEGER   A(1)
      DO 100 I = 1, N
 100     A(I ) = 1
      return
      END subroutine ione
!-----------------------------------------------------------------------
      subroutine rone(a,n)
      DIMENSION  A(1)
      DO 100 I = 1, N
 100     A(I ) = 1.0
      return
      END subroutine rone
!-----------------------------------------------------------------------
      subroutine cfill(a,b,n)
      DIMENSION  A(1)

      DO 100 I = 1, N
 100     A(I) = B
      return
      END subroutine cfill
!-----------------------------------------------------------------------
      subroutine ifill(ia,ib,n)
      DIMENSION IA(1)

      DO 100 I = 1, N
 100     IA(I) = IB
      return
      END subroutine ifill
!-----------------------------------------------------------------------
      subroutine copy(a,b,n)
      real a(1),b(1)

      do i=1,n
         a(i)=b(i)
      enddo

      return
      end subroutine copy
!-----------------------------------------------------------------------
      subroutine chcopy(a,b,n)
      CHARACTER*1 A(1), B(1)

      DO 100 I = 1, N
 100     A(I) = B(I)
      return
      END subroutine chcopy

      subroutine icopy(a,b,n)
      INTEGER A(1), B(1)

      DO 100 I = 1, N
 100     A(I) = B(I)
      return
      END subroutine icopy
!-----------------------------------------------------------------------
      subroutine i8copy(a,b,n)
      INTEGER*8 A(1), B(1)

      DO 100 I = 1, N
 100     A(I) = B(I)
      return
      END subroutine i8copy
!-----------------------------------------------------------------------
      subroutine chsign(a,n)
      REAL A(1)

      DO 100 I=1,N
         A(I) = -A(I)
 100  CONTINUE
      return
      END subroutine chsign

!-----------------------------------------------------------------------
      subroutine cmult(a,const,n)
      REAL A(1)


      DO 100 I=1,N
         A(I)=A(I)*CONST
 100  CONTINUE
      return
      END subroutine cmult
!-----------------------------------------------------------------------
      subroutine cadd(a,const,n)
      REAL A(1)


      DO 100 I=1,N
         A(I)=A(I)+CONST
 100  CONTINUE
      return
      END subroutine cadd
!-----------------------------------------------------------------------
      subroutine iadd(i1,iscal,n)
      DIMENSION I1(1)

      DO 10 I=1,N
         I1(I)=I1(I)+ISCAL
   10 CONTINUE
      return
      END subroutine iadd
!-----------------------------------------------------------------------
      subroutine cadd2(a,b,const,n)
      REAL A(1),B(1)

      DO 100 I=1,N
         A(I)=B(I)+CONST
 100  CONTINUE
      return
      END subroutine cadd2
!-----------------------------------------------------------------------
      real function vlmin(vec,n)
      REAL VEC(1)
      TMIN = 99.0E20

      DO 100 I=1,N
         TMIN = MIN(TMIN,VEC(I))
 100  CONTINUE
      VLMIN = TMIN
      return
      END function vlmin
!-----------------------------------------------------------------------
      integer function ivlmin(vec,n)
      integer vec(1),tmin
      if (n.eq.0) then
         ivlmin=0
         return
      endif
      tmin = 8888888
      do i=1,n
         tmin = min(tmin,vec(i))
      enddo
      ivlmin = tmin
      return
      end function ivlmin
!-----------------------------------------------------------------------
      integer function ivlmax(vec,n)
      integer vec(1),tmax
      if (n.eq.0) then
         ivlmax=0
         return
      endif
      TMAX =-8888888
      do i=1,n
         TMAX = MAX(TMAX,VEC(I))
      enddo
      Ivlmax = tmax
      return
      end function ivlmax
!-----------------------------------------------------------------------
      real function vlmax(vec,n)
      REAL VEC(1)
      TMAX =-99.0E20
      do i=1,n
         TMAX = MAX(TMAX,VEC(I))
      enddo
      VLMAX = TMAX
      return
      END function vlmax
!-----------------------------------------------------------------------
      real function vlamax(vec,n)
      REAL VEC(1)
      TAMAX = 0.0

      DO 100 I=1,N
         TAMAX = MAX(TAMAX,ABS(VEC(I)))
 100  CONTINUE
      VLAMAX = TAMAX
      return
      END function vlamax
!-----------------------------------------------------------------------
      real function vlsum(vec,n)
      REAL VEC(1)

      SUM = 0.

      DO 100 I=1,N
         SUM=SUM+VEC(I)
 100  CONTINUE
      VLSUM = SUM
      return
      END function vlsum
!-----------------------------------------------------------------------
      subroutine vcross (u1,u2,u3,v1,v2,v3,w1,w2,w3,n)

!     Compute a Cartesian vector cross product.

      DIMENSION U1(1),U2(1),U3(1)
      DIMENSION V1(1),V2(1),V3(1)
      DIMENSION W1(1),W2(1),W3(1)

      DO 100 I=1,N
         U1(I) = V2(I)*W3(I) - V3(I)*W2(I)
         U2(I) = V3(I)*W1(I) - V1(I)*W3(I)
         U3(I) = V1(I)*W2(I) - V2(I)*W1(I)
  100 CONTINUE
      return
      END subroutine vcross
!-----------------------------------------------------------------------
      subroutine vdot2 (dot,u1,u2,v1,v2,n)

!    Compute a Cartesian vector dot product. 2-d version

      DIMENSION DOT(1)
      DIMENSION U1(1),U2(1)
      DIMENSION V1(1),V2(1)


      DO 100 I=1,N
         DOT(I) = U1(I)*V1(I) + U2(I)*V2(I) 
  100 CONTINUE
      return
      END subroutine vdot2
!-----------------------------------------------------------------------
      subroutine vdot3 (dot,u1,u2,u3,v1,v2,v3,n)

!     Compute a Cartesian vector dot product. 3-d version
      real:: dot
      !DIMENSION DOT(1)
      DIMENSION U1(1),U2(1),U3(1)
      DIMENSION V1(1),V2(1),V3(1)


      DO 100 I=1,N
         dot = U1(I)*V1(I) + U2(I)*V2(I) + U3(I)*V3(I)
         !DOT(I) = U1(I)*V1(I) + U2(I)*V2(I) + U3(I)*V3(I)
  100 CONTINUE
      return
      END subroutine vdot3
!-----------------------------------------------------------------------
      subroutine addtnsr(s,h1,h2,h3,nx,ny,nz)

!     Map and add to S a tensor product form of the three functions H1,H2,H3.
!     This is a single element routine used for deforming geometry.

      DIMENSION H1(1),H2(1),H3(1)
      DIMENSION S(NX,NY,NZ)

      DO 200 IZ=1,NZ
      DO 200 IY=1,NY
         HH = H2(IY)*H3(IZ)
         DO 100 IX=1,NX
            S(IX,IY,IZ)=S(IX,IY,IZ)+HH*H1(IX)
  100    CONTINUE
  200 CONTINUE
      return
      END subroutine addtnsr
      function ltrunc(string,l)
      CHARACTER*1 STRING(L)
      CHARACTER*1   BLNK
      DATA BLNK/' '/

      DO 100 I=L,1,-1
         L1=I
         IF (STRING(I).NE.BLNK) GOTO 200
  100 CONTINUE
      L1=0
  200 CONTINUE
      LTRUNC=L1
      return
      END function ltrunc
!-----------------------------------------------------------------------
      function mod1(i,n)

!    Yields MOD(I,N) with the exception that if I=K*N, result is N.

      MOD1=0
      IF (I.EQ.0) THEN
         return
      ENDIF
      IF (N.EQ.0) THEN
         WRITE(6,*) &
       'WARNING:  Attempt to take MOD(I,0) in function mod1.'
         return
      ENDIF
      II = I+N-1
      MOD1 = MOD(II,N)+1
      return
      END function mod1
!-----------------------------------------------------------------------
      integer function log2(k)
      RK=(K)
      RLOG=LOG10(RK)
      RLOG2=LOG10(2.0)
      RLOG=RLOG/RLOG2+0.5
      LOG2=INT(RLOG)
      return
      END function log2
!-----------------------------------------------------------------------
      subroutine iflip(i1,n)
      DIMENSION I1(1)
      N1=N+1
      N2=N/2
      DO 10 I=1,N2
         ILAST=N1-I
         ITMP=I1(ILAST)
         I1(ILAST)=I1(I)
         I1(I)=ITMP
   10 CONTINUE
      return
      END subroutine iflip
!-----------------------------------------------------------------------
      subroutine iswap(b,ind,n,temp)
      INTEGER B(1),IND(1),TEMP(1)

!  SORT ASSOCIATED ELEMENTS BY PUTTING ITEM(JJ)
!  INTO ITEM(I), WHERE JJ=IND(I).

      DO 20 I=1,N
         JJ=IND(I)
         TEMP(I)=B(JJ)
   20 CONTINUE
      DO 30 I=1,N
   30 B(I)=TEMP(I)
      return
      END subroutine iswap
!-----------------------------------------------------------------------
      subroutine col2(a,b,n)
      real a(1),b(1)
      do i=1,n
         a(i)=a(i)*b(i)
      enddo

      return
      end subroutine col2
!-----------------------------------------------------------------------
      subroutine col2c(a,b,c,n)
      real a(1),b(1),c

      do i=1,n
         a(i)=a(i)*b(i)*c
      enddo

      return
      end subroutine col2c
!-----------------------------------------------------------------------
      subroutine col3(a,b,c,n)
      real a(1),b(1),c(1)
      do i=1,n
         a(i)=b(i)*c(i)
      enddo
      return
      end subroutine col3
!-----------------------------------------------------------------------
      subroutine add2(a,b,n)
      real a(1),b(1)
!xbm* unroll (10)
      do i=1,n
         a(i)=a(i)+b(i)
      enddo
      return
      end subroutine add2
!-----------------------------------------------------------------------
      subroutine add3(a,b,c,n)
      real a(1),b(1),c(1)
      do i=1,n
         a(i)=b(i)+c(i)
      enddo
      return
      end subroutine add3
!-----------------------------------------------------------------------
      subroutine addcol3(a,b,c,n)
      real a(1),b(1),c(1)
      do i=1,n
         a(i)=a(i)+b(i)*c(i)
      enddo
      return
      end subroutine addcol3
!-----------------------------------------------------------------------
      subroutine add2s1(a,b,c1,n)
      real a(1),b(1)

      DO 100 I=1,N
        A(I)=C1*A(I)+B(I)
  100 CONTINUE
      return
      END subroutine add2s1

!-----------------------------------------------------------------------
      subroutine add2s2(a,b,c1,n)
      real a(1),b(1)

      DO 100 I=1,N
        A(I)=A(I)+C1*B(I)
  100 CONTINUE
      return
      END subroutine add2s2

!-----------------------------------------------------------------------
      subroutine add3s2(a,b,c,c1,c2,n)
      real a(1),b(1),c(1)

      DO 100 I=1,N
        A(I)=C1*B(I)+C2*C(I)
  100 CONTINUE
      return
      END subroutine add3s2

!-----------------------------------------------------------------------
      subroutine add4(a,b,c,d,n)
      REAL A(1),B(1),C(1),D(1)

      DO 100 I=1,N
         A(I)=B(I)+C(I)+D(I)
 100  CONTINUE
      return
      END subroutine add4

      real function vlsc2(x,y,n)
      REAL X(1),Y(1)
      s = 0.
      do i=1,n
         s = s + x(i)*y(i)
      enddo
      vlsc2=s
      return
      end function vlsc2
!-----------------------------------------------------------------------
      real function vlsc21(x,y,n)
      real x(1),y(1)

      s = 0.
      do i=1,n
         s = s + x(i)*x(i)*y(i)
      enddo
      vlsc21=s
      return
      end function vlsc21


!----------------------------------------------------------------------------

!     Vector reduction routines which require communication 
!     on a parallel machine. These routines must be substituted with
!     appropriate routines which take into account the specific architecture.

!----------------------------------------------------------------------------




      function rnorm2(x,y,z)
      real rnorm2,x,y,z
      rnorm2=sqrt(x**2+y**2+z**2)
      return
      end function rnorm2

      subroutine ind2to1(ind,irow,jcol,rowmax,colmax)
      integer rowmax,colmax
      ind = (jcol-1)*rowmax + irow
      return
      end subroutine ind2to1

      subroutine ind1to2(irow,jcol,ind,rowmax,colmax)
      integer rowmax,colmax
      irow = mod(ind,rowmax)
      jcol = (ind-irow)/rowmax + 1
      if(irow.eq.0) then
      jcol = jcol - 1
      irow = rowmax
      end if
      return
      end subroutine ind1to2

      subroutine cal_dr(dx,dy,dz,x1,y1,z1,x2,y2,z2)
      dx = x1 - x2
      dy = y1 - y2
      dz = z1 - z2
      return
      end subroutine cal_dr



      function vec2ave(a,m,n)
      real vec2ave
      real a(m,n)
      vec2ave = 0
      vec2ave = sum(a)/(m*n)
      return
      end function vec2ave


      end module mod_math
