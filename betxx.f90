!       PROGRAM solving DU/Dt +  DE/Dx +DF/Dy = 0 using Runge Kutta...
!	variable RK accuracy...in fact anything!
!	plus hu/hussaini/manthey LDDRK coeffs put in 
!	1d working at moment....
!	2d with solid wall almost there... ie. checked major bits...10/3/98
!	2d  working too (26/11/98)

!	note buffer treatment not very elegant yet... but 
!	will do.....for moment..ie should have NPINTS specified
!	with DX DY calculated rather than as present.....
!	put in split equations and buffering treatment.. all well
!	 now working on variallble buffer regiin extent.... (2/12/98)

! basic array setup module     
        Module param_info
	integer, parameter :: MX=65, MT=25, LEQ=8
	integer, parameter :: MY=120
	integer, parameter :: NCBUF = 25
     end module param_info 
         
        Module EQUATON_GEQAT
	use param_info, only : LEQ
	DOUBLE PRECISION, DIMENSION(LEQ, LEQ) ::  E, F
	DOUBLE PRECISION, DIMENSION(LEQ) :: GB, HB, CB
     end module EQUATON_GEQAT 
!	
 
	PROGRAM TEST
	use param_info
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	integer, parameter :: MXMAX = 2*MX+NCBUF, MYMAX = MY+NCBUF
	DOUBLE PRECISION FSOL(-NCBUF:MXMAX,0:MYMAX,MT,LEQ),TPOS(MXMAX)
	DOUBLE PRECISION BCONDS(-NCBUF:MXMAX,0:MYMAX,LEQ)
	EXTERNAL UZERO
        DOUBLE PRECISION EX(MT)
	DOUBLE PRECISION XPOS(-NCBUF:MXMAX),YPOS(0:MYMAX)
	COMMON /INFLOW/XPOS,YPOS
	COMMON /OUTRET/EX
	COMMON /BIODATA/BCONDS
	COMMON /SIMPL/DT,TO,TSTART
	COMMON /ACCURCY/NORDER,ISIDE
	COMMON /OUTPUT/NTIME
	COMMON /SIGVALU/SIG0,BETA,DRAT
	COMMON /BUFFER/NC
	
	FSOL=0.d0 ; BCONDS=0.d0
	
!	this is really what matters....
!	final time step
	OPEN(9,FILE='3D.DAT')
	READ(9,*)TFINAL,DT,NTIME
	READ(9,*)XMIN,XMAX,DX,NORDER
	READ(9,*)YMAX,DY,NORDERY
	READ(9,*)FRACTX,FRACTY
	READ(9,*)ISTAGE,IEX,ITAM ! <<<
	READ(9,*)SIG0,BETA,NC
! ************************************************************
!	TFINAL=Final time step
!	DT=time interval
!	XMAX= maximum xvalue
!	XMIN = minimum xvalue
!	DX= Xstep size
!	NORDER = Order of finite differencing in X direction
!                up to 10 , multiples of 2
!	NORDERY = Order of finite differencing in y direction
!                up to 10 , multiples of 2
!	but you will get wrong answers if set to 10 (in SUBY)
!	though not much point with NORDER >6ish...probs
!	NTIME = number of time levels to be saved (max 25)
!	FRACTX,FRACTY = Buffer region as fraction of XMAX and YMAX...
!	ISTAGE = Runge kutta accuracy, ie. 4 order ISTAGE =4 etc.
!	IEX =1 Low dissipation dispersion active .. though doesn't realy
!	do that much.. just as effective to use Hi RK + smaller DX...
!	IEX =0 Convetional RK.....
!	ITAM=1 TAM WEBB LDR 4thorder scheme for fd's only..
!	ITAM= 0  6th order fd 

! ************************************************************

!	setting buffer region size	
!	NC=10

	TOLGRID=0.0000001D0
	TSTART=0.D0
	NDIV=TFINAL/DT

	H=DX

!	buffer absorption coeef == 1/NC effectively...
	DRAT=1.D0/NC

        CALL DIFFER(DX,1,ISTAGE,IEX,ITAM)
        CALL DIFFER(DY,2,ISTAGE,IEX,ITAM)

	IRK=ISTAGE
	ISIDE=NORDER/2

	PRINT*,'input NEQ'
	NEQ=8

	CALL SETMX(NEQ)

! PML attempt....
	BUFMIN=DX*NC
	BUFYIN=DY*NC
	XBUF=FRACTX*XMAX
	YBUF=FRACTY*YMAX

	IF(XBUF > BUFMIN)THEN
	BUF=XBUF
	ELSE
	BUF=BUFMIN
	END IF

        IF(YBUF > BUFYIN)THEN
        YUF=YBUF
        ELSE
        YUF=BUFYIN
        END IF

	YMAX=YMAX+YUF
	YZERO=0.D0

	
	XMAX=XMAX+BUF
	XZERO=XMIN-BUF

	IYCHEK=(YMAX-YZERO)/DY
	IXCHEK=(XMAX-XZERO)/DX

	IERR=0
	IF(IYCHEK < MYMAX)GOTO 11022
	WRITE(*,77077)IYCHEK,MYMAX
	IERR=1

11022	IF(IXCHEK < 2*MXMAX)GOTO 11011
	WRITE(*,77076)IXCHEK,MYMAX
	IERR=1

11011	CONTINUE
	IF(IERR == 1)GOTO 99099

	IX=-NC
	XO=XZERO
	PRINT*,XZERO,XMAX,'XZERO XMAX'
1111	CONTINUE
	XPOS(IX)=XO
	PRINT*,IX,' ix is....',XO
	IX=IX+1
	XO=XO+DX
	IF(XO <= XMAX+TOLGRID)GOTO 1111

!	hence inflow is at IX=0
!	outflow is at IX-1-NC	
	INFLO=-NC
	INOUT=IX-1-NC
!	but begin calculations at INFLO=1
	NPOS=INOUT-1

!	total number of stations is NPOS	
	N=NPOS
	NX=NPOS
	PRINT*,'IX,NX,N',IX,NX,N
	IFLAG=0
	PRINT*,'read out'

	IY=0
	YO=0.D0
2222	CONTINUE
	YPOS(IY)=YO
	PRINT*,' IY is..',IY,YPOS(IY)
	IY=IY+1
	YO=YO+DY
	IF(YO <= YMAX+TOLGRID)GOTO 2222
	NY=IY-1-NC
	PRINT*,NX,IX,XO-DX,'  NX  '
	PRINT*,NY,IY,YO-DY,'  NY   '

	PRINT*,NY,' NY is'

	DO I=-NC,NX+NC+1
	DO J=0,NY+NC
	BCONDS(I,J,1:2)=UZERO(XPOS(I),YPOS(J))	
	BCONDS(I,J,3:6)=0.D0
	BCONDS(I,J,7:8)=BCONDS(I,J,1)
	END DO
	END DO


	IEQ=1
        WRITE(*,13)(XPOS(I),I=-NC,NX+NC+1)
	DO J=NY+NC,0,-1
        WRITE(*,11)YPOS(J),(BCONDS(I,J,IEQ),I=-NC,NX+NC+1)
	END DO
	WRITE(*,*)

!	

!	xfinal value is at NX+NC+1
!	yfinal value is at NY+NC


	TO=TSTART
	CALL INFULL(NDIV,N,NY,IRK,TPOS,FSOL,ITI,NEQ)


	WRITE(*,*)'J =',J
	DO JT=1,ITI
        DO I=-NC,N+NC+1
	DO J=NY+NC,0,-1
        WRITE(JT+19,11)XPOS(I),YPOS(J), &
     &  ((FSOL(I,J,JT,IEQ)+FSOL(I,J,JT,IEQ+1)),IEQ=1,7,2), &
     &  (FSOL(I,J,JT,IEQ),IEQ=1,NEQ)
	END DO
	WRITE(JT+19,*)
	END DO
	END DO
	
	
	
	JT=ITI/8
	PRINT*,' ITT is ..',ITI
	IEQ=1
	DO JT=1,ITI
	DO J=NY+NC,0,-1
        WRITE(121+JT-1,'(2001E19.9)') &
     & (FSOL(I,J,JT,IEQ)+FSOL(I,J,JT,IEQ+1), I=-NC,N+NC+1)
	END DO
	END DO


	STOP

99099	PRINT*,'Some error occured'

77077   FORMAT(4X,I7,' Elements Exceeded in Y-dir:', &
     & I7,' is max allowed')
77076   FORMAT(4X,I7,' Elements Exceeded in X-dir:', &
     & I7,' is max allowed')

11	FORMAT(27g13.5)
17	FORMAT(27F12.7)
13	FORMAT(8X,27F8.4)
12	FORMAT(49X,21F7.3)
	END
	
	DOUBLE PRECISION FUNCTION UZERO(X,Y)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	Xa=-25.D0
	Ya=30.D0
	BIT=DLOG(2.D0)*((X-Xa)**2+(Y-Ya)**2)/9.D0
	IF(DABS(BIT) < 75.0)THEN
        UZERO=0.5D0*DEXP(-BIT)
	ELSE
	UZERO=0.D0
	END IF
!        UZERO=0.5D0*DEXP(-DLOG(2.D0)*X*X/9.D0)
	RETURN
	END

	SUBROUTINE SETMX(NEQ)
	use param_info, only : LEQ
	use EQUATON_GEQAT
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	COMMON /SIGVALU/SIG0,BETA,DRAT
	COMMON /SPEED/XMX,XMY
	COMMON /BUFFER/NC
	XMX=1.0025D0
	XMY=0.D0
	
	E=0.D0 ; F=0.D0

	SIGX=SIG0
	SIGY=SIG0

        GB(1)=SIGX
        GB(2)=SIGY*0.D0
        GB(3)=SIGX
        GB(4)=SIGX
        GB(5)=SIGY*0.D0
        GB(6)=SIGX
        GB(7)=SIGX
        GB(8)=SIGY*0.D0

        HB(1)=SIGX*0.D0
        HB(2)=SIGY
        HB(3)=SIGX*0.D0
        HB(4)=SIGX*0.D0
        HB(5)=SIGY
        HB(6)=SIGX*0.D0
        HB(7)=SIGX*0.D0
        HB(8)=SIGY

        CB(1)=SIGX
        CB(2)=SIGY
        CB(3)=SIGX
        CB(4)=SIGX
        CB(5)=SIGY
        CB(6)=SIGX
        CB(7)=SIGX
        CB(8)=SIGY

        E(1,1)=XMX
        E(1,2)=XMX
        E(1,3)=1.D0
        E(1,4)=1.D0

        E(3,7)=1.D0
        E(3,8)=1.D0
        E(4,3)=XMX
        E(4,4)=XMX

        E(6,5)=XMX
        E(6,6)=XMX

        E(7,3)=1.D0
        E(7,4)=1.D0
        E(7,7)=XMX
        E(7,8)=XMX

        F(2,5)=1.D0
        F(2,6)=1.D0
        F(5,7)=1.D0
        F(5,8)=1.D0
        F(8,5)=1.D0
        F(8,6)=1.D0


	RETURN
	END
 
	SUBROUTINE INFULL(NDIV,N,NY,IRK,TPOS,S1,ITI,NEQ)
	use param_info
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        integer, parameter :: MXMAX = 2*MX+NCBUF, MYMAX = MY+NCBUF
        DOUBLE PRECISION S1(-NCBUF:MXMAX,0:MYMAX,MT,LEQ)
	DOUBLE PRECISION TPOS(MXMAX)
	DOUBLE PRECISION BC1(-NCBUF:MXMAX,0:MYMAX,LEQ)
	DOUBLE PRECISION Z1(-NCBUF:MXMAX,0:MYMAX,LEQ)
	DOUBLE PRECISION EX(MT),TCAP(MT)
	COMMON /SIMPL/DT,TO,TSTART
	COMMON /OUTPUT/NTIME
	COMMON /OUTRET/EX
	COMMON /BUFFER/NC

	NBITS=29
	DTOLT=1.D-7

	TCAP(1)=1.D0-DTOLT
	TCAP(2)=5.D0-DTOLT
	TCAP(3)=10.D0-DTOLT
	TCAP(4)=40.D0-DTOLT
	TCAP(5)=60.D0-DTOLT
	TCAP(6)=90.D0-DTOLT
	TCAP(7)=91.D0-DTOLT
	TCAP(8)=92.D0-DTOLT
	TCAP(9)=95.D0-DTOLT
	TCAP(10)=100.D0-DTOLT
	TCAP(11)=110.D0-DTOLT
	TCAP(12)=120.D0-DTOLT
	TCAP(13)=130.D0-DTOLT
	TCAP(14)=134.D0-DTOLT

	
	TO=TSTART
	CALL BOUNDRY(N,NY,TO,BC1,NEQ)
	DO IEQ=1,NEQ
	DO JY=0,NY+NC
	DO I=-NC,N+NC+1
	S1(I,JY,1,IEQ)=BC1(I,JY,IEQ)
	END DO
	END DO
	END DO
	EX(1)=TO
	TPOS(1)=TO

	TDIV=NDIV/NTIME
 
	IHI=0
	ITI=1
	DO 10 J=1,NDIV
	PRINT*,J
	CALL INTEG(N,NY,IRK,DT,BC1,Z1,NEQ)
 
	TO=TO+DT
!	TPOS(1+J)=TO
	CALL TREAT(N,NY,TO,Z1,NEQ)

	DO IEQ=1,NEQ
	DO JY=0,NY+NC
	DO I=-NC,N+NC+1
	BC1(I,JY,IEQ)=Z1(I,JY,IEQ)
	END DO
	END DO
	END DO

	IF(TO >= TCAP(ITI)  .AND.  ITI <= NBITS)THEN
	ITI=ITI+1
	PRINT*,'ITI == ', ITI,TO 
	DO IEQ=1,NEQ
	DO JY=0,NY+NC
	DO I=-NC,N+NC+1
	S1(I,JY,ITI,IEQ)=Z1(I,JY,IEQ)
	END DO
	END DO
	END DO
	EX(ITI)=TO
	END IF

3434	CONTINUE

10	CONTINUE

	RETURN
	END
 
	SUBROUTINE TREAT(N,NY,T,Z1,NEQ)
	use param_info, ONLY : MX, LEQ, MY, NCBUF
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        integer, parameter :: MXMAX = 2*MX+NCBUF, MYMAX = MY+NCBUF
        DOUBLE PRECISION XPOS(-NCBUF:MXMAX),YPOS(0:MYMAX)
	DOUBLE PRECISION Z1(-NCBUF:MXMAX,0:MYMAX,LEQ)
        COMMON /INFLOW/XPOS,YPOS
	COMMON /BUFFER/NC

!       keep on treating the left boundary with eaxcat
!       state....

	DO IEQ=1,NEQ
	DO J=0,NY+NC
!	Z1(-NC,J,IEQ)=UZERO(XPOS(-NC)-T,YPOS(J)-T)
!	Z1(N+NC+1,J,IEQ)=UZERO(XPOS(N+NC+1)-T,YPOS(J)-T)
	END DO
	END DO

!	treat upper y boundarty....
	DO I=-NC+1,N+NC
!	Z1(I,0,1)=UZERO(XPOS(I)-T,YPOS(0)-T)
!	Z1(I,NY+NC,1)=UZERO(XPOS(I)-T,YPOS(NY+NC)-T)
	END DO

        J=0
        IEQ=5
        DO I=-NC,N+NC+1
        Z1(I,J,IEQ)=0.D0 ; Z1(I,J,IEQ+1)=0.D0
	END DO

        DO J=NY+NC,NY+NC+1
        DO IEQ=1,NEQ
        DO I=-NC,N+NC+1
        Z1(I,J,IEQ)=0.D0
	END DO
	END DO
	END DO

        I=-NC
        DO J=0,NY+NC+1
        DO IEQ=1,NEQ
        Z1(I,J,IEQ)=0.D0
	END DO
	END DO

	RETURN
	END
 
	
	SUBROUTINE BOUNDRY(N,NY,YSTART,Z,NEQ)
	use param_info, ONLY : MX, LEQ, MY, NCBUF
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        integer, parameter :: MXMAX = 2*MX+NCBUF, MYMAX = MY+NCBUF
	DOUBLE PRECISION Z(-NCBUF:MXMAX,0:MYMAX,LEQ)
        DOUBLE PRECISION BCONDS(-NCBUF:MXMAX,0:MYMAX,LEQ)
	COMMON /BIODATA/BCONDS
	COMMON /BUFFER/NC
	DO IEQ=1,NEQ
	DO I=-NC,N+NC+1
	DO J=0,NY+NC
	Z(I,J,IEQ)=BCONDS(I,J,IEQ)
	END DO
	END DO
	END DO
	J=0
	IEQ=5
	DO I=-NC,N+NC+1
	Z(I,J,IEQ)=0.D0
	Z(I,J,IEQ+1)=0.D0
	END DO

	RETURN
	END

 
	SUBROUTINE INTEG(N,NY,IRK,DT,Z,ZN,NEQ)
	use param_info, ONLY : MX, LEQ, MY, NCBUF
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	integer, parameter :: KLEVEL = 12
        integer, parameter :: MXMAX = 2*MX+NCBUF, MYMAX = MY+NCBUF
	DOUBLE PRECISION, DIMENSION(-NCBUF:MXMAX,0:MYMAX,LEQ) :: Z, &
     & ZN, XK
	DOUBLE PRECISION BETA(KLEVEL)
	COMMON /LDDRK/BETA
	COMMON /BUFFER/NC

	CALL SUBFN(N,NY,DT,Z,XK,NEQ) 
	DO 10 IP=2,IRK
	DO 20 IEQ=1,NEQ
	DO 20 J=0,NY
	DO 20 I=0,N+1
	ZN(I,J,IEQ)=Z(I,J,IEQ)+BETA(IP)*XK(I,J,IEQ)
	IF(J == 0)THEN
	ZN(I,J,5)=0.D0
	ZN(I,J,6)=0.D0
	END IF
	IF(J == NY+NC) ZN(I,J,IEQ) = 0.D0
	IF(I == -NC) ZN(I,J,IEQ) =0.D0
20	CONTINUE

        DO 204 IEQ=1,NEQ
        DO 204 J=0,NY+NC
        DO 204 I=-NC,-1
        ZN(I,J,IEQ)=Z(I,J,IEQ)+BETA(IP)*XK(I,J,IEQ)
	IF(J == 0)THEN
	ZN(I,J,5)=0.D0
	ZN(I,J,6)=0.D0
	END IF
        IF(J == NY+NC) ZN(I,J,IEQ) = 0.D0
	IF(I == -NC) ZN(I,J,IEQ) =0.D0
204      CONTINUE


        DO 203 IEQ=1,NEQ
        DO 203 J=0,NY+NC
        DO 203 I=N+2,N+1+NC
        ZN(I,J,IEQ)=Z(I,J,IEQ)+BETA(IP)*XK(I,J,IEQ)
	IF(J == 0)THEN
	ZN(I,J,5)=0.D0
	ZN(I,J,6)=0.D0
	END IF
        IF(J == NY+NC) ZN(I,J,IEQ) = 0.D0
	IF(I == -NC) ZN(I,J,IEQ) =0.D0
203      CONTINUE

        DO 202 IEQ=1,NEQ
        DO 202 J=NY+1,NY+NC
        DO 202 I=0,N+1
        ZN(I,J,IEQ)=Z(I,J,IEQ)+BETA(IP)*XK(I,J,IEQ)
	IF(J == 0)THEN
	ZN(I,J,5)=0.D0
	ZN(I,J,6)=0.D0
	END IF
        IF(J == NY+NC) ZN(I,J,IEQ) = 0.D0
	IF(I == -NC) ZN(I,J,IEQ) =0.D0
202     CONTINUE

	CALL SUBFN(N,NY,DT,ZN,XK,NEQ)
10	CONTINUE

	DO 30 IEQ=1,NEQ
	DO 30 J=0,NY
	DO 30 I=0,N+1
	ZN(I,J,IEQ)=Z(I,J,IEQ)+XK(I,J,IEQ)
	IF(J == 0)THEN
	ZN(I,J,5)=0.D0
	ZN(I,J,6)=0.D0
	END IF
	IF(J == NY+NC) ZN(I,J,IEQ) = 0.D0
	IF(I == -NC) ZN(I,J,IEQ) =0.D0
30	CONTINUE


        DO 302 IEQ=1,NEQ
        DO 302 J=0,NY+NC
        DO 302 I=-NC,-1
        ZN(I,J,IEQ)=Z(I,J,IEQ)+XK(I,J,IEQ)
	IF(J == 0)THEN
	ZN(I,J,5)=0.D0
	ZN(I,J,6)=0.D0
	END IF
        IF(J == NY+NC) ZN(I,J,IEQ) = 0.D0
	IF(I == -NC) ZN(I,J,IEQ) =0.D0
302     CONTINUE

        DO 303 IEQ=1,NEQ
        DO 303 J=0,NY+NC
        DO 303 I=N+2,N+1+NC
        ZN(I,J,IEQ)=Z(I,J,IEQ)+XK(I,J,IEQ)
	IF(J == 0)THEN
	ZN(I,J,5)=0.D0
	ZN(I,J,6)=0.D0
	END IF
        IF(J == NY+NC) ZN(I,J,IEQ) = 0.D0
	IF(I == -NC) ZN(I,J,IEQ) =0.D0
303     CONTINUE

        DO 304 IEQ=1,NEQ
        DO 304 J=NY+1,NY+NC
        DO 304 I=0,N+1
        ZN(I,J,IEQ)=Z(I,J,IEQ)+XK(I,J,IEQ)
	IF(J == 0)THEN
	ZN(I,J,5)=0.D0
	ZN(I,J,6)=0.D0
	END IF
        IF(J == NY+NC) ZN(I,J,IEQ) = 0.D0
	IF(I == -NC) ZN(I,J,IEQ) =0.D0
304     CONTINUE

	RETURN
	END
 
	SUBROUTINE SUBFN(N,NY,DT,Z,XK,NEQ)
	use param_info, only : MX, LEQ, MY, NCBUF
	use EQUATON_GEQAT
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        integer, parameter :: NA = 12, KLEVEL = 12
        integer, parameter :: MXMAX = 2*MX+NCBUF, MYMAX = MY+NCBUF
	DOUBLE PRECISION Z(-NCBUF:MXMAX,0:MYMAX,LEQ)
        DOUBLE PRECISION COEFF(-NA:NA,NA,2)
	DOUBLE PRECISION XK(-NCBUF:MXMAX,0:MYMAX,LEQ)
	COMMON /SIGVALU/SIG0,BETA,DRAT
        COMMON /FINITE/COEFF
	COMMON /ACCURCY/NORDER,ISIDE
	COMMON /BUFFER/NC
!	note cf with lddrk.f ok if same orders used at inflow...

	DO 1001 J=0,NY

	DO 10 I=0,N+1
	DO 10 LX=1,NEQ
	SUM=0.D0
	DO 20 LY=1,NEQ
	DO 20 IK=-ISIDE,ISIDE
	SUM=SUM+COEFF(IK,NORDER,1)*Z(I+IK,J,LY)*E(LX,LY)
20	CONTINUE
	CALL SUBY(I,J,LX,NY,Z,SUM,NEQ)
	XK(I,J,LX)=-SUM*DT
10	CONTINUE

!	inflow boundary....
	CALL INFLOX(J,N,NY,DT,Z,XK,NEQ)

!	outflow boundary....
	CALL SOUTX(J,N,NY,DT,Z,XK,NEQ)

1001	CONTINUE

!	outflow in y direction
	CALL SOUTY(N,NY,DT,Z,XK,NEQ)

	RETURN              
	END

 
	SUBROUTINE SOUTY(N,NY,DT,Z,XK,NEQ)
	use param_info, only : MX, LEQ, MY, NCBUF
	use EQUATON_GEQAT
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        integer, parameter :: NA = 12, KLEVEL = 12
        integer, parameter :: MXMAX = 2*MX+NCBUF, MYMAX = MY+NCBUF
	DOUBLE PRECISION Z(-NCBUF:MXMAX,0:MYMAX,LEQ)
        DOUBLE PRECISION COEFF(-NA:NA,NA,2)
	DOUBLE PRECISION XK(-NCBUF:MXMAX,0:MYMAX,LEQ)
	COMMON /SIGVALU/SIG0,BETA,DRAT
        COMMON /FINITE/COEFF
	COMMON /ACCURCY/NORDER,ISIDE
	COMMON /BUFFER/NC
!	note cf with lddrk.f ok if same orders used at inflow...

	DO 1001 J=NY+1,NY+NC

	JO=J-NY

	DO 10 I=0,N+1
	DO 10 LX=1,NEQ
	SUM=0.D0
	DO 20 LY=1,NEQ
	DO 20 IK=-ISIDE,ISIDE
	SUM=SUM+COEFF(IK,NORDER,1)*Z(I+IK,J,LY)*E(LX,LY)
20	CONTINUE
	CALL SUBY(I,J,LX,NY,Z,SUM,NEQ)
	SUM=SUM+Z(I,J,LX)*HB(LX)*(DRAT*JO)**BETA
	XK(I,J,LX)=-SUM*DT
10	CONTINUE

!	inflow boundary....
	CALL INUFLOX(J,N,NY,DT,Z,XK,NEQ)

!	outflow boundary....
	CALL SUOUTX(J,N,NY,DT,Z,XK,NEQ)

1001	CONTINUE

	RETURN              
	END


        SUBROUTINE SOUTX(J,N,NY,DT,Z,XK,NEQ)
	use param_info, only : MX, LEQ, MY, NCBUF
	use EQUATON_GEQAT
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        integer, parameter :: NA = 12, KLEVEL = 12
        integer, parameter :: MXMAX = 2*MX+NCBUF, MYMAX = MY+NCBUF
        DOUBLE PRECISION Z(-NCBUF:MXMAX,0:MYMAX,LEQ)
        DOUBLE PRECISION COEFF(-NA:NA,NA,2)
        DOUBLE PRECISION XK(-NCBUF:MXMAX,0:MYMAX,LEQ)
	COMMON /SIGVALU/SIG0,BETA,DRAT
        COMMON /FINITE/COEFF
        COMMON /ACCURCY/NORDER,ISIDE
	COMMON /BUFFER/NC

	DO 3081 I=N+2,N+NC-3
        DO 2081 LX=1,NEQ
        SUM=0.D0
        DO 208 LY=1,NEQ
        DO 208 IK=-4,4
        SUM=SUM+COEFF(IK,8,1)*Z(I+IK,J,LY)*E(LX,LY)
208     CONTINUE
        CALL SUBY(I,J,LX,NY,Z,SUM,NEQ)
        XK(I,J,LX)=-SUM*DT
2081    CONTINUE
3081	CONTINUE


        I=N+NC-2
        DO 7021 LX=1,NEQ
        SUM=0.D0
        DO 702 LY=1,NEQ
        DO 702 IK=-3,3
        SUM=SUM+COEFF(IK,6,1)*Z(I+IK,J,LY)*E(LX,LY)
702     CONTINUE
        CALL SUBY(I,J,LX,NY,Z,SUM,NEQ)
        XK(I,J,LX)=-SUM*DT
7021    CONTINUE

        I=N+NC-1
        DO 701 LX=1,NEQ
        SUM=0.D0
        DO 70 LY=1,NEQ
        DO 70 IK=-2,2
        SUM=SUM+COEFF(IK,4,1)*Z(I+IK,J,LY)*E(LX,LY)
70      CONTINUE
        CALL SUBY(I,J,LX,NY,Z,SUM,NEQ)
        XK(I,J,LX)=-SUM*DT
701     CONTINUE

        I=N+NC
        DO 801 LX=1,NEQ
        SUM=0.D0
        DO 80 LY=1,NEQ
        DO 80 IK=-1,1
        SUM=SUM+COEFF(IK,2,1)*Z(I+IK,J,LY)*E(LX,LY)
80      CONTINUE
        CALL SUBY(I,J,LX,NY,Z,SUM,NEQ)
        XK(I,J,LX)=-SUM*DT
801     CONTINUE

        I=N+NC+1
        DO 10 LX=1,NEQ
        SUM=0.D0
        DO 20 LY=1,NEQ
        DO 20 IK=0,-6,-1
        SUM=SUM-COEFF(-IK,7,1)*Z(I+IK,J,LY)*E(LX,LY)
20      CONTINUE
        CALL SUBY(I,J,LX,NY,Z,SUM,NEQ)
        XK(I,J,LX)=-SUM*DT
10      CONTINUE

	RETURN
	END



        SUBROUTINE SUOUTX(J,N,NY,DT,Z,XK,NEQ)
	use param_info, only : MX, LEQ, MY, NCBUF
	use EQUATON_GEQAT
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        integer, parameter :: NA = 12, KLEVEL = 12
        integer, parameter :: MXMAX = 2*MX+NCBUF, MYMAX = MY+NCBUF
        DOUBLE PRECISION Z(-NCBUF:MXMAX,0:MYMAX,LEQ)
        DOUBLE PRECISION COEFF(-NA:NA,NA,2)
        DOUBLE PRECISION XK(-NCBUF:MXMAX,0:MYMAX,LEQ)
	COMMON /SIGVALU/SIG0,BETA,DRAT
        COMMON /FINITE/COEFF
        COMMON /ACCURCY/NORDER,ISIDE
	COMMON /BUFFER/NC

	JO=J-NY

	DO 3081 I=N+2,N+NC-3
        DO 2081 LX=1,NEQ
        SUM=0.D0
        DO 208 LY=1,NEQ
        DO 208 IK=-4,4
        SUM=SUM+COEFF(IK,8,1)*Z(I+IK,J,LY)*E(LX,LY)
208     CONTINUE
        CALL SUBY(I,J,LX,NY,Z,SUM,NEQ)
        SUM=SUM+Z(I,J,LX)*HB(LX)*(DRAT*JO)**BETA
        XK(I,J,LX)=-SUM*DT
2081    CONTINUE
3081	CONTINUE

        I=N+NC-2
        DO 7021 LX=1,NEQ
        SUM=0.D0
        DO 702 LY=1,NEQ
        DO 702 IK=-3,3
        SUM=SUM+COEFF(IK,6,1)*Z(I+IK,J,LY)*E(LX,LY)
702     CONTINUE
        CALL SUBY(I,J,LX,NY,Z,SUM,NEQ)
	SUM=SUM+Z(I,J,LX)*HB(LX)*(DRAT*JO)**BETA
        XK(I,J,LX)=-SUM*DT
7021    CONTINUE

        I=N+NC-1
        DO 701 LX=1,NEQ
        SUM=0.D0
        DO 70 LY=1,NEQ
        DO 70 IK=-2,2
        SUM=SUM+COEFF(IK,4,1)*Z(I+IK,J,LY)*E(LX,LY)
70      CONTINUE
        CALL SUBY(I,J,LX,NY,Z,SUM,NEQ)
	SUM=SUM+Z(I,J,LX)*HB(LX)*(DRAT*JO)**BETA
        XK(I,J,LX)=-SUM*DT
701     CONTINUE

        I=N+NC
        DO 801 LX=1,NEQ
        SUM=0.D0
        DO 80 LY=1,NEQ
        DO 80 IK=-1,1
        SUM=SUM+COEFF(IK,2,1)*Z(I+IK,J,LY)*E(LX,LY)
80      CONTINUE
        CALL SUBY(I,J,LX,NY,Z,SUM,NEQ)
	SUM=SUM+Z(I,J,LX)*HB(LX)*(DRAT*JO)**BETA
        XK(I,J,LX)=-SUM*DT
801     CONTINUE

        I=N+NC+1
        DO 10 LX=1,NEQ
        SUM=0.D0
        DO 20 LY=1,NEQ
        DO 20 IK=0,-6,-1
        SUM=SUM-COEFF(-IK,7,1)*Z(I+IK,J,LY)*E(LX,LY)
20      CONTINUE
        CALL SUBY(I,J,LX,NY,Z,SUM,NEQ)
	SUM=SUM+Z(I,J,LX)*HB(LX)*(DRAT*JO)**BETA
        XK(I,J,LX)=-SUM*DT
10      CONTINUE

	RETURN
	END




        SUBROUTINE INFLOX(J,N,NY,DT,Z,XK,NEQ)
	use param_info, only : MX, LEQ, MY, NCBUF
	use EQUATON_GEQAT
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        integer, parameter :: NA = 12, KLEVEL = 12
        integer, parameter :: MXMAX = 2*MX+NCBUF, MYMAX = MY+NCBUF
        DOUBLE PRECISION Z(-NCBUF:MXMAX,0:MYMAX,LEQ)
        DOUBLE PRECISION COEFF(-NA:NA,NA,2)
        DOUBLE PRECISION XK(-NCBUF:MXMAX,0:MYMAX,LEQ)
	COMMON /SIGVALU/SIG0,BETA,DRAT
        COMMON /FINITE/COEFF
        COMMON /ACCURCY/NORDER,ISIDE
	COMMON /BUFFER/NC

!       inflow boundary....
!       2nd-order

        I=-NC
        DO 3001 LX=1,NEQ
        SUM=0.D0
        DO 1320 LY=1,NEQ
        DO 1320 IK=0,6
        SUM=SUM+COEFF(IK,7,1)*Z(I+IK,J,LY)*E(LX,LY)
1320     CONTINUE
        CALL SUBY(I,J,LX,NY,Z,SUM,NEQ)
	SUM=SUM+Z(I,J,LX)*GB(LX)*(DRAT*NC)**BETA
        XK(I,J,LX)=-SUM*DT
3001    CONTINUE

        I=-NC+1
        DO 301 LX=1,NEQ
        SUM=0.D0
        DO 30 LY=1,NEQ
        DO 30 IK=-1,1
        SUM=SUM+COEFF(IK,2,1)*Z(I+IK,J,LY)*E(LX,LY)
30      CONTINUE
        CALL SUBY(I,J,LX,NY,Z,SUM,NEQ)
	SUM=SUM+Z(I,J,LX)*GB(LX)*(DRAT*(NC-1))**BETA
        XK(I,J,LX)=-SUM*DT
301     CONTINUE

!       4th order
        I=-NC+2
        DO 501 LX=1,NEQ
        SUM=0.D0
        DO 50 LY=1,NEQ
        DO 50 IK=-2,2
        SUM=SUM+COEFF(IK,4,1)*Z(I+IK,J,LY)*E(LX,LY)
50      CONTINUE
        CALL SUBY(I,J,LX,NY,Z,SUM,NEQ)
	SUM=SUM+Z(I,J,LX)*GB(LX)*(DRAT*(NC-2))**BETA
        XK(I,J,LX)=-SUM*DT
501     CONTINUE

!       6th order
        I=-NC+3
        DO 401 LX=1,NEQ
        SUM=0.D0
        DO 40 LY=1,NEQ
        DO 40 IK=-3,3
        SUM=SUM+COEFF(IK,6,1)*Z(I+IK,J,LY)*E(LX,LY)
40      CONTINUE
        CALL SUBY(I,J,LX,NY,Z,SUM,NEQ)
	SUM=SUM+Z(I,J,LX)*GB(LX)*(DRAT*(NC-3))**BETA
        XK(I,J,LX)=-SUM*DT
401     CONTINUE

	DO 601 I=-NC+4,-1
	JO=-I
        DO 2071 LX=1,NEQ
        SUM=0.D0
        DO 207 LY=1,NEQ
        DO 207 IK=-4,4
        SUM=SUM+COEFF(IK,8,1)*Z(I+IK,J,LY)*E(LX,LY)
207     CONTINUE
        CALL SUBY(I,J,LX,NY,Z,SUM,NEQ)
        SUM=SUM+Z(I,J,LX)*GB(LX)*(DRAT*JO)**BETA
        XK(I,J,LX)=-SUM*DT
2071    CONTINUE

601	CONTINUE


        RETURN
        END




        SUBROUTINE INUFLOX(J,N,NY,DT,Z,XK,NEQ)
	use param_info, only : MX, LEQ, MY, NCBUF
	use EQUATON_GEQAT
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        integer, parameter :: NA = 12, KLEVEL = 12
        integer, parameter :: MXMAX = 2*MX+NCBUF, MYMAX = MY+NCBUF
        DOUBLE PRECISION Z(-NCBUF:MXMAX,0:MYMAX,LEQ)
        DOUBLE PRECISION COEFF(-NA:NA,NA,2)
        DOUBLE PRECISION XK(-NCBUF:MXMAX,0:MYMAX,LEQ)
	COMMON /SIGVALU/SIG0,BETA,DRAT
        COMMON /FINITE/COEFF
        COMMON /ACCURCY/NORDER,ISIDE
	COMMON /BUFFER/NC

!       inflow boundary....
!       2nd-order

        I=-NC
        DO 3001 LX=1,NEQ
        SUM=0.D0
        DO 1320 LY=1,NEQ
        DO 1320 IK=0,6
        SUM=SUM+COEFF(IK,7,1)*Z(I+IK,J,LY)*E(LX,LY)
1320     CONTINUE
        CALL SUBY(I,J,LX,NY,Z,SUM,NEQ)
	SUM=SUM+Z(I,J,LX)*GB(LX)*(DRAT*NC)**BETA
        XK(I,J,LX)=-SUM*DT
3001    CONTINUE

        I=-NC+1
        DO 301 LX=1,NEQ
        SUM=0.D0
        DO 30 LY=1,NEQ
        DO 30 IK=-1,1
        SUM=SUM+COEFF(IK,2,1)*Z(I+IK,J,LY)*E(LX,LY)
30      CONTINUE
        CALL SUBY(I,J,LX,NY,Z,SUM,NEQ)
	SUM=SUM+Z(I,J,LX)*GB(LX)*(DRAT*(NC-1))**BETA
        XK(I,J,LX)=-SUM*DT
301     CONTINUE

!       4th order
        I=-NC+2
        DO 501 LX=1,NEQ
        SUM=0.D0
        DO 50 LY=1,NEQ
        DO 50 IK=-2,2
        SUM=SUM+COEFF(IK,4,1)*Z(I+IK,J,LY)*E(LX,LY)
50      CONTINUE
        CALL SUBY(I,J,LX,NY,Z,SUM,NEQ)
	SUM=SUM+Z(I,J,LX)*GB(LX)*(DRAT*(NC-2))**BETA
        XK(I,J,LX)=-SUM*DT
501     CONTINUE

!       6th order
        I=-NC+3
        DO 401 LX=1,NEQ
        SUM=0.D0
        DO 40 LY=1,NEQ
        DO 40 IK=-3,3
        SUM=SUM+COEFF(IK,6,1)*Z(I+IK,J,LY)*E(LX,LY)
40      CONTINUE
        CALL SUBY(I,J,LX,NY,Z,SUM,NEQ)
	SUM=SUM+Z(I,J,LX)*GB(LX)*(DRAT*(NC-3))**BETA
        XK(I,J,LX)=-SUM*DT
401     CONTINUE

        DO 601 I=-NC+4,-1
        JO=-I
        DO 2071 LX=1,NEQ
        SUM=0.D0
        DO 207 LY=1,NEQ
        DO 207 IK=-4,4
        SUM=SUM+COEFF(IK,8,1)*Z(I+IK,J,LY)*E(LX,LY)
207     CONTINUE
        CALL SUBY(I,J,LX,NY,Z,SUM,NEQ)
        SUM=SUM+Z(I,J,LX)*HB(LX)*(DRAT*JO)**BETA
        XK(I,J,LX)=-SUM*DT
2071    CONTINUE
601	CONTINUE

        RETURN
        END




	SUBROUTINE SUBY(I,J,LX,NY,Z,SUM,NEQ)
	use param_info, only : MX, LEQ, MY, NCBUF
	use EQUATON_GEQAT
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        integer, parameter :: NA = 12, KLEVEL = 12
        integer, parameter :: MXMAX = 2*MX+NCBUF, MYMAX = MY+NCBUF
	DOUBLE PRECISION Z(-NCBUF:MXMAX,0:MYMAX,LEQ)
        DOUBLE PRECISION COEFF(-NA:NA,NA,2)
        COMMON /FINITE/COEFF
	COMMON /ACCURCY/NORDER,ISIDE
	COMMON /BUFFER/NC

        IF(J >= 5 .AND. J <= NY)THEN
        DO 20 LY=1,NEQ
        DO 20 IK=-ISIDE,ISIDE
        SUM=SUM+COEFF(IK,NORDER,2)*Z(I,J+IK,LY)*F(LX,LY)
20      CONTINUE


        ELSE IF(J > NY)THEN
	CALL BUFFY(I,J,LX,NY,Z,SUM,NEQ)

!	doing the wall
	ELSE IF(J == 0)THEN
        DO 21 LY=1,NEQ
        DO 21 IK=0,6
        SUM=SUM+COEFF(IK,7,2)*Z(I,J+IK,LY)*F(LX,LY)
21      CONTINUE

	ELSE IF(J == 1)THEN
        DO 23 LY=1,NEQ
        DO 23 IK=-1,1
        SUM=SUM+COEFF(IK,2,2)*Z(I,J+IK,LY)*F(LX,LY)
23      CONTINUE

	ELSE IF(J == 2)THEN

        DO 25 LY=1,NEQ
        DO 25 IK=-2,2
        SUM=SUM+COEFF(IK,4,2)*Z(I,J+IK,LY)*F(LX,LY)
25      CONTINUE

	ELSE IF(J == 3)THEN

        DO 27 LY=1,NEQ
        DO 27 IK=-3,3
        SUM=SUM+COEFF(IK,6,2)*Z(I,J+IK,LY)*F(LX,LY)
27      CONTINUE

        ELSE IF(J == 4)THEN
        DO 22 LY=1,NEQ
        DO 22 IK=-4,4
        SUM=SUM+COEFF(IK,8,2)*Z(I,J+IK,LY)*F(LX,LY)
22      CONTINUE
	
	END IF

	RETURN              
	END


        SUBROUTINE BUFFY(I,J,LX,NY,Z,SUM,NEQ)
	use param_info, only : MX, LEQ, MY, NCBUF
	use EQUATON_GEQAT
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        integer, parameter :: NA = 12, KLEVEL = 12
        integer, parameter :: MXMAX = 2*MX+NCBUF, MYMAX = MY+NCBUF
        DOUBLE PRECISION Z(-NCBUF:MXMAX,0:MYMAX,LEQ)
        DOUBLE PRECISION COEFF(-NA:NA,NA,2)
        COMMON /FINITE/COEFF
        COMMON /ACCURCY/NORDER,ISIDE
	COMMON /BUFFER/NC

	IF(J >= NY+1 .AND. J < NY+NC-3)THEN
        DO LY=1,NEQ
        DO IK=-4,4
        SUM=SUM+COEFF(IK,8,2)*Z(I,J+IK,LY)*F(LX,LY)
	END DO
	END DO

        ELSE IF(J == NY+NC-3)THEN

        DO LY=1,NEQ
        DO IK=-4,4
        SUM=SUM+COEFF(IK,8,2)*Z(I,J+IK,LY)*F(LX,LY)
	END DO
	END DO

        ELSE IF(J == NY+NC-2)THEN

        DO LY=1,NEQ
        DO IK=-3,3
        SUM=SUM+COEFF(IK,6,2)*Z(I,J+IK,LY)*F(LX,LY)
	END DO
	END DO

        ELSE IF(J == NY+NC-1)THEN

        DO LY=1,NEQ
        DO IK=-2,2
        SUM=SUM+COEFF(IK,4,2)*Z(I,J+IK,LY)*F(LX,LY)
	END DO
	END DO

        ELSE IF(J == NY+NC)THEN

        DO LY=1,NEQ
        DO IK=1,-7,-1
        SUM=SUM-COEFF(-IK,1,2)*Z(I,J+IK,LY)*F(LX,LY)
	END DO
	END DO
	END IF

	RETURN
	END


 
        SUBROUTINE DIFFER(H,IX,ISTAGE,IEX,ITAM)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        integer, parameter :: NA = 12, KLEVEL = 12
        DOUBLE PRECISION AICH(NA)
        DOUBLE PRECISION COEFF(-NA:NA,NA,2)
        DOUBLE PRECISION, DIMENSION(KLEVEL) :: BETA, CEES
        COMMON /FINITE/COEFF
        COMMON /LDDRK/BETA
	COMMON /BUFFER/NC

	IF(IX /= 1)GOTO 999

	BETA(1:KLEVEL)=0.D0
	CEES(1:KLEVEL)=0.D0

!	exact runge kutta......
	CEES(1)=1.D0
	CEES(2)=0.5D0

	IF(IEX == 0 .AND. ISTAGE == 4)THEN
	CEES(3)=1.D0/6.D0
	CEES(4)=1.D0/24.D0

	ELSE IF(IEX == 0 .AND. ISTAGE == 5)THEN
	CEES(3)=1.D0/6.D0
	CEES(4)=1.D0/24.D0
	CEES(5)=1.D0/120.D0
	ELSE IF(IEX == 0 .AND. ISTAGE == 6)THEN
	CEES(3)=1.D0/6.D0
	CEES(4)=1.D0/24.D0
	CEES(5)=1.D0/120.D0
	CEES(6)=1.D0/720.D0
	ELSE IF(IEX == 0 .AND. ISTAGE == 7)THEN
	CEES(3)=1.D0/6.D0
	CEES(4)=1.D0/24.D0
	CEES(5)=1.D0/120.D0
	CEES(6)=1.D0/720.D0
	CEES(7)=1.D0/5040.D0
	ELSE IF(IEX == 0 .AND. ISTAGE == 8)THEN
	CEES(3)=1.D0/6.D0
	CEES(4)=1.D0/24.D0
	CEES(5)=1.D0/120.D0
	CEES(6)=1.D0/720.D0
	CEES(7)=1.D0/5040.D0
	CEES(8)=1.D0/40320.D0
        ELSE IF(IEX == 0 .AND. ISTAGE == 9)THEN
        CEES(3)=1.D0/6.D0
        CEES(4)=1.D0/24.D0
        CEES(5)=1.D0/120.D0
        CEES(6)=1.D0/720.D0
        CEES(7)=1.D0/5040.D0
        CEES(8)=1.D0/40320.D0
	ELSE IF(IEX == 0 .AND. ISTAGE == 10)THEN
	CEES(3)=1.D0/6.D0
	CEES(4)=1.D0/24.D0
	CEES(5)=1.D0/120.D0
	CEES(6)=1.D0/720.D0
	CEES(7)=1.D0/5040.D0
	CEES(8)=1.D0/40320.D0
	CEES(9)=1.D0/362880.D0
	CEES(10)=1.D0/3628800.D0

	ELSE IF(IEX == 1 .AND. ISTAGE == 4)THEN
!	4 stage LDDRK
	CEES(3)=0.162997D0
	CEES(4)=0.0407574D0

	ELSE IF(IEX == 1 .AND. ISTAGE == 5)THEN
!	5 stage LDDRK
	CEES(3)=0.166558D0
	CEES(4)=0.0395041D0
	CEES(5)=0.00781071D0

	ELSE IF(IEX == 1 .AND. ISTAGE == 6)THEN
!	6 stage LDDRK
	CEES(3)=1.D0/6.D0
	CEES(4)=1.D0/24.D0
	CEES(5)=0.00781005D0
	CEES(6)=0.00132141D0
	ELSE
        CEES(3)=1.D0/6.D0
        CEES(4)=1.D0/24.D0
	ISTAGE=4
	END IF
	

	KORD=ISTAGE

	BETA(1)=0.D0
	BETA(KORD)=CEES(2)
	IB=3
	DO 20 I=KORD-1,2,-1
	BETA(I)=CEES(IB)
	DO IK=I+1,KORD
	BETA(I)=BETA(I)/BETA(IK)
	END DO
	IB=IB+1
20	CONTINUE

999	CONTINUE

        AICH(2)=2.D0*H
        AICH(4)=12.D0*H
        AICH(6)=60.D0*H
        AICH(8)=840.D0*H
	AICH(10)=20160.D0*H
	
	HIN=1.D0/H

        AICH(1)=H
        AICH(3)=H
        AICH(5)=H
        AICH(7)=H
        AICH(9)=H

        DO I=-NA,NA
        COEFF(I,1:NA,IX)=0.D0
	END DO

!	diferencing or NORDER 2 4 6 8 in main part
        COEFF(1,2,IX)=1.D0/AICH(2)
        COEFF(-1,2,IX)=-1.D0/AICH(2)

        COEFF(1,4,IX)=8.D0/AICH(4)
        COEFF(-1,4,IX)=-8.D0/AICH(4)
        COEFF(2,4,IX)=-1.D0/AICH(4)
        COEFF(-2,4,IX)=1.D0/AICH(4)

	IF(ITAM == 0)THEN
        COEFF(1,6,IX)=45.D0/AICH(6)
        COEFF(-1,6,IX)=-45.D0/AICH(6)
        COEFF(2,6,IX)=-9.D0/AICH(6)
        COEFF(-2,6,IX)=9.D0/AICH(6)
        COEFF(3,6,IX)=1.D0/AICH(6)
        COEFF(-3,6,IX)=-1.D0/AICH(6)
	ELSE
!	TAM WEBB scheme
!	not right...!
	COEFF(1,6,IX)=0.79926643D0
	COEFF(-1,6,IX)=-0.79926643D0
	COEFF(2,6,IX)=-0.18941314D0
	COEFF(-2,6,IX)=0.18941314D0
	COEFF(3,6,IX)=0.02651995D0
	COEFF(-3,6,IX)=-0.02651995D0
	END IF

        COEFF(1,8,IX)=672.D0/AICH(8)
        COEFF(-1,8,IX)=-672.D0/AICH(8)
        COEFF(2,8,IX)=-168.D0/AICH(8)
        COEFF(-2,8,IX)=168.D0/AICH(8)
        COEFF(3,8,IX)=32.D0/AICH(8)
        COEFF(-3,8,IX)=-32.D0/AICH(8)
        COEFF(4,8,IX)=-3.D0/AICH(8)
        COEFF(-4,8,IX)=3.D0/AICH(8)

	COEFF(1,10,IX)=17178.D0/AICH(10)
	COEFF(-1,10,IX)=-17178.D0/AICH(10)
	COEFF(2,10,IX)=-5232.D0/AICH(10)
	COEFF(-2,10,IX)=5232.D0/AICH(10)
	COEFF(3,10,IX)=1443.D0/AICH(10)
	COEFF(-3,10,IX)=-1443.D0/AICH(10)
	COEFF(4,10,IX)=-272.D0/AICH(10)
	COEFF(-4,10,IX)=272.D0/AICH(10)
	COEFF(5,10,IX)=25.D0/AICH(10)
	COEFF(-5,10,IX)=-25.D0/AICH(10)


!	forward differencin ogf 8 orfder at inflow

        COEFF(-1,1,IX)=-1.D0*HIN/8.D0
        COEFF(0,1,IX)=-223.D0*HIN/140.D0
        COEFF(1,1,IX)=3.5D0*HIN
        COEFF(2,1,IX)=-3.5D0*HIN
        COEFF(3,1,IX)=35.D0*HIN/12.D0
        COEFF(4,1,IX)=-7.D0*HIN/4.D0
        COEFF(5,1,IX)=0.7D0*HIN
        COEFF(6,1,IX)=-1.D0*HIN/6.D0
        COEFF(7,1,IX)=1.D0*HIN/56.D0


        COEFF(-2,3,IX)=1.D0*HIN/56.D0
        COEFF(-1,3,IX)=-2.D0*HIN/7.D0
        COEFF(0,3,IX)=-19.D0*HIN/20.D0
        COEFF(1,3,IX)=2.D0*HIN
        COEFF(2,3,IX)=-5.D0*HIN/4.D0
        COEFF(3,3,IX)=2.D0*HIN/3.D0
        COEFF(4,3,IX)=-0.25D0*HIN
        COEFF(5,3,IX)=2.D0*HIN/35.D0
        COEFF(6,3,IX)=-1.D0*HIN/168.D0

        COEFF(-3,5,IX)=-1.D0*HIN/168.D0
        COEFF(-2,5,IX)=1.D0*HIN/14.D0
        COEFF(-1,5,IX)=-0.5D0*HIN
        COEFF(0,5,IX)=-9.D0*HIN/20.D0
        COEFF(1,5,IX)=1.25D0*HIN
        COEFF(2,5,IX)=-0.5D0*HIN
        COEFF(3,5,IX)=1.D0*HIN/6.D0
        COEFF(4,5,IX)=-1.D0*HIN/28.D0
        COEFF(5,5,IX)=1.D0*HIN/280.D0

!	forward differencin ogf 6th order at wall...
        COEFF(0,7,IX)=-147.D0*HIN/60.D0
        COEFF(1,7,IX)=360.D0*HIN/60.D0
        COEFF(2,7,IX)=-450.D0*HIN/60.D0
        COEFF(3,7,IX)=400.D0*HIN/60.D0
        COEFF(4,7,IX)=-225.D0*HIN/60.D0
        COEFF(5,7,IX)=72.D0*HIN/60.D0
        COEFF(6,7,IX)=-10.D0*HIN/60.D0

        RETURN
        END

