!!!*************************************************************
! 文件/File: frxx16.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:17:01
!*************************************************************

!***********************************************************************
! 
!    Copyright 2018, I.J. Thompson
!
!    This file is part of FRESCOX.
!
!    FRESCOX is free software; you can redistribute it and/or
!    modify it under the terms of the GNU General Public License
!    as published by the Free Software Foundation; either version 2
!    of the License, or (at your option) any later version.
!    
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!    
!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, 
!    MA  02110-1301, USA.
!
!    OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC
!    LICENSE
!
!    The precise terms and conditions for copying, distribution and
!    modification are contained in the file COPYING.
!
!***********************************************************************
***readcp**************************************************************
      subroutine readcp(ki,icto,icfrom,kind,ip1,ip2,ip3,ip4,ip5,
     X                  p1,p2,jmax,rmax,infile)
	real*8 jmax,p1,p2,rmax
	integer infile
	namelist/coupling/icto,icfrom,kind,ip1,ip2,ip3,ip4,ip5,
     X			  p1,p2,jmax,rmax,kfrag,kcore,nforms,infile
      	icto=0;icfrom=0;kind=0;ip1=0;ip2=0;ip3=0;ip4=-1;ip5=-1
	p1=0;p2=0;jmax=0;rmax=0;kfrag=0;kcore=0;nforms=0;infile=4
	read(KI,nml=coupling)
	if(kfrag>0) p1 = kfrag ; if(kcore>0) p2 = kcore
	return
	end
      SUBROUTINE SUMX(SMAT,NCH,JIN,EL,JCOEF,XS,SIGT,SIGJ,SIGR,PART,
     X           EXCIT,AJUMP,FUSL,CSIG,JEX,ITC,K,RMASS,PEL,EXL,LVAL,
     X         NSA,NJA,JVAL,JPROJ,JTARG,JTOTAL,TOTFUS,FUSUM,
     X         CORFUS,SIGEL,SIGTOT,SMATS,ITCM,LXS,IF1,IF2,FRACFUS)
	use io
	use parameters
	use searchpar, only: final
	use parallel, only: mpisets
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(LXSEC=41)
      INTEGER PEL,EXL,C,ITC(MXP,MXX),EL,PART(NCH),EXCIT(NCH),LVAL(NCH),
     X        SMATS
      COMPLEX*16 SMAT(NCH)
      REAL*8 XS(MAXCH),SIGJ(0:MXPEX),SIGR(MXP,MXX),
     X       K(MXP,MXX),CSIG(LMAX1,MXPEX),SIGTOT(3),SIGEL(3)
      REAL*8 JCOEF,TOTFUS(3),CORFUS(3,NFUS1),FUSL(1+NFUS),SIGT(3)
      REAL*8 RMASS(MXP),JEX(6,MXP,MXX),JTOTAL,MSA,MJA
      REAL*8 JVAL(MAXCH),JPROJ(MAXCH),JTARG(MAXCH)
      COMPLEX*16 FUSUM(MAXCH,NSA*NJA),C6,S
C
      DATA  Z / 0D0 /
C
      SIGINEL = 0.0
      SIGTRAN = 0.0
      IFT = 1
      FUSL(1) = 0.0
      FRACFUS = 0.0
      SMT = 2.*(1. - REAL(SMAT(EL)))
      SMEL = abs(1. - SMAT(EL))**2
      DO 690 C=1,NCH
         IC = PART(C)
         IA = EXCIT(C)
         IT = ITC(IC,IA)
        IF(DBLE(K(IC,IA)**2) .LT. 0.0) GOTO 690
      AMDSQS = ABS(SMAT(C))**2
      IF(C.EQ.EL) AMDSQS = 1 - AMDSQS
!@@      
         IF(RMASS(IC).lt.1e-5) then
             S =          RMASS(PEL)*amu/ (HBC*K(PEL,EXL))
!@@@         S = K(IC,IA)*RMASS(PEL)*amu/ (HBC*K(PEL,EXL))
	 ELSE IF(RMASS(PEL).lt.1e-5) then
          S = HBC*K(IC,IA)/(RMASS(IC)*amu)
         ELSE
             S = K(IC,IA)/RMASS(IC)/(K(PEL,EXL)/RMASS(PEL))
         ENDIF
!      write(6,*) 'C,JCOEF,S =',C,real(JCOEF),real(S)
!        exactly the same changes as frxx3.f
!
!        S = K(IC,IA)/RMASS(IC)/(K(PEL,EXL)/RMASS(PEL))
!@@
      XS(C) = JCOEF* AMDSQS * S
      X = XS(C)*AJUMP
      IF(C.EQ.EL) THEN
         SIGT(1) = SIGT(1) + X
         SIGT(2) = SIGT(2) + X*JTOTAL
         SIGT(3) = SIGT(3) + X*JTOTAL**2
           XT = JCOEF * SMT * AJUMP
         SIGTOT(1) = SIGTOT(1) + XT
         SIGTOT(2) = SIGTOT(2) + XT*JTOTAL
         SIGTOT(3) = SIGTOT(3) + XT*JTOTAL**2
           XE = JCOEF * SMEL * AJUMP
         SIGEL(1) = SIGEL(1) + XE
         SIGEL(2) = SIGEL(2) + XE*JTOTAL
         SIGEL(3) = SIGEL(3) + XE*JTOTAL**2
         SIGJ(0) = SIGJ(0) + XS(C)
         FUSL(1) = FUSL(1) + XS(C)
        ELSE
         SIGR(IC,IA) = SIGR(IC,IA) + X
         SIGJ(IT) = SIGJ(IT) + XS(C)
         FUSL(1) = FUSL(1) - XS(C)
         IF(IC.EQ.PEL) SIGINEL = SIGINEL + XS(C)
         IF(IC.NE.PEL) SIGTRAN = SIGTRAN + XS(C)
        ENDIF
      if( PART(C)==PEL.and.EXCIT(C)==EXL) then
         if(C==EL) FRACFUS = FRACFUS + AMDSQS * S
         if(C/=EL) FRACFUS = FRACFUS - AMDSQS * S
      else
         FRACFUS = FRACFUS - AMDSQS * S   ! CN production excluding direct products
      endif
      IFT = IF2
      if(final) then
C      ................CALCULATE SUMS FOR FUSION POLARISATIONS
        S = SMAT(C)*SQRT(S)*EXP((0.,1.)*CSIG(LVAL(EL)+1,ITC(PEL,EXL)))
        IAM = 0
        DO 685 I=1,NSA
           MSA = I-1 - JEX(1,PEL,EXL)
        DO 685 J=1,NJA
           MJA = J-1 - JEX(2,PEL,EXL)
           IAM = IAM + 1
        C6 = S  * CLEB6(JVAL(EL),-MSA,JPROJ(EL),MSA,LVAL(EL)+Z,Z)
     X          * CLEB6(JTOTAL,-(MSA+MJA),JTARG(EL),MJA,JVAL(EL),-MSA)
         FUSUM(C,IAM) = FUSUM(C,IAM) + C6
  685    continue
        endif
  690 CONTINUE
C
      X = FUSL(1)*AJUMP
      TOTFUS(1) = TOTFUS(1) + X
      TOTFUS(2) = TOTFUS(2) + X*JTOTAL
      TOTFUS(3) = TOTFUS(3) + X*JTOTAL**2
      DO IT=1,NFUS
      X = FUSL(1+IT)*AJUMP
      CORFUS(1,IT) = CORFUS(1,IT) + X
      CORFUS(2,IT) = CORFUS(2,IT) + X * JTOTAL
      CORFUS(3,IT) = CORFUS(3,IT) + X * JTOTAL**2
      ENDDO
      IF(LXS.LT.0) WRITE(LXSEC,100) JTOTAL,LVAL(EL),SIGJ(0),SIGINEL,
     X                              SIGTRAN,FUSL(1)
      IF(LXS.LT.0) written(LXSEC) = .true.
 100  FORMAT(F6.1,I6,1P,4E10.3,' (SIG FOR J/L)')
	if(mpisets.gt.1) 
     X WRITE(48,1450) JTOTAL,' ',MOD(LVAL(EL),100),JIN
     X       , (SIGJ(IT),' ',IT=0,ITCM),(FUSL(I),' ',I=IF1,IFT)
 1450 FORMAT(' Reaction Xsec',F8.1,A1,'/',I2,' @',I2,' =',F10.3,A1,',',
     X   ' Out:', 9(F8.3,A1),:,/(11X,'Xsec',23X,10(F8.3,A1)))
C
      IF(SMATS==2.and.final)
     X  WRITE(7,719) SMAT(EL),LVAL(EL),JVAL(EL),JTOTAL
      IF(SMATS==3.and.final) then
	do C=1,NCH
	if(PART(C)==PEL.and.EXCIT(C)==EXL) 
     X    WRITE(7,750) SMAT(C),LVAL(C),JVAL(C),JTOTAL,PART(C),EXCIT(C),
     X    LVAL(EL),JVAL(EL)
	enddo
	endif
      IF(SMATS.GE.4.and.final) THEN
      DO 720 C=1,NCH
        IF(ABS(SMAT(C)).GT.0E-10)
     X WRITE(7,750) SMAT(C),LVAL(C),JVAL(C),JTOTAL,PART(C),EXCIT(C),
     X    LVAL(EL),JVAL(EL)
     *  ,jproj(c),jtarg(c),jproj(el),jtarg(el)
     *  ,C,el   ! added by neil summers  Nov 2009
     *  ,K(PART(C),EXCIT(C))/RMASS(PART(C))/(K(PEL,EXL)/RMASS(PEL))

!	write(79,750) SMAT(C),LVAL(C),JVAL(C),JTOTAL,PART(C),EXCIT(C),
!     X    LVAL(EL),JVAL(EL) ,jproj(c),jtarg(c),jproj(el),jtarg(el)
!	call flush(79)
!	write(80,750) SMAT(C),LVAL(C),JVAL(C),JTOTAL,PART(C),EXCIT(C),
!     X    LVAL(EL),JVAL(EL) ,jproj(c),jtarg(c),jproj(el),jtarg(el)
!      IC = PART(C)
!      IA = EXCIT(C)
!      if(JPROJ(C)/= JEX(1,IC,IA).or.JTARG(C)/= JEX(2,IC,IA)) then
!	write(78,750) SMAT(C),LVAL(C),JVAL(C),JTOTAL,PART(C),EXCIT(C),
!     X    LVAL(EL),JVAL(EL) ,jproj(c),jtarg(c)
!      endif
720    CONTINUE
        ENDIF
	if(SMATS.ge.2.and.final) written(7) = .true.
!	if(SMATS.ge.2) call flush(7)
 719  FORMAT(2F15.10,I6,2F6.1,' : S(L,J,JT)')
C750  FORMAT(2F15.10,I6,2F6.1,I6,F6.1,' S(L,J,JT,L-IN,J-IN)')
C750  FORMAT(2F15.10,I6,2F6.1,2I3,I6,F6.1,' S-MAT')
C750  FORMAT(2F15.10,I6,2F6.1,2I3,I6,F6.1,4F4.1)
 750  FORMAT(2F15.10,I6,2F6.1,2I3,I6,F6.1,4F5.1
     x       ,2I4,f12.8        ! Neil Summers
     x      )   

C		Print out cross section information in file 13
      RETURN
      END
	
****READIN**************************************************************
	SUBROUTINE READIN(WOUT,WOUT0,NF0,NF,STREN,ISNONO,LPMAX,MIXPOT,
     X     NAME,MASS,NEX,QVAL,RMASS,HP,COPY,EXTRA,IOFAM,GIVEXS,
     X     JEX,CPOT,ENEX,BAND,ITC,FLIP,NCHAN,ITCM,PSIGN,BHEAD,
     X     PWFLAG,FORMF,FORML,FORMC,PTYPE,CPSO,EPS,MMXQRN,FILE,ffreal,
     X  HPOT,LAMBDA,MATRIX,MEK,NIX ,NSP, QNF,D0,BE,AFRAC,CCFRAC,BPHASE,
     x  ICTO,ICFROM,KIND,QQ,IREM,KPCORE,BETAR,BETAI,QCM,JMAX,RMAX,REV,
     x  NFI,TWOW,OPN,NPRIOR,NPOST,REW14,LOCAL,NIB,LTRANS,FORMDEF,
     x  MLCALL,BPROJ,XA,XB,XP,XQ,SPINTR,FPT,NKP,VARYL,LSHAPE,LDEP,
     x  NLL,LOCF,KLT,COUPLE,ICOM,ICOR,NCP,IEXCH,GPT,POTCAP,INFILE,
     x  NSA, NJA, XCOEF, DISCIT, DISC8, LAMAX,NBINS,EMID,KMINX,NKBIN,
     x  NORBIN,STOP19,WID,CENTR,NPWCOUP,PWCOUP,JPWCOUP,BSMAT,RMAS)
	use parameters
	use factorials
	use drier
	use kcom
	use trace
	use io
	use fresco1
	use searchdata, only: energy_list,num_energies,pel_list
	use searchpar, only: final
	IMPLICIT NONE
C
C    FORM FACTORS AND THEIR PARAMETERS
C    ---------------------------------
      COMPLEX*16 FORMF(MAXM,MLOC),FORMC(MAXNLC,MSP,2),
     x           BSMAT(NCHBINMAX,NCHBINMAX,max(1,NKMAX),MSP)
      REAL*8 FORML(MAXNLR,MSP,2),VIMX,VARYL(2,MLOC),EMID(MSP),
     X		KMINX(2,MSP),VRMX,FORMDEF(12,MLOC),DELTAE(MSP),
     x          PWCOUP(3,MPWCOUP),RMAS(MSP)

C
      INTEGER NEX(MXP+1),BAND(2,MXP,MXX),COPY(2,MXP,MXX,2),
     X       ITC(MXP,MXX),CPOT(MXP,MXX),PTYPE(12,MLOC),NORBIN(MSP),
     x       LPMAX(MXP)
      INTEGER QNF(19,MSP),BHEAD(2,2,0:9),LSHAPE(MLOC),NKBIN(2,MSP)
      INTEGER CP,ICTO(MAXCPL+1),ICFROM(MAXCPL+1),KIND(MAXCPL+1),
     X        LOCF(MAXCPL),FPT(7,MAXQRN,MAXCPL),NKP(2,MAXCPL),
     X        IREM(MAXCPL+1),KPCORE(MAXCPL+1),ICOR(MAXCPL,2),
     X        GPT(2,MAXQRN,MAXCPL),FILE(MAXCPL),NFI(3),
     X        NLL(MAXCPL),KLT(MAXCPL),NIB(MAXCPL),NPWCOUP,
     X        QQ(MAXCPL+1),ICOM(MAXCPL,2),NBINS,JPWCOUP(9,MPWCOUP),
     X        QCM(MAXCPL+1),LAM(MAXCPL+1),INFILE(MAXCPL+1),
     x        MIXPOT(MXP)
C
C    CONTROL VARIABLES
C    -----------------
      LOGICAL GIVEXS(MXP),EXTRA(2,MXP,MXX),FLIP
      LOGICAL LOCAL(MAXCPL),REV(MAXCPL),USED(2,MSP),COUPLE(MAXCPL),
     X        NPRIOR,NPOST,OPN(2),ffreal(MAXCPL),REW14,WOUT,WOUT0
      LOGICAL LTRANS(MAXQRN,MAXCPL),MLCALL(2,MAXCPL),HASO(KPMAX)
      LOGICAL DISC8,DISCIT,PWFLAG(MXP),CPSO(MAXCPL),STOP19
	integer NF,NF0,LAMAX,KN,NSP,NCP,IEXCH,NSA,NJA,LDEP(MLOC)
        real*8 XCOEF,JAP,R
	integer NIX,NCHAN,IC,ITCM,L,JF,J,IN,I,IITER,IOFAM(2,MXP,MXX)
C
C    DEFINING THE MASS PARTITIONS AND THEIR EXCITED STATES
C    -----------------------------------------------------
      REAL*8 MASS(4,MXP+1),RMASS(MXP),HP(MXP),JEX(6,MXP,MXX)
      REAL*8 ENEX(2,MXP,MXX),QVAL(MXP+1),CONV
      REAL*8 BE(MSP,5),D0(MSP,2),BPROJ(MSP,MXP),HPOT(MLOC)
      INTEGER NCHPMAX(MXP)
C
C    COUPLINGS
C    ---------
      REAL*8 BETAR(MAXCPL+1),BETAI(MAXCPL+1),JMAX(MAXCPL+1),
     X	RMAX(MAXCPL+1),   BPHASE(2,2,max(1,NKMAX),MSP),
     x  WID(MAXCPL+1),CENTR(MAXCPL+1),
     X  XA(MAXCPL),XB(MAXCPL),XP(MAXCPL),XQ(MAXCPL),STREN(MLOC)
      REAL*8 AFRAC(MXPEX,MXPEX,2,MSP),MEK(MPAIR),SPINTR(2,MPAIR)
      REAL*8 CCFRAC(MXPEX,2,MXPEX)
      INTEGER MATRIX(6,MPAIR),LAMBDA(MLOC),MMXQRN,POTCAP(MLOC,2,MAXCPL)
	LOGICAL TWOW,ISNONO
C
    	character*70 TMP,TMPN
     	common/cfs/ TMP,TMPN  
      real*8 EPS
C
      CHARACTER*8 NAME(2,MXP+1)
      CHARACTER*3 WAYS(2),W
      CHARACTER*1 PSIGN(3)
      DATA WAYS/'ONE','TWO'/
C
      CALL CHECK(N,MAXN,7)
      CALL CHECK(MINT,MAXM,7)
      CALL CHECK(NLN,MAXNLN,8)
      CALL CHECK(INT(JTMAX),LMAX,11)

C
	if(melfil.eq.1) then
            open (53,access = 'sequential',
     x         form='unformatted',status = 'unknown')
            open (54,access = 'sequential',
     x         form='formatted',status = 'unknown')
	endif

      WOUT0 = VEFF.NE.0 .OR. MOD(WAVES,2).NE.0 .OR. WDISK.NE.0
      WOUT = WOUT0 .OR. NFUS.NE.0 
	MLCALL(:,:) = .false.

	if(final) then
	if(SMATS>=2) then
	   open(45,form='formatted',access='sequential')
	   rewind 45
	   endif
	if(SMATS>=2) then
	   open(38,form='formatted',access='sequential')
	   rewind 38
	   endif
	open(39,form='formatted',access='sequential'); rewind 39
	open(56,form='formatted',access='sequential'); rewind 56
	open(40,form='formatted',access='sequential',recl=535);rewind 40
	open(13,form='formatted',access='sequential'); rewind 13
	endif
	HASO(:) = .false.

C	NOW READ PRE-DIGESTED INPUT IN FILE 3:
C
C    DEFINING THE MASS PARTITIONS AND THEIR EXCITED STATES
C    -----------------------------------------------------
      CALL PARTEX(NAME,MASS,NEX,QVAL,RMASS,HCM,HP,COPY,EXTRA,LPMAX,
     X     GIVEXS,JEX,CPOT,ENEX,BAND,ITC,FLIP,NCHAN,ITCM,PSIGN,BHEAD,
     X     NCHPMAX,PWFLAG,FCWFN,IOFAM,MIXPOT)
C
C    ENTER THE VARIOUS TYPES OF POTENTIALS
C    -------------------------------------
      NF0 = 0
	STREN(:)=0d0
   

      CALL POTENT(FORMF,NF0,PTYPE,MAXM,HCM,TRENEG,MR,M,HASO,
     X            HP,HPOT,STREN,LAMBDA,NCHPMAX,VARYL,LDEP,LSHAPE,
     X            MATRIX,MEK,NIX, JEX,MASS,CPOT,NEX,NCHAN,COPY,STOP19,
     x 		  FORMDEF,lanecoup)

!  debugging only!:
!	STREN(:)=0d0
!	print *,' ALL ASYMPTOTIC STRENGTHS == ZERO!!!!'
C
C
C    ONE AND TWO - PARTICLE FORM FACTORS AND THEIR PARAMETERS
C    --------------------------------------------------------
C
      CALL INFORM(FORML,FORMC,NSP,FORMF,NF0,PTYPE,TRENEG,QVAL,BPHASE,
     X  QNF,D0,BE,NCHAN,NEX,ENEX,JEX,BAND,MASS,AFRAC,CCFRAC,ITC,COPY,N,
     X  NLN,MINT,NNN,RNN,RMIN,HP(1),RIN,MR,NNU,ERANGE,DK,PI,NAME,
     X	NBINS,EMID,DELTAE,KMINX,NKBIN,NORBIN,BSMAT,RMAS)
        do i=43,44
        if(written(i)) close(i)
        enddo
C
C
C    COUPLINGS
C    ---------
      CP = 1
      NF = NF0
      NFI(1) = 1
      NFI(3) = 1
      TWOW = .false.
      ISNONO = .false.
      LREC66 = 0
      OPN(1) = .FALSE.
      OPN(2) = .FALSE.
      NPRIOR =.FALSE.
      NPOST = .FALSE.
      REW14 = .FALSE.
      BPROJ(:,:) = 0.
      WID(:) = 0.; CENTR(:) = 0.
	IITER=ITER
	if(IBLOCK<0) IITER=1
	MCALLS = .false.
	ALLPART = .false.
   90 continue
!      READ(KI,1220) ICTO(CP),ICFROM(CP),KIND(CP),QQ(CP),IREM(CP),
!     X         KPCORE(CP),BETAR(CP),BETAI(CP),JMAX(CP),RMAX(CP)
      CALL readcp(KI,ICTO(CP),ICFROM(CP),KIND(CP),
     X         QQ(CP),IREM(CP),KPCORE(CP),QCM(CP),LAM(CP),
     X         BETAR(CP),BETAI(CP),JMAX(CP),RMAX(CP),INFILE(CP))
!      WRITE(KO,1220) ICTO(CP),ICFROM(CP),KIND(CP),QQ(CP),IREM(CP),
!     X         KPCORE(CP),BETAR(CP),BETAI(CP),JMAX(CP),RMAX(CP)
!1220 FORMAT(3I4,3I2,2F8.2,2F4.1)
  100 IF(ICTO(CP).EQ.0.OR.ICFROM(CP).EQ.0) GO TO 120
!!      if(KIND(CP)>8) KIND(CP) = KIND(CP)-8   ! temporary compatibility

      CALL CHECK(CP,MAXCPL,4)
      LOCAL(CP) = KIND(CP).eq.9.and.QQ(CP).eq.0 .or.    ! non-local partial-wave couplings
     x            KIND(CP).LE.6.and..not.
     x                (KIND(CP)==1.and.QQ(CP)==1 .or.   ! non-local general transitions
     x                 KIND(CP)==2.and.QCM(CP)>0 )  	! semi-direct captures
      REV(CP)   = ICTO(CP).GT.0
	TWOW = TWOW .or. (REV(CP).and.KIND(CP).ge.5)
	ISNONO = ISNONO .or. KIND(CP).eq.8 .or. KIND(CP)>=11
      W = WAYS(1)
      IF(REV(CP)) W = WAYS(2)
      IF(JMAX(CP).lt..01) JMAX(CP) = JTMAX + 0.51
      IF(RMAX(CP).lt..01) RMAX(CP) = abs(RMATCH) - 3*HCM
!      IF(RMAX(CP).lt..01) RMAX(CP) = abs(RMATCH) 
      NLL(CP) = MIN(NLN, NINT(RMAX(CP)/RINTP)+1 )
      LOCF(CP) = 0
      KLT(CP) = 1
      WRITE(KO,1230) W,CP,ICTO(CP),ICFROM(CP),KIND(CP),
     X           QQ(CP),IREM(CP),KPCORE(CP),QCM(CP),LAM(CP),
     X           BETAR(CP),BETAI(CP),JMAX(CP),RMAX(CP)
 1230 FORMAT(/' ',132('*')//' ',A3,'-way',
     X        ' COUPLING #',I2,' for partitions ',I2,' <- ',I2,
     X' of KIND',I3,', '               ,5I3,
     X' & P1,P2 =',F10.4,F9.4,' : for J <=',F8.1,' & R <',
     X F6.1,' fm.'/)
       CALL FLUSH(KO)
      ICTO(CP)  = ABS(ICTO(CP))
      IC = ICTO(CP)
      COUPLE(CP) = IC.LE.NCHAN .AND. ICFROM(CP).LE.NCHAN
	ALLPART = ALLPART .or. 
     x      (LOCAL(CP).and.ICTO(CP).ne.ICFROM(CP))
!     x      (LOCAL(CP).and.REV(CP).and.ICTO(CP).ne.ICFROM(CP))
      DO 110 IN=1,2
        ICOM(CP,IN) = IC
        ICOR(CP,IN) = ICFROM(CP)
        IF(MASS(IN,ICOR(CP,IN)).LT.MASS(IN,ICOM(CP,IN))) GO TO 110
        ICOM(CP,IN) = ICFROM(CP)
        ICOR(CP,IN) = IC
C       NOW, ICOR(CP,IN) = CORE PARTITION & ICOM(CP,IN) = COMPOSITE.
  110   CONTINUE
C
C    DEAL WITH PARTICULAR KIND OF COUPLING
C    .....................................

      CALL INTER(ICTO(CP),ICFROM(CP),IC,KIND(CP),
     X           QQ(CP),IREM(CP),KPCORE(CP),QCM(CP),LAM(CP),
     X           BETAR(CP),BETAI(CP), KLT(CP),NIB(CP),HPOT,
     X        REV(CP),NLL(CP),LOCF(CP),COUPLE(CP),ICOR,ICOM,CP,IITER,
     X        XA(CP),XB(CP),XP(CP),XQ(CP),RSP,ffreal(CP),FILE(CP),
     X        MASS,HP,JEX,ENEX,QVAL,NEX,BAND,COPY,ITC,CPOT,
     X        PTYPE,NAME,FORMF,FORML,FORMC,BE,D0,BPROJ,POTCAP(1,1,CP),
     X        QNF,AFRAC,CCFRAC,FPT(1,1,CP),NKP(1,CP),GPT(1,1,CP),
     X        NFI,USED,NPRIOR,NPOST,OPN,M,HCM,EPS,MMXQRN,
     X        NF,NF0,NSP,NLCN,RIN,WID(CP),CENTR(CP),INFILE(CP),
     X        MATRIX,MEK,NIX,SPINTR, STREN,LAMBDA,HASO,CPSO(CP),
     X        LTRANS(1,CP),REW14,MTMIN,MLCALL(1,CP),MLCALL(2,CP),MCALLS,
     X        NBINS,PSIGN,EMID,DELTAE,NPWCOUP,PWCOUP,JPWCOUP)
C
      CP = CP + 1
      IF(KIND(CP-1).NE.7.OR.ABS(IREM(CP-1)).NE.2) GO TO 90
C          CONSTRUCT A NON-ORTHOGONALITY SUPPLEMENT FOR ABOVE COUPLING:
      ICTO(CP) = ICTO(CP-1)
          IF(.NOT.REV(CP-1)) ICTO(CP) = - ICTO(CP)
      ICFROM(CP) = ICFROM(CP-1)
      KIND(CP) = 8
      QQ(CP) = QQ(CP-1)
      IREM(CP) = 0
      KPCORE(CP) = 0
      BETAR(CP) = 0.0
      RMAX(CP) = RMAX(CP-1)
      JMAX(CP) = JMAX(CP-1)
      GO TO 100
C
  120 NCP= CP -1
      DO 121 JF=NF0+1,NF
  121 PTYPE(6,JF) = 0
      DO 125 KN=1,NSP
       L = QNF(9,KN)
       IF(L.EQ.0) GO TO 125
       DO 1245 J=1,2
       if(.not.cxwf) then
       DO 124 I=1,NLN
        R = (I-1)*RINTP
  124   FORML(I,KN,J) = FORML(I,KN,J) * R**L
	else
       DO 1241 I=1,NLN
        R = (I-1)*RINTP
 1241   FORMC(I,KN,J) = FORMC(I,KN,J) * R**L
	endif
 1245  CONTINUE
  125  CONTINUE
      IF(LREC66.GT.0) open(66,access='sequential',status='scratch',
     X                    form='unformatted')
      IF(LREC66.GT.0) REWIND 66

C
C    END OF COUPLING TYPES.    NOW FIND INCOMING BEAM
C    ................................................
 
!1240 FORMAT('0File ',I3,' needs NR =',I6)
         IEXCH = 0
       IF(COPY(2,PEL,EXL,2).EQ.PEL) THEN
         IEXCH = 1
         IF(JEX(1,PEL,EXL).GT. 0.25) IEXCH = -1
       ENDIF
!      if(final) WRITE(KO,1250) PEL,EXL,LAB,LIN,LEX
                WRITE(KO,1250) PEL,EXL,LAB,LIN,LEX
 1250 FORMAT(/' Incoming partition',I3,' in excitation state #',I2,
     X ',   Laboratory Energy given for partition',I3,' Nucleus',I2,
     X ' in Excitation pair',I2/)
      IF(IEXCH.NE.0) WRITE(KO,*) 'Scattering of identical particles with
     X SYMMETRISATION FACTOR = ',IEXCH
      DO 140 J=1,NJ+1
         JAP = JBORD(J) + JEX(1,PEL,EXL) + JEX(2,PEL,EXL)
  140    IF(JAP.GT.AINT(JAP+.1)+.1) JBORD(J) = JBORD(J) + 0.5
C
      NSA = 2*JEX(1,PEL,EXL) + 1.1
      NJA = 2*JEX(2,PEL,EXL) + 1.1
      XCOEF = 10.0/(NSA * NJA)
      DISCIT= NPOST .OR. NPRIOR .or. INITWF/=0
      DISC8 = DISCIT
	LAMAX = maxval(LAMBDA(1:MLOC))
	LAMAX = max(LAMAX,1)
C
C    INCOMING ENERGIES
C    -----------------
!      READ(KI,1255) (ELAB(I),NLAB(I),I=1,3),ELAB(4)
!      READ(KI,nml=incident)
!1255 FORMAT(3(F8.4,I8),F8.4)
      if(num_energies==0) then
      IF(abs(NLAB(1))+abs(NLAB(2))+abs(NLAB(3)).GT.-1)
     X WRITE(KO,1256) (ELAB(I),ELAB(I+1),NLAB(I),I=1,3)
 1256 FORMAT('0Lab. ENERGY ranges :',/,
     X       3(/'   from ',G14.6,' to',G14.6,' in',i9,' intervals'))
      else
       if(final) WRITE(KO,1257) 
     X        (energy_list(I),pel_list(I),I=1,num_energies)
 1257 FORMAT('0Lab. ENERGY ranges :',/(1x,10(f8.3,',',i1)))
      endif
        ELAB(5) = ELAB(4)
        NLAB(4) = 1

 	I = RMATR/HCM+1 - 10
	I = min(I,MINT-5)
	VRMX = 0d0
	VIMX = 0d0
	do JF=1,NF
        VRMX = max(VRMX,abs(real(FORMF(I,JF))))
        VIMX = max(VIMX,abs(aimag(FORMF(I,JF))))
        enddo
	if(final) write(ko,1280) (I-1)*HCM,VRMX,VIMX
1280	format(/' Largest real,imaginary parts of any form factor at R='
     X   ,f8.2,' are ',1p,2E10.2,' MeV')

        RETURN
	END SUBROUTINE READIN

	subroutine fkind(written,KO)
	character*40 flkind(299)
	integer writf(299),nwrit
	logical written(299)
	flkind(:) = ' '
	flkind( 1) = 'two-particle multipoles'
	flkind( 2) = 'KIND=9 nonlocal formfactor'
	flkind( 3) = 'local copy of User input'
	flkind( 4) = 'external form factors & potentials'
	flkind( 5) = 'standard input'
	flkind( 6) = 'standard output'
	flkind( 7) = 'elastic  S-matrix elements'
	flkind( 8) = 'single-particle & scattering wfs'
	flkind( 9) = 'complex transfer multipoles'
	flkind(10) = 'S-matrix elements'
	flkind(11) = 'real transfer multipole'
	flkind(12) = 'transfer kernels'
	flkind(13) = 'total cross sections/state'
	flkind(14) = 'interaction potentials'
	flkind(15) = 'local equivalent potentials'
	flkind(16) = 'tables of cross sections'
	flkind(17) = 'output scattering waves'
	flkind(18) = 'wfns of ''best'' iterate'
!	flkind(19) = 'temporary output amplitudes'
	flkind(20:33) = 'user files'
	flkind(34) = 'input optical potentials'
 	flkind(35) = 'Astrophysics S-factors  / Ecm'
	flkind(36) = 'output scattering AMPL amplitudes'
	flkind(37) = 'output scattering FAM amplitudes'
	flkind(38) = 'cross sections for each J/pi'
	flkind(39) = 'cross sections for each Ecm'
	flkind(40) = 'all cross sectns. for each Elab'
	flkind(41) = 'source terms at each iteration'
	flkind(42) = 'bin wavefunctions for each E'
	flkind(43) = 'bin phase shifts as k functions'
	flkind(44) = 'bin phase shifts as E functions'
	flkind(45) = 'scat phase shift as E functions'
	flkind(46) = 'bs wave functions & ANC ratios'
	flkind(47) = 'reduced matrix elements'
	flkind(48) = 'concurrency log file'
	flkind(49) = 'S-matrix elements concurrent->10'
	flkind(50) = 'standard output concurrent->6'
	flkind(51) = 'standard output for node N>0'
	flkind(52) = 'standard output for node N>0'
	flkind(53) = 'spec file for sturmx(x)'
	flkind(54) = 'mel file for sturmx(x)'
	flkind(55) = 'Single-particle wfns'
	flkind(56) = 'Fusion for each Jtotal'
	flkind(57) = 'CDCC amplitudes'
	flkind(58) = 'Single-particle wfns'
	flkind(59) = 'Single-particle vertex functions'
	flkind(60) = 'trace R-matrix calculations'
	flkind(61) = 'single-channel eigenenergies'
	flkind(63) = 'Asymptotic couplings'
	flkind(70) = 'Fusion for each J/pi'
	flkind(75) = 'S-factors for lab energies'
	flkind(88) = 'R-matrix eigenstate wfs'
	flkind(89) = 'All local couplings'
	flkind(90) = 'average polarisation potential'
	flkind(91) = 'polarisation potl for rereading'
	flkind(92) = 'partial fusion cross sections'
	flkind(99) = 'input J/L rescaling factors'
	flkind(201:210) = 'Separate cross sections'
	nwrit = 0
	do i=1,299
	 if(written(i)) then
		flkind(i) = trim(flkind(i))//'.'
	   nwrit = nwrit+1
	   writf(nwrit) = i
	 endif
	enddo
	write(KO,990) (writf(i),flkind(writf(i)),i=1,nwrit)
990	format(/'  The following files have been created:',
     X	/(1x,2(i3,':',a35)))
	return
	end

****WRITESPEC***********************************************************
	SUBROUTINE WRITESPEC(ipmax,HEADNG,RMASS,PEL,NLN,NAME,symm,
     X	  RINTP,NEX,COPY,MASS,JEX,BAND,ENEX,ntarg,nproj,levelp,levelt,
     X    ECM,MAL1)
	use parameters
	use io
	implicit none
	logical symm,CCMP,TEXT
	integer i,IC,IA,ipmax,nproj,ntarg,PEL,NLN,MAL1
      	INTEGER NEX(MXP+1),BAND(2,MXP,MXX),COPY(2,MXP,MXX,2)
       	integer  levelp(*),levelt(*)		
	real*8 RINTP,ECM
	REAL*8 MASS(4,MXP+1),RMASS(MXP),JEX(6,MXP,MXX),ENEX(2,MXP,MXX)
    	character*120 headng
      	CHARACTER*8 NAME(2,MXP+1)
	CCMP = abs(MELFIL)>=2
	TEXT = MELFIL<0
	IC = PEL
	written(54)=.true.
!				WRITE spec file for 1st partition ONLY
	write(54,'(A120)')   HEADNG
	write(54,*) 1./real(FMSCAL*RMASS(IC)),real(COULCN),real(ECM)
	write(54,102) 0.,RINTP,NLN,0.,RINTP,NLN,(NLN-1)*RINTP,CCMP,TEXT
102	format(2f8.4,i5,2f8.4,i5,f10.5,1x,2L2)
	 nproj = 0
	 ntarg = 0
	 do IA=1,NEX(IC)
	   if(COPY(1,IC,IA,1).eq.0) then
	    nproj = nproj + 1
	    levelp(nproj) = IA
 	   endif
	   if(COPY(2,IC,IA,1).eq.0) then
	    ntarg = ntarg + 1
	    levelt(ntarg) = IA
 	   endif
	 enddo
	write(6,*) nproj,' projectile states: ',(levelp(i),i=1,nproj)
	write(6,*) ntarg,' target     states: ',(levelt(i),i=1,ntarg)
	write(54,*) 2,0,4,2,ipmax,1,symm,.true.,1,max(nproj,ntarg),0,
     X			MAL1
	write(54,'(''  L Ip  J It JT'')')
* multipliers, jgt
	write(54,'(''  1  2  2  2  2  5'')')
* 4 coupling orders
	write(54,'(''  1  2  3  3  4  5'')')
* masses, charges, sizes
	write(54,*) MASS(1,IC),MASS(2,IC), MASS(3,IC),MASS(4,IC),0,0

        write(54,1045) NAME(1,IC),nproj,JEX(1,IC,levelp(1)),0d0,0d0,
     X	(ENEX(1,IC,levelp(i)),sign(1,BAND(1,IC,levelp(i))),i=1,nproj)
        write(54,1045) NAME(2,IC),ntarg,JEX(2,IC,levelt(1)),0d0,0d0,
     X	(ENEX(2,IC,levelt(i)),sign(1,BAND(2,IC,levelt(i))),i=1,ntarg)
1045      format(a8,i3,f4.1,2f10.4,10(f8.4,i4))
	write(54,1046) .false.,0.
1046 	  format(3(L2,f4.1))
	call flush(54)
	RETURN
	END SUBROUTINE WRITESPEC
****WRITEMEL************************************************************
	SUBROUTINE WRITEMEL(PEL,jtotal,psign,parity,nch,part,MINTL,
     X	  ipmax,ETOTAL,LVAL,JVAL,JPROJ,JTARG,CI,NLN,FORMF,STREN,
     X    LAMBDA,NF,MR,EXCIT,symm,nproj,ntarg,levelp,levelt,MASS,INITL,
     X    CLIST,NFLIST,NCLIST,ET,ryev,CF)
	use parameters
	use io
      	implicit none
	integer PEL,parity,nch,i,MINTL,levelp(*),levelt(*),ntarg,nproj,
     X   part(MAXCH,3),LVAL(nch),INITL(MINTL),IF,ia,JIN,EXCIT(MAXCH,3),
     X	 ipmax,NLN,NF,LAMBDA(NF),IC,IMAX,C,C2,MR,NCLIST(MAXCH,MAXCH),
     X   NFLIST(MAXCH,MAXCH,MCLIST),NC
	character psign(3)
	real*8 jtotal,ETOTAL,JVAL(NCH),JPROJ(NCH),JTARG(NCH),ET(NCH)
	REAL*8 MASS(4,MXP+1),cfp(ipmax),STREN(NF),T,
     x		ryev,CF(MAXCH,MAXCH,ipmax)
  	integer,allocatable:: ichl(:)
	complex*16 CI,FORMF(MAXM,MLOC),CLIST(MAXCH,MAXCH,MCLIST),S,
     X     ww(NLN)
	logical symm,CCMP,TEXT
C
	written(53) = .true.
	IC = PEL
	CCMP = abs(MELFIL)>=2
	TEXT = MELFIL<0
c---------------------------------- write files 
   	write(54,1361) jtotal,psign(parity+2)
1361   	format(' Coupled Channels Set ',f5.1,1x,a1)
	 imax = 0
	 do i=1,NCH
	  if(PART(i,1)==IC) imax=i
	 enddo

	write(54,'(4i5,f10.5,g12.4)') 
     X 		nint(2*jtotal),parity,imax,MINTL,0d0,ETOTAL

       allocate (ichl(imax))
	write(54,*) (1,i=1,imax)
*					In order: projectile,target
	 do i=1,imax
	  ichl(i) = 0
	  do ia=1,nproj
	   if(EXCIT(i,2).eq.levelp(ia)) ichl(i) = ia
	  enddo
	 enddo
	write(54,*) (ichl(i),i=1,imax)
!	write(6,*) 'proj:',(ichl(i),i=1,imax)
	 do i=1,imax
	  ichl(i) = 0
	  do ia=1,ntarg
	   if(EXCIT(i,3).eq.levelt(ia)) ichl(i) = ia
	  enddo
	 enddo
	write(54,*) (ichl(i),i=1,imax)
!	write(6,*) 'targ:',(ichl(i),i=1,imax)
        deallocate (ichl)

c quantum numbers written in order  L,Ip,J,It,JTOT
 	write(54,*) (LVAL(i),i=1,imax)
 	write(54,*) (nint(2*JPROJ(i)),i=1,imax)
 	write(54,*) (nint(2*JVAL(i)),i=1,imax)
 	write(54,*) (nint(2*JTARG(i)),i=1,imax)
	write(54,*) (nint(2*JTOTAL),i=1,imax)
	if(MINTL>0) write(54,*) (INITL(JIN),JIN=1,MINTL)
	write(54,*) 0,0
	call flush(54)
!!!!			Now some extra output for  Gailitis tests
!@	write(63,*) imax,ipmax,ETOTAL,ryev,jtotal,parity,symm

	if(.not.symm) then
        do 342 c=1,imax
!@		write(63,*) LVAL(c),ET(c)
        do 342 c2=1,imax
          S = CI**(-LVAL(C2)+LVAL(C))
	  ww(1:NLN) = 0d0
	  cfp(1:ipmax) = CF(c,c2,1:ipmax)*ryev
	  do NC=1,NCLIST(c,c2)
	    IF = NFLIST(c,c2,NC)
	    t = S*CLIST(c,c2,NC)
	    if(abs(t)>1d-20) then
	    do i=1,NLN
	    ww(i) = ww(i) + t*FORMF((i-1)*MR+1,IF)
 	    enddo
	    endif
	  enddo
	if(.not.TEXT) then
        if(.not.CCMP) write(53) (real(ww(i)),i=1,NLN),cfp(1:ipmax)
        if(     CCMP) write(53) (ww(i),i=1,NLN),cfp(1:ipmax)
	else
        if(.not.CCMP) write(53,34) (real(ww(i)),i=1,NLN),cfp(1:ipmax)
        if(     CCMP) write(53,34) (ww(i),i=1,NLN),cfp(1:ipmax)
34	format(1p,6e12.4)
	endif
!@	write(63,*) cfp(1:ipmax)
342	continue
        do 343 i=1,imax
	if(.not.TEXT) write(53) (0d0,c=1,i-1),1d0,(0d0,c=i+1,imax)
	if(     TEXT) write(53,34) (0d0,c=1,i-1),1d0,(0d0,c=i+1,imax)
343     continue
        if(.not.TEXT) write(53) (1,i=1,imax)
        if(     TEXT) write(53,'(40i2)') (1,i=1,imax)
 	else
        do 344 c=1,imax
!@		write(63,*) LVAL(c),ET(c)
        do 344 c2=c,imax
          S = CI**(-LVAL(C2)+LVAL(C))
	  ww(1:NLN) = 0d0
	  cfp(1:ipmax) = CF(c,c2,1:ipmax)*ryev
	  do NC=1,NCLIST(c,c2)
	    IF = NFLIST(c,c2,NC)
	    t = S*CLIST(c,c2,NC)
	    if(abs(t)>1d-20) then
	    do i=1,NLN
	    ww(i) = ww(i) + t*FORMF((i-1)*MR+1,IF)
 	    enddo
	    endif
! 	if(i>1)     write(105,3431) c,c2,IF,i,cfp(i),t,STREN(IF)
!3431	  format('CFF:',4i4,6g12.3)
	  enddo
	if(.not.TEXT) then
           if(.not.CCMP) write(53) (real(ww(i)),i=1,NLN),cfp(1:ipmax)
           if(     CCMP) write(53) (ww(i),i=1,NLN),cfp(1:ipmax)
	   else
           if(.not.CCMP) write(53,34) (real(ww(i)),i=1,NLN),cfp(1:ipmax)
           if(     CCMP) write(53,34) (ww(i),i=1,NLN),cfp(1:ipmax)
	   endif
!@	write(63,*) cfp(1:ipmax)
344	continue
	endif
	call flush(53)
!@	call flush(63)
        RETURN
	END SUBROUTINE WRITEMEL
