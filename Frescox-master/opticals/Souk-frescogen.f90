!!!*************************************************************
! 文件/File: Souk-frescogen.f90
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
!######### OPTICAL MODEL PARAMETERS for FLAP2.2 ##########

!          neutron  on 238U 

! Energy    V     rv    av      W     rw    aw      Vd   rvd   avd      Wd   rwd   awd     Vso   rvso  avso    Wso   rwso  awso  rc

!f7.3,f7.3,f5.3,f5.3,f7.3,f5.3,f5.3,f7.3,f5.3,f5.3,f7.3,f5.3,f5.3,f7.3,f5.3,f5.3,f7.3,f5.3,f5.3,f5.3f
	character*4 NAME
	character*8 POTL
	character*100 fname
	real elevels(20), jlevels(20)
	
	Z = 92
	A = 238.052163
	NAME = '238U'
	POTL = 'MASLOV03'
	POTL = 'KoningDR'
	POTL = 'Soukhovi'
	kp = 1
	nex = 5
	elevels(1:10) =(/ 0., .044916, .14838, 0.30718, 0.5181, 0.7759, 1.0767, 1.4155, 1.7884, 2.1911/)
	jlevels(1:10) =(/ 0., 2., 4., 6., 8., 10., 12., 14., 16., 18. /)
	
	EMIN = 2.5
	EMAX = 2.5
	NE = 1
	
	write(6,1) POTL,NAME
1	format('####### OPTICAL PARAMETERS for ',A8,' ########'//'      neutron on ',a4,// &
     & ' Energy    V     rv    av      W     rw    aw     ', & 
     &           ' Vd   rvd   avd      Wd   rwd   awd     ', &
     &           'Vso   rvso  avso    Wso   rwso  awso  rc')
	
	DE = 1.
	IF(NE>1) DE = (EMAX-EMIN)/(NE-1)
      NA = nint(A); NZ = nint(Z)
	AC = real(NA)**(1./3.)

	DO IE=1,NE
	E = EMIN + (IE-1)*DE
	RC = 1.2 ! default
	
	if(POTL=='MASLOV03') then
	 call MASLOV03(E,VR,RR,AR, W,RW,AW, VD,RVD,AVD, WD,RD,AD,VSO,RSO,ASO, WSO,WRSO,WASO, RC)
	else if(POTL=='KoningDR') then
           NTYPE = 1 ! neutron!!
	 call KDParam(NA,NZ,NTYPE,E,VR,RR,AR,W,RW,AW,WD,RD,AD,VSO,RSO,ASO,WSO,WRSO,WASO,RC,D3)
		VD = 0.0; RVD=0.; AVD = 0.0
	else if(POTL=='Soukhovi') then
           NTYPE = 1 ! neutron!!
	 call soukhovitskii(A,Z,NTYPE,E,VR,RR,AR,VD,RVD,AVD,W,RW,AW,WD,RD,AD,VSO,RSO,ASO,WSO,WRSO,WASO,RC)
		VD = 0.0; RVD=0.; AVD = 0.0
	else 
	  write(0,*) ' Potential <'//trim(POTL)//'> not recognized'
	  stop 
	endif
	
	
	
	RVOLV = AC * RR
	RVOLW = AC * RW
	RSURF = AC * RD
	
	BETA2 = 0.195
	BETA4 = 0.078	
	
	write(6,10) E,VR,RR,AR, W,RW,AW, VD,RVD,AVD, WD,RD,AD,VSO,RSO,ASO, WSO,WRSO,WASO, RC
10	format(f7.3, 6(f8.3,2f6.3),f6.3)
      fname='fresco-00-'//POTL//'-s'//CHAR(ICHAR('0')+nex)//'-E0000000.in'
        write(fname(8:8),'(i1)') mod(nint(z),10)
        write(fname(9:9),'(i1)') mod(nint(a),10)
        write(fname(24:30),'(f7.3)') e
        write(fname(24:26),'(i3.3)') int(e)
        write(0,*) 'Create file <'//trim(fname)//'>'
        open(1,form='formatted',file=trim(fname))
	write(1,'(a)') 'n+'//trim(NAME)//' with '//trim(POTL)//', s='//CHAR(ICHAR('0')+nex)//' at E ='//fname(24:30)
	write(1,'(a)') 'NAMELIST'
	write(1,'(a)') ' &Fresco  hcm= 0.1 rmatch=  20.000'
	write(1,'(a)') '    jtmin=   0.0 jtmax= 20 absend= 0.000001 '
	write(1,14) nex
14	format('    thmin=0.0 thinc=2 thmax=180. iblock=',i3)
	write(1,'(a)') '    chans= 1 smats= 2 xstabl= 1'
	write(1,15) E
15	format('    elab=',f10.5,'  pel=1 exl=1 lab=1 lin=1 lex=1 /')
	write(1,*) 
 	write(1,16) nex
16	format('&Partition namep=''n       '' massp=  1.008665 zp=  0 nex=',i3)
	write(1,17) NAME,A,Z
17	format('            namet=''',a8,''' masst=',f10.6,' zt=',f5.1,' qval=  0.000/')
 	write(1,'(a)') '&States jp= 0.5 ptyp= 1 ep=  0.000000  cpot=  1 jt= 0.0 ptyt= 1 et= 0.000000/'
 	do 20 il = 2,nex
20 	write(1,21) kp,jlevels(il),elevels(il)
21	format('&States copyp= 1                       cpot=',i3,' jt=',f4.1,' ptyt= 1 et=',f8.4,'/')
 	write(1,'(a)') '&Partition /'
	write(1,*) 

	write(1,30) kp,0,0,real(NA),0.0,RC
30	format('&POT kp=',i3,' type =',i2,' shape=',i2,' p(1:3)=',3f9.4,'/')
31	format('&POT kp=',i3,' type =',i2,' shape=',i2,' p(1:6)=',6f9.4,'/')
32	format('&POT /'/)
	write(1,30) kp,1,0,VR,RR,AR
	DELTA2 = BETA2 * RVOLV
	DELTA4 = BETA4 * RVOLV
	write(1,31) kp,11,13, 0., DELTA2, 0., DELTA4, 0., 0.
	if(abs(W)>1e-10) then
	write(1,31) kp,1,0,0.,0.,0., W,RW,AW
	DELTA2 = BETA2 * RVOLW
	DELTA4 = BETA4 * RVOLW
	write(1,31) kp,11,13, 0., DELTA2, 0., DELTA4, 0., 0.
	endif
	write(1,31) kp,2,0,VD,RVD,AVD, WD,RD,AD
	DELTA2 = BETA2 * RSURF  ! assume VD=0
	DELTA4 = BETA4 * RSURF
	write(1,31) kp,11,13, 0., DELTA2, 0., DELTA4, 0., 0.
	write(1,31) kp,3,0,VSO,RSO,ASO, WSO,WRSO,WASO
	write(1,32)
	
	write(1,*) '&Overlap /'
	write(1,*) '&Coupling /'

	close(1)
	enddo
	end
	subroutine maslov03(E,VR,RR,AR, W,RW,AW, VD,RVD,AVD, WD,RD,AD,VSO,RSO,ASO, WSO,WRSO,WASO, RC)
	VR = 45.93 - 0.28*E + 0.000573*E*E
	RR = 1.26; 	AR = 0.63
	if(E<8.) then
		WD = 3.14 + 0.436*E
	   else
	    	WD = 6.628
	   endif
	RD = 1.26; 	AD = 0.52
	W = 0.0; RW = RR; AW = AR	
	RC = 0.0
	VSO = 6.2; 	RSO =1.120; ASO = 0.47; WSO = 0.0; 	WRSO=RSO; WASO=ASO
	return
	end