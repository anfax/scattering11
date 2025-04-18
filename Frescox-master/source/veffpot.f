!!!*************************************************************
! 文件/File: veffpot.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:17:02
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
****VEFFPOT**************************************************************
!	SUBROUTINE VEFFPOT
C
C    CALCULATE COUPLED-CHANNELS EFFECTIVE POTENTIAL FOR ELASTIC CHANNELS
C    -------------------------------------------------------------------
C      (AVERAGED OVER ALL INCOMING PARTIAL WAVES)
C
      REWIND 15
         LAP = 0.6 * JTMAX + 0.4 * ABS(JTMIN)
         JN =JEX(1,PEL,EXL)
      DO 790 IMA=1,NSA
      READ(15) (FNC(I,IMA),I=1,NLN)
         DO 790 I=1,NLN
         VOPT(I) = 0.0
  790    WNM(I,IMA,1) = 0.0
      DO 800 IMA=1,NSA
      READ(15) (WNM(I,IMA,2),I=1,NLN)
      IF(VEFF.GT.0) GO TO 800
        JAP = LAP + IMA - JN - 1.
        IF(JAP.LT.0.0) GO TO 800
         DO 792 I=1,NLN
  792    WNM(I,IMA,1) = 0.0
         DO 795 IF=1,NF
          IF(EL.LE.IEX .AND.PTYPE(3,IF).GT.0) GO TO 795
         DO 793 I=1,NLN
  793   WNM(I,IMA,1) = WNM(I,IMA,1)+ELCOEF(IF)*FORMF((I-1)*MR+1,IF)
  795   CONTINUE
  800 CONTINUE
      DONE = 0
      REPEAT = .FALSE.
      DO 820 I=1,NLN
      DO 810 IMA=1,NSA
          IF(ABS(WNM(I,IMA,2)).LT.1E-35) GO TO 810
            REPEAT = .TRUE.
      WNM(I,IMA,2) = WNM(I,IMA,1) + FNC(I,IMA) / WNM(I,IMA,2)
  810 CONTINUE
      T = (I-1)*HP(PEL)*MR
      IF(REPEAT) DONE = DONE + 1
      IF(DONE.EQ.1) WRITE(KO,1480) VEFF
 1480 FORMAT(///' Mean components of Elastic Optical Potential generated
     X by channel couplings of type',I3/)
      IF(REPEAT) WRITE(KO,1358) T,(WNM(I,IMA,2),IMA=1,NSA)
*** 1358 FORMAT(' ',F9.3,11F11.4/(11X,11F11.4) )
  820 CONTINUE
      REWIND 15
         DO 840 IMA=1,NSA
         IF(VEFF.GT.0) GO TO 835
            DO 830 I=1,NLN
  830       WNM(I,IMA,2) = WNM(I,IMA,2) - WNM(I,IMA,1)
  835     WRITE(15) (WNM(I,IMA,2),I=1,NLN)
         DO 838 I=1,NLN
C                                VOPT is spin-averaged potential.
  838    VOPT(I) = VOPT(I) + WNM(I,IMA,2)/NSA
  840    CONTINUE
      REWIND 15
	   call rewop(90)
	   call rewop(91)
         WRITE(91,*) 'Mean Elastic Optical Potential generated ',VEFF
         WRITE(90,845) NLN,RINTP,ENLAB
845	format('#',i4,f8.4,f10.4,' = NLN,RINTP,ELAB')
         WRITE(91,*) NLN,RINTP,0.
	 DO I=1,NLN
         WRITE(90,'(f8.3,2f12.6)') (I-1)*RINTP,VOPT(I)
	 enddo
         WRITE(90,*) 'END'
         WRITE(91,'(1p,6e12.4)') VOPT(1:NLN)
         call flush(90)
         call flush(91)
         written(90) = .true.
         written(91) = .true.


!	END SUBROUTINE VEFFPOT


