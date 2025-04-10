!!!*************************************************************
! 文件/File: wfpr.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:17:01
!*************************************************************

       implicit real*8(a-h,o-z)
       parameter (MAXN=100,MAXCH=10)
       real*8 JTOTAL,HCM,ENLAB,JVAL(MAXCH)
       integer PARITY,LVAL(MAXCH),ITC(1,1),PART(MAXCH),EXCIT(MAXCH),
     x		WDISK
       complex*16 psi(MAXN,MAXCH),smat(MAXCH)


       WRITE(17,657) N,HCM,ENLAB,JTOTAL,PARITY
657   FORMAT(I4,2F8.4,F8.1,I3)
      DO 690 IT=1,ITCM
          IF(ABS(MOD(WDISK,2)).EQ.1.AND.IT.NE.ITC(PEL,EXL)) GOTO 690
          DO 670 C=1,NCH
            IF(ITC(PART(C),EXCIT(C)).NE.IT) GOTO 670
            IF(ABS(MOD(WDISK,2)).EQ.1.AND.C.ne.JIN) GOTO 670
          WRITE(17,660) IT,LVAL(C),JVAL(C),JTOTAL,LVAL(JIN),JVAL(JIN),
     X                   SMAT(C)
  660         FORMAT(2I4,2F6.1,I4,F6.1,2F15.10)
              WRITE(17,680) (PSI(I,C),I=1,M)
  670      CONTINUE
  680      FORMAT(1P,6E12.4)
  690  CONTINUE
         WRITE(17,660) -1
       end
