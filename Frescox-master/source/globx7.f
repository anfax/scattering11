!!!*************************************************************
! 文件/File: globx7.f
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
	module parameters
	integer mxx,mxp,mxpex,maxit,maxcpl,maxnnu,
     X		maxn,maxnln,maxnlo,msp,lmax1,mpair,nkmax,maxm,mpwcoup
        integer   mloc,maxch,maxqrn,maxqrn2,mfnl,mclist,maxpw,maxb
	integer inh,maxich,maxcch,nfdec
!					derivative parameters:
	integer maxmul
	real*8 unitmass,finec,fmscal,coulcn,amu,hbc,etacns
	real*8 pmass,mun
	integer melfil,nfus,nfus1,kpmax,buttle
	integer maxnlr,maxnlc,nchbinmax,dgam
	end module parameters

	module io
      	integer KI,KO,KOI
	logical written(399),grace
        character*14 pre_is,post_is
	end module io

	module factorials
	integer mfact
	real*8 pi,r4pi
	real*8,save,allocatable:: fact(:) 
!  	real*8 fact(1:2000) 
	end module factorials

	module drier
	integer mach,eigens
	logical dry,psiren
	real*8 fpmax,acc8
	end module drier
	
	module kcom
    	real*8 rintp,epc
    	integer nl,maxl,mlt,nln,nlo,nlc,locpw,nlocpw,
     X         n,mr,nbb(3),minl,jtest,lrec66,nlt,nlm,nnu,mtau(0:3)
	data mtau /  0,0, 3, 3 /
	end module kcom
	module kcom2
      	real*8,save,allocatable::     FNL(:,:),QRLN(:,:,:)
      	complex*16,save,allocatable:: FNC(:,:),QERN(:,:,:)
      	integer,save,allocatable::    WHOL(:),WHOI(:),WHERE(:,:)
      	real*8 z,eps,hfht,rinto
      	integer maxl1,minl1,minlrq,maxlrq,ib
	end module kcom2
	
	module trace
	integer chans,listcc,treneg,cdetr,smats,smatl,xstabl,nlpl,waves
     X         ,lampl,veff,kfus,wdisk,bpm,cdcc,tcfile
	logical say
	end module trace

	module fresco1
        character*120 headng,tcfilename
      	character*1 iso,btype
      	character*2 rela
	real*8 hcm,rmatch,hnl,rnl,centre,hnn,rnn,rmin,rasym,accrcy,rsp,
     X   switch,ajswtch,sinjmax,jtmin,jtmax,absend,erange,dk,hktarg,
     X   thmin,thmax,thinc,cutl,cutr,cutc,ips,jbord(8),elab(5),phase,
     x   jleast,time0,sock(5,2),dist0,elpmax
 	integer nearfa,jump(7,3),koords,kqmax,pp,isocen,nnn,ngail,
     X		it0,iter,iblock,pade,mtmin,nlab(4),m,md,lmax,meigs,
     X		pel,exl,lin,lab,lex,nrbases,nrbmin,pcon,mrm,plane,
     x          vsearch,echan,enodes,nlagcc,KRM,relref
	logical fatal,nosol,fcwfn,pralpha,symm,locfil,mcalls,rterms,
     x		allpart,cxwf,ccbins,nn2wf,sumccbins,sumkql,bloch,
     x          ccreal,gscouplonly,lanecoup,eobs
	real*8 rin,bndx(2),rmatr,smallchan,smallcoup,gailacc,weak
	integer mlm,nlcn,mint,mintm2,nj,jset,pset,icutc,itmin,nsol,
     x          nforml1,nforml2,sumform,ompform,initwf,npluto,pluto(10)
	character*100 ppots(9,20)
	integer nextpot(20)
        character*10 ppkind(0:4)
        data ppkind / 'projectile','target','ejectile','residual',
     x		      'proj+eject' /
	character*4 coords(4)
      	data coords / 'Mads','Tran','Recl','H-J ' /
 	end module fresco1


	module gails
	implicit real*8(a-h,o-z)
	real*8 ryev
	logical cfalloc
     	real*8, save,allocatable:: cf(:,:,:)
	integer numax
!	real*8, save,allocatable:: gaa(:,:,:),gba(:,:,:)
!	common /machs/ scalx,rscalx
!       real*8 scalx,rscalx
	end module gails

	module parallel
#ifdef MPI
        use mpi
#endif /* MPI */
!      	integer numt,msgtype,info,mytid
!       integer msgtag,atag,alen,atid
        integer numthread,iame,iams,mpinodes,mpisets,
     x          mpihelp,ihelp,master,commgroup
	integer,allocatable:: ICNTXT(:) 
#ifdef MPI
	integer ierr,status(MPI_STATUS_SIZE)
#endif /* MPI */
        logical isparallel
        logical MPIC,STDINALL
	end module parallel

	module searchpar
 	integer nvars
	parameter (mvars=4000)
	integer srch_kind(mvars)
!  Kind of search parameter: 0=ignore, 1=potential, 2=afrac, 3=R-mat energy,
!                            4=R-mat width, 5=dataset normalisation
!
	character*15 srch_name(0:mvars)
	real*8 srch_value(mvars),srch_step(mvars),srch_damp(mvars),
     x      srch_minvalue(mvars),srch_maxvalue(mvars),nul,srch_B(mvars),
     x      srch_leff(mvars)
	integer srch_kp(mvars),srch_pline(mvars),srch_col(mvars)
	integer srch_pline2(mvars),srch_col2(mvars)
	integer srch_nafrac(mvars),numafrac
	character*50 srch_afrac_overlap(mvars)
	real*8 srch_jtot(mvars),penalty,fine,srch_ratio2(mvars),
     x         srch_r_jch(mvars),srch_r_sch(mvars),srch_r_shift(mvars)
	integer srch_par(mvars),srch_r_ch(mvars),srch_r_ic(mvars),
     x          srch_r_ia(mvars),srch_r_lch(mvars)
	integer srch_datanorm(2,mvars),srch_rterm(mvars)
        character*80 srch_reffile(mvars)
	logical srch_nopot(mvars),srch_rwa(mvars),srch_Brune(mvars)
	logical interactive,final,boldplot,hasEshift,firstE
        logical,save,allocatable::  rm_Brune(:)
        real*8,save,allocatable::  E_Brune(:),W_Brune(:)
	integer number_calls,maxleg,mterms
	end module searchpar
	
	module searchdata
	!parameter(mds=160, mdl=1000, maxen=mdl*2)
	parameter(mds=260, mdl=1000, maxen=200000)
	integer datasets,datalen(mds),data_ien(mdl,mds)
	real*8 datangles(mdl,mds),datavals(mdl,mds),dataerr(mdl,mds)
	real*8 theoryvals(mdl,mds),data_energies(mdl,mds),data_jtot(mds)
        real*8,save,allocatable::  theoryplot(:,:,:)
	character*11 data_info(mdl,mds),dat_info(mds)
	logical data_lab(mds),data_matched(mds),gettheoryplot
	logical dat_eplots(mds),dat_logy(mds),data_Mflip(mds)
	integer data_type(mds),bs_no(mdl,mds),data_kind(mds),ntheory_pts
	integer data_pel(mds),data_exl(mds),data_labe(mds),
     x 	        data_lin(mds),data_lex(mds),data_ib(mds)
	integer data_ic(mds),data_ia(mds),data_thfile(mds),ndof
	integer data_idir(mds),data_idir1(mds),data_rank_k(mds)
	integer data_ch(mds),data_par(mds),data_rank_q(mds),
     x          data_normvar(mds),data_shiftvar(mds),data_term(mds)
        character*80 data_reffile(mds),dat_file(mds)
        real*8 data_chisq(mds),ecmrat,etarat
	real*8 energy_list(maxen),esmin,esmax
	integer num_energies,energy_count(maxen),pel_list(maxen)
	end module searchdata


