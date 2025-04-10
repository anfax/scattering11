!!!*************************************************************
! 文件/File: ion_atom.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:17:02
!*************************************************************

!!      ion_atom.f90
!!-------------------------------------------------------------------------!!
!!                                                                         !!
!!                                                                         !!
!!                    I O N   A T O M   W A V E                            !!
!!                                                                         !!
!!                                                                         !!
!!                                                                         !!
!!           A Fortran 90 PROGRAM TO CALCULATE CROSS SECTIONS              !! 
              
!!           FOR ION - He COLLISIONS AT INTERMEDIATE ENERGIES              !!
!!                     USING THE WAVE TREATMENT.                           !!
!!                                                                         !!
!!                                BY                                       !!
!!                                                                         !!
!!                          Brian Nesbitt *,                               !!
!!                         Francesca O'Rourke                              !!
!!                                 &                                       !!
!!                        Derrick SF Crothers                              !!
!!                                                                         !!
!!             Theoretical & Computational Research Division,              !!
!!        Department of Applied Mathematics & Theoretical Physics,         !!
!!                                                                         !!
!!                   The Queen's University of Belfast,                    !!
!!                           Belfast Bt7 1NN,                              !!
!!                          Northern  Ireland.                             !!
!!                                                                         !!
!!                                                                         !!
!!    *  e-mail:  b.nesbitt@am.qub.ac.uk                                   !!
!!                                                                         !!
!!                                                                         !!
!!    Version 1.0:  1st April 1998                                         !!
!!                                                                         !!
!!                                                                         !!
!!      NOTE                                                               !!
!!     ------                                                              !!
!!                                                                         !!
!!  The quadrature in this program is performed by the QUADPACK [1]        !!
!!  routines QAGE using  the 15-21 and 61 point Gauss-Kronrod rules which  !!
!!  have been updated to Fortran 90 by the author.                         !!
!!                                                                         !!
!!                                                                         !!
!!      DISCLAIMER                                                         !!
!!     ------------                                                        !!
!!                                                                         !!
!!  Despite the fact that the modified Fortran 90 QUADPACK subroutines     !!
!!  have been tested and appear to be bug free, the authors would like to  !!
!!  point out that if any errors are found then these failings should not  !!
!!  in any way be reflected on the original code by R. Piessens et al..    !!
!!  The original code [1] has been subjected to  numerous tests and has    !!
!!  been found to be extremely reliable and accurate.                      !! 
!!                                                                         !!
!!      References                                                         !!
!!     ------------                                                        !!
!!                                                                         !!
!!  [1] R. Piessens, E. de Doncker-Kapenga, C.W. Uberhuber & D.K. Kahane,  !!
!!      QUADPACK, A Subroutine Package for Automatic Integration, Springer !!
!!      Series in Comp. Math. 1., Springer-Verlag, New York, 1983.         !!
!!                                                                         !!
!!      QUADPACK is freely available at:                                   !!
!!             www.netlib.org                                              !!
!!      or at the UK mirror:                                               !!
!!          www.hensa.ac.uk/cgi-bin/cftp/unix.hensa.ac.uk/mirrors/netlib   !!
!!                                                                         !!
!!-------------------------------------------------------------------------!!




Module Precision

       !!---------------------------------------------------------------!!
       !!                                                               !!
       !!  Module Precision defines the complie-time constants  mp and  !!
       !!  c_mp. The user need only redefine these constants to change  !!
       !!  from one precision to another. This feature enhances the     !!
       !!  portability of the program.                                  !!
       !!                                                               !!
       !!---------------------------------------------------------------!!     

       implicit none

       !!  Define precision for real values.  !!

       integer, parameter :: mp = selected_real_kind(10,30) 

       !!  Define precision for complex values.  !!

       integer, parameter :: c_mp = mp 

       real(kind=mp) :: epmach, uflow, tolerance  ! Machine dependent constants.

Contains

      Subroutine Machine_Accuracy

             !!----------------------------------------------------------------!!
             !!                                                                !!
             !!      Description                                               !!
             !!     -------------                                              !!
             !!                                                                !!
             !!  This subroutine defines the machine dependent constants:      !!
             !!                                                                !!
             !!      epmach (real)     -     Largest relative spacing.         !!
             !!                                                                !!
             !!      uflow (real)      -     Smallest postive spacing.         !!
             !!                                                                !!
             !!      tolerance (real)  -     Machine dependent tolerance.      !!
             !!                                                                !!
             !!----------------------------------------------------------------!!

             implicit none

             real(kind=mp) :: dummy_variable

             epmach = epsilon(dummy_variable)
             uflow =  tiny(dummy_variable)
             tolerance  = 1.e+6_mp*epmach

      end Subroutine Machine_Accuracy

end Module Precision



Module Global_Data
        
       !!-------------------------------------------------------------------!!
       !!                                                                   !!
       !!              Description                                          !!
       !!             -------------                                         !!
       !!                                                                   !!
       !!   This module reads in the input data from files when the         !! 
       !!   subroutine 'Get_Data' is called by the main program. This data  !!
       !!   may then be accessed by other program units using the 'use      !!
       !!   Global_Data' statement.                                         !!
       !!                                                                   !!
       !!                                                                   !! 
       !!              Details                                              !1
       !!             ---------                                             !!
       !!                                                                   !!
       !!   This module contains one subroutine which is called to read     !!
       !!   input data from units 'get1' and 'get2'. The file associated    !!
       !!   with unit 'get1' contains the following information:            !!
       !!                                                                   !!
       !!       target (integer)  -  0 for hydrogen-like atom.              !!
       !!                            1 for helium-like atom.                !!
       !!                            2 for Molecular hydrogen.              !!
       !!                                                                   !!
       !!       no (integer)  -   1 or 2 depending on the chosen theory.    !!
       !!                                                                   !!
       !!       my_section (integer)  -  0 for total cross section.         !!
       !!                                1 for single differential          !!
       !!                                2 for double differential          !!
       !!                                                                   !! 
       !!       ZP (real)     -   The projectile charge.                    !!
       !!                                                                   !!
       !!                                                                   !!
       !!       Projectile_Energy_Start  (real)  -  Starting point of       !!
       !!                                         projectile energy in keV  !!
       !!                                                                   !!
       !!       Projectile_Energy_Finish (real)  -  End point.              !!
       !!                                                                   !!
       !!       Projectile_Energy_Increment (real) - Determines number of   !!
       !!                                     points between start and      !!
       !!                                     finish.                       !!
       !!                                                                   !!  
       !!                                                                   !!
       !!       Electron_Energy_Start (real)  -  Energy of ejected electron !!
       !!                                        measured in  eV.           !!
       !!                                        Required for DDCS and      !!
       !!                                        SDCS.                      !!
       !!                                                                   !!
       !!       Electron_Energy_Finish (real) -  End point.                 !!
       !!                                                                   !!
       !!       Electron_Energy_Step (real) -  Determines number of points  !!
       !!                                      between start and finish     !!
       !!                                                                   !!
       !!       acc (real)    -   Tolerance used for the outer integral.    !!
       !!                                                                   !!
       !!       tol (real)    -   Global tolerance used throughout code.    !!
       !!                                                                   !!
       !!       angle (real)  -   The polar angle (wrt the velocity vector  !!
       !!                         v of the incident ion) of emission for    !!
       !!                         the ejected electron, measured in radians.!!
       !!                         This value  is only  required for  the    !!
       !!                         double differential cross section since   !!
       !!                         the single differential cross section is  !!
       !!                         calculated with respect to ejected        !!
       !!                         electron energy.                          !!
       !!                                                                   !!
       !!                                                                   !!
       !!   The file associated with unit 'get2' contains the following     !!
       !!   information :                                                   !!
       !!                                                                   !!
       !!       Z_orbital (real)                                            !!
       !!                            -   He orbitals.                       !!
       !!       N_orbital (real)                                            !!
       !!                                                                   !!
       !!       num_orb (real)       -   Size of vectors Z_orbital and      !!
       !!                                N_orbital.                         !!
       !!                                                                   !!
       !!                                                                   !!
       !!-------------------------------------------------------------------!! 

       use Precision, only : mp
       implicit none

       real(kind=mp) :: ZP, tol, acc, pi, angle, Projectile_Energy_Start,           &
                        Projectile_Energy_Finish, Projectile_Energy_Increment,      &
                        Electron_Energy_Start, Electron_Energy_Finish,              &
                        Electron_Energy_Step

       real(kind=mp), parameter :: zero = 0.0_mp
       real(kind=mp), parameter :: Rydberg = 27.2116_mp   !  The Rydberg constant.  

       integer, parameter :: global_limit = 500    !  Upper bound on the number of subintervals
                                                   !  in the partion (a,b) of the integral to be
                                                   !  evaluated. 

       integer, parameter :: iout = 6   ! Output channel.
       integer, parameter :: get1 = 8   ! Input channel for general data.
       integer, parameter :: get2 = 9   ! Input channel for atomic orbitals.
        
       integer :: no, num_orb, my_section, target

       real(kind=mp), allocatable, dimension(:) :: Z_orbital, N_orbital
     
       character(len=8) :: Theory(2)  ! CDW or CDW-EIS.

Contains

       Subroutine Get_Data

            implicit none

            integer :: m, loop
            integer :: stat  ! Status of every open and read.

            pi = 4.0_mp*atan(1.0_mp)

            !!  Cross sections are evaluated using the CDW or CDW-EIS theory.  !!

            Theory(1) = '   CDW  '
            Theory(2) = ' CDW-EIS'

            !!  Open files for input and output. The file he_data.dat contains the  !!
            !!  detailed formatted output.                                          !!

            open(unit=get1,file='info_total.dat',status='unknown',iostat=stat)
            if (stat .ne. 0) then
                write(iout,*)
                write(iout,*) "Error: could not open file info.dat."
                write(iout,*) "Either the file does not exist or you need to check"
                write(iout,*) "the read/write permissions on your directory."
                write(iout,*) 'Code has terminated.'
                stop
            end if

            !!  Read input data from file.  !!
 
            read(get1,*) target, no, my_section, ZP, Projectile_Energy_Start,     & 
                         Projectile_Energy_Finish, & 
                         Projectile_Energy_Increment, Electron_Energy_Start, &
                         Electron_Energy_Finish, Electron_Energy_Step, angle, tol, acc   
            
            if (target .eq. 1) then  !  Helium orbitals
                open(unit=get2,file='helium.dat',status='old',iostat=stat)
                if (stat .ne. 0) then
                    write(iout,*)
                    write(iout,*) "Error: could not open file helium.dat."
                    write(iout,*) "Either  the file does not exist or you need to check"
                    write(iout,*) "the read/write permissions on your directory."
                    write(iout,*) 'Code has terminated.'
                    stop
                end if
            else
                open(unit=get2,file='hydrogen.dat',status='old',iostat=stat)
                if (stat .ne. 0) then
                    write(iout,*)
                    write(iout,*) "Error: could not open file hydrogen.dat."
                    write(iout,*) "Either  the file does not exist or you need to check"
                    write(iout,*) "the read/write permissions on your directory."
                    write(iout,*) 'Code has terminated.'
                    stop
                end if
            end if

            read(get2,*) m
            num_orb = m  ! Determine the size of the arrays Z_orbital and N_orbital.

            !!  Allocate space for the arrays Z_orbital and N_orbital.  !!

            allocate (Z_orbital(m), N_orbital(m))

            !!  Read He orbitals from file.  !!

            do loop = 1,m
               read(get2,*) Z_orbital(loop), N_orbital(loop)
            end do

            close(unit=get1)  ! Finished reading in all the input data.
            close(unit=get2)
 
       end Subroutine Get_Data

end Module Global_Data





Module Complex_Functions

       !!-----------------------------------------------------------------!!
       !!                                                                 !!
       !!            Description                                          !!
       !!           -------------                                         !!
       !!                                                                 !!
       !!                                                                 !!
       !!  This module contains the following real and complex functions  !!
       !!                                                                 !!
       !!       Sqmodn                                                    !!
       !!                                                                 !!
       !!       Cgamma                                                    !!
       !!                                                                 !!
       !!       Feq                                                       !!
       !!                                                                 !!
       !!       F_n                                                       !!
       !!                                                                 !!
       !!       Cdigam                                                    !!
       !!                                                                 !!
       !!-----------------------------------------------------------------!!

       use Precision, only : mp, c_mp, epmach, tolerance 
       implicit none
      
       real(kind=mp), parameter, private :: pi = 3.141592653589793_mp

       integer, parameter, private :: iout = 6  ! Output channel.


Contains

       real(kind=mp) Function Sqmodn(x)

           !!---------------------------------------------------------------!!
           !!                                                               !!
           !!           Decsription                                         !!
           !!          -------------                                        !!
           !!                                                               !!
           !!  The functiom 'Sqmodn' evaluates the square modulus of N(x).  !!
           !!  Equation (6.1.31) of reference 1 is employed.                !!
           !!                                                               !!
           !!                                                               !!
           !!   Reference                                                   !!
           !!  -----------                                                  !!
           !!                                                               !!
           !!  1)  Abramowitz M & Stegun IA - The Handbook of Mathematical  !!
           !!      Function (Dover: NY 1970).                               !!
           !!                                                               !!
           !!---------------------------------------------------------------!!

           implicit none

           real(kind=mp), intent(in) :: x


           !!  First executable statement of Sqmodn.  !!

           if (x .eq. 0.0_mp) then 
               Sqmodn = 1.0_mp
           else if (x .ge. 7.0_mp) then
               Sqmodn = 2.0_mp*pi*x
           else
               Sqmodn = 2.0_mp*pi*x/(1.0_mp - exp(-2.0_mp*pi*x))
           end if

           return

       end Function Sqmodn




       complex(kind=c_mp) Function Cgamma(z)

           !!---------------------------------------------------------------------!!
           !!                                                                     !!
           !!      Description                                                    !!
           !!     -------------                                                   !!
           !!                                                                     !!
           !!  CGAMMA(Z) calculates the gamma function for complex argument Z     !!
           !!  using  the Lanczos Approximation ( ref 2 , chap. 2 ) .             !!
           !!                                                                     !!
           !!      Details                                                        !!
           !!      -------                                                        !!
           !!                                                                     !!
           !!  For 1>Re(Z)>=0 ( i.e. L=2 ) the factorial properties of the gamma  !!
           !!  function are used ( ref. 1 eqn. 6.1.15 ) and for Re(Z) < 0         !!
           !!                                                                     !!
           !!  ( i.e. L=1 ) analytic continuation is used ( ref.1 eqn. 6.1.17 ).  !!
           !!                                                                     !!
           !!      References                                                     !!
           !!      ----------                                                     !!
           !!                                                                     !!
           !!       1) Abramowitz M and Stegun IA     Handbook of Mathematical    !!
           !!          Functions (Dover NY 1970 ).                                !!
           !!       2) Luke Y L   The Special Functions and Their Approximations  !!
           !!          ( Academic Press NY 1969 ).                                !!
           !!                                                                     !!
           !!---------------------------------------------------------------------!!


           implicit none

           real(kind=mp), parameter :: root_pi = 2.506628274631001_mp  !  The square root of pi.
           real(kind=mp) :: x, y, rem, fk, fk1

           complex(kind=c_mp), intent(in) :: z
           complex(kind=c_mp) :: a, v, h, s, Cgammav

           integer :: l, k

           real(kind=mp), dimension(16) :: G = (/&
               41.624436916439068_mp, -51.224241022374774_mp,+11.338755813488977_mp,   &
               -0.747732687772388_mp, +0.008782877493061_mp, -0.000001899030264_mp,    &
               +0.000000001946335_mp, -0.000000000199345_mp, +0.000000000008433_mp,    &
               +0.000000000001486_mp, -0.000000000000806_mp, +0.000000000000293_mp,    &
               -0.000000000000102_mp, +0.000000000000037_mp, -0.000000000000014_mp,    &
               +0.000000000000006_mp /)

           !!  First executable statement of Cgamma.  !!

           x = real(z,mp)
           y = aimag(z)
           rem = x - nint(x)

           !!  The gamma function is analytic throughout the whole complex   !!
           !!  plane except at the points z=0,-1,-2,... where it has simple  !!
           !!  poles. If an attempt is made to evaluate the gamma function   !!
           !!  at one of these points the program is terminated.             !!

           if (x .le. 0.0_mp .and. abs(rem) .lt. tolerance .and. abs(y) .lt. tolerance) then
              write(iout,*) 'Argument of gamma function is zero'
              write(iout,*) 'or a negative number, therfore the '
              write(iout,*) 'program is terminated.'

              stop  !  Well there's no point in continuing is there?  

           end if

           if (x .ge. 1.0_mp) then
              v = z
              l = 3
           else if (x .ge. 0.0_mp) then
              v = z + 1.0_mp
              l = 2
           else
              v = 1.0_mp - z
              l = 1
           end if

          !!  Cgammav (the gamma function for the argument V) is  !!
          !!  evaluated using the Lancozos approximation.         !!

          h = 1.0_mp
          s = G(1)

          do k = 2,16
             fk = real(k,mp) - 2.0_mp
             fk1 = real(k,mp) - 1.0_mp
             h = ((v - fk1)/(v + fk))*h
             s = s + G(k)*h
          end do

          a = v + 4.5_mp

          Cgammav = root_pi*exp((v - 0.5_mp)*log(a) - a)*s

          !!  Factorial property or analytic continuation is used as appropriate.  !!

          if (l .eq. 3) then
              Cgamma = Cgammav
          else if (l .eq. 2) then
              Cgamma = Cgammav/z
          else
              Cgamma = pi/(sin(pi*z)*Cgammav)
          end if

          return

       end Function Cgamma



       complex(kind=c_mp) Function Feq(Ip,Ihv,zeta,c,x)

          !!--------------------------------------------------------------------------!!
          !!                                                                          !!
          !!     Description                                                          !!
          !!     -----------                                                          !!
          !!                                                                          !!
          !!    Feq(Ip,Ihv,zeta,c,x) evaluates the hypergeometric function            !!
          !!    F(C-1+iZETA,C-1+iZETA,C,X) by using the expansion given by            !!
          !!    eqn.29 of write-up.                                                   !!
          !!    C,X,and zeta  are all real , C must be equal to 1 or 2 and            !!
          !!    mod(x)>1 .                                                            !!
          !!                                                                          !! 
          !!     Other Functions Called                                               !!
          !!     ----------------------                                               !!
          !!                                                                          !!
          !!    Cdigam(Z)  evaluates the  Digamma Function  for complex               !!
          !!    argument Z .                                                          !!
          !!                                                                          !!
          !!     Details                                                              !!
          !!     -------                                                              !!
          !!                                                                          !!
          !!    Eqn.29 of write-up is only valid for mod(arg(-x))<pi but it           !!
          !!    can easily be generalized povided we know whether  - x                !! 
          !!    approaches the real negative axis from above or below ;               !!
          !!    ie whether we should write -1 as EXP(i*pi) or EXP(-i*pi).             !!
          !!    The arguments Ip and Ihv deal specifically with this ;                !! 
          !!                                                                          !!
          !!    (I)   if x>1 { ie mod(arg(-x)) = pi } then                            !! 
          !!          Ihv=1 and Ip=-1 if -x approaches the -ve real axis from above   !! 
          !!          Ihv=1 and Ip=-1 "   "     "       "   "    "    "    "  below   !!
          !!                                                                          !!
          !!   (II)   if x<-1 { ie mod(arg(-x)) < pi } then                           !!
          !!          Ihv=-1 and Ip=0  ( ie the standard formula )                    !!
          !!                                                                          !!
          !!    Also note that Euler is Eulers constant .                             !!
          !!                                                                          !!
          !!     References                                                           !!
          !!     ----------                                                           !!
          !!    Abramowitz M and Stegun IA  Handbook of Mathematical                  !!
          !!    Functions (Dover NY 1970)  .                                          !!
          !!                                                                          !!
          !!--------------------------------------------------------------------------!!

           implicit none

           real(kind=mp), parameter :: Euler = 0.57721566490153_mp  !  Eulers constant.
           real(kind=mp), intent(in) :: x, c, zeta
           real(kind=mp) :: r

           complex(kind=c_mp), parameter :: i = cmplx(0.0_mp,1.0_mp,c_mp)  !  The complex number i.
           complex(kind=c_mp) :: a, t1, t2, t3, fac, pfac, my_sum, term

           integer, intent(in) :: Ip, Ihv
           integer :: n

           !!  First executable statement of function Feq.  !!

           a = c - 1.0_mp + i*zeta
           t1 = log(Ihv*x) - pi*i*Ip
           t2 = -2.0_mp*Euler - Cdigam(a) - Cdigam(c - a)
           t3 = cmplx(1.0_mp,0.0_mp,c_mp)

           !!  Abramowitz and Stegun eqn. (6.1.31) can be used to express  !!
           !!  the product of gamma functions in terms of sinhs.           !!

           if (c .eq. 1.0_mp) then
               pfac = i*sinh(pi*zeta)/pi
           else

               !! ie. c = 2.  !!

               pfac = sinh(pi*zeta)/(pi*zeta)
           end if

           fac = pfac*exp(i*pi*a*Ip)*(Ihv*x)**(-a)
           my_sum = t1 + t2
 
           n = 0  !  Initialise counter. 

           do
              n = n + 1
              r = real(n,mp)
              t3 = t3*(a + r - 1.0_mp)*(r + a - c)/(x*r**2)
              t2 = t2 + 2.0_mp/r - 1.0_mp/(a + r - 1.0_mp) - 1.0_mp/(a - c + r)
              term = t3*(t2 + t1)
              my_sum = my_sum + term

              if (abs(term) .lt. tolerance) then
                  exit  !  Convergence achieved, exit from do loop. 
              else if (n .gt. 10000) then

                  !!  Safeguard against infinite loop. !!

                  write(iout,*) 'Convergence not achieved in function Feq'
                  write(iout,*) 'program has been terminated.'
                  stop
              end if

           end do

           Feq = my_sum*fac

           return

       end Function Feq


       complex(kind=c_mp) Function F_n(A,B,C,x)

           !!------------------------------------------------------------------!!          
           !!                                                                  !!
           !!        Description                                               !!
           !!        -----------                                               !!
           !!                                                                  !!
           !!  F_n(A,B,C,X) calculates the  hypergeometric series for complex  !! 
           !!  A,B,C,and real X (see eqn.25 of write-up) .                     !! 
           !!  The series is truncated when the last term in the expansion is  !! 
           !!  less than the tolerence TOL .                                   !!
           !!  The series is absolutely convergent provided, Re(C-A-B)>0 ,     !!
           !!  and -1<X<1 .                                                    !!
           !!                                                                  !!
           !!                                                                  !!
           !!------------------------------------------------------------------!!

           implicit none

           real(kind=mp), intent(in) :: x 
           real(kind=mp) :: n

           complex(kind=c_mp), intent(in) :: A, B, C
           complex(kind=c_mp) :: s1, t1

           !!  First executable statement of function F_n.  !!

           s1 = (0.0_mp,0.0_mp) 
           t1 = (1.0_mp,0.0_mp)
   
           n = 1.0_mp

           do
              s1 = s1 + t1
              t1 = t1*x*(A + n - 1.0_mp)*(B + n - 1.0_mp)/((C + n - 1.0_mp)*n)

              if (abs(t1) .le. tolerance) then
                  exit  !  Convergence achieved, exit from loop.
              else if (n .gt. 10000.0_mp) then

                  !!  Safeguard against infinite loop.  !!

                  write(iout,*) 'Converence not achieved in function F_n'
                  write(iout,*) 'program is terminating.'
                  stop
              end if
      
              n = n + 1.0_mp

           end do

           F_n = s1

           return

       end Function F_n



       complex(kind=c_mp) Function Cdigam(z)

           !!---------------------------------------------------------------------!!
           !!                                                                     !!
           !!        Description                                                  !!    
           !!        -----------                                                  !!    
           !!                                                                     !!
           !!    Cdigam(Z)   calculates the  Digamma Function  for complex        !!
           !!    argument Z  using an asymptotic expansion ( ref.1 eqn 6.3.18 ).  !! 
           !!                                                                     !!
           !!        Details                                                      !!    
           !!        -------                                                      !!    
           !!                                                                     !! 
           !!    If abs(Re(Z))<15 the recurrence property of the Digamma          !!
           !!    Function ( ref.1 eqn 6.3.6 ) is used to retain accuracy  when    !!
           !!    using the asymptotic expansion.If Re(Z)<0 analytic continuation  !!
           !!    is  employed ( ref.1 eqn 6.3.7 ).                                !!
           !!                                                                     !! 
           !!        References                                                   !!    
           !!        ----------                                                   !!   
           !!                                                                     !!
           !!    1) Abramowitz M and Stegun IA    Handbook of Mathematical        !!
           !!       Functions (Dover NY 1970) .                                   !!
           !!                                                                     !!
           !!---------------------------------------------------------------------!! 

           implicit none

           real(kind=mp) :: x, y, a, rem, eps

           complex(kind=c_mp), intent(in) :: z
           complex(kind=c_mp) :: v, h, r, Cdigamv, piz

           integer :: n, i

           real(kind=mp), dimension(6) :: B = (/ &
                   +8.333333333333333e-2_mp, -8.333333333333333e-3_mp,  &
                   +3.968253968253968e-3_mp, -4.166666666666667e-3_mp,  &
                   +7.575757575757576e-3_mp, -2.109279609279609e-2_mp /)

           !!  First executable statement of function Cdigam.  !!

           x = real(z,mp)
           y = aimag(z)
           a = abs(x)
           rem = x - nint(x)

           eps = tolerance ! or what value we require.

           !!  The diagamma function is analytic throughout the whole complex    !!
           !!  plane except at the points z = 0, -1, -2,... where it has poles.  !!
           !!  If an attempt is made to evaluate the diagamma function at one    !!
           !!  of these points the program is terminated.                        !!

           if (x .le. 0.0_mp .and. abs(rem) .lt. eps .and. abs(y) .lt. eps) then
               write(iout,*) 'Arguement of diagamma function is zero'
               write(iout,*) 'or a negative integer ie. there is a pole.'
               stop
           end if

           !!  Cdigamv (the diagamma function for complex arguement v) is  !!
           !!  evaluated using the asymptotic expansion. If Re(v)<15 the   !!
           !!  recurrenve property of the diagamma function is used.       !!

           if (x .lt. 0.0_mp) then 
               v = -z
           else 
               v = z
           end if

           h = 0.0_mp 

           if (a .lt. 15.0_mp) then 
              n = 14 - nint(a)
              h = 1.0_mp/v

              if (n .gt. 0) then
                  do i = 1,n,2
                     v = v + 1.0_mp
                     h = h + 1.0_mp/v
                     v = v + 1.0_mp
                     h = h + 1.0_mp/v
                  end do
              end if

              v = v + 1.0_mp
           end if
              
           r = 1.0_mp/v**2.0_mp  

           cdigamv = log(v) - 0.5_mp/v - r*(B(1) + r*(B(2) + r*(B(3) +  &
                     r*(B(4) + r*(B(5) + r*(B(6) + r*B(1)))))))- h

           !!  If Re(z) < 0 then analytic continuation is used.  !!

           if (x .lt. 0.0_mp) then 
              piz = -pi*z
              Cdigam = Cdigamv - 1.0_mp/z + pi*cos(piz)/sin(piz)
           else
              Cdigam = Cdigamv
           end if

           return

       end Function Cdigam 

end Module Complex_Functions



Module Cross_Section
        
       use Precision
       use Global_Data
       use Complex_Functions

       real(kind=mp) :: ZT, velocity, global_k, global_p 

       real(kind=mp), private :: nu, smnnu, smnxi, delta_epsilon, q_min,    &
                                 q_dot_v, snt, cst, k_dot_v, p_dot_v,       &
                                 zeta, AB, xi, csp

       complex(kind=c_mp), private :: g1, g2, g3, g4, g5, g6  ! Complex gamma functions.

Contains

       Subroutine Total_Differential(Model_Energy,TDCS)
     
           !!---------------------------------------------------------------------!!
           !!                                                                     !!
           !!     Description                                                     !!
           !!    -------------                                                    !!
           !!                                                                     !!
           !!  The subroutine Total_Differential returns the value TDCS when      !!
           !!  called where;                                                      !!
           !!                                                                     !!
           !!      TDCS =  total cross section if 'my_section' = 0.               !!
           !!           =  single differential cross section if 'my_section' = 1. !!
           !!           =  double differential cross section if 'my_section' = 2. !!
           !!                                                                     !!
           !!     Details                                                         !!
           !!    ---------                                                        !! 
           !!                                                                     !!
           !!  Since the single differential cross section is calculated as a     !!
           !!  function of ejected electron energy then we set                    !!
           !!                                                                     !!
           !!            Ejected_Electron_Energy = Model_Energy                   !!
           !!                                                                     !!
           !!  and by deafult we set                                              !!
           !!                                                                     !!    
           !!            Projectile_Energy = Projectile_Energy_Start.             !!
           !!                                                                     !! 
           !!  This is also the case for the double differential cross section;   !!
           !!  additionally, the DDCS is also a function of theta therefore we    !!
           !!  use the value 'angle' read from file.                              !!
           !!        Others values of importance are                              !!
           !!                                                                     !!
           !!    global_k (real)  -  Global variable; the ejected electron        !!
           !!                        velocity in au.                              !!
           !!                                                                     !!
           !!    velocity (real)  -  Module variable; the velocity of the         !!
           !!                        projectile ion in au.                        !!
           !!                                                                     !!
           !!    nu (real)        -  Module variable; defined by eqn. ? of the    !!
           !!                        write up.                                    !!
           !!                                                                     !!
           !!    aa (real)        -  An amalgamation of constants that can be     !!
           !!                        taken outside the integration.               !!
           !!                                                                     !!
           !!    k_max (real)     -  Upper limit for the outer integral.          !!
           !!                                                                     !!
           !!---------------------------------------------------------------------!!


           implicit none 

           real(kind=mp), intent(out) :: TDCS 
           real(kind=mp), intent(in) :: Model_Energy
           
           real(kind=mp), allocatable, dimension(:) :: work

           real(kind=mp) :: k_max, ss, aa, Q_Total, epsabs, epsrel,       &
                            start_1, finish_1, abserr, answer, k_dummy,   &
                            Projectile_Energy, Ejected_Electron_Energy
                            
           real(kind=mp), parameter :: Bohr = 0.529177_mp,        &   !  The Bohr radius
                                       delta_v = 0.0000001_mp

           integer, allocatable, dimension(:) :: iwork
           integer :: key, neval, ier, last, limit, lenw, all_status, loop

           !!  First executable statement of Total_Differential.  !!
           
           if (my_section .gt. 0) then  ! single or double differential cross sections.
               Projectile_Energy = Projectile_Energy_Start
               Ejected_Electron_Energy = Model_Energy
               global_k = (2.0_mp*Ejected_Electron_Energy/Rydberg)**0.5_mp
           else if (my_section .eq. 0) then  ! Total cross sections.
               Projectile_Energy = Model_Energy
           end if

           select case(target)
               case(0)
                 ZT = 1.0_mp  ! Hydrogen-like target.
               case(1)
                 ss = 2.0_mp*(24.58_mp/Rydberg)
                 ZT = ss**0.5_mp   ! The effectivc charge of a helium-like target.
               case(2)
                 ZT = 1.064_mp  !  H2-like target.
           end select

           velocity = (Projectile_Energy/24.8018418_mp)**0.5_mp  ! Velocity of incoming ion.

           !!  We must avoid the singularity that occurs when then velocity of  !!
           !!  the projectile ion = momentum of the ejected electron in the     !!
           !!  forward direction (theta = 0).                                   !!

           if (my_section .eq. 2 .and. abs(angle) .lt. tol) then
               if (abs(velocity - global_k) .lt. tol) then
                   if (velocity .lt. global_k) then
                       velocity = velocity - delta_v
                   else
                       velocity = velocity + delta_v
                  end if
               end if
           end if                 
           
           aa = 2.0_mp*Bohr*Bohr*pi

           nu = ZP/velocity
           smnnu = Sqmodn(nu)
           g4 = Cgamma(cmplx(1.0_mp,-nu,c_mp))
           g6 = conjg(g4)/cmplx(0.0_mp,nu,c_mp)
       
           if (my_section .eq. 0) then 
               k_max = (Projectile_Energy*0.5_mp)**0.5_mp   ! Upper limit in the integral 
                                                            ! over k. Only needed for the
                                                            ! TDCS since the SDCS is 
                                                            ! evaluated with respect to 
                                                            ! ejected electron energy.
           end if

           select case(my_section)        
               case(2)  ! double differential cross sections.
                  k_dummy = global_k
                  answer = DD1(k_dummy)
                  TDCS = answer*aa/(2.0_mp*pi*Rydberg)
               case(1)  ! single differential cross sections.
                  k_dummy = global_k
                  answer = DD1(k_dummy)
                  TDCS = answer*aa/Rydberg
               case(0)  ! Total differential cross sections.

                  !!  Prepare for call to Quad4.  !!

                  key = 6
                  limit = global_limit
                  lenw = 4*limit

                  !!  Allocate space for the array's work and iwork.  !!
 
                  allocate (work(lenw), stat=all_status)
                  if (all_status .ne. 0) then
                      write(iout,*)
                      write(iout,*) "Error: Can't allocate enough memory for the array work."
                      stop
                  end if

                  allocate (iwork(limit), stat=all_status)
                  if (all_status .ne. 0) then
                      write(iout,*)
                      write(iout,*) "Error: Can't allocate enough memory for the array iwork."
                      stop
                  end if

                  epsabs = 0.0_mp
                  epsrel = acc
		   
                  start_1 = 0.0_mp
                  finish_1 = k_max
                  
                  call Quad4(DD1,start_1,finish_1,epsabs,epsrel,key,Q_Total,   &
                             abserr,neval,ier,limit,lenw,last,iwork,work) 

                  deallocate (work,iwork)

                  TDCS = Q_Total*aa*2.0_mp   ! The total differential cross section.
           end select
           
           return

       end Subroutine Total_Differential



       real(kind=mp) Function DD1(dummy_k)

           !!---------------------------------------------------------------------!!
           !!                                                                     !!
           !!       Description                                                   !!
           !!      -------------                                                  !!
           !!                                                                     !!
           !!  If 'my_section' = 0 or 1 (ie. evaluation of total and single       !!
           !!  differential cross sections) DD1 integrates the function DD2       !!
           !!  over a range 0 to pi for a given ejected electron velocity.        !!
           !!  If 'my_section' = 2 (double differential cross section) then DD1   !!
           !!  evaluates DD2 at a particular value 'angle' (a global parameter)   !!
           !!  and 'dummy_k' (= ejected electron velocity).                       !!
           !!                                                                     !!
           !!                                                                     !!
           !!       Details                                                       !!
           !!      ---------                                                      !!
           !!                                                                     !!
           !!  In particular, we note that                                        !!
           !!                                                                     !!
           !!     delta_epsilon (real)  -  Module variable; refer to eqn. ? of    !!
           !!                              the write up.                          !!
           !!                                                                     !!
           !!     q_dot_v (real)        -  Module variable; the dot product of    !!
           !!                              vectors q and v where q is defined     !!
           !!                              by eqn. ? of write up and v=velocity   !!
           !!                              is the velocity of the projectile ion. !!
           !!                                                                     !!
           !!     q_min (real)          -  Module variable; the lower limit of    !!
           !!                              integral given by eqn ? of write up.   !!
           !!                                                                     !!
           !!     xi (real)             -  Module variable; refer to eqn ? of     !!
           !!                              write up.                              !!
           !!                                                                     !!
           !!---------------------------------------------------------------------!!
 

           implicit none

           real(kind=mp), intent(in) :: dummy_k
           real(kind=mp), allocatable, dimension(:) :: work
           
           real(kind=mp) :: answer, epsabs, epsrel, start_2, finish_2, abserr, &
                            dummy_angle

           integer, allocatable, dimension(:) :: iwork
           integer :: key, neval, ier, last, limit, lenw, all_status

           if (my_section .eq. 0) then  ! For TDCS only.
               global_k = dummy_k
           end if
        
           if (dummy_k .eq. 0.0_mp) then 
              DD1 = 0.0_mp
              return
           end if

           xi = ZT/dummy_k
           smnxi = Sqmodn(xi)

           delta_epsilon = (ZT*ZT + dummy_k*dummy_k)*0.5_mp
           q_min = delta_epsilon/velocity
           q_dot_v = -delta_epsilon 
        
           select case(my_section)
               case(2)  ! double differential cross sections.
                   dummy_angle = angle
                   DD1 = DD2(dummy_angle)*dummy_k
               case default  ! Total and single differential cross sections. 

                   !!  Prepare for call to Quad3.  !!

                   !!  The function we are about to integrate contains a possible  !!
                   !!  large peak (EEC peak) and therefore it is recommended that  !!
                   !!  the Gauss-Kronrod 15-21 point rule is used.                 !!
                   
                   key = 1  !   Gauss-Kronrod 15-21 point rule.
                   limit = global_limit - 100
                   lenw = 4*limit

                   !!  Allocate space for the array's work and iwork.  !!

                   allocate (work(lenw), stat=all_status)
                   if (all_status .ne. 0) then 
                       write(iout,*)
                       write(iout,*) "Error: Can't allocate enough memory for the array work."
                       stop
                   end if

                   allocate (iwork(limit), stat=all_status)
                   if (all_status .ne. 0) then
                       write(iout,*)
                       write(iout,*) "Error: Can't allocate enough memory for the array iwork."
                       stop
                   end if

                   epsabs = 0.0_mp
                   epsrel = acc*0.1_mp

                   start_2 = 0.0_mp
                   finish_2 = pi

                   call Quad3(DD2,start_2,finish_2,epsabs,epsrel,key,answer,      &
                              abserr,neval,ier,limit,lenw,last,iwork,work)

                   deallocate (work,iwork)

                   if (my_section .eq. 1) then
                       DD1 = answer*dummy_k
                   else if (my_section .eq. 0) then 
                       DD1 = answer*dummy_k*dummy_k
                   end if

           end select

           return        

       end Function DD1 




       real(kind=mp) Function DD2(polar)

           !!-----------------------------------------------------------------!!
           !!                                                                 !!
           !!         Description                                             !!
           !!        -------------                                            !!
           !!                                                                 !! 
           !!  The function DD2(polar) integrates the function DD3(phi) over  !!
           !!  the range 0 to pi for a given angle 'polar' (polar is the      !!
           !!  polar angle of the wave vector k wrt the polar axis vector v). !!
           !!                                                                 !!
           !!                                                                 !!
           !!         Details                                                 !!
           !!        ---------                                                !!
           !!                                                                 !!
           !!  In particular the global variables declared here are :         !!
           !!                                                                 !!
           !!    k_dot_v (real)  -  vector k dot vector v.                    !!
           !!                                                                 !!
           !!    p (real)        -  the magnitude of the vector p ie the      !!
           !!                       momentum of the electron relative to the  !!
           !!                       projectile.                               !!
           !!                                                                 !!
           !!    p_dot_v (real)  -  vector p dot vector v.                    !!
           !!                                                                 !!
           !!    zeta (real)     -  ZP/p.  Projectile charge divided by the   !!
           !!                       momentum p.                               !!
           !!                                                                 !!
           !!  It should be noted that the integral over the azimuthal angle  !!
           !!  'phi' is reduced from o to 2*pi to 0 to pi by using the fact   !!
           !!  that the integrand is a function of cos(phi) only. Hence the   !!
           !!  integral from 0 to 2*pi of DD3(cos(phi)) can be transformed to !!
           !!  0 to pi of 2*DD3(-cos(phi)); this is why csp is defined to be  !!
           !!  -cos(phi) in the function DD3.                                 !!
           !!                                                                 !!
           !!-----------------------------------------------------------------!!

           implicit none

           real(kind=mp), intent(in) :: polar      
           real(kind=mp), allocatable, dimension(:) :: work
           real(kind=mp) :: epsabs, epsrel, start, finish, result,   &
                            abserr, smnze, p

           integer :: key, neval, ier, last, limit, lenw, status
           integer, allocatable, dimension(:) :: iwork
           
           !!  First executable statement of DD2.  !!

           if (polar .eq. pi) then
               DD2 = 0.0_mp
               return
           end if

           cst = cos(polar)
           snt = sin(polar)
           k_dot_v = global_k*velocity*cst
           p_dot_v = (k_dot_v - velocity*velocity)
           p = (global_k*global_k + velocity*velocity - 2.0_mp*k_dot_v)**0.5_mp
           global_p = p

           !!  There is a simple pole in DDCS which is often referred to as the     !!
           !!  electron capture to the continuum (EEC) peak. This pole only occurs  !!
           !!  in the forward direction (theta=0) when then velocity of the         !!
           !!  incoming ion = momentum of the ejected electron. In such a case p is !!
           !!  zero and the singularity is given by zeta below. Obviously we cannot !!
           !!  allow this to happen so we must give a lower bound for p.            !! 
           
           if (p .lt. tolerance) then
               p = 1.e-7_mp     !  The user is free to change this value, but bear in 
           end if               !  mind that extremely small values of p can lead to
                                !  round-off.

           zeta = ZP/p
           smnze = Sqmodn(zeta)

           g1 = Cgamma(cmplx(0.0_mp,nu - zeta,c_mp))
           g2 = Cgamma(cmplx(0.0_mp,nu + zeta,c_mp))
           g3 = Cgamma(cmplx(1.0_mp,-zeta,c_mp))
           g5 = conjg(g3)/cmplx(0.0_mp,zeta,c_mp)

           AB = smnnu*smnze*smnxi

           !!  Prepare for call to QUAD2.  !!

           key = 6
           limit = global_limit - 200
           lenw = 4*limit

           !!  Allocate space for the array's work and iwork.  !!

           allocate (work(lenw),stat=status)
           if (status .ne. 0) then
               write(iout,*)
               write(iout,*) "Error: Can't allocate memory for the array work."
               stop
           end if 

           allocate (iwork(limit),stat=status)
           if (status .ne. 0) then
               write(iout,*)
               write(iout,*) "Error: Can't allocate memory for the array iwork."
               stop
           end if

           epsabs = 0.0_mp
           epsrel = acc*0.01_mp 
  
           start = 0.0_mp
           finish = pi

           call QUAD2(DD3,start,finish,epsabs,epsrel,key,result,           &
                      abserr,neval,ier,limit,lenw,last,iwork,work)

           select case(my_section)
               case(2)  ! double differential cross sections.  
                   DD2 = 2.0_mp*result
               case default  ! Total differential cross sections.
                   DD2 = 2.0_mp*result*snt
           end select

	   deallocate (work, iwork)

           return

      end Function DD2 




      real(kind=mp) Function DD3(phi)

           !!-----------------------------------------------------------!!
           !!                                                           !!
           !!      Decription                                           !!
           !!     ------------                                          !!
           !!                                                           !!
           !!  DD3 integrates the function CS over the range -1 to +1   !!
           !!  for a given phi (phi is the azimuthal angle of the wave  !!
           !!  vector k with respect to the collision plane).           !!
           !!                                                           !!
           !!-----------------------------------------------------------!!

           implicit none
       
           real(kind=mp), intent(in) :: phi
           real(kind=mp) :: start, finish, epsabs, epsrel, result, abserr
           real(kind=mp), allocatable, dimension(:) :: work

           integer :: key, neval, ier, limit, lenw, last, status
           integer, allocatable, dimension(:) :: iwork

           csp = -cos(phi)

           !!  Prepare for call to QUAD1.  !!

           limit = global_limit - 300
           lenw = 4*limit

           !!  Allocate space for the array's work and iwork.  !!

           allocate (work(lenw),stat=status)
           if (status .ne. 0) then
               write(iout,*)
               write(iout,*) "Error: Can't allocate memory for the array work."
               stop
           end if

           allocate (iwork(limit),stat=status)
           if (status .ne. 0) then
               write(iout,*)
               write(iout,*) "Error: Can't allocate memory for the array iwork."
               stop
           end if

           epsabs = 0.0_mp
           epsrel = acc*0.001_mp

           key = 6

           start = -1.0_mp
           finish = 1.0_mp

           call QUAD1(CS,start,finish,epsabs,epsrel,key,result,            &
                      abserr,neval,ier,limit,lenw,last,iwork,work)
  
           DD3 = result

           deallocate (work, iwork)

           return

      end Function DD3          




      real(kind=mp) Function CS(x)

           !!----------------------------------------------------------------------!!
           !!                                                                      !!
           !!      Description                                                     !!
           !!     -------------                                                    !!
           !!                                                                      !!
           !!  CS(x) is the central integrand.                                     !!
           !!                                                                      !!
           !!                                                                      !!
           !!      Details                                                         !!
           !!     ---------                                                        !!
           !!                                                                      !!
           !!  1)  q_dot_k (real) -  The dot product of vectors q and k.           !!
           !!                                                                      !!
           !!  2)  q_dot_p (real) -  The dor product of vectors q and p.           !!
           !!                                                                      !!
           !!  3)  alpha, beta, gamma and delta are given by eqn ? of write up.    !!
           !!                                                                      !!
           !!  4)  tau (real)     -  Given by eqn. ? of write up if no = 1 (CDW).  !!
           !!                  Given by eqn. ? of write up if no = 2 (CDW-EIS).    !!
           !!                                                                      !! 
           !!  5)  test (real)    -  There are restrictions on the four analytic   !!
           !!                  continuations (eqns. ? to ? of write up) which are  !!
           !!                  used in evaluating the Gauss Hypergeometric         !!
           !!                  function in the functions V1, W1 and Feq.           !!
           !!                      Equations ? and ? are restricted to             !!
           !!                  mod(arg(1-z)) < pi and ? and ? restricted to        !!
           !!                  mod(arg(-z)) < pi. These restrictions can be        !!
           !!                  removed provided we know whether (1-z) and (-z)     !!
           !!                  approach the -ve real axis from above or below (ie. !!
           !!                  whether we should write -1 as exp(i*pi) or          !!
           !!                  exp(-i*pi)). In these cases the sign of test        !!
           !!                  indicates which way we should write -1; if test > 0 !!
           !!                  then z (in this case tau) is approaching the -ve    !!
           !!                  real axis from below and if test < 0 then z (tau)   !!
           !!                  is approaching the -ve real axis from above.        !!
           !!                                                                      !!
           !!  6)  omega (complex)   -  Given by eqn. ? or ? of write depending on !!
           !!                           whether no = 1 (CDW) or n = 2 (CDW-EIS).   !!
           !!                                                                      !!
           !!----------------------------------------------------------------------!!

           implicit none
           
           real(kind=mp), intent(in) :: x
           real(kind=mp) :: QSQ, q, csb, snb, q_dot_k, q_dot_p,                 &
                            ba, Aqk, alpha, beta, gamma, alga1,                 &
                            delta, pi_v, pam,                                   &
                            XX(num_orb), Alp(num_orb), alpbe, tau, test,        &
                            alpha_beta, Y(num_orb) 

           complex(kind=c_mp) :: cc, ZZ, pp, A1, AA2, A2, A3(num_orb),               &
                                 T(num_orb), Bet(num_orb), F(num_orb), G(num_orb),   &
                                 D(num_orb), Omega(num_orb)

           integer :: lambda, I2, Ip, Ihv, Ipw, Ihw

           if ( (1.0_mp + X) .lt. tolerance) then
               CS = 0.0_mp   !  If (1+X) < tolerance then we run the risk of dividing by 
                             !  a number very close to zero which can lead to huge 
                             !  round-off problems.
               return
           end if

           QSQ = q_min**2 + (1.0_mp - X)/(1.0_mp + X)
           q = QSQ**0.5_mp

           csb = q_min/q
           snb = ((1.0_mp - csb)*(1.0_mp + csb))**0.5_mp
           q_dot_k = -q*global_k*(cst*csb + snt*snb*csp)
           q_dot_p = q_dot_k - q_dot_v

           ba = QSQ - 2.0_mp*delta_epsilon

           Aqk = 1.0_mp

           if (no .eq. 1 .and. ba .lt. 0.0_mp) then
               Aqk = exp(-2.0_mp*pi*nu)
           end if

           if (no .eq. 2) then
               Aqk = exp(-2.0_mp*pi*nu)
           end if 

           ZZ = cmplx(1.0_mp,xi,c_mp)
           pp = cmplx(0.0_mp,xi,c_mp)

           alpha = 0.5_mp*QSQ
           beta = q_dot_v
           gamma = q_dot_p + alpha
           alga1 = (alpha + beta)*gamma
           delta = p_dot_v - global_p*velocity + q_dot_v
           pam = 1.0_mp/(alpha*alpha*gamma*gamma)
           pi_v = 1.0_mp/(2.0_mp*pi*pi*velocity*velocity)

           A1 = cmplx(0.0_mp,0.0_mp,c_mp)
           AA2 = cmplx(0.0_mp,0.0_mp,c_mp)

           XX = Z_orbital/global_k
           T = cmplx(1.0_mp,XX,c_mp)
           Alp = 0.5*(QSQ + global_k**2 + 2*q_dot_k + &
                            Z_orbital**2)
           Bet = -(q_dot_k + global_k*global_k*T)
           Y = QSQ + Z_orbital*Z_orbital - global_k*global_k
           F = (Alp/(Alp + Bet))**pp

           where (Y .lt. 0.0_mp)
               F = F*exp(-pi*XX)
           end where

           G = QSQ + (q_dot_k*T)
           D = delta_epsilon - T*k_dot_v +    &
               (delta_epsilon + q_dot_k + T*  &
               (global_k*global_k - k_dot_v))*velocity/global_p 
           Omega = gamma*D + delta*G       

           A3 = (N_orbital*Z_orbital**1.5_mp)*F/(Alp*(Alp + Bet))

           do lambda = 1,num_orb
              A1 = A1 + A3(lambda)*G(lambda)
              AA2 = AA2 + A3(lambda)*Omega(lambda)
           end do


           if (no .eq. 1) then
               A2 = (alpha/alga1)*AA2
           else
               A2 = (alpha/(beta*gamma))*AA2
           end if

           if (no .eq. 1) then
               alpha_beta = alpha + beta

               if ( abs(alpha_beta) .lt. tolerance) then
                   alpha_beta = tolerance
               end if

               tau = (beta*gamma - alpha*delta)/(gamma*alpha_beta)
               test = global_p*(alpha_beta)*delta +    &
                      velocity*gamma*(gamma + delta - alpha_beta)
           else
               tau = (beta*gamma - alpha*delta)/(gamma*beta)
               test = global_p*beta*delta + velocity*gamma*(delta - beta)
           end if

           I2 = 1

           if (tau .lt. 0.0_mp) then
               Ip = 0
               Ihv = -1
           else

              if (test .gt. 0.0_mp) then
                 I2 = -1
              end if

             Ip = I2
             Ihv = 1
           end if

           if (tau .lt. 1.0_mp) then
              Ipw = 0
              Ihw = -1
           else
              Ihw = 1
              Ipw = Ip
           end if

           if ( (tau .ge. 1.5_mp .or. tau .lt. -2.0_mp) .and.             &
                 (abs(nu - zeta) .lt. 1.e-3_mp) ) then
                      cc = A1*Feq(Ip,Ihv,zeta,1.0_mp,tau) -                &
                      cmplx(0.0_mp,nu,c_mp)*A2*Feq(Ip,Ihv,zeta,2.0_mp,tau)
           else
                      cc = A1*V1(Ip,Ipw,Ihv,Ihw,zeta,nu,tau) - &
                      cmplx(0.0_mp,nu,c_mp)*A2*W1(Ip,Ipw,Ihv,Ihw,zeta,nu,tau)

           end if

           CS = ZP**2*ZT**2*AB*pi_v*pam*Aqk*(abs(cc)/(1.0_mp + X))**2

           return

       end Function CS


        
       complex(kind=c_mp) Function V1(Ip,Ipw,Ihv,Ihw,zeta,nu,z)

           implicit none

           real(kind=mp) :: zinv, w
           real(kind=mp), intent(in) :: z, zeta, nu

           complex(kind=c_mp) :: A, B, C, U
           complex(kind=c_mp), parameter :: i = cmplx(0.0_mp, 1.0_mp,c_mp)  ! The complex number i.

           integer, intent(in) :: Ip, Ipw, Ihv, Ihw

           !!  First executable statement of V1.  !!

           A = cmplx(0.0_mp,zeta,c_mp)
           B = cmplx(0.0_mp,nu,c_mp)
           C = cmplx(1.0_mp,0.0_mp,c_mp)

           if (z .ge. 1.5_mp .or. z .lt. -2.0_mp) then

               !!  Use analytic continuation.  !!

               zinv = 1.0_mp/z
               V1 = (g1/(g6*g3))*(Ihv*z)**(-A)*exp(i*A*pi*Ip)            &
                    *F_n(A,1.0_mp - C + A,1.0_mp - B + A,zinv)           &
                    + (conjg(g1)/(g5*g4))*(Ihv*z)**(-B)*exp(i*B*pi*Ip)   &
                    *F_n(B,1.0_mp - C + B,1.0_mp - A + B,zinv)

           else if (z .ge. -2.0_mp .and. z .lt. -0.666_mp) then

               !!  Use analytic continuation.  !!
 
               V1 = (1.0_mp - z)**(-A)*F_n(A,C-B,C,Z/(Z - 1.0_mp))

           else if (z .ge. 0.5_mp .and. z .lt. 1.5_mp) then
               W = z - 1.0_mp

               if (abs(w) .lt. tolerance) then

                  !!  Use Abramowitz and Stegun eqn. (15.1.20) to  !!
                  !!  evaluate F(A,B;C;1).                         !!

                  V1 = conjg(g2)*cmplx(0.0_mp,-nu-zeta,c_mp)/(g3*g4)
               else
                  U = C - A - B
                  V1 = conjg(g2)*cmplx(0.0_mp,-nu-zeta,c_mp)/(g3*g4)     &
                       *F_n(A,B,A+B-C+1.0_mp,-w)                         &
                       + (Ihw*w)**U*exp(-Ipw*pi*i*U)*g2*(-nu*zeta)/      &
                       (conjg(g3*g4)*cmplx(-1.0_mp,nu+zeta,c_mp))        &
                       *F_n(C-A,C-B,C-A-B+1.0_mp,-w)
               end if

           else

               !!  If z > -0.66 and z < 0.5 then                    !!
               !!  the function F_n(A,B,C,Z) can be used directly.  !!

               V1 = F_n(A,B,C,z)
           end if
        
           return

       end Function V1


       complex(kind=c_mp) Function W1(Ip,Ipw,Ihv,Ihw,zeta,nu,z)

           !!----------------------------------------------------------------------------!!
           !!                                                                            !!
           !!      DESCRIPTION                                                           !!
           !!      -----------                                                           !!
           !!                                                                            !!
           !!      W1 evaluates the hypergeometric function F(1+iNU,1+iZETA,2,Z)         !!
           !!      for REAL values of NU, ZETA, and , Z. ( i = SQRT(-1) )                !!
           !!                                                                            !!
           !!                                                                            !!
           !!                  OTHER FUNCTIONS CALLED                                    !!
           !!                  ----------------------                                    !!
           !!                                                                            !! 
           !!      W1  calls the function F(A,B,C,X) . F(A,B,C,X) evaluates the          !!
           !!      hypergeometric function for COMPLEX A, B, and C , and REAL X          !!
           !!      provided Re(C-A-B)>0 and -1<X<+1.                                     !!
           !!                                                                            !! 
           !!                                                                            !!
           !!                       DETAILS                                              !!
           !!                       -------                                              !!
           !!                                                                            !!
           !!                                                                            !!
           !!         The appropriate analytical continuations of the                    !! 
           !!      hypergeometric functions ( eqns.26-29 of write-up )                   !!
           !!      are used in conjunction with F(A,B,C,X) to ensure that                !!
           !!      the hypergeometric function can be evaluated for all Z .              !!
           !!      Equations 26 and 27 are restricted to mod(arg(1-Z))<pi                !!
           !!      and 28 and 29 are restricted to mod(arg(-Z))<pi.                      !!
           !!      These restrictions can be removed provided we know whether            !!
           !!      (1-Z) and -Z approach the -ve real axis form above or below           !!
           !!      ie whether we should write -1 as CEXP(i*PI) or CEXP(-i*PI).           !!
           !!      This information is provided thru the      arguments IP, IHV,         !!
           !!      IPW , AND IHW. Their values are set as follows ;                      !!
           !!                                                                            !!
           !!                                                                            !!
           !!      (Ia) if Z.GE.0 ( ie mod(arg(-Z))= pi ) then                           !!
           !!             IHV=1 and IP=-1 if -Z approaches the -ve real axis from below  !!
           !!         OR  IHV=1 and IP=+1 "  "      "       "   "    "    "    "  above  !!
           !!
           !!      (Ib) if Z.LT.0 ( ie mod(arg(-Z)) < pi ) then
           !!             IHV= -1 and IP = 0 ( ie the standard formula)
           !!
           !!      (IIa) if Z.GE.1 ( ie mod(arg(1-Z)) = pi ) then
           !!             IHW=1 and IPW=-1 if 1-Z approaches the -ve real axis from below
           !!          OR IHW=-1and IPW=+1 "   "      "       "   "    "   "    "   above
           !!
           !!      (IIb) if Z.LT.1 ( ie mod(arg(1-Z)) = pi ) then
           !!             IHW = -1 and  IPW = 0 ( ie the standard formula )
           !!
           !!
           !!            REFERENCE
           !!            ---------
           !!            Abramowitz M and Stegun IA    Handbook of Mathematical
           !!            Functions (Dover NY 1970).
           !!
           !!--------------------------------------------------------------------------------

           implicit none

           real(kind=mp), intent(in) :: zeta, nu, z
           real(kind=mp) :: zinv, w

           complex(kind=c_mp) :: A, B, C, U
           complex(kind=c_mp), parameter :: i = cmplx(0.0_mp,1.0_mp,c_mp)  ! The complex number i.

           integer, intent(in) :: Ip, Ipw, Ihv, Ihw

           A = cmplx(1.0_mp,zeta)
           B = cmplx(1.0_mp,nu)
           C = cmplx(2.0_mp,0.0_mp) 

           if (z .ge. 1.5_mp .or. z .lt. -2.0_mp) then

              !! Use analytical continuation.  !!

              zinv = 1.0_mp/z

              W1 = (g1/(conjg(g4)*g3))*(Ihv*z)**(-A)*exp(i*pi*A*Ip)          &
                   *F_n(A,1.0_mp - C + A,1.0_mp - B + A,zinv)                &
                   + (conjg(g1)/(conjg(g3)*g4))*(Ihv*z)**(-B)                &
                   *exp(i*Ip*pi*B)*F_n(B,1.0_mp - C + B,1.0_mp - A + B,zinv)

           else if (z .ge. -2.0_mp .and. z .lt. -0.666_mp) then
              
              !! Use analytic continuation eqn 26.  !!

              W1 = (1.0_mp - z)**(-A)*F_n(A,C - B,C,z/(z - 1.0_mp))  
              
           else if (z .ge. 0.5_mp .and. z .lt. 1.5_mp) then 

              W = z - 1.0_mp

              if (abs(W) .lt. tolerance) then 

                 !!  Use Abramowitz & Stegun (15.1.20).  !!

                 W1 = conjg(g2)/(g3*g4)

              else

                 !!  Use analytic continuation 27.  !!

                 U = C - A - B

                 W1 = (conjg(g2)/(g3*g4))*F_n(A,B,A + B - C + 1.0_mp,-W)    &
                      + (Ihw*W)**U*exp(-Ipw*pi*i*U)*(g2/conjg(g3*g4))       &
                      *F_n(C - A,C - B,C - A - B + 1.0_mp,-W)

              end if

           else

              !!  ie if z > -0.66 and z < 0.5.                        !!
              !!  Thus the function F(A,B,C,Z) can be used directly.  !! 

              W1 = F_n(A,B,C,z)

           end if

           return

       end function W1


        
      Subroutine Quad1(f,a,b,epsabs,epsrel,key,result,abserr,neval,ier, &
                        limit,lenw,last,iwork,work)
       !!----------------------------------------------------------------------!!
       !!                                                                      !!
       !!   begin prologue  Quad                                               !!
       !!                                                                      !!
       !!  date written   800101   (yymmdd) (original code)                    !!
       !!  revision date  830518   (yymmdd)                                    !!
       !!  Date modified to f90 971001 (yymmdd)                                !!
       !!                                                                      !!
       !!  category no.  h2a1a1                                                !!
       !!                                                                      !!
       !!  Keywords    -                                                       !!
       !!                   Automatic integrator, general-purpose,             !!
       !!                   integrand examinator, globally adaptive,           !!
       !!                   Gauss-Kronrod                                      !!
       !!                                                                      !!
       !!  Author -                                                            !!
       !!           Piessens,Robert,appl. math. & progr. div - k.u.leuven      !!
       !!           de doncker,elise,appl. math. & progr. div. - k.u.leuven    !!
       !!                                                                      !! 
       !!  Modified to f90 -                                                   !! 
       !!                   Brian Nesbitt                                      !!
       !!                   Dept. of Applied Maths & Theoretical Physics       !!
       !!                   The Queen's University of Belfast.                 !!
       !!                                                                      !!
       !!                                                                      !!
       !!                    Purpose                                           !!
       !!                   ---------                                          !!
       !!                                                                      !!
       !!    the routine calculates an approximation result to a given         !!
       !!    definite integral i = integral of f over (a,b),                   !!
       !!    hopefully satisfying following claim for accuracy                 !!
       !!    abs(i-result)le.max(epsabs,epsrel*abs(i)).                        !!
       !!                                                                      !!
       !!                                                                      !!
       !!                                                                      !!
       !!                    Description                                       !!
       !!                   -------------                                      !!
       !!                                                                      !!
       !!    Computation of a definite integral                                !!
       !!    standard fortran subroutine                                       !!
       !!    Fortran90 precision version                                       !!
       !!                                                                      !!
       !!     f      -    Real                                                 !!
       !!                 Function subprogam defining the integrand            !! 
       !!                 function f(x). the actual name for f needs to be     !! 
       !!                 declared e x t e r n a l in the driver program.      !!
       !!                                                                      !!
       !!     a      -    Real                                                 !! 
       !!                 Lower limit of integration.                          !!
       !!                                                                      !!
       !!     b      -    Real                                                 !! 
       !!                 Upper limit of integration.                          !!
       !!                                                                      !!
       !!     epsabs -    Real                                                 !!
       !!                 Absolute accuracy requested.                         !!
       !!                                                                      !!
       !!     epsrel -    Real                                                 !! 
       !!                 Relative accuracy requested                          !!
       !!                 if  epsabs.le.0                                      !!
       !!                 and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),         !!
       !!                 the routine will end with ier = 6.                   !!
       !!                                                                      !!
       !!     key    -    Integer                                              !!
       !!                 Key for choice of local integration rule             !!
       !!                 a gauss-kronrod pair is used with                    !!
       !!                   7 - 15 points if key.lt.2,                         !!
       !!                  10 - 21 points if key = 2,                          !!
       !!                  15 - 31 points if key = 3,                          !!
       !!                  20 - 41 points if key = 4,                          !!
       !!                  25 - 51 points if key = 5,                          !!
       !!                  30 - 61 points if key .gt. 5.                       !!
       !!                                                                      !!
       !!                    On return                                         !!
       !!                   -----------                                        !! 
       !!                                                                      !!
       !!     result -    Real                                                 !!
       !!                 approximation to the integral.                       !!
       !!                                                                      !!
       !!     abserr -    Real                                                 !!
       !!                 estimate of the modulus of the absolute error,       !!
       !!                 which should equal or exceed abs(i-result).          !!
       !!                                                                      !!
       !!     neval  -    Integer                                              !!
       !!                 number of integrand evaluations.                     !!
       !!                                                                      !!
       !!     ier    -    Integer                                              !!
       !!                 ier = 0 normal and reliable termination of the       !!
       !!                      routine. it is assumed that the requested       !!
       !!                      accuracy has been achieved.                     !!
       !!                                                                      !!
       !!                 ier.gt.0 abnormal termination of the routine         !!
       !!                      the estimates for result and error are          !!
       !!                      less reliable. it is assumed that the           !!
       !!                      requested accuracy has not been achieved.       !!
       !!                                                                      !!
       !!                                                                      !!
       !!                    Error Messages                                    !!
       !!                   ----------------                                   !!
       !!                                                                      !!
       !!              ier = 1 maximum number of subdivisions allowed          !!
       !!                      has been achieved. one can allow more           !!
       !!                      subdivisions by increasing the value of         !!
       !!                      limit (and taking the according dimension       !!
       !!                      adjustments into account). however, if          !!
       !!                      this yield no improvement it is advised         !!
       !!                      to analyze the integrand in order to            !!
       !!                      determine the integration difficulaties.        !!
       !!                      if the position of a local difficulty can       !!
       !!                      be determined (i.e.singularity,                 !!
       !!                      discontinuity within the interval) one          !!
       !!                      will probably gain from splitting up the        !!
       !!                      interval at this point and calling the          !!
       !!                      integrator on the subranges. if possible,       !!
       !!                      an appropriate special-purpose integrator       !!
       !!                      should be used which is designed for            !!
       !!                      handling the type of difficulty involved.       !!
       !!                                                                      !!
       !!                  = 2 the occurrence of roundoff error is             !!
       !!                      detected, which prevents the requested          !!
       !!                      tolerance from being achieved.                  !!
       !!                                                                      !!
       !!                  = 3 extremely bad integrand behaviour occurs        !!
       !!                      at some points of the integration               !!
       !!                      interval.                                       !!
       !!                                                                      !!
       !!                  = 6 the input is invalid, because                   !!
       !!                      (epsabs.le.0 and                                !!
       !!                       epsrel.lt.max(50*rel.mach.acc.,0.5d-28))       !!
       !!                      or limit.lt.1 or lenw.lt.limit*4.               !!
       !!                      result, abserr, neval, last are set             !!
       !!                      to zero.                                        !!
       !!                      except when lenw is invalid, iwork(1),          !!
       !!                      work(limit*2+1) and work(limit*3+1) are         !!
       !!                      set to zero, work(1) is set to a and            !!
       !!                      work(limit+1) to b.                             !!
       !!                                                                      !!
       !!                                                                      !!
       !!                    Dimensioning Parameters                           !!
       !!                   -------------------------                          !!
       !!                                                                      !!
       !!                                                                      !!
       !!     limit -   Integer                                                !!
       !!               dimensioning parameter for iwork                       !!
       !!               limit determines the maximum number of subintervals    !!
       !!               in the partition of the given integration interval     !!
       !!               (a,b), limit.ge.1.                                     !!
       !!               if limit.lt.1, the routine will end with ier = 6.      !!
       !!                                                                      !!
       !!     lenw  -   Integer                                                !!
       !!               dimensioning parameter for work                        !!
       !!               lenw must be at least limit*4.                         !!
       !!               if lenw.lt.limit*4, the routine will end with          !!
       !!               ier = 6.                                               !!
       !!                                                                      !!
       !!     last  -   Integer                                                !!
       !!               on return, last equals the number of subintervals      !!
       !!               produced in the subdiviosion process, which            !!
       !!               determines the number of significant elements          !!
       !!               actually in the work arrays.                           !!
       !!                                                                      !!
       !!                                                                      !!
       !!                    Work Arrays                                       !!
       !!                   -------------                                      !!
       !!                                                                      !!
       !!                                                                      !!
       !!     iwork -   Integer                                                !!
       !!               vector of dimension at least limit, the first k        !!
       !!               elements of which contain pointers to the error        !!
       !!               estimates over the subintervals, such that             !!
       !!               work(limit*3+iwork(1)),... , work(limit*3+iwork(k))    !!
       !!               form a decreasing sequence with k = last if            !!
       !!               last.le.(limit/2+2), and k = limit+1-last otherwise    !!
       !!                                                                      !!
       !!     work  -   Real                                                   !!
       !!               vector of dimension at least lenw                      !!
       !!               on return                                              !!
       !!               work(1), ..., work(last) contain the left end          !!
       !!               points of the subintervals in the partition of         !!
       !!               (a,b),                                                 !!
       !!               work(limit+1), ..., work(limit+last) contain the       !!
       !!               right end points,                                      !!
       !!               work(limit*2+1), ..., work(limit*2+last) contain       !!
       !!               the integral approximations over the subintervals,     !!
       !!               work(limit*3+1), ..., work(limit*3+last) contain       !!
       !!               the error estimates.                                   !!
       !!                                                                      !!
       !!             References  (none)                                       !!
       !!                                                                      !!
       !!          Routines called  Dqage,Xerror                               !!
       !!                                                                      !! 
       !! end prologue  dqag                                                   !!
       !!                                                                      !!
       !!----------------------------------------------------------------------!!

           implicit none

           real(kind=mp), intent(in) :: a, b, epsabs, epsrel
           real(kind=mp), intent(out) :: abserr, result, work(*)

           integer, intent(in) :: lenw, limit, key
           integer, intent(out) :: ier, iwork(*), last, neval
           integer :: lvl, l1, l2, l3

           real(kind=mp), external :: f  ! external function defining the integrand.
       
           !!  Check validity of lenw. !!

           ier = 6
           neval = 0
           last = neval 
           result = 0.0_mp 
           abserr = result 

           if (limit .lt. 1 .or. lenw .lt. 4*limit) then
 
               if (ier .eq. 6) then
                  lvl = 1
               end if

           else
         
              !!  Prepare for call to DQAGE  !!

              l1 = limit + 1
              l2 = limit + l1
              l3 = limit + l2

              call DQAGE1(f,a,b,epsabs,epsrel,key,limit,result,abserr,neval,ier, &
                         work(1),work(l1),work(l2),work(l3),iwork,last)

              !!  Call error handler if necessary.  !!

              lvl = 0

           end if 
 
           if (ier .eq. 6) then
              lvl = 1
           else if (ier .ne. 0) then 
              call Xerror(ier)  !  Quadrature has failed.
           end if

     
          return

       end Subroutine Quad1   



       Subroutine DQAGE1(f,a,b,epsabs,epsrel,key,limit,result,abserr,neval,ier, &
                        alist,blist,rlist,elist,iord,last)

       !!-------------------------------------------------------------------------!!
       !!                                                                         !!
       !!    Begin prologue  dqage                                                !!
       !!                                                                         !!
       !!  Date written   800101   (yymmdd) (original code)                       !!
       !!  Revision date  830518   (yymmdd)                                       !!
       !!  Date modifies to f90 971001 (yymmdd)                                   !!
       !!                                                                         !!
       !!  category no.  h2a1a1                                                   !!
       !!                                                                         !!
       !!  Keywords   -                                                           !!
       !!                     automatic integrator, general-purpose,              !!
       !!                     integrand examinator, globally adaptive,            !!
       !!                     Gauss-Kronrod.                                      !!
       !!                                                                         !!
       !!  Author -                                                               !!
       !!           Piessens,Robert,appl. math. & progr. div - k.u.leuven         !!
       !!           de Doncker,Elise,appl. math. & progr. div. - k.u.leuven       !!
       !!                                                                         !!
       !!  Modified to f90 -                                                      !!
       !!                   Brian Nesbitt                                         !!
       !!                   Dept. of Applied Maths & Theoretical Physics          !!
       !!                   The Queen's University of Belfast.                    !!
       !!                                                                         !!
       !!                                                                         !!
       !!                    Purpose                                              !! 
       !!                   ---------                                             !!
       !!                                                                         !!
       !!   The routine calculates an approximation result to a given             !!
       !!   definite integral i = integral of f over (a,b),                       !!
       !!   hopefully satisfying following claim for accuracy                     !!
       !!   abs(i-result)le.max(epsabs,epsrel*abs(i)).                            !!
       !!                                                                         !!
       !!                                                                         !!
       !!                                                                         !!
       !!                    Description                                          !!
       !!                   -------------                                         !!
       !!                                                                         !!        
       !!                                                                         !!
       !!    Computation of a definite integral                                   !!
       !!    standard fortran subroutine.                                         !!
       !!                                                                         !!
       !!                                                                         !!        
       !!                    Parameters                                           !!
       !!                   ------------                                          !!
       !!                                                                         !!
       !!                                                                         !! 
       !!                    On Entry                                             !!
       !!                   ----------                                            !!
       !!                                                                         !!
       !!     f      -    Real                                                    !!
       !!                 function subprogram defining the integrand              !!
       !!                 function f(x). the actual name for f needs to be        !!
       !!                 declared e x t e r n a l in the driver program.         !!
       !!                                                                         !!
       !!     a      -    Real                                                    !! 
       !!                 lower limit of integration.                             !!         
       !!                                                                         !!
       !!     b      -    Real                                                    !!
       !!                 upper limit of integration.                             !!
       !!                                                                         !!
       !!     epsabs -    Real                                                    !!
       !!                 absolute accuracy requested.                            !!
       !!                                                                         !!
       !!     epsrel -    Real                                                    !!
       !!                 relative accuracy requested.                            !! 
       !!                 if  epsabs.le.0                                         !!
       !!                 and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),            !!
       !!                 the routine will end with ier = 6.                      !!
       !!                                                                         !!
       !!     key    -    Integer                                                 !!
       !!                 key for choice of local integration rule                !! 
       !!                 a Gauss-Kronrod pair is used with                       !!
       !!                     7 - 15 points if key.lt.2,                          !!
       !!                    10 - 21 points if key = 2,                           !!
       !!                    15 - 31 points if key = 3,                           !!
       !!                    20 - 41 points if key = 4,                           !!
       !!                    25 - 51 points if key = 5,                           !!
       !!                    30 - 61 points if key.gt.5.                          !!
       !!                                                                         !!
       !!     limit  -    Integer                                                 !!
       !!                 gives an upperbound on the number of subintervals       !!
       !!                 in the partition of (a,b), limit.ge.1.                  !!
       !!                                                                         !!
       !!                    On Return                                            !!
       !!                   -----------                                           !!
       !!                                                                         !!
       !!     result -    Real                                                    !!
       !!                 approximation to the integral                           !!
       !!                                                                         !!
       !!     abserr -    Real                                                    !!  
       !!                 estimate of the modulus of the absolute error,          !!
       !!                 which should equal or exceed abs(i-result).             !!
       !!                                                                         !!
       !!     neval  -    Integer                                                 !!
       !!                 number of integrand evaluations.                        !!
       !!                                                                         !!
       !!     ier    -    Integer                                                 !!
       !!                 ier = 0     normal and reliable termination of the      !!
       !!                             routine. it is assumed that the requested   !!
       !!                             accuracy has been achieved.                 !!
       !!                                                                         !!
       !!                 ier.gt.0    abnormal termination of the routine         !!
       !!                             the estimates for result and error are      !!
       !!                             less reliable. it is assumed that the       !!
       !!                             requested accuracy has not been achieved.   !!
       !!                                                                         !!
       !!                    Error Messages                                       !!
       !!                   ----------------                                      !!
       !!                                                                         !!
       !!              ier = 1 maximum number of subdivisions allowed             !!
       !!                      has been achieved. one can allow more              !!
       !!                      subdivisions by increasing the value               !!
       !!                      of limit.                                          !!
       !!                      however, if this yields no improvement it          !!
       !!                      is rather advised to analyze the integrand         !!
       !!                      in order to determine the integration              !!
       !!                      difficulties. if the position of a local           !!
       !!                      difficulty can be determined(e.g.                  !!
       !!                      singularity, discontinuity within the              !!
       !!                      interval) one will probably gain from              !!
       !!                      splitting up the interval at this point            !!
       !!                      and calling the integrator on the                  !!
       !!                      subranges. if possible, an appropriate             !!
       !!                      special-purpose integrator should be used          !!
       !!                      which is designed for handling the type of         !!
       !!                      difficulty involved.                               !!
       !!                                                                         !!
       !!                 =  2 the occurrence of roundoff error is                !!
       !!                      detected, which prevents the requested             !!
       !!                      tolerance from being achieved.                     !!
       !!                                                                         !!
       !!                  = 3 extremely bad integrand behaviour occurs           !!
       !!                      at some points of the integration                  !!
       !!                      interval.                                          !!
       !!                                                                         !!
       !!                  = 6 the input is invalid, because                      !!
       !!                      (epsabs.le.0 and                                   !!
       !!                      epsrel.lt.max(50*rel.mach.acc.,0.5d-28),           !!
       !!                      result, abserr, neval, last, rlist(1) ,            !!
       !!                      elist(1) and iord(1) are set to zero.              !!
       !!                      alist(1) and blist(1) are set to a and b           !!
       !!                      respectively.                                      !!
       !!                                                                         !!
       !!     alist   -   Real                                                    !!
       !!                 vector of dimension at least limit, the first           !!
       !!                 last  elements of which are the left                    !!
       !!                 end points of the subintervals in the partition         !!
       !!                 of the given integration range (a,b)                    !!
       !!                                                                         !!
       !!     blist   -   Real                                                    !!
       !!                 vector of dimension at least limit, the first           !!
       !!                 last  elements of which are the right                   !!
       !!                 end points of the subintervals in the partition         !!
       !!                 of the given integration range (a,b)                    !!
       !!                                                                         !!
       !!     rlist   -   Real                                                    !!
       !!                 vector of dimension at least limit, the first           !!         
       !!                 last  elements of which are the                         !!
       !!                 integral approximations on the subintervals             !! 
       !!                                                                         !!
       !!     elist   -   Real                                                    !! 
       !!                 vector of dimension at least limit, the first           !!
       !!                 last  elements of which are the moduli of the           !!
       !!                 absolute error estimates on the subintervals            !!
       !!                                                                         !!
       !!     iord    -   Integer                                                 !!
       !!                 vector of dimension at least limit, the first k         !!
       !!                 elements of which are pointers to the                   !!
       !!                 error estimates over the subintervals,                  !!
       !!                 such that elist(iord(1)), ...,                          !!
       !!                 elist(iord(k)) form a decreasing sequence,              !!
       !!                 with k = last if last.le.(limit/2+2), and               !!
       !!                 k = limit+1-last otherwise.                             !!
       !!                                                                         !!
       !!     last    -   Integer                                                 !!         
       !!                 number of subintervals actually produced in the         !!
       !!                 subdivision process.                                    !!
       !!                                                                         !!
       !!               references  (none)                                        !!
       !!            routines called  d1mach,dqk15,dqk21,dqk31,                   !!
       !!                             dqk41,dqk51,dqk61,dqpsrt                    !!
       !!                                                                         !!
       !! end prologue  dqage                                                     !!
       !!                                                                         !!
       !!-------------------------------------------------------------------------!!

           implicit none

           real(kind=mp), intent(in) :: a, b, epsrel, epsabs
           real(kind=mp), intent(out) :: abserr, result
           real(kind=mp), intent(out) :: alist(*), blist(*), elist(*), rlist(*)

           real(kind=mp) :: area, area1, area12, area2, a1, a2,        &
                            b1, b2, defabs, defab1, defab2,            &
                            errbnd, errmax, error1, erro12, errsum,    &
                            resabs, error2

           integer, intent(in) :: limit, key 
           integer, intent(out) :: ier, iord(*), last, neval

           integer :: iroff1, iroff2, k, keyf, maxerr, nrmax

           real(kind=mp), external :: f ! external function defining the integrand.

           !!---------------------------------------------------------------------!!
           !!                                                                     !!  
           !!    list of major variables                                          !!  
           !!    -----------------------                                          !! 
           !!                                                                     !!  
           !!          alist     -   Real array                                   !! 
           !!                        list of left end points of all subintervals  !!  
           !!                        considered up to now.                        !!
           !!                                                                     !!
           !!          blist     -   Real array                                   !! 
           !!                        list of right end points of all subinterval  !!
           !!                        considered up to now.                        !!
           !!                                                                     !!
           !!          rlist(i)  -   Real                                         !!
           !!                        approximation to the integral over           !!
           !!                        (alist(i),blist(i)).                         !!
           !!                                                                     !!
           !!          elist(i)  -   Real                                         !!
           !!                        error estimate applying to rlist(i).         !!
           !!                                                                     !!
           !!          maxerr    -   Integer                                      !!
           !!                        pointer to the interval with largest.        !!
           !!                                                                     !!
           !!                                                                     !!
           !!            Error Estimate                                           !!
           !!           ----------------                                          !!
           !!                                                                     !!
           !!          errmax    -   Real                                         !! 
           !!                        elist(maxerr).                               !!
           !!                                                                     !!
           !!          area      -   Real                                         !!
           !!                        sum of the integrals over the subintervals.  !!
           !!                                                                     !!
           !!          errsum    -   Real                                         !!
           !!                        sum of the errors over the subintervals.     !!
           !!                                                                     !!
           !!          errbnd    -   Real                                         !!
           !!                        requested accuracy max(epsabs,epsrel*        !!
           !!                        abs(result)).                                !!
           !!                                                                     !!
           !!          *****1    -   variable for the left subinterval.           !!
           !!          *****2    -   variable for the right subinterval.          !!
           !!                                                                     !!
           !!          last      -   Integer                                      !!
           !!                        index for subdivision.                       !!
           !!                                                                     !!
           !!                                                                     !!
           !!              Machine Dependent Constants                            !!
           !!             -----------------------------                           !!
           !!                                                                     !!
           !!      epmach  -   is the largest relative spacing.                   !!
           !!      uflow   -   is the smallest positive magnitude.                !!
           !!                                                                     !!
           !!---------------------------------------------------------------------!!
           
           !!  First executable statement of DQAGE. !!


           !!  Test on validity of parameters.      !!

           ier = 0
           neval = ier 
           last = ier 
           result = 0.0_mp 
           abserr = result 
           alist(1) = a
           blist(1) = b
           rlist(1) = result 
           elist(1) = result 
           iord(1) = ier 

           if (epsabs .le. 0.0_mp .and. epsrel .lt.     &
               max(0.5e+2_mp*epmach,0.5e-28_mp)) then
               ier = 6
           end if

           if (ier .eq. 6) then
              call Xerror(ier)  ! Quadrature has failed.
           end if

           keyf = key

           if (key .le. 0) then 
               keyf = 1
           else if (key .ge. 7) then 
               keyf = 6
           end if

           neval = 0

      !     call DQK61_1(f,a,b,result,abserr,defabs,resabs)

           select case(keyf)
                  case(1)
                    call DQK15_1(f,a,b,result,abserr,defabs,resabs)
                  case(6)
                    call DQK61_1(f,a,b,result,abserr,defabs,resabs)
           end select 

           last = 1
           rlist(1) = result
           elist(1) = abserr
           iord(1) = 1

           !!  Test on accuracy.  !!

           errbnd = max(epsabs,epsrel*abs(result))

           if (abserr .le. 0.5e+2_mp*epmach*defabs .and. abserr .gt. errbnd) then 
               ier = 2 
           end if

           if (limit .eq. 1) then 
               ier = 1
           end if 

           if (ier .ne. 0 .or. (abserr .le. errbnd .and. abserr    &
              .ne. resabs) .or. abserr .eq. 0.0_mp) then
              
              if (keyf .ne. 1) then
                  neval = (10*keyf + 1)*(2*neval + 1)
              else if (keyf .eq. 1) then
                  neval = 30*neval + 15
              end if

              return
           end if

           !!  Initialization  !!

           errmax = abserr
           maxerr = 1
           area = result
           errsum = abserr
           nrmax = 1
           iroff1 = 0
           iroff2 = 0

           !!  Main do loop  !!
    
           do last = 2,limit

               !!  Bisect the subinterval with the largest error estimate.  !!

               a1 = alist(maxerr)
               b1 = 0.5_mp*(alist(maxerr) + blist(maxerr))
               a2 = b1
               b2 = blist(maxerr)

               !!  Use 61-point rule only at this stage.  !!

               select case(keyf)
                      case(1)
                        call DQK15_1(f,a1,b1,area1,error1,resabs,defab1)
                        call DQK15_1(f,a2,b2,area2,error2,resabs,defab2)
                      case(6)
                        call DQK61_1(f,a1,b1,area1,error1,resabs,defab1)
                        call DQK61_1(f,a2,b2,area2,error2,resabs,defab2)
               end select 

          !     call DQK61_1(f,a1,b1,area1,error1,resabs,defab1)
          !     call DQK61_1(f,a2,b2,area2,error2,resabs,defab2)

               !!  Improve the previous approximations to the integral  !!
               !!  and error and test for accuracy.                     !!

               neval = neval + 1
               area12 = area1 + area2
               erro12 = error1 + error2
               errsum = errsum + erro12 - errmax
               area = area + area12 - rlist(maxerr)

               if (defab1 .ne. error1 .and. defab2 .ne. error2) then

                   if (abs(rlist(maxerr) - area12) .le. 0.1e-4_mp*abs(area12) &
                      .and. erro12 .ge. 0.99_mp*errmax) then
                      iroff1 = iroff1 + 1
                   end if

                   if (last .gt. 10 .and. erro12 .gt. errmax) then
                      iroff2 = iroff2 + 1
                   end if

               end if

               rlist(maxerr) = area1
               rlist(last)  = area2

               errbnd = max(epsabs,epsrel*abs(area))

               if (errsum .gt. errbnd) then 

                  !!  Test for roundoff error and eventually set  !!
                  !!  error flag.                                 !!
                  
                  if (iroff1 .ge. 6 .or. iroff2 .ge. 20) then 
                      ier = 2
                  end if

                  !!  Set error flag in the case that the number  !!
                  !!  of subintervals equals limit.               !!

                  if (last .eq. limit) then 
                     ier = 1
                  end if

                  !!  Set error flag in the case of bad integrand     !!   
                  !!  behaviour at a point of the integration range.  !!

                  if (max(abs(a1),abs(b2)) .le.  & 
                     (0.1e+1_mp + 0.1e+3_mp*epmach)*(abs(a2) + 0.1e+4_mp*uflow)) then 
                      ier = 3
                  end if
        
               end if 

               !!  Apend the newly-created intervals to the list.  !!

               if (error2 .gt. error1) then 
                   alist(maxerr) = a2
                   alist(last) = a1
                   blist(last) = b1
                   rlist(maxerr) = area2
                   rlist(last) = area1
                   elist(maxerr) = error2
                   elist(last) = error1
               else
                   alist(last) = a2
                   blist(maxerr) = b1
                   blist(last) = b2
                   elist(maxerr) = error1
                   elist(last) = error2
               end if

               !!  Call subroutine DQPSRT to maintain the descending ordering  !!
               !!  in the list of error estimates and select the subinterval   !!
               !!  with the largest error estimate (to be bisected next).      !!

               call DQPSRT_1(limit,last,maxerr,errmax,elist,iord,nrmax)

               !!  Exit do loop.  !!
 
               if (ier .ne. 0 .or. errsum .le. errbnd) then 
                  exit
               end if

           end do

           !!  Compare the final result. !!

           result = 0.0_mp 

           do k = 1,last
              result = result + rlist(k)
           end do
 
           abserr = errsum

           if (keyf .ne. 1) then
               neval = (10*keyf + 1)*(2*neval + 1)
           end if

           if (keyf .eq. 1) then 
               neval = 30*neval + 15
           end if

           return

       end Subroutine DQAGE1



       Subroutine DQK15_1(f,a,b,result,abserr,resabs,resasc)

           implicit none

           real(kind=mp), intent(in) :: a, b
           real(kind=mp), intent(out) :: abserr, resabs, resasc, result            

           real(kind=mp) :: absc, centr, dhlgth, fc, fsum, fval1,     &
                            fval2, hlgth, resg, resk, reskh

           real(kind=mp) :: fv1(7), fv2(7), wg(4), wgk(8), xgk(8)

           integer :: j, jtw, jtwm1

           real(kind=mp), external :: f

           !!-------------------------------------------------------------------------!!
           !!                                                                         !!
           !!  The abscissae and weights are given for the interval (-1,1).           !!
           !!  Because of symmetry only the positive abscissae and their              !!
           !!  corresponding weights are given.                                       !!
           !!                                                                         !!
           !!  xgk    - abscissae of the 15-point Kronrod rule                        !!
           !!         xgk(2), xgk(4), ...  abscissae of the 7-point                   !!
           !!         Gauss rule                                                      !!
           !!         xgk(1), xgk(3), ...  abscissae which are optimally              !!
           !!         added to the 7-point Gauss rule                                 !!
           !!                                                                         !!
           !!  wgk    - weights of the 15-point Kronrod rule                          !!
           !!                                                                         !!
           !!  wg     - weights of the 7-point Gauss rule                             !!
           !!                                                                         !!
           !!                                                                         !!
           !!  Gauss quadrature weights and Kronrod quadrature abscissae and weights  !!
           !!  as evaluated with 80 decimal digit arithmetic by L. W. Fullerton,      !!
           !!  Bell labs, Nov. 1981.                                                  !!
           !!                                                                         !!
           !!-------------------------------------------------------------------------!!

           data wg  (1) / 0.129484966168869693270611432679082_mp  /
           data wg  (2) / 0.279705391489276667901467771423780_mp  /
           data wg  (3) / 0.381830050505118944950369775488975_mp  /
           data wg  (4) / 0.417959183673469387755102040816327_mp  /

           data xgk (1) / 0.991455371120812639206854697526329_mp  /
           data xgk (2) / 0.949107912342758524526189684047851_mp  /
           data xgk (3) / 0.864864423359769072789712788640926_mp  /
           data xgk (4) / 0.741531185599394439863864773280788_mp  /
           data xgk (5) / 0.586087235467691130294144838258730_mp  /
           data xgk (6) / 0.405845151377397166906606412076961_mp  /
           data xgk (7) / 0.207784955007898467600689403773245_mp  /
           data xgk (8) / 0.000000000000000000000000000000000_mp  /

           data wgk (1) / 0.022935322010529224963732008058970_mp  /
           data wgk (2) / 0.063092092629978553290700663189204_mp  /
           data wgk (3) / 0.104790010322250183839876322541518_mp  /
           data wgk (4) / 0.140653259715525918745189590510238_mp  /
           data wgk (5) / 0.169004726639267902826583426598550_mp  /
           data wgk (6) / 0.190350578064785409913256402421014_mp  /
           data wgk (7) / 0.204432940075298892414161999234649_mp  /
           data wgk (8) / 0.209482141084727828012999174891714_mp  /

           !!  First executable statements of DQK15.  !!

           centr = 0.5_mp*(a + b)
           hlgth = 0.5_mp*(b - a)
           dhlgth = abs(hlgth)

           !!  Compute the 15-point Kronrod approximation to the  !!
           !!  integral, and estimate the absolute error.         !!

           fc = f(centr)
           resg = fc*wg(4)
           resk = fc*wgk(8)
           resabs = abs(resk)

           do j = 1,3
              jtw = j*2
              absc = hlgth*xgk(jtw)
              fval1 = f(centr - absc)
              fval2 = f(centr + absc)
              fv1(jtw) = fval1
              fv2(jtw) = fval2
              fsum = fval1 + fval2
              resg = resg + wg(j)*fsum
              resk  = resk + wgk(jtw)*fsum
              resabs = resabs + wgk(jtw)*(abs(fval1) + abs(fval2))
           end do

           do j = 1,4
              jtwm1 = j*2 - 1
              absc = hlgth*xgk(jtwm1)
              fval1 = f(centr - absc)
              fval2 = f(centr + absc)
              fv1(jtwm1) = fval1
              fv2(jtwm1) = fval2
              fsum = fval1 + fval2
              resk = resk + wgk(jtwm1)*fsum
              resabs = resabs + wgk(jtwm1)*(abs(fval1) + abs(fval2))
           end do

           reskh = resk*0.5_mp
           resasc = wgk(8)*abs(fc - reskh)

           do j = 1,7
              resasc = resasc + wgk(j)*(abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
           end do

           result = resk*hlgth
           resabs = resabs*dhlgth
           resasc = resasc*dhlgth
           abserr = abs((resk - resg)*hlgth)

           if (resasc .ne. 0.0_mp .and. abserr .ne. 0.0_mp) then
               abserr = resasc*min(0.1e+1_mp,(0.2e+3_mp*abserr/resasc)**1.5_mp)
           end if

           if (resabs .gt. uflow/(0.5e+2_mp*epmach)) then
               abserr = max((epmach*0.5e+2_mp)*resabs,abserr)
           end if

           return

       end Subroutine DQK15_1 



       Subroutine DQK61_1(f,a,b,result,abserr,resabs,resasc)
            
           implicit none

           real(kind=mp), intent(in) :: a, b
           real(kind=mp), intent(out) :: abserr, resabs, resasc, result

           real(kind=mp) :: dabsc, centr, dhlgth,            &        
                            fc, fsum, fval1, fval2, hlgth,   &
                            resg, resk, reskh

           real(kind=mp) :: fv1(30), fv2(30), wg(15), wgk(31), xgk(31)

           real(kind=mp), external :: f

           integer :: j, jtw, jtwm1
       
           !!---------------------------------------------------------!!
           !!                                                         !!
           !!  The abscissae and weights are given for the interval   !!
           !!  (-1,+1). Because of symmetry only the positive         !!  
           !!  abscissae and their corresponding weights are given.   !! 
           !!                                                         !!
           !!                                                         !!
           !!   xgk  -  abscissae of the 51-point Kronrod rule:       !!
           !!             xgk(2), xgk(4), ...                         !!
           !!                                                         !! 
           !!           abscissae of the 25-point Gauss rule:         !! 
           !!             xgk(1), xgk(3), ...                         !!
           !!                                                         !!
           !!           abscissae which are optimally added to        !!
           !!           the 25-point Gauss rule.                      !!
           !!                                                         !!
           !!   wgk  -  weights of the 51-point Kronrod rule.         !!
           !!                                                         !!
           !!   wg   -  weights of the 25-point Gauss rule.           !!
           !!                                                         !!
           !!   Gauss quadrature weights and Kronrod Quadrature       !!
           !!   abscissae and weights as evaluated with 80 decimal    !!
           !!   arithmetic by L.W. Fullerton, Bell Labs., Nov 1981.   !!
           !!                                                         !!
           !!---------------------------------------------------------!!
  
           data wg  (  1) / 0.007968192496166605615465883474674_mp  /
           data wg  (  2) / 0.018466468311090959142302131912047_mp  /
           data wg  (  3) / 0.028784707883323369349719179611292_mp  /
           data wg  (  4) / 0.038799192569627049596801936446348_mp  /
           data wg  (  5) / 0.048402672830594052902938140422808_mp  /
           data wg  (  6) / 0.057493156217619066481721689402056_mp  /
           data wg  (  7) / 0.065974229882180495128128515115962_mp  /
           data wg  (  8) / 0.073755974737705206268243850022191_mp  /
           data wg  (  9) / 0.080755895229420215354694938460530_mp  /
           data wg  ( 10) / 0.086899787201082979802387530715126_mp  /
           data wg  ( 11) / 0.092122522237786128717632707087619_mp  /
           data wg  ( 12) / 0.096368737174644259639468626351810_mp  /
           data wg  ( 13) / 0.099593420586795267062780282103569_mp  /
           data wg  ( 14) / 0.101762389748405504596428952168554_mp  /
           data wg  ( 15) / 0.102852652893558840341285636705415_mp  /

           data xgk (  1) / 0.999484410050490637571325895705811_mp  /
           data xgk (  2) / 0.996893484074649540271630050918695_mp  /
           data xgk (  3) / 0.991630996870404594858628366109486_mp  /
           data xgk (  4) / 0.983668123279747209970032581605663_mp  /
           data xgk (  5) / 0.973116322501126268374693868423707_mp  /
           data xgk (  6) / 0.960021864968307512216871025581798_mp  /
           data xgk (  7) / 0.944374444748559979415831324037439_mp  /
           data xgk (  8) / 0.926200047429274325879324277080474_mp  /
           data xgk (  9) / 0.905573307699907798546522558925958_mp  /
           data xgk ( 10) / 0.882560535792052681543116462530226_mp  /
           data xgk ( 11) / 0.857205233546061098958658510658944_mp  /
           data xgk ( 12) / 0.829565762382768397442898119732502_mp  /
           data xgk ( 13) / 0.799727835821839083013668942322683_mp  /
           data xgk ( 14) / 0.767777432104826194917977340974503_mp  /
           data xgk ( 15) / 0.733790062453226804726171131369528_mp  /
           data xgk ( 16) / 0.697850494793315796932292388026640_mp  /
           data xgk ( 17) / 0.660061064126626961370053668149271_mp  /
           data xgk ( 18) / 0.620526182989242861140477556431189_mp  /
           data xgk ( 19) / 0.579345235826361691756024932172540_mp  /
           data xgk ( 20) / 0.536624148142019899264169793311073_mp  /
           data xgk ( 21) / 0.492480467861778574993693061207709_mp  /
           data xgk ( 22) / 0.447033769538089176780609900322854_mp  /
           data xgk ( 23) / 0.400401254830394392535476211542661_mp  /
           data xgk ( 24) / 0.352704725530878113471037207089374_mp  /
           data xgk ( 25) / 0.304073202273625077372677107199257_mp  /
           data xgk ( 26) / 0.254636926167889846439805129817805_mp  /
           data xgk ( 27) / 0.204525116682309891438957671002025_mp  /
           data xgk ( 28) / 0.153869913608583546963794672743256_mp  /
           data xgk ( 29) / 0.102806937966737030147096751318001_mp  /
           data xgk ( 30) / 0.051471842555317695833025213166723_mp  /
           data xgk ( 31) / 0.000000000000000000000000000000000_mp  /

           data wgk (  1) / 0.001389013698677007624551591226760_mp  /
           data wgk (  2) / 0.003890461127099884051267201844516_mp  /
           data wgk (  3) / 0.006630703915931292173319826369750_mp  /
           data wgk (  4) / 0.009273279659517763428441146892024_mp  /
           data wgk (  5) / 0.011823015253496341742232898853251_mp  /
           data wgk (  6) / 0.014369729507045804812451432443580_mp  /
           data wgk (  7) / 0.016920889189053272627572289420322_mp  /
           data wgk (  8) / 0.019414141193942381173408951050128_mp  /
           data wgk (  9) / 0.021828035821609192297167485738339_mp  /
           data wgk ( 10) / 0.024191162078080601365686370725232_mp  /
           data wgk ( 11) / 0.026509954882333101610601709335075_mp  /
           data wgk ( 12) / 0.028754048765041292843978785354334_mp  /
           data wgk ( 13) / 0.030907257562387762472884252943092_mp  /
           data wgk ( 14) / 0.032981447057483726031814191016854_mp  /
           data wgk ( 15) / 0.034979338028060024137499670731468_mp  /
           data wgk ( 16) / 0.036882364651821229223911065617136_mp  /
           data wgk ( 17) / 0.038678945624727592950348651532281_mp  /
           data wgk ( 18) / 0.040374538951535959111995279752468_mp  /
           data wgk ( 19) / 0.041969810215164246147147541285970_mp  /
           data wgk ( 20) / 0.043452539701356069316831728117073_mp  /
           data wgk ( 21) / 0.044814800133162663192355551616723_mp  /
           data wgk ( 22) / 0.046059238271006988116271735559374_mp  /
           data wgk ( 23) / 0.047185546569299153945261478181099_mp  /
           data wgk ( 24) / 0.048185861757087129140779492298305_mp  /
           data wgk ( 25) / 0.049055434555029778887528165367238_mp  /
           data wgk ( 26) / 0.049795683427074206357811569379942_mp  /
           data wgk ( 27) / 0.050405921402782346840893085653585_mp  /
           data wgk ( 28) / 0.050881795898749606492297473049805_mp  /
           data wgk ( 29) / 0.051221547849258772170656282604944_mp  /
           data wgk ( 30) / 0.051426128537459025933862879215781_mp  /
           data wgk ( 31) / 0.051494729429451567558340433647099_mp  /
           
           !!-------------------------------------------------------!!
           !!                                                       !!
           !!     List of major variables                           !!          
           !!     -----------------------                           !!
           !!                                                       !!
           !!     centr  -  mid point of the interval               !! 
           !!     hlgth  -  half-length of the interval             !!
           !!     dabsc  -  abscissa                                !!
           !!     fval*  -  function value                          !!
           !!     resg   -  result of the 30-point gauss rule       !!
           !!     resk   -  result of the 61-point kronrod rule     !!
           !!     reskh  -  approximation to the mean value of f    !! 
           !!               over (a,b), i.e. to i/(b-a)             !! 
           !!                                                       !! 
           !!     Machine Dependent Constants                       !!
           !!     ---------------------------                       !!
           !!                                                       !!
           !!     epmach is the largest relative spacing.           !!
           !!     uflow is the smallest positive magnitude.         !! 
           !!                                                       !!
           !!-------------------------------------------------------!!
      
           centr = 0.5_mp*(b + a)
           hlgth = 0.5_mp*(b - a)
           dhlgth = abs(hlgth)

           !!  Compute the 61-point Kronrod approximation to the  !!
           !!  integral, and estimate the absolute error.         !!

           !!  First executable statement of DQK61.  !!

           resg = 0.0_mp 
           fc = f(centr)
           resk = wgk(31)*fc
           resabs = abs(resk)

           do j = 1,15
               jtw = j*2
               dabsc = hlgth*xgk(jtw)
               fval1 = f(centr - dabsc)
               fval2 = f(centr + dabsc)
               fv1(jtw) = fval1
               fv2(jtw) = fval2
               fsum = fval1 + fval2
               resg = resg + wg(j)*fsum
               resk = resk + wgk(jtw)*fsum
               resabs = resabs + wgk(jtw)*(abs(fval1) + abs(fval2))
           end do    

           do j = 1,15
               jtwm1 = j*2 - 1
               dabsc = hlgth*xgk(jtwm1)
               fval1 = f(centr - dabsc)
               fval2 = f(centr + dabsc)
               fv1(jtwm1) = fval1
               fv2(jtwm1) = fval2
               fsum = fval1 + fval2
               resk = resk + wgk(jtwm1)*fsum                    
               resabs = resabs + wgk(jtwm1)*(abs(fval1) + abs(fval2))
           end do 

           reskh = resk*0.5_mp
           resasc = wgk(31)*abs(fc - reskh)

           do j = 1,30,2
               resasc = resasc + wgk(j)*(abs(fv1(j) - reskh) + &
                        abs(fv2(j) - reskh))
               resasc = resasc + wgk(j+1)*(abs(fv1(j+1) - reskh) + &
                        abs(fv2(j+1) - reskh))
           end do

           result = resk*hlgth
           resabs = resabs*dhlgth
           resasc = resasc*dhlgth
           abserr = abs((resk - resg)*hlgth)

           if (resasc .ne. 0.0_mp .and. abserr .ne. 0.0_mp) then 
               abserr = resasc*min(0.1e+1_mp,(0.2e+3_mp*abserr/resasc)**1.5_mp)
           end if

           if (resabs .gt. uflow/(0.5e+2_mp*epmach)) then
              abserr = max((epmach*0.5e+2_mp)*resabs,abserr)
           end if

           return

       end Subroutine DQK61_1



       Subroutine DQPSRT_1(limit,last,maxerr,ermax,elist,iord,nrmax)
   
           !!----------------------------------------------------------------------!!
           !!                                                                      !!
           !!   Begin prologue DQPSRT                                              !!       
           !!                                                                      !!
           !!   Refer to  dqage,dqagie,dqagpe,dqawse                               !!
           !!   Routines called  (none)                                            !!
           !!   Revision date  810101   (yymmdd)                                   !!
           !!                                                                      !!
           !!   keywords  -   sequential sorting.                                  !!
           !!                                                                      !!
           !!   Author  Piessens,Robert,appl. math. & progr. div. - k.u.leuven     !!
           !!   de doncker,elise,appl. math. & progr. div. - k.u.leuven.           !!
           !!                                                                      !!
           !!                                                                      !!
           !!      Purpose                                                         !!
           !!     ---------                                                        !!
           !!                                                                      !!
           !!      This routine maintains the descending ordering in the           !!
           !!      list of the local error estimated resulting from the            !!
           !!      interval subdivision process. at each call two error            !!
           !!      estimates are inserted using the sequential search              !!
           !!      method, top-down for the largest error estimate and             !!
           !!      bottom-up for the smallest error estimate.                      !!
           !!                                                                      !!
           !!                                                                      !!
           !!      Description                                                     !!
           !!     -------------                                                    !!
           !!                                                                      !!
           !!      ordering routine                                                !!
           !!      standard fortran subroutine                                     !!
           !!                                                                      !!
           !!      Parameters:                                                     !!
           !!     -------------                                                    !!
           !!                                                                      !! 
           !!                  Meaning At Output                                   !!
           !!                 -------------------                                  !!
           !!                                                                      !!
           !!         limit  -   Integer                                           !! 
           !!                    Maximum number of error estimates the list        !!
           !!                    can contain                                       !!
           !!                                                                      !! 
           !!         last   -   Integer                                           !!
           !!                    Number of error estimates currently in the list.  !!
           !!                                                                      !!
           !!         maxerr -   Integer                                           !!
           !!                    maxerr points to the nrmax-th largest error       !! 
           !!                    estimate currently in the list                    !!
           !!                                                                      !!
           !!         ermax  -   Real                                              !!
           !!                    nrmax-th largest error estimate.                  !!
           !!                    ermax = elist(maxerr)                             !!
           !!                                                                      !!
           !!         elist  -   Real                                              !!
           !!                    Vector of dimension last containing               !!
           !!                    the error estimates.                              !!
           !!                                                                      !!
           !!         iord   -   Integer                                           !!
           !!                    Vector of dimension last, the first k elements    !!
           !!                    of which contain pointers to the error            !!
           !!                    estimates, such that                              !!
           !!                    elist(iord(1)),...,  elist(iord(k))               !!
           !!                    form a decreasing sequence, with                  !!
           !!                    k = last if last.le.(limit/2+2), and              !!
           !!                    k = limit+1-last otherwise                        !!
           !!                                                                      !!
           !!         nrmax  -   Integer                                           !!
           !!                                                                      !!
           !!         maxerr =   iord(nrmax)                                       !!
           !!                                                                      !!
           !! end prologue  dqpsrt                                                 !!
           !!                                                                      !!
           !!----------------------------------------------------------------------!!
           

           implicit none

           real(kind=mp), intent(inout) :: ermax, elist(*)      
      
           real(kind=mp) :: errmax, errmin

           integer, intent(inout) :: iord(*), maxerr, nrmax
           integer, intent(in) :: limit

           integer :: i, ibeg, last, isucc, j, jbnd, jupbn, k, ido 

           !!  Check whether the list contains more than  !!
           !!  two error estimates.                       !!
            
           if (last .le. 2) then
               iord(1) = 1
               iord(2) = 2
               maxerr = iord(nrmax)
               ermax = elist(maxerr)
               return
           end if

           !!  This part of the routine is only executed if, due to a    !!
           !!  difficult integrand, subdivision increased the error      !!
           !!  estimate. In the normal case the insert procedure should  !!
           !!  start after the nrmax-th largest error estimate.          !!

           errmax = elist(maxerr)

           if (nrmax .ne. 1) then
               ido = nrmax - 1

               do i = 1, ido
                  isucc = iord(nrmax - 1)
                  
                  !!  Exit from do  loop.  !!
                  
                  if (errmax .le. elist(isucc)) then
                      exit
                  end if
            
                  iord(nrmax) = isucc
                  nrmax = nrmax - 1
               end do
           end if

           !!  Compute the number of elements in the list to be maintained  !!
           !!  in descending order. This number depends on the number of    !!
           !!  subdivisions still allowed.                                  !!

           jupbn = last

           if (last .gt. (limit/2 + 2)) then 
               jupbn = limit + 3 - last
           end if

           errmin = elist(last)

           !!  Insert errmax by traversing the list top-down, starting  !!
           !!  comparison from the element list elist(iord(nrmax + 1)). !!

           jbnd = jupbn - 1
           ibeg = nrmax + 1

           if (ibeg .le. jbnd) then 

               do i = ibeg, jbnd
                  isucc = iord(i)

                  !!  Exit from do loop.  !!

                  if (errmax .ge. elist(isucc)) then 
                      exit
                  end if

                  iord(i-1) = isucc
               end do

               if (errmax .ge. elist(isucc)) then 

                   !! Insert errmin by traversing the list bottom up.  !!

                   iord(i-1) = maxerr
                   k = jbnd

                   do j = i, jbnd
                      isucc = iord(k)

                      !!  Exit from do loop.  !!

                      if (errmin .lt. elist(isucc)) then 
                          exit
                      end if
           
                      iord(k+1) = isucc
                      k = k - 1
                   end do

                   if (errmin .lt. elist(isucc)) then
                       iord(k+1) = last
                   else
                       iord(i) = last
                   end if
          
               else
                   iord(jbnd) = maxerr
                   iord(jupbn) = last
               end if

           else
               iord(jbnd) = maxerr
               iord(jupbn) = last
           end if

           !!  Set maxerr and ermax.  !!

           maxerr = iord(nrmax)
           ermax = elist(maxerr)  
           
           return

       end Subroutine DQPSRT_1 


        
         
      Subroutine Quad2(f,a,b,epsabs,epsrel,key,result,abserr,neval,ier, &
                        limit,lenw,last,iwork,work)
       !!----------------------------------------------------------------------!!
       !!                                                                      !!
       !!   begin prologue  Quad                                               !!
       !!                                                                      !!
       !!  date written   800101   (yymmdd) (original code)                    !!
       !!  revision date  830518   (yymmdd)                                    !!
       !!  Date modified to f90 971001 (yymmdd)                                !!
       !!                                                                      !!
       !!  category no.  h2a1a1                                                !!
       !!                                                                      !!
       !!  Keywords    -                                                       !!
       !!                   Automatic integrator, general-purpose,             !!
       !!                   integrand examinator, globally adaptive,           !!
       !!                   Gauss-Kronrod                                      !!
       !!                                                                      !!
       !!  Author -                                                            !!
       !!           Piessens,Robert,appl. math. & progr. div - k.u.leuven      !!
       !!           de doncker,elise,appl. math. & progr. div. - k.u.leuven    !!
       !!                                                                      !! 
       !!  Modified to f90 -                                                   !! 
       !!                   Brian Nesbitt                                      !!
       !!                   Dept. of Applied Maths & Theoretical Physics       !!
       !!                   The Queen's University of Belfast.                 !!
       !!                                                                      !!
       !!                                                                      !!
       !!                    Purpose                                           !!
       !!                   ---------                                          !!
       !!                                                                      !!
       !!    the routine calculates an approximation result to a given         !!
       !!    definite integral i = integral of f over (a,b),                   !!
       !!    hopefully satisfying following claim for accuracy                 !!
       !!    abs(i-result)le.max(epsabs,epsrel*abs(i)).                        !!
       !!                                                                      !!
       !!                                                                      !!
       !!                                                                      !!
       !!                    Description                                       !!
       !!                   -------------                                      !!
       !!                                                                      !!
       !!    Computation of a definite integral                                !!
       !!    standard fortran subroutine                                       !!
       !!    Fortran90 precision version                                       !!
       !!                                                                      !!
       !!     f      -    Real                                                 !!
       !!                 Function subprogam defining the integrand            !! 
       !!                 function f(x). the actual name for f needs to be     !! 
       !!                 declared e x t e r n a l in the driver program.      !!
       !!                                                                      !!
       !!     a      -    Real                                                 !! 
       !!                 Lower limit of integration.                          !!
       !!                                                                      !!
       !!     b      -    Real                                                 !! 
       !!                 Upper limit of integration.                          !!
       !!                                                                      !!
       !!     epsabs -    Real                                                 !!
       !!                 Absolute accuracy requested.                         !!
       !!                                                                      !!
       !!     epsrel -    Real                                                 !! 
       !!                 Relative accuracy requested                          !!
       !!                 if  epsabs.le.0                                      !!
       !!                 and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),         !!
       !!                 the routine will end with ier = 6.                   !!
       !!                                                                      !!
       !!     key    -    Integer                                              !!
       !!                 Key for choice of local integration rule             !!
       !!                 a gauss-kronrod pair is used with                    !!
       !!                   7 - 15 points if key.lt.2,                         !!
       !!                  10 - 21 points if key = 2,                          !!
       !!                  15 - 31 points if key = 3,                          !!
       !!                  20 - 41 points if key = 4,                          !!
       !!                  25 - 51 points if key = 5,                          !!
       !!                  30 - 61 points if key .gt. 5.                       !!
       !!                                                                      !!
       !!                    On return                                         !!
       !!                   -----------                                        !! 
       !!                                                                      !!
       !!     result -    Real                                                 !!
       !!                 approximation to the integral.                       !!
       !!                                                                      !!
       !!     abserr -    Real                                                 !!
       !!                 estimate of the modulus of the absolute error,       !!
       !!                 which should equal or exceed abs(i-result).          !!
       !!                                                                      !!
       !!     neval  -    Integer                                              !!
       !!                 number of integrand evaluations.                     !!
       !!                                                                      !!
       !!     ier    -    Integer                                              !!
       !!                 ier = 0 normal and reliable termination of the       !!
       !!                      routine. it is assumed that the requested       !!
       !!                      accuracy has been achieved.                     !!
       !!                                                                      !!
       !!                 ier.gt.0 abnormal termination of the routine         !!
       !!                      the estimates for result and error are          !!
       !!                      less reliable. it is assumed that the           !!
       !!                      requested accuracy has not been achieved.       !!
       !!                                                                      !!
       !!                                                                      !!
       !!                    Error Messages                                    !!
       !!                   ----------------                                   !!
       !!                                                                      !!
       !!              ier = 1 maximum number of subdivisions allowed          !!
       !!                      has been achieved. one can allow more           !!
       !!                      subdivisions by increasing the value of         !!
       !!                      limit (and taking the according dimension       !!
       !!                      adjustments into account). however, if          !!
       !!                      this yield no improvement it is advised         !!
       !!                      to analyze the integrand in order to            !!
       !!                      determine the integration difficulaties.        !!
       !!                      if the position of a local difficulty can       !!
       !!                      be determined (i.e.singularity,                 !!
       !!                      discontinuity within the interval) one          !!
       !!                      will probably gain from splitting up the        !!
       !!                      interval at this point and calling the          !!
       !!                      integrator on the subranges. if possible,       !!
       !!                      an appropriate special-purpose integrator       !!
       !!                      should be used which is designed for            !!
       !!                      handling the type of difficulty involved.       !!
       !!                                                                      !!
       !!                  = 2 the occurrence of roundoff error is             !!
       !!                      detected, which prevents the requested          !!
       !!                      tolerance from being achieved.                  !!
       !!                                                                      !!
       !!                  = 3 extremely bad integrand behaviour occurs        !!
       !!                      at some points of the integration               !!
       !!                      interval.                                       !!
       !!                                                                      !!
       !!                  = 6 the input is invalid, because                   !!
       !!                      (epsabs.le.0 and                                !!
       !!                       epsrel.lt.max(50*rel.mach.acc.,0.5d-28))       !!
       !!                      or limit.lt.1 or lenw.lt.limit*4.               !!
       !!                      result, abserr, neval, last are set             !!
       !!                      to zero.                                        !!
       !!                      except when lenw is invalid, iwork(1),          !!
       !!                      work(limit*2+1) and work(limit*3+1) are         !!
       !!                      set to zero, work(1) is set to a and            !!
       !!                      work(limit+1) to b.                             !!
       !!                                                                      !!
       !!                                                                      !!
       !!                    Dimensioning Parameters                           !!
       !!                   -------------------------                          !!
       !!                                                                      !!
       !!                                                                      !!
       !!     limit -   Integer                                                !!
       !!               dimensioning parameter for iwork                       !!
       !!               limit determines the maximum number of subintervals    !!
       !!               in the partition of the given integration interval     !!
       !!               (a,b), limit.ge.1.                                     !!
       !!               if limit.lt.1, the routine will end with ier = 6.      !!
       !!                                                                      !!
       !!     lenw  -   Integer                                                !!
       !!               dimensioning parameter for work                        !!
       !!               lenw must be at least limit*4.                         !!
       !!               if lenw.lt.limit*4, the routine will end with          !!
       !!               ier = 6.                                               !!
       !!                                                                      !!
       !!     last  -   Integer                                                !!
       !!               on return, last equals the number of subintervals      !!
       !!               produced in the subdiviosion process, which            !!
       !!               determines the number of significant elements          !!
       !!               actually in the work arrays.                           !!
       !!                                                                      !!
       !!                                                                      !!
       !!                    Work Arrays                                       !!
       !!                   -------------                                      !!
       !!                                                                      !!
       !!                                                                      !!
       !!     iwork -   Integer                                                !!
       !!               vector of dimension at least limit, the first k        !!
       !!               elements of which contain pointers to the error        !!
       !!               estimates over the subintervals, such that             !!
       !!               work(limit*3+iwork(1)),... , work(limit*3+iwork(k))    !!
       !!               form a decreasing sequence with k = last if            !!
       !!               last.le.(limit/2+2), and k = limit+1-last otherwise    !!
       !!                                                                      !!
       !!     work  -   Real                                                   !!
       !!               vector of dimension at least lenw                      !!
       !!               on return                                              !!
       !!               work(1), ..., work(last) contain the left end          !!
       !!               points of the subintervals in the partition of         !!
       !!               (a,b),                                                 !!
       !!               work(limit+1), ..., work(limit+last) contain the       !!
       !!               right end points,                                      !!
       !!               work(limit*2+1), ..., work(limit*2+last) contain       !!
       !!               the integral approximations over the subintervals,     !!
       !!               work(limit*3+1), ..., work(limit*3+last) contain       !!
       !!               the error estimates.                                   !!
       !!                                                                      !!
       !!             References  (none)                                       !!
       !!                                                                      !!
       !!          Routines called  Dqage,Xerror                               !!
       !!                                                                      !! 
       !! end prologue  dqag                                                   !!
       !!                                                                      !!
       !!----------------------------------------------------------------------!!

           implicit none

           real(kind=mp), intent(in) :: a, b, epsabs, epsrel
           real(kind=mp), intent(out) :: abserr, result, work(*)

           integer, intent(in) :: lenw, limit, key
           integer, intent(out) :: ier, iwork(*), last, neval
           integer :: lvl, l1, l2, l3

           real(kind=mp), external :: f  ! external function defining the integrand.
       
           !!  Check validity of lenw. !!

           ier = 6
           neval = 0
           last = neval 
           result = 0.0_mp 
           abserr = result 

           if (limit .lt. 1 .or. lenw .lt. 4*limit) then
 
               if (ier .eq. 6) then
                  lvl = 1
               end if

           else
         
              !!  Prepare for call to DQAGE  !!

              l1 = limit + 1
              l2 = limit + l1
              l3 = limit + l2

              call DQAGE2(f,a,b,epsabs,epsrel,key,limit,result,abserr,neval,ier, &
                         work(1),work(l1),work(l2),work(l3),iwork,last)

              !!  Call error handler if necessary.  !!

              lvl = 0

           end if 
 
           if (ier .eq. 6) then
              lvl = 1
           else if (ier .ne. 0) then 
              call Xerror(ier)  !  Quadrature has failed.
           end if

     
          return

       end Subroutine Quad2



       Subroutine DQAGE2(f,a,b,epsabs,epsrel,key,limit,result,abserr,neval,ier, &
                        alist,blist,rlist,elist,iord,last)

       !!-------------------------------------------------------------------------!!
       !!                                                                         !!
       !!    Begin prologue  dqage                                                !!
       !!                                                                         !!
       !!  Date written   800101   (yymmdd) (original code)                       !!
       !!  Revision date  830518   (yymmdd)                                       !!
       !!  Date modifies to f90 971001 (yymmdd)                                   !!
       !!                                                                         !!
       !!  category no.  h2a1a1                                                   !!
       !!                                                                         !!
       !!  Keywords   -                                                           !!
       !!                     automatic integrator, general-purpose,              !!
       !!                     integrand examinator, globally adaptive,            !!
       !!                     Gauss-Kronrod.                                      !!
       !!                                                                         !!
       !!  Author -                                                               !!
       !!           Piessens,Robert,appl. math. & progr. div - k.u.leuven         !!
       !!           de Doncker,Elise,appl. math. & progr. div. - k.u.leuven       !!
       !!                                                                         !!
       !!  Modified to f90 -                                                      !!
       !!                   Brian Nesbitt                                         !!
       !!                   Dept. of Applied Maths & Theoretical Physics          !!
       !!                   The Queen's University of Belfast.                    !!
       !!                                                                         !!
       !!                                                                         !!
       !!                    Purpose                                              !! 
       !!                   ---------                                             !!
       !!                                                                         !!
       !!   The routine calculates an approximation result to a given             !!
       !!   definite integral i = integral of f over (a,b),                       !!
       !!   hopefully satisfying following claim for accuracy                     !!
       !!   abs(i-result)le.max(epsabs,epsrel*abs(i)).                            !!
       !!                                                                         !!
       !!                                                                         !!
       !!                                                                         !!
       !!                    Description                                          !!
       !!                   -------------                                         !!
       !!                                                                         !!        
       !!                                                                         !!
       !!    Computation of a definite integral                                   !!
       !!    standard fortran subroutine.                                         !!
       !!                                                                         !!
       !!                                                                         !!        
       !!                    Parameters                                           !!
       !!                   ------------                                          !!
       !!                                                                         !!
       !!                                                                         !! 
       !!                    On Entry                                             !!
       !!                   ----------                                            !!
       !!                                                                         !!
       !!     f      -    Real                                                    !!
       !!                 function subprogram defining the integrand              !!
       !!                 function f(x). the actual name for f needs to be        !!
       !!                 declared e x t e r n a l in the driver program.         !!
       !!                                                                         !!
       !!     a      -    Real                                                    !! 
       !!                 lower limit of integration.                             !!         
       !!                                                                         !!
       !!     b      -    Real                                                    !!
       !!                 upper limit of integration.                             !!
       !!                                                                         !!
       !!     epsabs -    Real                                                    !!
       !!                 absolute accuracy requested.                            !!
       !!                                                                         !!
       !!     epsrel -    Real                                                    !!
       !!                 relative accuracy requested.                            !! 
       !!                 if  epsabs.le.0                                         !!
       !!                 and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),            !!
       !!                 the routine will end with ier = 6.                      !!
       !!                                                                         !!
       !!     key    -    Integer                                                 !!
       !!                 key for choice of local integration rule                !! 
       !!                 a Gauss-Kronrod pair is used with                       !!
       !!                     7 - 15 points if key.lt.2,                          !!
       !!                    10 - 21 points if key = 2,                           !!
       !!                    15 - 31 points if key = 3,                           !!
       !!                    20 - 41 points if key = 4,                           !!
       !!                    25 - 51 points if key = 5,                           !!
       !!                    30 - 61 points if key.gt.5.                          !!
       !!                                                                         !!
       !!     limit  -    Integer                                                 !!
       !!                 gives an upperbound on the number of subintervals       !!
       !!                 in the partition of (a,b), limit.ge.1.                  !!
       !!                                                                         !!
       !!                    On Return                                            !!
       !!                   -----------                                           !!
       !!                                                                         !!
       !!     result -    Real                                                    !!
       !!                 approximation to the integral                           !!
       !!                                                                         !!
       !!     abserr -    Real                                                    !!  
       !!                 estimate of the modulus of the absolute error,          !!
       !!                 which should equal or exceed abs(i-result).             !!
       !!                                                                         !!
       !!     neval  -    Integer                                                 !!
       !!                 number of integrand evaluations.                        !!
       !!                                                                         !!
       !!     ier    -    Integer                                                 !!
       !!                 ier = 0     normal and reliable termination of the      !!
       !!                             routine. it is assumed that the requested   !!
       !!                             accuracy has been achieved.                 !!
       !!                                                                         !!
       !!                 ier.gt.0    abnormal termination of the routine         !!
       !!                             the estimates for result and error are      !!
       !!                             less reliable. it is assumed that the       !!
       !!                             requested accuracy has not been achieved.   !!
       !!                                                                         !!
       !!                    Error Messages                                       !!
       !!                   ----------------                                      !!
       !!                                                                         !!
       !!              ier = 1 maximum number of subdivisions allowed             !!
       !!                      has been achieved. one can allow more              !!
       !!                      subdivisions by increasing the value               !!
       !!                      of limit.                                          !!
       !!                      however, if this yields no improvement it          !!
       !!                      is rather advised to analyze the integrand         !!
       !!                      in order to determine the integration              !!
       !!                      difficulties. if the position of a local           !!
       !!                      difficulty can be determined(e.g.                  !!
       !!                      singularity, discontinuity within the              !!
       !!                      interval) one will probably gain from              !!
       !!                      splitting up the interval at this point            !!
       !!                      and calling the integrator on the                  !!
       !!                      subranges. if possible, an appropriate             !!
       !!                      special-purpose integrator should be used          !!
       !!                      which is designed for handling the type of         !!
       !!                      difficulty involved.                               !!
       !!                                                                         !!
       !!                 =  2 the occurrence of roundoff error is                !!
       !!                      detected, which prevents the requested             !!
       !!                      tolerance from being achieved.                     !!
       !!                                                                         !!
       !!                  = 3 extremely bad integrand behaviour occurs           !!
       !!                      at some points of the integration                  !!
       !!                      interval.                                          !!
       !!                                                                         !!
       !!                  = 6 the input is invalid, because                      !!
       !!                      (epsabs.le.0 and                                   !!
       !!                      epsrel.lt.max(50*rel.mach.acc.,0.5d-28),           !!
       !!                      result, abserr, neval, last, rlist(1) ,            !!
       !!                      elist(1) and iord(1) are set to zero.              !!
       !!                      alist(1) and blist(1) are set to a and b           !!
       !!                      respectively.                                      !!
       !!                                                                         !!
       !!     alist   -   Real                                                    !!
       !!                 vector of dimension at least limit, the first           !!
       !!                 last  elements of which are the left                    !!
       !!                 end points of the subintervals in the partition         !!
       !!                 of the given integration range (a,b)                    !!
       !!                                                                         !!
       !!     blist   -   Real                                                    !!
       !!                 vector of dimension at least limit, the first           !!
       !!                 last  elements of which are the right                   !!
       !!                 end points of the subintervals in the partition         !!
       !!                 of the given integration range (a,b)                    !!
       !!                                                                         !!
       !!     rlist   -   Real                                                    !!
       !!                 vector of dimension at least limit, the first           !!         
       !!                 last  elements of which are the                         !!
       !!                 integral approximations on the subintervals             !! 
       !!                                                                         !!
       !!     elist   -   Real                                                    !! 
       !!                 vector of dimension at least limit, the first           !!
       !!                 last  elements of which are the moduli of the           !!
       !!                 absolute error estimates on the subintervals            !!
       !!                                                                         !!
       !!     iord    -   Integer                                                 !!
       !!                 vector of dimension at least limit, the first k         !!
       !!                 elements of which are pointers to the                   !!
       !!                 error estimates over the subintervals,                  !!
       !!                 such that elist(iord(1)), ...,                          !!
       !!                 elist(iord(k)) form a decreasing sequence,              !!
       !!                 with k = last if last.le.(limit/2+2), and               !!
       !!                 k = limit+1-last otherwise.                             !!
       !!                                                                         !!
       !!     last    -   Integer                                                 !!         
       !!                 number of subintervals actually produced in the         !!
       !!                 subdivision process.                                    !!
       !!                                                                         !!
       !!               references  (none)                                        !!
       !!            routines called  d1mach,dqk15,dqk21,dqk31,                   !!
       !!                             dqk41,dqk51,dqk61,dqpsrt                    !!
       !!                                                                         !!
       !! end prologue  dqage                                                     !!
       !!                                                                         !!
       !!-------------------------------------------------------------------------!!

           implicit none

           real(kind=mp), intent(in) :: a, b, epsrel, epsabs
           real(kind=mp), intent(out) :: abserr, result
           real(kind=mp), intent(out) :: alist(*), blist(*), elist(*), rlist(*)

           real(kind=mp) :: area, area1, area12, area2, a1, a2,        &
                            b1, b2, defabs, defab1, defab2,            &
                            errbnd, errmax, error1, erro12, errsum,    &
                            resabs, error2

           integer, intent(in) :: limit, key 
           integer, intent(out) :: ier, iord(*), last, neval

           integer :: iroff1, iroff2, k, keyf, maxerr, nrmax

           real(kind=mp), external :: f ! external function defining the integrand.

           !!---------------------------------------------------------------------!!
           !!                                                                     !!  
           !!    list of major variables                                          !!  
           !!    -----------------------                                          !! 
           !!                                                                     !!  
           !!          alist     -   Real array                                   !! 
           !!                        list of left end points of all subintervals  !!  
           !!                        considered up to now.                        !!
           !!                                                                     !!
           !!          blist     -   Real array                                   !! 
           !!                        list of right end points of all subinterval  !!
           !!                        considered up to now.                        !!
           !!                                                                     !!
           !!          rlist(i)  -   Real                                         !!
           !!                        approximation to the integral over           !!
           !!                        (alist(i),blist(i)).                         !!
           !!                                                                     !!
           !!          elist(i)  -   Real                                         !!
           !!                        error estimate applying to rlist(i).         !!
           !!                                                                     !!
           !!          maxerr    -   Integer                                      !!
           !!                        pointer to the interval with largest.        !!
           !!                                                                     !!
           !!                                                                     !!
           !!            Error Estimate                                           !!
           !!           ----------------                                          !!
           !!                                                                     !!
           !!          errmax    -   Real                                         !! 
           !!                        elist(maxerr).                               !!
           !!                                                                     !!
           !!          area      -   Real                                         !!
           !!                        sum of the integrals over the subintervals.  !!
           !!                                                                     !!
           !!          errsum    -   Real                                         !!
           !!                        sum of the errors over the subintervals.     !!
           !!                                                                     !!
           !!          errbnd    -   Real                                         !!
           !!                        requested accuracy max(epsabs,epsrel*        !!
           !!                        abs(result)).                                !!
           !!                                                                     !!
           !!          *****1    -   variable for the left subinterval.           !!
           !!          *****2    -   variable for the right subinterval.          !!
           !!                                                                     !!
           !!          last      -   Integer                                      !!
           !!                        index for subdivision.                       !!
           !!                                                                     !!
           !!                                                                     !!
           !!              Machine Dependent Constants                            !!
           !!             -----------------------------                           !!
           !!                                                                     !!
           !!      epmach  -   is the largest relative spacing.                   !!
           !!      uflow   -   is the smallest positive magnitude.                !!
           !!                                                                     !!
           !!---------------------------------------------------------------------!!
           
           !!  First executable statement of DQAGE. !!


           !!  Test on validity of parameters.      !!

           ier = 0
           neval = ier 
           last = ier 
           result = 0.0_mp 
           abserr = result 
           alist(1) = a
           blist(1) = b
           rlist(1) = result 
           elist(1) = result 
           iord(1) = ier 

           if (epsabs .le. 0.0_mp .and. epsrel .lt.     &
               max(0.5e+2_mp*epmach,0.5e-28_mp)) then
               ier = 6
           end if

           if (ier .eq. 6) then
              call Xerror(ier)  ! Quadrature has failed.
           end if

           keyf = key

           if (key .le. 0) then 
               keyf = 1
           else if (key .ge. 7) then 
               keyf = 6
           end if

           neval = 0

      !     call DQK61_2(f,a,b,result,abserr,defabs,resabs)

           select case(keyf)
                  case(1)
                    call DQK15_2(f,a,b,result,abserr,defabs,resabs)
                  case(6)
                    call DQK61_2(f,a,b,result,abserr,defabs,resabs)
           end select 

           last = 1
           rlist(1) = result
           elist(1) = abserr
           iord(1) = 1

           !!  Test on accuracy.  !!

           errbnd = max(epsabs,epsrel*abs(result))

           if (abserr .le. 0.5e+2_mp*epmach*defabs .and. abserr .gt. errbnd) then 
               ier = 2 
           end if

           if (limit .eq. 1) then 
               ier = 1
           end if 

           if (ier .ne. 0 .or. (abserr .le. errbnd .and. abserr    &
              .ne. resabs) .or. abserr .eq. 0.0_mp) then
              
              if (keyf .ne. 1) then
                  neval = (10*keyf + 1)*(2*neval + 1)
              else if (keyf .eq. 1) then
                  neval = 30*neval + 15
              end if

              return
           end if

           !!  Initialization  !!

           errmax = abserr
           maxerr = 1
           area = result
           errsum = abserr
           nrmax = 1
           iroff1 = 0
           iroff2 = 0

           !!  Main do loop  !!
    
           do last = 2,limit

               !!  Bisect the subinterval with the largest error estimate.  !!

               a1 = alist(maxerr)
               b1 = 0.5_mp*(alist(maxerr) + blist(maxerr))
               a2 = b1
               b2 = blist(maxerr)

               !!  Use 61-point rule only at this stage.  !!

               select case(keyf)
                      case(1)
                        call DQK15_2(f,a1,b1,area1,error1,resabs,defab1)
                        call DQK15_2(f,a2,b2,area2,error2,resabs,defab2)
                      case(6)
                        call DQK61_2(f,a1,b1,area1,error1,resabs,defab1)
                        call DQK61_2(f,a2,b2,area2,error2,resabs,defab2)
               end select 

               !!  Improve the previous approximations to the integral  !!
               !!  and error and test for accuracy.                     !!

               neval = neval + 1
               area12 = area1 + area2
               erro12 = error1 + error2
               errsum = errsum + erro12 - errmax
               area = area + area12 - rlist(maxerr)

               if (defab1 .ne. error1 .and. defab2 .ne. error2) then

                   if (abs(rlist(maxerr) - area12) .le. 0.1e-4_mp*abs(area12) &
                      .and. erro12 .ge. 0.99_mp*errmax) then
                      iroff1 = iroff1 + 1
                   end if

                   if (last .gt. 10 .and. erro12 .gt. errmax) then
                      iroff2 = iroff2 + 1
                   end if

               end if

               rlist(maxerr) = area1
               rlist(last)  = area2

               errbnd = max(epsabs,epsrel*abs(area))

               if (errsum .gt. errbnd) then 

                  !!  Test for roundoff error and eventually set  !!
                  !!  error flag.                                 !!
                  
                  if (iroff1 .ge. 6 .or. iroff2 .ge. 20) then 
                      ier = 2
                  end if

                  !!  Set error flag in the case that the number  !!
                  !!  of subintervals equals limit.               !!

                  if (last .eq. limit) then 
                     ier = 1
                  end if

                  !!  Set error flag in the case of bad integrand     !!   
                  !!  behaviour at a point of the integration range.  !!

                  if (max(abs(a1),abs(b2)) .le.  & 
                     (0.1e+1_mp + 0.1e+3_mp*epmach)*(abs(a2) + 0.1e+4_mp*uflow)) then 
                      ier = 3
                  end if
        
               end if 

               !!  Apend the newly-created intervals to the list.  !!

               if (error2 .gt. error1) then 
                   alist(maxerr) = a2
                   alist(last) = a1
                   blist(last) = b1
                   rlist(maxerr) = area2
                   rlist(last) = area1
                   elist(maxerr) = error2
                   elist(last) = error1
               else
                   alist(last) = a2
                   blist(maxerr) = b1
                   blist(last) = b2
                   elist(maxerr) = error1
                   elist(last) = error2
               end if

               !!  Call subroutine DQPSRT to maintain the descending ordering  !!
               !!  in the list of error estimates and select the subinterval   !!
               !!  with the largest error estimate (to be bisected next).      !!

               call DQPSRT_2(limit,last,maxerr,errmax,elist,iord,nrmax)

               !!  Exit do loop.  !!
 
               if (ier .ne. 0 .or. errsum .le. errbnd) then 
                  exit
               end if

           end do

           !!  Compare the final result. !!

           result = 0.0_mp 

           do k = 1,last
              result = result + rlist(k)
           end do
 
           abserr = errsum

           if (keyf .ne. 1) then
               neval = (10*keyf + 1)*(2*neval + 1)
           end if

           if (keyf .eq. 1) then 
               neval = 30*neval + 15
           end if

           return

       end Subroutine DQAGE2



       Subroutine DQK15_2(f,a,b,result,abserr,resabs,resasc)

           implicit none

           real(kind=mp), intent(in) :: a, b
           real(kind=mp), intent(out) :: abserr, resabs, resasc, result            

           real(kind=mp) :: absc, centr, dhlgth, fc, fsum, fval1,     &
                            fval2, hlgth, resg, resk, reskh

           real(kind=mp) :: fv1(7), fv2(7), wg(4), wgk(8), xgk(8)

           integer :: j, jtw, jtwm1

           real(kind=mp), external :: f

           !!-------------------------------------------------------------------------!!
           !!                                                                         !!
           !!  The abscissae and weights are given for the interval (-1,1).           !!
           !!  Because of symmetry only the positive abscissae and their              !!
           !!  corresponding weights are given.                                       !!
           !!                                                                         !!
           !!  xgk    - abscissae of the 15-point Kronrod rule                        !!
           !!         xgk(2), xgk(4), ...  abscissae of the 7-point                   !!
           !!         Gauss rule                                                      !!
           !!         xgk(1), xgk(3), ...  abscissae which are optimally              !!
           !!         added to the 7-point Gauss rule                                 !!
           !!                                                                         !!
           !!  wgk    - weights of the 15-point Kronrod rule                          !!
           !!                                                                         !!
           !!  wg     - weights of the 7-point Gauss rule                             !!
           !!                                                                         !!
           !!                                                                         !!
           !!  Gauss quadrature weights and Kronrod quadrature abscissae and weights  !!
           !!  as evaluated with 80 decimal digit arithmetic by L. W. Fullerton,      !!
           !!  Bell labs, Nov. 1981.                                                  !!
           !!                                                                         !!
           !!-------------------------------------------------------------------------!!

           data wg  (1) / 0.129484966168869693270611432679082_mp  /
           data wg  (2) / 0.279705391489276667901467771423780_mp  /
           data wg  (3) / 0.381830050505118944950369775488975_mp  /
           data wg  (4) / 0.417959183673469387755102040816327_mp  /

           data xgk (1) / 0.991455371120812639206854697526329_mp  /
           data xgk (2) / 0.949107912342758524526189684047851_mp  /
           data xgk (3) / 0.864864423359769072789712788640926_mp  /
           data xgk (4) / 0.741531185599394439863864773280788_mp  /
           data xgk (5) / 0.586087235467691130294144838258730_mp  /
           data xgk (6) / 0.405845151377397166906606412076961_mp  /
           data xgk (7) / 0.207784955007898467600689403773245_mp  /
           data xgk (8) / 0.000000000000000000000000000000000_mp  /

           data wgk (1) / 0.022935322010529224963732008058970_mp  /
           data wgk (2) / 0.063092092629978553290700663189204_mp  /
           data wgk (3) / 0.104790010322250183839876322541518_mp  /
           data wgk (4) / 0.140653259715525918745189590510238_mp  /
           data wgk (5) / 0.169004726639267902826583426598550_mp  /
           data wgk (6) / 0.190350578064785409913256402421014_mp  /
           data wgk (7) / 0.204432940075298892414161999234649_mp  /
           data wgk (8) / 0.209482141084727828012999174891714_mp  /

           !!  First executable statements of DQK15.  !!

           centr = 0.5_mp*(a + b)
           hlgth = 0.5_mp*(b - a)
           dhlgth = abs(hlgth)

           !!  Compute the 15-point Kronrod approximation to the  !!
           !!  integral, and estimate the absolute error.         !!

           fc = f(centr)
           resg = fc*wg(4)
           resk = fc*wgk(8)
           resabs = abs(resk)

           do j = 1,3
              jtw = j*2
              absc = hlgth*xgk(jtw)
              fval1 = f(centr - absc)
              fval2 = f(centr + absc)
              fv1(jtw) = fval1
              fv2(jtw) = fval2
              fsum = fval1 + fval2
              resg = resg + wg(j)*fsum
              resk  = resk + wgk(jtw)*fsum
              resabs = resabs + wgk(jtw)*(abs(fval1) + abs(fval2))
           end do

           do j = 1,4
              jtwm1 = j*2 - 1
              absc = hlgth*xgk(jtwm1)
              fval1 = f(centr - absc)
              fval2 = f(centr + absc)
              fv1(jtwm1) = fval1
              fv2(jtwm1) = fval2
              fsum = fval1 + fval2
              resk = resk + wgk(jtwm1)*fsum
              resabs = resabs + wgk(jtwm1)*(abs(fval1) + abs(fval2))
           end do

           reskh = resk*0.5_mp
           resasc = wgk(8)*abs(fc - reskh)

           do j = 1,7
              resasc = resasc + wgk(j)*(abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
           end do

           result = resk*hlgth
           resabs = resabs*dhlgth
           resasc = resasc*dhlgth
           abserr = abs((resk - resg)*hlgth)

           if (resasc .ne. 0.0_mp .and. abserr .ne. 0.0_mp) then
               abserr = resasc*min(0.1e+1_mp,(0.2e+3_mp*abserr/resasc)**1.5_mp)
           end if

           if (resabs .gt. uflow/(0.5e+2_mp*epmach)) then
               abserr = max((epmach*0.5e+2_mp)*resabs,abserr)
           end if

           return

       end Subroutine DQK15_2 



       Subroutine DQK61_2(f,a,b,result,abserr,resabs,resasc)
            
           implicit none

           real(kind=mp), intent(in) :: a, b
           real(kind=mp), intent(out) :: abserr, resabs, resasc, result

           real(kind=mp) :: dabsc, centr, dhlgth,            &        
                            fc, fsum, fval1, fval2, hlgth,   &
                            resg, resk, reskh

           real(kind=mp) :: fv1(30), fv2(30), wg(15), wgk(31), xgk(31)

           real(kind=mp), external :: f

           integer :: j, jtw, jtwm1
       
           !!---------------------------------------------------------!!
           !!                                                         !!
           !!  The abscissae and weights are given for the interval   !!
           !!  (-1,+1). Because of symmetry only the positive         !!  
           !!  abscissae and their corresponding weights are given.   !! 
           !!                                                         !!
           !!                                                         !!
           !!   xgk  -  abscissae of the 51-point Kronrod rule:       !!
           !!             xgk(2), xgk(4), ...                         !!
           !!                                                         !! 
           !!           abscissae of the 25-point Gauss rule:         !! 
           !!             xgk(1), xgk(3), ...                         !!
           !!                                                         !!
           !!           abscissae which are optimally added to        !!
           !!           the 25-point Gauss rule.                      !!
           !!                                                         !!
           !!   wgk  -  weights of the 51-point Kronrod rule.         !!
           !!                                                         !!
           !!   wg   -  weights of the 25-point Gauss rule.           !!
           !!                                                         !!
           !!   Gauss quadrature weights and Kronrod Quadrature       !!
           !!   abscissae and weights as evaluated with 80 decimal    !!
           !!   arithmetic by L.W. Fullerton, Bell Labs., Nov 1981.   !!
           !!                                                         !!
           !!---------------------------------------------------------!!
  
           data wg  (  1) / 0.007968192496166605615465883474674_mp  /
           data wg  (  2) / 0.018466468311090959142302131912047_mp  /
           data wg  (  3) / 0.028784707883323369349719179611292_mp  /
           data wg  (  4) / 0.038799192569627049596801936446348_mp  /
           data wg  (  5) / 0.048402672830594052902938140422808_mp  /
           data wg  (  6) / 0.057493156217619066481721689402056_mp  /
           data wg  (  7) / 0.065974229882180495128128515115962_mp  /
           data wg  (  8) / 0.073755974737705206268243850022191_mp  /
           data wg  (  9) / 0.080755895229420215354694938460530_mp  /
           data wg  ( 10) / 0.086899787201082979802387530715126_mp  /
           data wg  ( 11) / 0.092122522237786128717632707087619_mp  /
           data wg  ( 12) / 0.096368737174644259639468626351810_mp  /
           data wg  ( 13) / 0.099593420586795267062780282103569_mp  /
           data wg  ( 14) / 0.101762389748405504596428952168554_mp  /
           data wg  ( 15) / 0.102852652893558840341285636705415_mp  /

           data xgk (  1) / 0.999484410050490637571325895705811_mp  /
           data xgk (  2) / 0.996893484074649540271630050918695_mp  /
           data xgk (  3) / 0.991630996870404594858628366109486_mp  /
           data xgk (  4) / 0.983668123279747209970032581605663_mp  /
           data xgk (  5) / 0.973116322501126268374693868423707_mp  /
           data xgk (  6) / 0.960021864968307512216871025581798_mp  /
           data xgk (  7) / 0.944374444748559979415831324037439_mp  /
           data xgk (  8) / 0.926200047429274325879324277080474_mp  /
           data xgk (  9) / 0.905573307699907798546522558925958_mp  /
           data xgk ( 10) / 0.882560535792052681543116462530226_mp  /
           data xgk ( 11) / 0.857205233546061098958658510658944_mp  /
           data xgk ( 12) / 0.829565762382768397442898119732502_mp  /
           data xgk ( 13) / 0.799727835821839083013668942322683_mp  /
           data xgk ( 14) / 0.767777432104826194917977340974503_mp  /
           data xgk ( 15) / 0.733790062453226804726171131369528_mp  /
           data xgk ( 16) / 0.697850494793315796932292388026640_mp  /
           data xgk ( 17) / 0.660061064126626961370053668149271_mp  /
           data xgk ( 18) / 0.620526182989242861140477556431189_mp  /
           data xgk ( 19) / 0.579345235826361691756024932172540_mp  /
           data xgk ( 20) / 0.536624148142019899264169793311073_mp  /
           data xgk ( 21) / 0.492480467861778574993693061207709_mp  /
           data xgk ( 22) / 0.447033769538089176780609900322854_mp  /
           data xgk ( 23) / 0.400401254830394392535476211542661_mp  /
           data xgk ( 24) / 0.352704725530878113471037207089374_mp  /
           data xgk ( 25) / 0.304073202273625077372677107199257_mp  /
           data xgk ( 26) / 0.254636926167889846439805129817805_mp  /
           data xgk ( 27) / 0.204525116682309891438957671002025_mp  /
           data xgk ( 28) / 0.153869913608583546963794672743256_mp  /
           data xgk ( 29) / 0.102806937966737030147096751318001_mp  /
           data xgk ( 30) / 0.051471842555317695833025213166723_mp  /
           data xgk ( 31) / 0.000000000000000000000000000000000_mp  /

           data wgk (  1) / 0.001389013698677007624551591226760_mp  /
           data wgk (  2) / 0.003890461127099884051267201844516_mp  /
           data wgk (  3) / 0.006630703915931292173319826369750_mp  /
           data wgk (  4) / 0.009273279659517763428441146892024_mp  /
           data wgk (  5) / 0.011823015253496341742232898853251_mp  /
           data wgk (  6) / 0.014369729507045804812451432443580_mp  /
           data wgk (  7) / 0.016920889189053272627572289420322_mp  /
           data wgk (  8) / 0.019414141193942381173408951050128_mp  /
           data wgk (  9) / 0.021828035821609192297167485738339_mp  /
           data wgk ( 10) / 0.024191162078080601365686370725232_mp  /
           data wgk ( 11) / 0.026509954882333101610601709335075_mp  /
           data wgk ( 12) / 0.028754048765041292843978785354334_mp  /
           data wgk ( 13) / 0.030907257562387762472884252943092_mp  /
           data wgk ( 14) / 0.032981447057483726031814191016854_mp  /
           data wgk ( 15) / 0.034979338028060024137499670731468_mp  /
           data wgk ( 16) / 0.036882364651821229223911065617136_mp  /
           data wgk ( 17) / 0.038678945624727592950348651532281_mp  /
           data wgk ( 18) / 0.040374538951535959111995279752468_mp  /
           data wgk ( 19) / 0.041969810215164246147147541285970_mp  /
           data wgk ( 20) / 0.043452539701356069316831728117073_mp  /
           data wgk ( 21) / 0.044814800133162663192355551616723_mp  /
           data wgk ( 22) / 0.046059238271006988116271735559374_mp  /
           data wgk ( 23) / 0.047185546569299153945261478181099_mp  /
           data wgk ( 24) / 0.048185861757087129140779492298305_mp  /
           data wgk ( 25) / 0.049055434555029778887528165367238_mp  /
           data wgk ( 26) / 0.049795683427074206357811569379942_mp  /
           data wgk ( 27) / 0.050405921402782346840893085653585_mp  /
           data wgk ( 28) / 0.050881795898749606492297473049805_mp  /
           data wgk ( 29) / 0.051221547849258772170656282604944_mp  /
           data wgk ( 30) / 0.051426128537459025933862879215781_mp  /
           data wgk ( 31) / 0.051494729429451567558340433647099_mp  /
           
           !!-------------------------------------------------------!!
           !!                                                       !!
           !!     List of major variables                           !!          
           !!     -----------------------                           !!
           !!                                                       !!
           !!     centr  -  mid point of the interval               !! 
           !!     hlgth  -  half-length of the interval             !!
           !!     dabsc  -  abscissa                                !!
           !!     fval*  -  function value                          !!
           !!     resg   -  result of the 30-point gauss rule       !!
           !!     resk   -  result of the 61-point kronrod rule     !!
           !!     reskh  -  approximation to the mean value of f    !! 
           !!               over (a,b), i.e. to i/(b-a)             !! 
           !!                                                       !! 
           !!     Machine Dependent Constants                       !!
           !!     ---------------------------                       !!
           !!                                                       !!
           !!     epmach is the largest relative spacing.           !!
           !!     uflow is the smallest positive magnitude.         !! 
           !!                                                       !!
           !!-------------------------------------------------------!!
      
           centr = 0.5_mp*(b + a)
           hlgth = 0.5_mp*(b - a)
           dhlgth = abs(hlgth)

           !!  Compute the 61-point Kronrod approximation to the  !!
           !!  integral, and estimate the absolute error.         !!

           !!  First executable statement of DQK61.  !!

           resg = 0.0_mp 
           fc = f(centr)
           resk = wgk(31)*fc
           resabs = abs(resk)

           do j = 1,15
               jtw = j*2
               dabsc = hlgth*xgk(jtw)
               fval1 = f(centr - dabsc)
               fval2 = f(centr + dabsc)
               fv1(jtw) = fval1
               fv2(jtw) = fval2
               fsum = fval1 + fval2
               resg = resg + wg(j)*fsum
               resk = resk + wgk(jtw)*fsum
               resabs = resabs + wgk(jtw)*(abs(fval1) + abs(fval2))
           end do    

           do j = 1,15
               jtwm1 = j*2 - 1
               dabsc = hlgth*xgk(jtwm1)
               fval1 = f(centr - dabsc)
               fval2 = f(centr + dabsc)
               fv1(jtwm1) = fval1
               fv2(jtwm1) = fval2
               fsum = fval1 + fval2
               resk = resk + wgk(jtwm1)*fsum                    
               resabs = resabs + wgk(jtwm1)*(abs(fval1) + abs(fval2))
           end do 

           reskh = resk*0.5_mp
           resasc = wgk(31)*abs(fc - reskh)

           do j = 1,30,2
               resasc = resasc + wgk(j)*(abs(fv1(j) - reskh) + &
                        abs(fv2(j) - reskh))
               resasc = resasc + wgk(j+1)*(abs(fv1(j+1) - reskh) + &
                        abs(fv2(j+1) - reskh))
           end do

           result = resk*hlgth
           resabs = resabs*dhlgth
           resasc = resasc*dhlgth
           abserr = abs((resk - resg)*hlgth)

           if (resasc .ne. 0.0_mp .and. abserr .ne. 0.0_mp) then 
               abserr = resasc*min(0.1e+1_mp,(0.2e+3_mp*abserr/resasc)**1.5_mp)
           end if

           if (resabs .gt. uflow/(0.5e+2_mp*epmach)) then
              abserr = max((epmach*0.5e+2_mp)*resabs,abserr)
           end if

           return

       end Subroutine DQK61_2



       Subroutine DQPSRT_2(limit,last,maxerr,ermax,elist,iord,nrmax)
   
           !!----------------------------------------------------------------------!!
           !!                                                                      !!
           !!   Begin prologue DQPSRT                                              !!       
           !!                                                                      !!
           !!   Refer to  dqage,dqagie,dqagpe,dqawse                               !!
           !!   Routines called  (none)                                            !!
           !!   Revision date  810101   (yymmdd)                                   !!
           !!                                                                      !!
           !!   keywords  -   sequential sorting.                                  !!
           !!                                                                      !!
           !!   Author  Piessens,Robert,appl. math. & progr. div. - k.u.leuven     !!
           !!   de doncker,elise,appl. math. & progr. div. - k.u.leuven.           !!
           !!                                                                      !!
           !!                                                                      !!
           !!      Purpose                                                         !!
           !!     ---------                                                        !!
           !!                                                                      !!
           !!      This routine maintains the descending ordering in the           !!
           !!      list of the local error estimated resulting from the            !!
           !!      interval subdivision process. at each call two error            !!
           !!      estimates are inserted using the sequential search              !!
           !!      method, top-down for the largest error estimate and             !!
           !!      bottom-up for the smallest error estimate.                      !!
           !!                                                                      !!
           !!                                                                      !!
           !!      Description                                                     !!
           !!     -------------                                                    !!
           !!                                                                      !!
           !!      ordering routine                                                !!
           !!      standard fortran subroutine                                     !!
           !!                                                                      !!
           !!      Parameters:                                                     !!
           !!     -------------                                                    !!
           !!                                                                      !! 
           !!                  Meaning At Output                                   !!
           !!                 -------------------                                  !!
           !!                                                                      !!
           !!         limit  -   Integer                                           !! 
           !!                    Maximum number of error estimates the list        !!
           !!                    can contain                                       !!
           !!                                                                      !! 
           !!         last   -   Integer                                           !!
           !!                    Number of error estimates currently in the list.  !!
           !!                                                                      !!
           !!         maxerr -   Integer                                           !!
           !!                    maxerr points to the nrmax-th largest error       !! 
           !!                    estimate currently in the list                    !!
           !!                                                                      !!
           !!         ermax  -   Real                                              !!
           !!                    nrmax-th largest error estimate.                  !!
           !!                    ermax = elist(maxerr)                             !!
           !!                                                                      !!
           !!         elist  -   Real                                              !!
           !!                    Vector of dimension last containing               !!
           !!                    the error estimates.                              !!
           !!                                                                      !!
           !!         iord   -   Integer                                           !!
           !!                    Vector of dimension last, the first k elements    !!
           !!                    of which contain pointers to the error            !!
           !!                    estimates, such that                              !!
           !!                    elist(iord(1)),...,  elist(iord(k))               !!
           !!                    form a decreasing sequence, with                  !!
           !!                    k = last if last.le.(limit/2+2), and              !!
           !!                    k = limit+1-last otherwise                        !!
           !!                                                                      !!
           !!         nrmax  -   Integer                                           !!
           !!                                                                      !!
           !!         maxerr =   iord(nrmax)                                       !!
           !!                                                                      !!
           !! end prologue  dqpsrt                                                 !!
           !!                                                                      !!
           !!----------------------------------------------------------------------!!


           implicit none

           real(kind=mp), intent(inout) :: ermax, elist(*)      
      
           real(kind=mp) :: errmax, errmin

           integer, intent(inout) :: iord(*), maxerr, nrmax
           integer, intent(in) :: limit

           integer :: i, ibeg, last, isucc, j, jbnd, jupbn, k, ido 

           !!  Check whether the list contains more than  !!
           !!  two error estimates.                       !!
            
           if (last .le. 2) then
               iord(1) = 1
               iord(2) = 2
               maxerr = iord(nrmax)
               ermax = elist(maxerr)
               return
           end if

           !!  This part of the routine is only executed if, due to a    !!
           !!  difficult integrand, subdivision increased the error      !!
           !!  estimate. In the normal case the insert procedure should  !!
           !!  start after the nrmax-th largest error estimate.          !!

           errmax = elist(maxerr)

           if (nrmax .ne. 1) then
               ido = nrmax - 1

               do i = 1, ido
                  isucc = iord(nrmax - 1)
                  
                  !!  Exit from do  loop.  !!
                  
                  if (errmax .le. elist(isucc)) then
                      exit
                  end if
            
                  iord(nrmax) = isucc
                  nrmax = nrmax - 1
               end do
           end if

           !!  Compute the number of elements in the list to be maintained  !!
           !!  in descending order. This number depends on the number of    !!
           !!  subdivisions still allowed.                                  !!

           jupbn = last

           if (last .gt. (limit/2 + 2)) then 
               jupbn = limit + 3 - last
           end if

           errmin = elist(last)

           !!  Insert errmax by traversing the list top-down, starting  !!
           !!  comparison from the element list elist(iord(nrmax + 1)). !!

           jbnd = jupbn - 1
           ibeg = nrmax + 1

           if (ibeg .le. jbnd) then 

               do i = ibeg, jbnd
                  isucc = iord(i)

                  !!  Exit from do loop.  !!

                  if (errmax .ge. elist(isucc)) then 
                      exit
                  end if

                  iord(i-1) = isucc
               end do

               if (errmax .ge. elist(isucc)) then 

                   !! Insert errmin by traversing the list bottom up.  !!

                   iord(i-1) = maxerr
                   k = jbnd

                   do j = i, jbnd
                      isucc = iord(k)

                      !!  Exit from do loop.  !!

                      if (errmin .lt. elist(isucc)) then 
                          exit
                      end if
           
                      iord(k+1) = isucc
                      k = k - 1
                   end do

                   if (errmin .lt. elist(isucc)) then
                       iord(k+1) = last
                   else
                       iord(i) = last
                   end if
          
               else
                   iord(jbnd) = maxerr
                   iord(jupbn) = last
               end if

           else
               iord(jbnd) = maxerr
               iord(jupbn) = last
           end if

           !!  Set maxerr and ermax.  !!

           maxerr = iord(nrmax)
           ermax = elist(maxerr)  
           
           return

       end Subroutine DQPSRT_2


        
      Subroutine Quad3(f,a,b,epsabs,epsrel,key,result,abserr,neval,ier, &
                        limit,lenw,last,iwork,work)
       !!----------------------------------------------------------------------!!
       !!                                                                      !!
       !!   begin prologue  Quad                                               !!
       !!                                                                      !!
       !!  date written   800101   (yymmdd) (original code)                    !!
       !!  revision date  830518   (yymmdd)                                    !!
       !!  Date modified to f90 971001 (yymmdd)                                !!
       !!                                                                      !!
       !!  category no.  h2a1a1                                                !!
       !!                                                                      !!
       !!  Keywords    -                                                       !!
       !!                   Automatic integrator, general-purpose,             !!
       !!                   integrand examinator, globally adaptive,           !!
       !!                   Gauss-Kronrod                                      !!
       !!                                                                      !!
       !!  Author -                                                            !!
       !!           Piessens,Robert,appl. math. & progr. div - k.u.leuven      !!
       !!           de doncker,elise,appl. math. & progr. div. - k.u.leuven    !!
       !!                                                                      !! 
       !!  Modified to f90 -                                                   !! 
       !!                   Brian Nesbitt                                      !!
       !!                   Dept. of Applied Maths & Theoretical Physics       !!
       !!                   The Queen's University of Belfast.                 !!
       !!                                                                      !!
       !!                                                                      !!
       !!                    Purpose                                           !!
       !!                   ---------                                          !!
       !!                                                                      !!
       !!    the routine calculates an approximation result to a given         !!
       !!    definite integral i = integral of f over (a,b),                   !!
       !!    hopefully satisfying following claim for accuracy                 !!
       !!    abs(i-result)le.max(epsabs,epsrel*abs(i)).                        !!
       !!                                                                      !!
       !!                                                                      !!
       !!                                                                      !!
       !!                    Description                                       !!
       !!                   -------------                                      !!
       !!                                                                      !!
       !!    Computation of a definite integral                                !!
       !!    standard fortran subroutine                                       !!
       !!    Fortran90 precision version                                       !!
       !!                                                                      !!
       !!     f      -    Real                                                 !!
       !!                 Function subprogam defining the integrand            !! 
       !!                 function f(x). the actual name for f needs to be     !! 
       !!                 declared e x t e r n a l in the driver program.      !!
       !!                                                                      !!
       !!     a      -    Real                                                 !! 
       !!                 Lower limit of integration.                          !!
       !!                                                                      !!
       !!     b      -    Real                                                 !! 
       !!                 Upper limit of integration.                          !!
       !!                                                                      !!
       !!     epsabs -    Real                                                 !!
       !!                 Absolute accuracy requested.                         !!
       !!                                                                      !!
       !!     epsrel -    Real                                                 !! 
       !!                 Relative accuracy requested                          !!
       !!                 if  epsabs.le.0                                      !!
       !!                 and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),         !!
       !!                 the routine will end with ier = 6.                   !!
       !!                                                                      !!
       !!     key    -    Integer                                              !!
       !!                 Key for choice of local integration rule             !!
       !!                 a gauss-kronrod pair is used with                    !!
       !!                   7 - 15 points if key.lt.2,                         !!
       !!                  10 - 21 points if key = 2,                          !!
       !!                  15 - 31 points if key = 3,                          !!
       !!                  20 - 41 points if key = 4,                          !!
       !!                  25 - 51 points if key = 5,                          !!
       !!                  30 - 61 points if key .gt. 5.                       !!
       !!                                                                      !!
       !!                    On return                                         !!
       !!                   -----------                                        !! 
       !!                                                                      !!
       !!     result -    Real                                                 !!
       !!                 approximation to the integral.                       !!
       !!                                                                      !!
       !!     abserr -    Real                                                 !!
       !!                 estimate of the modulus of the absolute error,       !!
       !!                 which should equal or exceed abs(i-result).          !!
       !!                                                                      !!
       !!     neval  -    Integer                                              !!
       !!                 number of integrand evaluations.                     !!
       !!                                                                      !!
       !!     ier    -    Integer                                              !!
       !!                 ier = 0 normal and reliable termination of the       !!
       !!                      routine. it is assumed that the requested       !!
       !!                      accuracy has been achieved.                     !!
       !!                                                                      !!
       !!                 ier.gt.0 abnormal termination of the routine         !!
       !!                      the estimates for result and error are          !!
       !!                      less reliable. it is assumed that the           !!
       !!                      requested accuracy has not been achieved.       !!
       !!                                                                      !!
       !!                                                                      !!
       !!                    Error Messages                                    !!
       !!                   ----------------                                   !!
       !!                                                                      !!
       !!              ier = 1 maximum number of subdivisions allowed          !!
       !!                      has been achieved. one can allow more           !!
       !!                      subdivisions by increasing the value of         !!
       !!                      limit (and taking the according dimension       !!
       !!                      adjustments into account). however, if          !!
       !!                      this yield no improvement it is advised         !!
       !!                      to analyze the integrand in order to            !!
       !!                      determine the integration difficulaties.        !!
       !!                      if the position of a local difficulty can       !!
       !!                      be determined (i.e.singularity,                 !!
       !!                      discontinuity within the interval) one          !!
       !!                      will probably gain from splitting up the        !!
       !!                      interval at this point and calling the          !!
       !!                      integrator on the subranges. if possible,       !!
       !!                      an appropriate special-purpose integrator       !!
       !!                      should be used which is designed for            !!
       !!                      handling the type of difficulty involved.       !!
       !!                                                                      !!
       !!                  = 2 the occurrence of roundoff error is             !!
       !!                      detected, which prevents the requested          !!
       !!                      tolerance from being achieved.                  !!
       !!                                                                      !!
       !!                  = 3 extremely bad integrand behaviour occurs        !!
       !!                      at some points of the integration               !!
       !!                      interval.                                       !!
       !!                                                                      !!
       !!                  = 6 the input is invalid, because                   !!
       !!                      (epsabs.le.0 and                                !!
       !!                       epsrel.lt.max(50*rel.mach.acc.,0.5d-28))       !!
       !!                      or limit.lt.1 or lenw.lt.limit*4.               !!
       !!                      result, abserr, neval, last are set             !!
       !!                      to zero.                                        !!
       !!                      except when lenw is invalid, iwork(1),          !!
       !!                      work(limit*2+1) and work(limit*3+1) are         !!
       !!                      set to zero, work(1) is set to a and            !!
       !!                      work(limit+1) to b.                             !!
       !!                                                                      !!
       !!                                                                      !!
       !!                    Dimensioning Parameters                           !!
       !!                   -------------------------                          !!
       !!                                                                      !!
       !!                                                                      !!
       !!     limit -   Integer                                                !!
       !!               dimensioning parameter for iwork                       !!
       !!               limit determines the maximum number of subintervals    !!
       !!               in the partition of the given integration interval     !!
       !!               (a,b), limit.ge.1.                                     !!
       !!               if limit.lt.1, the routine will end with ier = 6.      !!
       !!                                                                      !!
       !!     lenw  -   Integer                                                !!
       !!               dimensioning parameter for work                        !!
       !!               lenw must be at least limit*4.                         !!
       !!               if lenw.lt.limit*4, the routine will end with          !!
       !!               ier = 6.                                               !!
       !!                                                                      !!
       !!     last  -   Integer                                                !!
       !!               on return, last equals the number of subintervals      !!
       !!               produced in the subdiviosion process, which            !!
       !!               determines the number of significant elements          !!
       !!               actually in the work arrays.                           !!
       !!                                                                      !!
       !!                                                                      !!
       !!                    Work Arrays                                       !!
       !!                   -------------                                      !!
       !!                                                                      !!
       !!                                                                      !!
       !!     iwork -   Integer                                                !!
       !!               vector of dimension at least limit, the first k        !!
       !!               elements of which contain pointers to the error        !!
       !!               estimates over the subintervals, such that             !!
       !!               work(limit*3+iwork(1)),... , work(limit*3+iwork(k))    !!
       !!               form a decreasing sequence with k = last if            !!
       !!               last.le.(limit/2+2), and k = limit+1-last otherwise    !!
       !!                                                                      !!
       !!     work  -   Real                                                   !!
       !!               vector of dimension at least lenw                      !!
       !!               on return                                              !!
       !!               work(1), ..., work(last) contain the left end          !!
       !!               points of the subintervals in the partition of         !!
       !!               (a,b),                                                 !!
       !!               work(limit+1), ..., work(limit+last) contain the       !!
       !!               right end points,                                      !!
       !!               work(limit*2+1), ..., work(limit*2+last) contain       !!
       !!               the integral approximations over the subintervals,     !!
       !!               work(limit*3+1), ..., work(limit*3+last) contain       !!
       !!               the error estimates.                                   !!
       !!                                                                      !!
       !!             References  (none)                                       !!
       !!                                                                      !!
       !!          Routines called  Dqage,Xerror                               !!
       !!                                                                      !! 
       !! end prologue  dqag                                                   !!
       !!                                                                      !!
       !!----------------------------------------------------------------------!!

           implicit none

           real(kind=mp), intent(in) :: a, b, epsabs, epsrel
           real(kind=mp), intent(out) :: abserr, result, work(*)

           integer, intent(in) :: lenw, limit, key
           integer, intent(out) :: ier, iwork(*), last, neval
           integer :: lvl, l1, l2, l3

           real(kind=mp), external :: f  ! external function defining the integrand.
       
           !!  Check validity of lenw. !!

           ier = 6
           neval = 0
           last = neval 
           result = 0.0_mp 
           abserr = result 

           if (limit .lt. 1 .or. lenw .lt. 4*limit) then
 
               if (ier .eq. 6) then
                  lvl = 1
               end if

           else
         
              !!  Prepare for call to DQAGE  !!

              l1 = limit + 1
              l2 = limit + l1
              l3 = limit + l2

              call DQAGE3(f,a,b,epsabs,epsrel,key,limit,result,abserr,neval,ier, &
                         work(1),work(l1),work(l2),work(l3),iwork,last)

              !!  Call error handler if necessary.  !!

              lvl = 0

           end if 
 
           if (ier .eq. 6) then
              lvl = 1
           else if (ier .ne. 0) then 
              call Xerror(ier)  !  Quadrature has failed.
           end if

     
          return

       end Subroutine Quad3



       Subroutine DQAGE3(f,a,b,epsabs,epsrel,key,limit,result,abserr,neval,ier, &
                        alist,blist,rlist,elist,iord,last)

       !!-------------------------------------------------------------------------!!
       !!                                                                         !!
       !!    Begin prologue  dqage                                                !!
       !!                                                                         !!
       !!  Date written   800101   (yymmdd) (original code)                       !!
       !!  Revision date  830518   (yymmdd)                                       !!
       !!  Date modifies to f90 971001 (yymmdd)                                   !!
       !!                                                                         !!
       !!  category no.  h2a1a1                                                   !!
       !!                                                                         !!
       !!  Keywords   -                                                           !!
       !!                     automatic integrator, general-purpose,              !!
       !!                     integrand examinator, globally adaptive,            !!
       !!                     Gauss-Kronrod.                                      !!
       !!                                                                         !!
       !!  Author -                                                               !!
       !!           Piessens,Robert,appl. math. & progr. div - k.u.leuven         !!
       !!           de Doncker,Elise,appl. math. & progr. div. - k.u.leuven       !!
       !!                                                                         !!
       !!  Modified to f90 -                                                      !!
       !!                   Brian Nesbitt                                         !!
       !!                   Dept. of Applied Maths & Theoretical Physics          !!
       !!                   The Queen's University of Belfast.                    !!
       !!                                                                         !!
       !!                                                                         !!
       !!                    Purpose                                              !! 
       !!                   ---------                                             !!
       !!                                                                         !!
       !!   The routine calculates an approximation result to a given             !!
       !!   definite integral i = integral of f over (a,b),                       !!
       !!   hopefully satisfying following claim for accuracy                     !!
       !!   abs(i-result)le.max(epsabs,epsrel*abs(i)).                            !!
       !!                                                                         !!
       !!                                                                         !!
       !!                                                                         !!
       !!                    Description                                          !!
       !!                   -------------                                         !!
       !!                                                                         !!        
       !!                                                                         !!
       !!    Computation of a definite integral                                   !!
       !!    standard fortran subroutine.                                         !!
       !!                                                                         !!
       !!                                                                         !!        
       !!                    Parameters                                           !!
       !!                   ------------                                          !!
       !!                                                                         !!
       !!                                                                         !! 
       !!                    On Entry                                             !!
       !!                   ----------                                            !!
       !!                                                                         !!
       !!     f      -    Real                                                    !!
       !!                 function subprogram defining the integrand              !!
       !!                 function f(x). the actual name for f needs to be        !!
       !!                 declared e x t e r n a l in the driver program.         !!
       !!                                                                         !!
       !!     a      -    Real                                                    !! 
       !!                 lower limit of integration.                             !!         
       !!                                                                         !!
       !!     b      -    Real                                                    !!
       !!                 upper limit of integration.                             !!
       !!                                                                         !!
       !!     epsabs -    Real                                                    !!
       !!                 absolute accuracy requested.                            !!
       !!                                                                         !!
       !!     epsrel -    Real                                                    !!
       !!                 relative accuracy requested.                            !! 
       !!                 if  epsabs.le.0                                         !!
       !!                 and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),            !!
       !!                 the routine will end with ier = 6.                      !!
       !!                                                                         !!
       !!     key    -    Integer                                                 !!
       !!                 key for choice of local integration rule                !! 
       !!                 a Gauss-Kronrod pair is used with                       !!
       !!                     7 - 15 points if key.lt.2,                          !!
       !!                    10 - 21 points if key = 2,                           !!
       !!                    15 - 31 points if key = 3,                           !!
       !!                    20 - 41 points if key = 4,                           !!
       !!                    25 - 51 points if key = 5,                           !!
       !!                    30 - 61 points if key.gt.5.                          !!
       !!                                                                         !!
       !!     limit  -    Integer                                                 !!
       !!                 gives an upperbound on the number of subintervals       !!
       !!                 in the partition of (a,b), limit.ge.1.                  !!
       !!                                                                         !!
       !!                    On Return                                            !!
       !!                   -----------                                           !!
       !!                                                                         !!
       !!     result -    Real                                                    !!
       !!                 approximation to the integral                           !!
       !!                                                                         !!
       !!     abserr -    Real                                                    !!  
       !!                 estimate of the modulus of the absolute error,          !!
       !!                 which should equal or exceed abs(i-result).             !!
       !!                                                                         !!
       !!     neval  -    Integer                                                 !!
       !!                 number of integrand evaluations.                        !!
       !!                                                                         !!
       !!     ier    -    Integer                                                 !!
       !!                 ier = 0     normal and reliable termination of the      !!
       !!                             routine. it is assumed that the requested   !!
       !!                             accuracy has been achieved.                 !!
       !!                                                                         !!
       !!                 ier.gt.0    abnormal termination of the routine         !!
       !!                             the estimates for result and error are      !!
       !!                             less reliable. it is assumed that the       !!
       !!                             requested accuracy has not been achieved.   !!
       !!                                                                         !!
       !!                    Error Messages                                       !!
       !!                   ----------------                                      !!
       !!                                                                         !!
       !!              ier = 1 maximum number of subdivisions allowed             !!
       !!                      has been achieved. one can allow more              !!
       !!                      subdivisions by increasing the value               !!
       !!                      of limit.                                          !!
       !!                      however, if this yields no improvement it          !!
       !!                      is rather advised to analyze the integrand         !!
       !!                      in order to determine the integration              !!
       !!                      difficulties. if the position of a local           !!
       !!                      difficulty can be determined(e.g.                  !!
       !!                      singularity, discontinuity within the              !!
       !!                      interval) one will probably gain from              !!
       !!                      splitting up the interval at this point            !!
       !!                      and calling the integrator on the                  !!
       !!                      subranges. if possible, an appropriate             !!
       !!                      special-purpose integrator should be used          !!
       !!                      which is designed for handling the type of         !!
       !!                      difficulty involved.                               !!
       !!                                                                         !!
       !!                 =  2 the occurrence of roundoff error is                !!
       !!                      detected, which prevents the requested             !!
       !!                      tolerance from being achieved.                     !!
       !!                                                                         !!
       !!                  = 3 extremely bad integrand behaviour occurs           !!
       !!                      at some points of the integration                  !!
       !!                      interval.                                          !!
       !!                                                                         !!
       !!                  = 6 the input is invalid, because                      !!
       !!                      (epsabs.le.0 and                                   !!
       !!                      epsrel.lt.max(50*rel.mach.acc.,0.5d-28),           !!
       !!                      result, abserr, neval, last, rlist(1) ,            !!
       !!                      elist(1) and iord(1) are set to zero.              !!
       !!                      alist(1) and blist(1) are set to a and b           !!
       !!                      respectively.                                      !!
       !!                                                                         !!
       !!     alist   -   Real                                                    !!
       !!                 vector of dimension at least limit, the first           !!
       !!                 last  elements of which are the left                    !!
       !!                 end points of the subintervals in the partition         !!
       !!                 of the given integration range (a,b)                    !!
       !!                                                                         !!
       !!     blist   -   Real                                                    !!
       !!                 vector of dimension at least limit, the first           !!
       !!                 last  elements of which are the right                   !!
       !!                 end points of the subintervals in the partition         !!
       !!                 of the given integration range (a,b)                    !!
       !!                                                                         !!
       !!     rlist   -   Real                                                    !!
       !!                 vector of dimension at least limit, the first           !!         
       !!                 last  elements of which are the                         !!
       !!                 integral approximations on the subintervals             !! 
       !!                                                                         !!
       !!     elist   -   Real                                                    !! 
       !!                 vector of dimension at least limit, the first           !!
       !!                 last  elements of which are the moduli of the           !!
       !!                 absolute error estimates on the subintervals            !!
       !!                                                                         !!
       !!     iord    -   Integer                                                 !!
       !!                 vector of dimension at least limit, the first k         !!
       !!                 elements of which are pointers to the                   !!
       !!                 error estimates over the subintervals,                  !!
       !!                 such that elist(iord(1)), ...,                          !!
       !!                 elist(iord(k)) form a decreasing sequence,              !!
       !!                 with k = last if last.le.(limit/2+2), and               !!
       !!                 k = limit+1-last otherwise.                             !!
       !!                                                                         !!
       !!     last    -   Integer                                                 !!         
       !!                 number of subintervals actually produced in the         !!
       !!                 subdivision process.                                    !!
       !!                                                                         !!
       !!               references  (none)                                        !!
       !!            routines called  d1mach,dqk15,dqk21,dqk31,                   !!
       !!                             dqk41,dqk51,dqk61,dqpsrt                    !!
       !!                                                                         !!
       !! end prologue  dqage                                                     !!
       !!                                                                         !!
       !!-------------------------------------------------------------------------!!

           implicit none

           real(kind=mp), intent(in) :: a, b, epsrel, epsabs
           real(kind=mp), intent(out) :: abserr, result
           real(kind=mp), intent(out) :: alist(*), blist(*), elist(*), rlist(*)

           real(kind=mp) :: area, area1, area12, area2, a1, a2,        &
                            b1, b2, defabs, defab1, defab2,            &
                            errbnd, errmax, error1, erro12, errsum,    &
                            resabs, error2

           integer, intent(in) :: limit, key 
           integer, intent(out) :: ier, iord(*), last, neval

           integer :: iroff1, iroff2, k, keyf, maxerr, nrmax

           real(kind=mp), external :: f ! external function defining the integrand.

           !!---------------------------------------------------------------------!!
           !!                                                                     !!  
           !!    list of major variables                                          !!  
           !!    -----------------------                                          !! 
           !!                                                                     !!  
           !!          alist     -   Real array                                   !! 
           !!                        list of left end points of all subintervals  !!  
           !!                        considered up to now.                        !!
           !!                                                                     !!
           !!          blist     -   Real array                                   !! 
           !!                        list of right end points of all subinterval  !!
           !!                        considered up to now.                        !!
           !!                                                                     !!
           !!          rlist(i)  -   Real                                         !!
           !!                        approximation to the integral over           !!
           !!                        (alist(i),blist(i)).                         !!
           !!                                                                     !!
           !!          elist(i)  -   Real                                         !!
           !!                        error estimate applying to rlist(i).         !!
           !!                                                                     !!
           !!          maxerr    -   Integer                                      !!
           !!                        pointer to the interval with largest.        !!
           !!                                                                     !!
           !!                                                                     !!
           !!            Error Estimate                                           !!
           !!           ----------------                                          !!
           !!                                                                     !!
           !!          errmax    -   Real                                         !! 
           !!                        elist(maxerr).                               !!
           !!                                                                     !!
           !!          area      -   Real                                         !!
           !!                        sum of the integrals over the subintervals.  !!
           !!                                                                     !!
           !!          errsum    -   Real                                         !!
           !!                        sum of the errors over the subintervals.     !!
           !!                                                                     !!
           !!          errbnd    -   Real                                         !!
           !!                        requested accuracy max(epsabs,epsrel*        !!
           !!                        abs(result)).                                !!
           !!                                                                     !!
           !!          *****1    -   variable for the left subinterval.           !!
           !!          *****2    -   variable for the right subinterval.          !!
           !!                                                                     !!
           !!          last      -   Integer                                      !!
           !!                        index for subdivision.                       !!
           !!                                                                     !!
           !!                                                                     !!
           !!              Machine Dependent Constants                            !!
           !!             -----------------------------                           !!
           !!                                                                     !!
           !!      epmach  -   is the largest relative spacing.                   !!
           !!      uflow   -   is the smallest positive magnitude.                !!
           !!                                                                     !!
           !!---------------------------------------------------------------------!!
           
           !!  First executable statement of DQAGE. !!


           !!  Test on validity of parameters.      !!

           ier = 0
           neval = ier 
           last = ier 
           result = 0.0_mp 
           abserr = result 
           alist(1) = a
           blist(1) = b
           rlist(1) = result 
           elist(1) = result 
           iord(1) = ier 

           if (epsabs .le. 0.0_mp .and. epsrel .lt.     &
               max(0.5e+2_mp*epmach,0.5e-28_mp)) then
               ier = 6
           end if

           if (ier .eq. 6) then
              call Xerror(ier)  ! Quadrature has failed.
           end if

           keyf = key

           if (key .le. 0) then 
               keyf = 1
           else if (key .ge. 7) then 
               keyf = 6
           end if

           neval = 0

           select case(keyf)
                  case(1)
                    call DQK15_3(f,a,b,result,abserr,defabs,resabs)
                  case(6)
                    call DQK61_3(f,a,b,result,abserr,defabs,resabs)
           end select 

           last = 1
           rlist(1) = result
           elist(1) = abserr
           iord(1) = 1

           !!  Test on accuracy.  !!

           errbnd = max(epsabs,epsrel*abs(result))

           if (abserr .le. 0.5e+2_mp*epmach*defabs .and. abserr .gt. errbnd) then 
               ier = 2 
           end if

           if (limit .eq. 1) then 
               ier = 1
           end if 

           if (ier .ne. 0 .or. (abserr .le. errbnd .and. abserr    &
              .ne. resabs) .or. abserr .eq. 0.0_mp) then
              
              if (keyf .ne. 1) then
                  neval = (10*keyf + 1)*(2*neval + 1)
              else if (keyf .eq. 1) then
                  neval = 30*neval + 15
              end if

              return
           end if

           !!  Initialization  !!

           errmax = abserr
           maxerr = 1
           area = result
           errsum = abserr
           nrmax = 1
           iroff1 = 0
           iroff2 = 0

           !!  Main do loop  !!
    
           do last = 2,limit

               !!  Bisect the subinterval with the largest error estimate.  !!

               a1 = alist(maxerr)
               b1 = 0.5_mp*(alist(maxerr) + blist(maxerr))
               a2 = b1
               b2 = blist(maxerr)

               !!  Use 61-point rule only at this stage.  !!

               select case(keyf)
                      case(1)
                        call DQK15_3(f,a1,b1,area1,error1,resabs,defab1)
                        call DQK15_3(f,a2,b2,area2,error2,resabs,defab2)
                      case(6)
                        call DQK61_3(f,a1,b1,area1,error1,resabs,defab1)
                        call DQK61_3(f,a2,b2,area2,error2,resabs,defab2)
               end select 


               !!  Improve the previous approximations to the integral  !!
               !!  and error and test for accuracy.                     !!

               neval = neval + 1
               area12 = area1 + area2
               erro12 = error1 + error2
               errsum = errsum + erro12 - errmax
               area = area + area12 - rlist(maxerr)

               if (defab1 .ne. error1 .and. defab2 .ne. error2) then

                   if (abs(rlist(maxerr) - area12) .le. 0.1e-4_mp*abs(area12) &
                      .and. erro12 .ge. 0.99_mp*errmax) then
                      iroff1 = iroff1 + 1
                   end if

                   if (last .gt. 10 .and. erro12 .gt. errmax) then
                      iroff2 = iroff2 + 1
                   end if

               end if

               rlist(maxerr) = area1
               rlist(last)  = area2

               errbnd = max(epsabs,epsrel*abs(area))

               if (errsum .gt. errbnd) then 

                  !!  Test for roundoff error and eventually set  !!
                  !!  error flag.                                 !!
                  
                  if (iroff1 .ge. 6 .or. iroff2 .ge. 20) then 
                      ier = 2
                  end if

                  !!  Set error flag in the case that the number  !!
                  !!  of subintervals equals limit.               !!

                  if (last .eq. limit) then 
                     ier = 1
                  end if

                  !!  Set error flag in the case of bad integrand     !!   
                  !!  behaviour at a point of the integration range.  !!

                  if (max(abs(a1),abs(b2)) .le.  & 
                     (0.1e+1_mp + 0.1e+3_mp*epmach)*(abs(a2) + 0.1e+4_mp*uflow)) then 
                      ier = 3
                  end if
        
               end if 

               !!  Apend the newly-created intervals to the list.  !!

               if (error2 .gt. error1) then 
                   alist(maxerr) = a2
                   alist(last) = a1
                   blist(last) = b1
                   rlist(maxerr) = area2
                   rlist(last) = area1
                   elist(maxerr) = error2
                   elist(last) = error1
               else
                   alist(last) = a2
                   blist(maxerr) = b1
                   blist(last) = b2
                   elist(maxerr) = error1
                   elist(last) = error2
               end if

               !!  Call subroutine DQPSRT to maintain the descending ordering  !!
               !!  in the list of error estimates and select the subinterval   !!
               !!  with the largest error estimate (to be bisected next).      !!

               call DQPSRT_3(limit,last,maxerr,errmax,elist,iord,nrmax)

               !!  Exit do loop.  !!
 
               if (ier .ne. 0 .or. errsum .le. errbnd) then 
                  exit
               end if

           end do

           !!  Compare the final result. !!

           result = 0.0_mp 

           do k = 1,last
              result = result + rlist(k)
           end do
 
           abserr = errsum

           if (keyf .ne. 1) then
               neval = (10*keyf + 1)*(2*neval + 1)
           end if

           if (keyf .eq. 1) then 
               neval = 30*neval + 15
           end if

           return

       end Subroutine DQAGE3



       Subroutine DQK15_3(f,a,b,result,abserr,resabs,resasc)

           implicit none

           real(kind=mp), intent(in) :: a, b
           real(kind=mp), intent(out) :: abserr, resabs, resasc, result            

           real(kind=mp) :: absc, centr, dhlgth, fc, fsum, fval1,     &
                            fval2, hlgth, resg, resk, reskh

           real(kind=mp) :: fv1(7), fv2(7), wg(4), wgk(8), xgk(8)

           integer :: j, jtw, jtwm1

           real(kind=mp), external :: f

           !!-------------------------------------------------------------------------!!
           !!                                                                         !!
           !!  The abscissae and weights are given for the interval (-1,1).           !!
           !!  Because of symmetry only the positive abscissae and their              !!
           !!  corresponding weights are given.                                       !!
           !!                                                                         !!
           !!  xgk    - abscissae of the 15-point Kronrod rule                        !!
           !!         xgk(2), xgk(4), ...  abscissae of the 7-point                   !!
           !!         Gauss rule                                                      !!
           !!         xgk(1), xgk(3), ...  abscissae which are optimally              !!
           !!         added to the 7-point Gauss rule                                 !!
           !!                                                                         !!
           !!  wgk    - weights of the 15-point Kronrod rule                          !!
           !!                                                                         !!
           !!  wg     - weights of the 7-point Gauss rule                             !!
           !!                                                                         !!
           !!                                                                         !!
           !!  Gauss quadrature weights and Kronrod quadrature abscissae and weights  !!
           !!  as evaluated with 80 decimal digit arithmetic by L. W. Fullerton,      !!
           !!  Bell labs, Nov. 1981.                                                  !!
           !!                                                                         !!
           !!-------------------------------------------------------------------------!!

           data wg  (1) / 0.129484966168869693270611432679082_mp  /
           data wg  (2) / 0.279705391489276667901467771423780_mp  /
           data wg  (3) / 0.381830050505118944950369775488975_mp  /
           data wg  (4) / 0.417959183673469387755102040816327_mp  /

           data xgk (1) / 0.991455371120812639206854697526329_mp  /
           data xgk (2) / 0.949107912342758524526189684047851_mp  /
           data xgk (3) / 0.864864423359769072789712788640926_mp  /
           data xgk (4) / 0.741531185599394439863864773280788_mp  /
           data xgk (5) / 0.586087235467691130294144838258730_mp  /
           data xgk (6) / 0.405845151377397166906606412076961_mp  /
           data xgk (7) / 0.207784955007898467600689403773245_mp  /
           data xgk (8) / 0.000000000000000000000000000000000_mp  /

           data wgk (1) / 0.022935322010529224963732008058970_mp  /
           data wgk (2) / 0.063092092629978553290700663189204_mp  /
           data wgk (3) / 0.104790010322250183839876322541518_mp  /
           data wgk (4) / 0.140653259715525918745189590510238_mp  /
           data wgk (5) / 0.169004726639267902826583426598550_mp  /
           data wgk (6) / 0.190350578064785409913256402421014_mp  /
           data wgk (7) / 0.204432940075298892414161999234649_mp  /
           data wgk (8) / 0.209482141084727828012999174891714_mp  /

           !!  First executable statements of DQK15.  !!

           centr = 0.5_mp*(a + b)
           hlgth = 0.5_mp*(b - a)
           dhlgth = abs(hlgth)

           !!  Compute the 15-point Kronrod approximation to the  !!
           !!  integral, and estimate the absolute error.         !!

           fc = f(centr)
           resg = fc*wg(4)
           resk = fc*wgk(8)
           resabs = abs(resk)

           do j = 1,3
              jtw = j*2
              absc = hlgth*xgk(jtw)
              fval1 = f(centr - absc)
              fval2 = f(centr + absc)
              fv1(jtw) = fval1
              fv2(jtw) = fval2
              fsum = fval1 + fval2
              resg = resg + wg(j)*fsum
              resk  = resk + wgk(jtw)*fsum
              resabs = resabs + wgk(jtw)*(abs(fval1) + abs(fval2))
           end do

           do j = 1,4
              jtwm1 = j*2 - 1
              absc = hlgth*xgk(jtwm1)
              fval1 = f(centr - absc)
              fval2 = f(centr + absc)
              fv1(jtwm1) = fval1
              fv2(jtwm1) = fval2
              fsum = fval1 + fval2
              resk = resk + wgk(jtwm1)*fsum
              resabs = resabs + wgk(jtwm1)*(abs(fval1) + abs(fval2))
           end do

           reskh = resk*0.5_mp
           resasc = wgk(8)*abs(fc - reskh)

           do j = 1,7
              resasc = resasc + wgk(j)*(abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
           end do

           result = resk*hlgth
           resabs = resabs*dhlgth
           resasc = resasc*dhlgth
           abserr = abs((resk - resg)*hlgth)

           if (resasc .ne. 0.0_mp .and. abserr .ne. 0.0_mp) then
               abserr = resasc*min(0.1e+1_mp,(0.2e+3_mp*abserr/resasc)**1.5_mp)
           end if

           if (resabs .gt. uflow/(0.5e+2_mp*epmach)) then
               abserr = max((epmach*0.5e+2_mp)*resabs,abserr)
           end if

           return

       end Subroutine DQK15_3 



       Subroutine DQK61_3(f,a,b,result,abserr,resabs,resasc)
            
           implicit none

           real(kind=mp), intent(in) :: a, b
           real(kind=mp), intent(out) :: abserr, resabs, resasc, result

           real(kind=mp) :: dabsc, centr, dhlgth,            &        
                            fc, fsum, fval1, fval2, hlgth,   &
                            resg, resk, reskh

           real(kind=mp) :: fv1(30), fv2(30), wg(15), wgk(31), xgk(31)

           real(kind=mp), external :: f

           integer :: j, jtw, jtwm1
       
           !!---------------------------------------------------------!!
           !!                                                         !!
           !!  The abscissae and weights are given for the interval   !!
           !!  (-1,+1). Because of symmetry only the positive         !!  
           !!  abscissae and their corresponding weights are given.   !! 
           !!                                                         !!
           !!                                                         !!
           !!   xgk  -  abscissae of the 51-point Kronrod rule:       !!
           !!             xgk(2), xgk(4), ...                         !!
           !!                                                         !! 
           !!           abscissae of the 25-point Gauss rule:         !! 
           !!             xgk(1), xgk(3), ...                         !!
           !!                                                         !!
           !!           abscissae which are optimally added to        !!
           !!           the 25-point Gauss rule.                      !!
           !!                                                         !!
           !!   wgk  -  weights of the 51-point Kronrod rule.         !!
           !!                                                         !!
           !!   wg   -  weights of the 25-point Gauss rule.           !!
           !!                                                         !!
           !!   Gauss quadrature weights and Kronrod Quadrature       !!
           !!   abscissae and weights as evaluated with 80 decimal    !!
           !!   arithmetic by L.W. Fullerton, Bell Labs., Nov 1981.   !!
           !!                                                         !!
           !!---------------------------------------------------------!!
  
           data wg  (  1) / 0.007968192496166605615465883474674_mp  /
           data wg  (  2) / 0.018466468311090959142302131912047_mp  /
           data wg  (  3) / 0.028784707883323369349719179611292_mp  /
           data wg  (  4) / 0.038799192569627049596801936446348_mp  /
           data wg  (  5) / 0.048402672830594052902938140422808_mp  /
           data wg  (  6) / 0.057493156217619066481721689402056_mp  /
           data wg  (  7) / 0.065974229882180495128128515115962_mp  /
           data wg  (  8) / 0.073755974737705206268243850022191_mp  /
           data wg  (  9) / 0.080755895229420215354694938460530_mp  /
           data wg  ( 10) / 0.086899787201082979802387530715126_mp  /
           data wg  ( 11) / 0.092122522237786128717632707087619_mp  /
           data wg  ( 12) / 0.096368737174644259639468626351810_mp  /
           data wg  ( 13) / 0.099593420586795267062780282103569_mp  /
           data wg  ( 14) / 0.101762389748405504596428952168554_mp  /
           data wg  ( 15) / 0.102852652893558840341285636705415_mp  /

           data xgk (  1) / 0.999484410050490637571325895705811_mp  /
           data xgk (  2) / 0.996893484074649540271630050918695_mp  /
           data xgk (  3) / 0.991630996870404594858628366109486_mp  /
           data xgk (  4) / 0.983668123279747209970032581605663_mp  /
           data xgk (  5) / 0.973116322501126268374693868423707_mp  /
           data xgk (  6) / 0.960021864968307512216871025581798_mp  /
           data xgk (  7) / 0.944374444748559979415831324037439_mp  /
           data xgk (  8) / 0.926200047429274325879324277080474_mp  /
           data xgk (  9) / 0.905573307699907798546522558925958_mp  /
           data xgk ( 10) / 0.882560535792052681543116462530226_mp  /
           data xgk ( 11) / 0.857205233546061098958658510658944_mp  /
           data xgk ( 12) / 0.829565762382768397442898119732502_mp  /
           data xgk ( 13) / 0.799727835821839083013668942322683_mp  /
           data xgk ( 14) / 0.767777432104826194917977340974503_mp  /
           data xgk ( 15) / 0.733790062453226804726171131369528_mp  /
           data xgk ( 16) / 0.697850494793315796932292388026640_mp  /
           data xgk ( 17) / 0.660061064126626961370053668149271_mp  /
           data xgk ( 18) / 0.620526182989242861140477556431189_mp  /
           data xgk ( 19) / 0.579345235826361691756024932172540_mp  /
           data xgk ( 20) / 0.536624148142019899264169793311073_mp  /
           data xgk ( 21) / 0.492480467861778574993693061207709_mp  /
           data xgk ( 22) / 0.447033769538089176780609900322854_mp  /
           data xgk ( 23) / 0.400401254830394392535476211542661_mp  /
           data xgk ( 24) / 0.352704725530878113471037207089374_mp  /
           data xgk ( 25) / 0.304073202273625077372677107199257_mp  /
           data xgk ( 26) / 0.254636926167889846439805129817805_mp  /
           data xgk ( 27) / 0.204525116682309891438957671002025_mp  /
           data xgk ( 28) / 0.153869913608583546963794672743256_mp  /
           data xgk ( 29) / 0.102806937966737030147096751318001_mp  /
           data xgk ( 30) / 0.051471842555317695833025213166723_mp  /
           data xgk ( 31) / 0.000000000000000000000000000000000_mp  /

           data wgk (  1) / 0.001389013698677007624551591226760_mp  /
           data wgk (  2) / 0.003890461127099884051267201844516_mp  /
           data wgk (  3) / 0.006630703915931292173319826369750_mp  /
           data wgk (  4) / 0.009273279659517763428441146892024_mp  /
           data wgk (  5) / 0.011823015253496341742232898853251_mp  /
           data wgk (  6) / 0.014369729507045804812451432443580_mp  /
           data wgk (  7) / 0.016920889189053272627572289420322_mp  /
           data wgk (  8) / 0.019414141193942381173408951050128_mp  /
           data wgk (  9) / 0.021828035821609192297167485738339_mp  /
           data wgk ( 10) / 0.024191162078080601365686370725232_mp  /
           data wgk ( 11) / 0.026509954882333101610601709335075_mp  /
           data wgk ( 12) / 0.028754048765041292843978785354334_mp  /
           data wgk ( 13) / 0.030907257562387762472884252943092_mp  /
           data wgk ( 14) / 0.032981447057483726031814191016854_mp  /
           data wgk ( 15) / 0.034979338028060024137499670731468_mp  /
           data wgk ( 16) / 0.036882364651821229223911065617136_mp  /
           data wgk ( 17) / 0.038678945624727592950348651532281_mp  /
           data wgk ( 18) / 0.040374538951535959111995279752468_mp  /
           data wgk ( 19) / 0.041969810215164246147147541285970_mp  /
           data wgk ( 20) / 0.043452539701356069316831728117073_mp  /
           data wgk ( 21) / 0.044814800133162663192355551616723_mp  /
           data wgk ( 22) / 0.046059238271006988116271735559374_mp  /
           data wgk ( 23) / 0.047185546569299153945261478181099_mp  /
           data wgk ( 24) / 0.048185861757087129140779492298305_mp  /
           data wgk ( 25) / 0.049055434555029778887528165367238_mp  /
           data wgk ( 26) / 0.049795683427074206357811569379942_mp  /
           data wgk ( 27) / 0.050405921402782346840893085653585_mp  /
           data wgk ( 28) / 0.050881795898749606492297473049805_mp  /
           data wgk ( 29) / 0.051221547849258772170656282604944_mp  /
           data wgk ( 30) / 0.051426128537459025933862879215781_mp  /
           data wgk ( 31) / 0.051494729429451567558340433647099_mp  /
           
           !!-------------------------------------------------------!!
           !!                                                       !!
           !!     List of major variables                           !!          
           !!     -----------------------                           !!
           !!                                                       !!
           !!     centr  -  mid point of the interval               !! 
           !!     hlgth  -  half-length of the interval             !!
           !!     dabsc  -  abscissa                                !!
           !!     fval*  -  function value                          !!
           !!     resg   -  result of the 30-point gauss rule       !!
           !!     resk   -  result of the 61-point kronrod rule     !!
           !!     reskh  -  approximation to the mean value of f    !! 
           !!               over (a,b), i.e. to i/(b-a)             !! 
           !!                                                       !! 
           !!     Machine Dependent Constants                       !!
           !!     ---------------------------                       !!
           !!                                                       !!
           !!     epmach is the largest relative spacing.           !!
           !!     uflow is the smallest positive magnitude.         !! 
           !!                                                       !!
           !!-------------------------------------------------------!!
      
           centr = 0.5_mp*(b + a)
           hlgth = 0.5_mp*(b - a)
           dhlgth = abs(hlgth)

           !!  Compute the 61-point Kronrod approximation to the  !!
           !!  integral, and estimate the absolute error.         !!

           !!  First executable statement of DQK61.  !!

           resg = 0.0_mp 
           fc = f(centr)
           resk = wgk(31)*fc
           resabs = abs(resk)

           do j = 1,15
               jtw = j*2
               dabsc = hlgth*xgk(jtw)
               fval1 = f(centr - dabsc)
               fval2 = f(centr + dabsc)
               fv1(jtw) = fval1
               fv2(jtw) = fval2
               fsum = fval1 + fval2
               resg = resg + wg(j)*fsum
               resk = resk + wgk(jtw)*fsum
               resabs = resabs + wgk(jtw)*(abs(fval1) + abs(fval2))
           end do    

           do j = 1,15
               jtwm1 = j*2 - 1
               dabsc = hlgth*xgk(jtwm1)
               fval1 = f(centr - dabsc)
               fval2 = f(centr + dabsc)
               fv1(jtwm1) = fval1
               fv2(jtwm1) = fval2
               fsum = fval1 + fval2
               resk = resk + wgk(jtwm1)*fsum                    
               resabs = resabs + wgk(jtwm1)*(abs(fval1) + abs(fval2))
           end do 

           reskh = resk*0.5_mp
           resasc = wgk(31)*abs(fc - reskh)

           do j = 1,30,2
               resasc = resasc + wgk(j)*(abs(fv1(j) - reskh) + &
                        abs(fv2(j) - reskh))
               resasc = resasc + wgk(j+1)*(abs(fv1(j+1) - reskh) + &
                        abs(fv2(j+1) - reskh))
           end do

           result = resk*hlgth
           resabs = resabs*dhlgth
           resasc = resasc*dhlgth
           abserr = abs((resk - resg)*hlgth)

           if (resasc .ne. 0.0_mp .and. abserr .ne. 0.0_mp) then 
               abserr = resasc*min(0.1e+1_mp,(0.2e+3_mp*abserr/resasc)**1.5_mp)
           end if

           if (resabs .gt. uflow/(0.5e+2_mp*epmach)) then
              abserr = max((epmach*0.5e+2_mp)*resabs,abserr)
           end if

           return

       end Subroutine DQK61_3



       Subroutine DQPSRT_3(limit,last,maxerr,ermax,elist,iord,nrmax)
   
           !!----------------------------------------------------------------------!!
           !!                                                                      !!
           !!   Begin prologue DQPSRT                                              !!       
           !!                                                                      !!
           !!   Refer to  dqage,dqagie,dqagpe,dqawse                               !!
           !!   Routines called  (none)                                            !!
           !!   Revision date  810101   (yymmdd)                                   !!
           !!                                                                      !!
           !!   keywords  -   sequential sorting.                                  !!
           !!                                                                      !!
           !!   Author  Piessens,Robert,appl. math. & progr. div. - k.u.leuven     !!
           !!   de doncker,elise,appl. math. & progr. div. - k.u.leuven.           !!
           !!                                                                      !!
           !!                                                                      !!
           !!      Purpose                                                         !!
           !!     ---------                                                        !!
           !!                                                                      !!
           !!      This routine maintains the descending ordering in the           !!
           !!      list of the local error estimated resulting from the            !!
           !!      interval subdivision process. at each call two error            !!
           !!      estimates are inserted using the sequential search              !!
           !!      method, top-down for the largest error estimate and             !!
           !!      bottom-up for the smallest error estimate.                      !!
           !!                                                                      !!
           !!                                                                      !!
           !!      Description                                                     !!
           !!     -------------                                                    !!
           !!                                                                      !!
           !!      ordering routine                                                !!
           !!      standard fortran subroutine                                     !!
           !!                                                                      !!
           !!      Parameters:                                                     !!
           !!     -------------                                                    !!
           !!                                                                      !! 
           !!                  Meaning At Output                                   !!
           !!                 -------------------                                  !!
           !!                                                                      !!
           !!         limit  -   Integer                                           !! 
           !!                    Maximum number of error estimates the list        !!
           !!                    can contain                                       !!
           !!                                                                      !! 
           !!         last   -   Integer                                           !!
           !!                    Number of error estimates currently in the list.  !!
           !!                                                                      !!
           !!         maxerr -   Integer                                           !!
           !!                    maxerr points to the nrmax-th largest error       !! 
           !!                    estimate currently in the list                    !!
           !!                                                                      !!
           !!         ermax  -   Real                                              !!
           !!                    nrmax-th largest error estimate.                  !!
           !!                    ermax = elist(maxerr)                             !!
           !!                                                                      !!
           !!         elist  -   Real                                              !!
           !!                    Vector of dimension last containing               !!
           !!                    the error estimates.                              !!
           !!                                                                      !!
           !!         iord   -   Integer                                           !!
           !!                    Vector of dimension last, the first k elements    !!
           !!                    of which contain pointers to the error            !!
           !!                    estimates, such that                              !!
           !!                    elist(iord(1)),...,  elist(iord(k))               !!
           !!                    form a decreasing sequence, with                  !!
           !!                    k = last if last.le.(limit/2+2), and              !!
           !!                    k = limit+1-last otherwise                        !!
           !!                                                                      !!
           !!         nrmax  -   Integer                                           !!
           !!                                                                      !!
           !!         maxerr =   iord(nrmax)                                       !!
           !!                                                                      !!
           !! end prologue  dqpsrt                                                 !!
           !!                                                                      !!
           !!----------------------------------------------------------------------!!

           implicit none

           real(kind=mp), intent(inout) :: ermax, elist(*)      
      
           real(kind=mp) :: errmax, errmin

           integer, intent(inout) :: iord(*), maxerr, nrmax
           integer, intent(in) :: limit

           integer :: i, ibeg, last, isucc, j, jbnd, jupbn, k, ido 

           !!  Check whether the list contains more than  !!
           !!  two error estimates.                       !!
            
           if (last .le. 2) then
               iord(1) = 1
               iord(2) = 2
               maxerr = iord(nrmax)
               ermax = elist(maxerr)
               return
           end if

           !!  This part of the routine is only executed if, due to a    !!
           !!  difficult integrand, subdivision increased the error      !!
           !!  estimate. In the normal case the insert procedure should  !!
           !!  start after the nrmax-th largest error estimate.          !!

           errmax = elist(maxerr)

           if (nrmax .ne. 1) then
               ido = nrmax - 1

               do i = 1, ido
                  isucc = iord(nrmax - 1)
                  
                  !!  Exit from do  loop.  !!
                  
                  if (errmax .le. elist(isucc)) then
                      exit
                  end if
            
                  iord(nrmax) = isucc
                  nrmax = nrmax - 1
               end do
           end if

           !!  Compute the number of elements in the list to be maintained  !!
           !!  in descending order. This number depends on the number of    !!
           !!  subdivisions still allowed.                                  !!

           jupbn = last

           if (last .gt. (limit/2 + 2)) then 
               jupbn = limit + 3 - last
           end if

           errmin = elist(last)

           !!  Insert errmax by traversing the list top-down, starting  !!
           !!  comparison from the element list elist(iord(nrmax + 1)). !!

           jbnd = jupbn - 1
           ibeg = nrmax + 1

           if (ibeg .le. jbnd) then 

               do i = ibeg, jbnd
                  isucc = iord(i)

                  !!  Exit from do loop.  !!

                  if (errmax .ge. elist(isucc)) then 
                      exit
                  end if

                  iord(i-1) = isucc
               end do

               if (errmax .ge. elist(isucc)) then 

                   !! Insert errmin by traversing the list bottom up.  !!

                   iord(i-1) = maxerr
                   k = jbnd

                   do j = i, jbnd
                      isucc = iord(k)

                      !!  Exit from do loop.  !!

                      if (errmin .lt. elist(isucc)) then 
                          exit
                      end if
           
                      iord(k+1) = isucc
                      k = k - 1
                   end do

                   if (errmin .lt. elist(isucc)) then
                       iord(k+1) = last
                   else
                       iord(i) = last
                   end if
          
               else
                   iord(jbnd) = maxerr
                   iord(jupbn) = last
               end if

           else
               iord(jbnd) = maxerr
               iord(jupbn) = last
           end if

           !!  Set maxerr and ermax.  !!

           maxerr = iord(nrmax)
           ermax = elist(maxerr)  
           
           return

       end Subroutine DQPSRT_3

        

        
      Subroutine Quad4(f,a,b,epsabs,epsrel,key,result,abserr,neval,ier, &
                        limit,lenw,last,iwork,work)
       !!----------------------------------------------------------------------!!
       !!                                                                      !!
       !!   begin prologue  Quad                                               !!
       !!                                                                      !!
       !!  date written   800101   (yymmdd) (original code)                    !!
       !!  revision date  830518   (yymmdd)                                    !!
       !!  Date modified to f90 971001 (yymmdd)                                !!
       !!                                                                      !!
       !!  category no.  h2a1a1                                                !!
       !!                                                                      !!
       !!  Keywords    -                                                       !!
       !!                   Automatic integrator, general-purpose,             !!
       !!                   integrand examinator, globally adaptive,           !!
       !!                   Gauss-Kronrod                                      !!
       !!                                                                      !!
       !!  Author -                                                            !!
       !!           Piessens,Robert,appl. math. & progr. div - k.u.leuven      !!
       !!           de doncker,elise,appl. math. & progr. div. - k.u.leuven    !!
       !!                                                                      !! 
       !!  Modified to f90 -                                                   !! 
       !!                   Brian Nesbitt                                      !!
       !!                   Dept. of Applied Maths & Theoretical Physics       !!
       !!                   The Queen's University of Belfast.                 !!
       !!                                                                      !!
       !!                                                                      !!
       !!                    Purpose                                           !!
       !!                   ---------                                          !!
       !!                                                                      !!
       !!    the routine calculates an approximation result to a given         !!
       !!    definite integral i = integral of f over (a,b),                   !!
       !!    hopefully satisfying following claim for accuracy                 !!
       !!    abs(i-result)le.max(epsabs,epsrel*abs(i)).                        !!
       !!                                                                      !!
       !!                                                                      !!
       !!                                                                      !!
       !!                    Description                                       !!
       !!                   -------------                                      !!
       !!                                                                      !!
       !!    Computation of a definite integral                                !!
       !!    standard fortran subroutine                                       !!
       !!    Fortran90 precision version                                       !!
       !!                                                                      !!
       !!     f      -    Real                                                 !!
       !!                 Function subprogam defining the integrand            !! 
       !!                 function f(x). the actual name for f needs to be     !! 
       !!                 declared e x t e r n a l in the driver program.      !!
       !!                                                                      !!
       !!     a      -    Real                                                 !! 
       !!                 Lower limit of integration.                          !!
       !!                                                                      !!
       !!     b      -    Real                                                 !! 
       !!                 Upper limit of integration.                          !!
       !!                                                                      !!
       !!     epsabs -    Real                                                 !!
       !!                 Absolute accuracy requested.                         !!
       !!                                                                      !!
       !!     epsrel -    Real                                                 !! 
       !!                 Relative accuracy requested                          !!
       !!                 if  epsabs.le.0                                      !!
       !!                 and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),         !!
       !!                 the routine will end with ier = 6.                   !!
       !!                                                                      !!
       !!     key    -    Integer                                              !!
       !!                 Key for choice of local integration rule             !!
       !!                 a gauss-kronrod pair is used with                    !!
       !!                   7 - 15 points if key.lt.2,                         !!
       !!                  10 - 21 points if key = 2,                          !!
       !!                  15 - 31 points if key = 3,                          !!
       !!                  20 - 41 points if key = 4,                          !!
       !!                  25 - 51 points if key = 5,                          !!
       !!                  30 - 61 points if key .gt. 5.                       !!
       !!                                                                      !!
       !!                    On return                                         !!
       !!                   -----------                                        !! 
       !!                                                                      !!
       !!     result -    Real                                                 !!
       !!                 approximation to the integral.                       !!
       !!                                                                      !!
       !!     abserr -    Real                                                 !!
       !!                 estimate of the modulus of the absolute error,       !!
       !!                 which should equal or exceed abs(i-result).          !!
       !!                                                                      !!
       !!     neval  -    Integer                                              !!
       !!                 number of integrand evaluations.                     !!
       !!                                                                      !!
       !!     ier    -    Integer                                              !!
       !!                 ier = 0 normal and reliable termination of the       !!
       !!                      routine. it is assumed that the requested       !!
       !!                      accuracy has been achieved.                     !!
       !!                                                                      !!
       !!                 ier.gt.0 abnormal termination of the routine         !!
       !!                      the estimates for result and error are          !!
       !!                      less reliable. it is assumed that the           !!
       !!                      requested accuracy has not been achieved.       !!
       !!                                                                      !!
       !!                                                                      !!
       !!                    Error Messages                                    !!
       !!                   ----------------                                   !!
       !!                                                                      !!
       !!              ier = 1 maximum number of subdivisions allowed          !!
       !!                      has been achieved. one can allow more           !!
       !!                      subdivisions by increasing the value of         !!
       !!                      limit (and taking the according dimension       !!
       !!                      adjustments into account). however, if          !!
       !!                      this yield no improvement it is advised         !!
       !!                      to analyze the integrand in order to            !!
       !!                      determine the integration difficulaties.        !!
       !!                      if the position of a local difficulty can       !!
       !!                      be determined (i.e.singularity,                 !!
       !!                      discontinuity within the interval) one          !!
       !!                      will probably gain from splitting up the        !!
       !!                      interval at this point and calling the          !!
       !!                      integrator on the subranges. if possible,       !!
       !!                      an appropriate special-purpose integrator       !!
       !!                      should be used which is designed for            !!
       !!                      handling the type of difficulty involved.       !!
       !!                                                                      !!
       !!                  = 2 the occurrence of roundoff error is             !!
       !!                      detected, which prevents the requested          !!
       !!                      tolerance from being achieved.                  !!
       !!                                                                      !!
       !!                  = 3 extremely bad integrand behaviour occurs        !!
       !!                      at some points of the integration               !!
       !!                      interval.                                       !!
       !!                                                                      !!
       !!                  = 6 the input is invalid, because                   !!
       !!                      (epsabs.le.0 and                                !!
       !!                       epsrel.lt.max(50*rel.mach.acc.,0.5d-28))       !!
       !!                      or limit.lt.1 or lenw.lt.limit*4.               !!
       !!                      result, abserr, neval, last are set             !!
       !!                      to zero.                                        !!
       !!                      except when lenw is invalid, iwork(1),          !!
       !!                      work(limit*2+1) and work(limit*3+1) are         !!
       !!                      set to zero, work(1) is set to a and            !!
       !!                      work(limit+1) to b.                             !!
       !!                                                                      !!
       !!                                                                      !!
       !!                    Dimensioning Parameters                           !!
       !!                   -------------------------                          !!
       !!                                                                      !!
       !!                                                                      !!
       !!     limit -   Integer                                                !!
       !!               dimensioning parameter for iwork                       !!
       !!               limit determines the maximum number of subintervals    !!
       !!               in the partition of the given integration interval     !!
       !!               (a,b), limit.ge.1.                                     !!
       !!               if limit.lt.1, the routine will end with ier = 6.      !!
       !!                                                                      !!
       !!     lenw  -   Integer                                                !!
       !!               dimensioning parameter for work                        !!
       !!               lenw must be at least limit*4.                         !!
       !!               if lenw.lt.limit*4, the routine will end with          !!
       !!               ier = 6.                                               !!
       !!                                                                      !!
       !!     last  -   Integer                                                !!
       !!               on return, last equals the number of subintervals      !!
       !!               produced in the subdiviosion process, which            !!
       !!               determines the number of significant elements          !!
       !!               actually in the work arrays.                           !!
       !!                                                                      !!
       !!                                                                      !!
       !!                    Work Arrays                                       !!
       !!                   -------------                                      !!
       !!                                                                      !!
       !!                                                                      !!
       !!     iwork -   Integer                                                !!
       !!               vector of dimension at least limit, the first k        !!
       !!               elements of which contain pointers to the error        !!
       !!               estimates over the subintervals, such that             !!
       !!               work(limit*3+iwork(1)),... , work(limit*3+iwork(k))    !!
       !!               form a decreasing sequence with k = last if            !!
       !!               last.le.(limit/2+2), and k = limit+1-last otherwise    !!
       !!                                                                      !!
       !!     work  -   Real                                                   !!
       !!               vector of dimension at least lenw                      !!
       !!               on return                                              !!
       !!               work(1), ..., work(last) contain the left end          !!
       !!               points of the subintervals in the partition of         !!
       !!               (a,b),                                                 !!
       !!               work(limit+1), ..., work(limit+last) contain the       !!
       !!               right end points,                                      !!
       !!               work(limit*2+1), ..., work(limit*2+last) contain       !!
       !!               the integral approximations over the subintervals,     !!
       !!               work(limit*3+1), ..., work(limit*3+last) contain       !!
       !!               the error estimates.                                   !!
       !!                                                                      !!
       !!             References  (none)                                       !!
       !!                                                                      !!
       !!          Routines called  Dqage,Xerror                               !!
       !!                                                                      !! 
       !! end prologue  dqag                                                   !!
       !!                                                                      !!
       !!----------------------------------------------------------------------!!

           implicit none

           real(kind=mp), intent(in) :: a, b, epsabs, epsrel
           real(kind=mp), intent(out) :: abserr, result, work(*)

           integer, intent(in) :: lenw, limit, key
           integer, intent(out) :: ier, iwork(*), last, neval
           integer :: lvl, l1, l2, l3

           real(kind=mp), external :: f  ! external function defining the integrand.
       
           !!  Check validity of lenw. !!

           ier = 6
           neval = 0
           last = neval 
           result = 0.0_mp 
           abserr = result 

           if (limit .lt. 1 .or. lenw .lt. 4*limit) then
 
               if (ier .eq. 6) then
                  lvl = 1
               end if

           else
         
              !!  Prepare for call to DQAGE  !!

              l1 = limit + 1
              l2 = limit + l1
              l3 = limit + l2

              call DQAGE4(f,a,b,epsabs,epsrel,key,limit,result,abserr,neval,ier, &
                         work(1),work(l1),work(l2),work(l3),iwork,last)

              !!  Call error handler if necessary.  !!

              lvl = 0

           end if 
 
           if (ier .eq. 6) then
              lvl = 1
           else if (ier .ne. 0) then 
              call Xerror(ier)  !  Quadrature has failed.
           end if

     
          return

       end Subroutine Quad4



       Subroutine DQAGE4(f,a,b,epsabs,epsrel,key,limit,result,abserr,neval,ier, &
                        alist,blist,rlist,elist,iord,last)

       !!-------------------------------------------------------------------------!!
       !!                                                                         !!
       !!    Begin prologue  dqage                                                !!
       !!                                                                         !!
       !!  Date written   800101   (yymmdd) (original code)                       !!
       !!  Revision date  830518   (yymmdd)                                       !!
       !!  Date modifies to f90 971001 (yymmdd)                                   !!
       !!                                                                         !!
       !!  category no.  h2a1a1                                                   !!
       !!                                                                         !!
       !!  Keywords   -                                                           !!
       !!                     automatic integrator, general-purpose,              !!
       !!                     integrand examinator, globally adaptive,            !!
       !!                     Gauss-Kronrod.                                      !!
       !!                                                                         !!
       !!  Author -                                                               !!
       !!           Piessens,Robert,appl. math. & progr. div - k.u.leuven         !!
       !!           de Doncker,Elise,appl. math. & progr. div. - k.u.leuven       !!
       !!                                                                         !!
       !!  Modified to f90 -                                                      !!
       !!                   Brian Nesbitt                                         !!
       !!                   Dept. of Applied Maths & Theoretical Physics          !!
       !!                   The Queen's University of Belfast.                    !!
       !!                                                                         !!
       !!                                                                         !!
       !!                    Purpose                                              !! 
       !!                   ---------                                             !!
       !!                                                                         !!
       !!   The routine calculates an approximation result to a given             !!
       !!   definite integral i = integral of f over (a,b),                       !!
       !!   hopefully satisfying following claim for accuracy                     !!
       !!   abs(i-result)le.max(epsabs,epsrel*abs(i)).                            !!
       !!                                                                         !!
       !!                                                                         !!
       !!                                                                         !!
       !!                    Description                                          !!
       !!                   -------------                                         !!
       !!                                                                         !!        
       !!                                                                         !!
       !!    Computation of a definite integral                                   !!
       !!    standard fortran subroutine.                                         !!
       !!                                                                         !!
       !!                                                                         !!        
       !!                    Parameters                                           !!
       !!                   ------------                                          !!
       !!                                                                         !!
       !!                                                                         !! 
       !!                    On Entry                                             !!
       !!                   ----------                                            !!
       !!                                                                         !!
       !!     f      -    Real                                                    !!
       !!                 function subprogram defining the integrand              !!
       !!                 function f(x). the actual name for f needs to be        !!
       !!                 declared e x t e r n a l in the driver program.         !!
       !!                                                                         !!
       !!     a      -    Real                                                    !! 
       !!                 lower limit of integration.                             !!         
       !!                                                                         !!
       !!     b      -    Real                                                    !!
       !!                 upper limit of integration.                             !!
       !!                                                                         !!
       !!     epsabs -    Real                                                    !!
       !!                 absolute accuracy requested.                            !!
       !!                                                                         !!
       !!     epsrel -    Real                                                    !!
       !!                 relative accuracy requested.                            !! 
       !!                 if  epsabs.le.0                                         !!
       !!                 and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),            !!
       !!                 the routine will end with ier = 6.                      !!
       !!                                                                         !!
       !!     key    -    Integer                                                 !!
       !!                 key for choice of local integration rule                !! 
       !!                 a Gauss-Kronrod pair is used with                       !!
       !!                     7 - 15 points if key.lt.2,                          !!
       !!                    10 - 21 points if key = 2,                           !!
       !!                    15 - 31 points if key = 3,                           !!
       !!                    20 - 41 points if key = 4,                           !!
       !!                    25 - 51 points if key = 5,                           !!
       !!                    30 - 61 points if key.gt.5.                          !!
       !!                                                                         !!
       !!     limit  -    Integer                                                 !!
       !!                 gives an upperbound on the number of subintervals       !!
       !!                 in the partition of (a,b), limit.ge.1.                  !!
       !!                                                                         !!
       !!                    On Return                                            !!
       !!                   -----------                                           !!
       !!                                                                         !!
       !!     result -    Real                                                    !!
       !!                 approximation to the integral                           !!
       !!                                                                         !!
       !!     abserr -    Real                                                    !!  
       !!                 estimate of the modulus of the absolute error,          !!
       !!                 which should equal or exceed abs(i-result).             !!
       !!                                                                         !!
       !!     neval  -    Integer                                                 !!
       !!                 number of integrand evaluations.                        !!
       !!                                                                         !!
       !!     ier    -    Integer                                                 !!
       !!                 ier = 0     normal and reliable termination of the      !!
       !!                             routine. it is assumed that the requested   !!
       !!                             accuracy has been achieved.                 !!
       !!                                                                         !!
       !!                 ier.gt.0    abnormal termination of the routine         !!
       !!                             the estimates for result and error are      !!
       !!                             less reliable. it is assumed that the       !!
       !!                             requested accuracy has not been achieved.   !!
       !!                                                                         !!
       !!                    Error Messages                                       !!
       !!                   ----------------                                      !!
       !!                                                                         !!
       !!              ier = 1 maximum number of subdivisions allowed             !!
       !!                      has been achieved. one can allow more              !!
       !!                      subdivisions by increasing the value               !!
       !!                      of limit.                                          !!
       !!                      however, if this yields no improvement it          !!
       !!                      is rather advised to analyze the integrand         !!
       !!                      in order to determine the integration              !!
       !!                      difficulties. if the position of a local           !!
       !!                      difficulty can be determined(e.g.                  !!
       !!                      singularity, discontinuity within the              !!
       !!                      interval) one will probably gain from              !!
       !!                      splitting up the interval at this point            !!
       !!                      and calling the integrator on the                  !!
       !!                      subranges. if possible, an appropriate             !!
       !!                      special-purpose integrator should be used          !!
       !!                      which is designed for handling the type of         !!
       !!                      difficulty involved.                               !!
       !!                                                                         !!
       !!                 =  2 the occurrence of roundoff error is                !!
       !!                      detected, which prevents the requested             !!
       !!                      tolerance from being achieved.                     !!
       !!                                                                         !!
       !!                  = 3 extremely bad integrand behaviour occurs           !!
       !!                      at some points of the integration                  !!
       !!                      interval.                                          !!
       !!                                                                         !!
       !!                  = 6 the input is invalid, because                      !!
       !!                      (epsabs.le.0 and                                   !!
       !!                      epsrel.lt.max(50*rel.mach.acc.,0.5d-28),           !!
       !!                      result, abserr, neval, last, rlist(1) ,            !!
       !!                      elist(1) and iord(1) are set to zero.              !!
       !!                      alist(1) and blist(1) are set to a and b           !!
       !!                      respectively.                                      !!
       !!                                                                         !!
       !!     alist   -   Real                                                    !!
       !!                 vector of dimension at least limit, the first           !!
       !!                 last  elements of which are the left                    !!
       !!                 end points of the subintervals in the partition         !!
       !!                 of the given integration range (a,b)                    !!
       !!                                                                         !!
       !!     blist   -   Real                                                    !!
       !!                 vector of dimension at least limit, the first           !!
       !!                 last  elements of which are the right                   !!
       !!                 end points of the subintervals in the partition         !!
       !!                 of the given integration range (a,b)                    !!
       !!                                                                         !!
       !!     rlist   -   Real                                                    !!
       !!                 vector of dimension at least limit, the first           !!         
       !!                 last  elements of which are the                         !!
       !!                 integral approximations on the subintervals             !! 
       !!                                                                         !!
       !!     elist   -   Real                                                    !! 
       !!                 vector of dimension at least limit, the first           !!
       !!                 last  elements of which are the moduli of the           !!
       !!                 absolute error estimates on the subintervals            !!
       !!                                                                         !!
       !!     iord    -   Integer                                                 !!
       !!                 vector of dimension at least limit, the first k         !!
       !!                 elements of which are pointers to the                   !!
       !!                 error estimates over the subintervals,                  !!
       !!                 such that elist(iord(1)), ...,                          !!
       !!                 elist(iord(k)) form a decreasing sequence,              !!
       !!                 with k = last if last.le.(limit/2+2), and               !!
       !!                 k = limit+1-last otherwise.                             !!
       !!                                                                         !!
       !!     last    -   Integer                                                 !!         
       !!                 number of subintervals actually produced in the         !!
       !!                 subdivision process.                                    !!
       !!                                                                         !!
       !!               references  (none)                                        !!
       !!            routines called  d1mach,dqk15,dqk21,dqk31,                   !!
       !!                             dqk41,dqk51,dqk61,dqpsrt                    !!
       !!                                                                         !!
       !! end prologue  dqage                                                     !!
       !!                                                                         !!
       !!-------------------------------------------------------------------------!!

           implicit none

           real(kind=mp), intent(in) :: a, b, epsrel, epsabs
           real(kind=mp), intent(out) :: abserr, result
           real(kind=mp), intent(out) :: alist(*), blist(*), elist(*), rlist(*)

           real(kind=mp) :: area, area1, area12, area2, a1, a2,        &
                            b1, b2, defabs, defab1, defab2,            &
                            errbnd, errmax, error1, erro12, errsum,    &
                            resabs, error2

           integer, intent(in) :: limit, key 
           integer, intent(out) :: ier, iord(*), last, neval

           integer :: iroff1, iroff2, k, keyf, maxerr, nrmax

           real(kind=mp), external :: f ! external function defining the integrand.

           !!---------------------------------------------------------------------!!
           !!                                                                     !!  
           !!    list of major variables                                          !!  
           !!    -----------------------                                          !! 
           !!                                                                     !!  
           !!          alist     -   Real array                                   !! 
           !!                        list of left end points of all subintervals  !!  
           !!                        considered up to now.                        !!
           !!                                                                     !!
           !!          blist     -   Real array                                   !! 
           !!                        list of right end points of all subinterval  !!
           !!                        considered up to now.                        !!
           !!                                                                     !!
           !!          rlist(i)  -   Real                                         !!
           !!                        approximation to the integral over           !!
           !!                        (alist(i),blist(i)).                         !!
           !!                                                                     !!
           !!          elist(i)  -   Real                                         !!
           !!                        error estimate applying to rlist(i).         !!
           !!                                                                     !!
           !!          maxerr    -   Integer                                      !!
           !!                        pointer to the interval with largest.        !!
           !!                                                                     !!
           !!                                                                     !!
           !!            Error Estimate                                           !!
           !!           ----------------                                          !!
           !!                                                                     !!
           !!          errmax    -   Real                                         !! 
           !!                        elist(maxerr).                               !!
           !!                                                                     !!
           !!          area      -   Real                                         !!
           !!                        sum of the integrals over the subintervals.  !!
           !!                                                                     !!
           !!          errsum    -   Real                                         !!
           !!                        sum of the errors over the subintervals.     !!
           !!                                                                     !!
           !!          errbnd    -   Real                                         !!
           !!                        requested accuracy max(epsabs,epsrel*        !!
           !!                        abs(result)).                                !!
           !!                                                                     !!
           !!          *****1    -   variable for the left subinterval.           !!
           !!          *****2    -   variable for the right subinterval.          !!
           !!                                                                     !!
           !!          last      -   Integer                                      !!
           !!                        index for subdivision.                       !!
           !!                                                                     !!
           !!                                                                     !!
           !!              Machine Dependent Constants                            !!
           !!             -----------------------------                           !!
           !!                                                                     !!
           !!      epmach  -   is the largest relative spacing.                   !!
           !!      uflow   -   is the smallest positive magnitude.                !!
           !!                                                                     !!
           !!---------------------------------------------------------------------!!
           
           !!  First executable statement of DQAGE. !!


           !!  Test on validity of parameters.      !!

           ier = 0
           neval = ier 
           last = ier 
           result = 0.0_mp 
           abserr = result 
           alist(1) = a
           blist(1) = b
           rlist(1) = result 
           elist(1) = result 
           iord(1) = ier 

           if (epsabs .le. 0.0_mp .and. epsrel .lt.     &
               max(0.5e+2_mp*epmach,0.5e-28_mp)) then
               ier = 6
           end if

           if (ier .eq. 6) then
              call Xerror(ier)  ! Quadrature has failed.
           end if

           keyf = key

           if (key .le. 0) then 
               keyf = 1
           else if (key .ge. 7) then 
               keyf = 6
           end if

           neval = 0

           select case(keyf)
                  case(1)
                    call DQK15_4(f,a,b,result,abserr,defabs,resabs)
                  case(6)
                    call DQK61_4(f,a,b,result,abserr,defabs,resabs)
           end select 

           last = 1
           rlist(1) = result
           elist(1) = abserr
           iord(1) = 1

           !!  Test on accuracy.  !!

           errbnd = max(epsabs,epsrel*abs(result))

           if (abserr .le. 0.5e+2_mp*epmach*defabs .and. abserr .gt. errbnd) then 
               ier = 2 
           end if

           if (limit .eq. 1) then 
               ier = 1
           end if 

           if (ier .ne. 0 .or. (abserr .le. errbnd .and. abserr    &
              .ne. resabs) .or. abserr .eq. 0.0_mp) then
              
              if (keyf .ne. 1) then
                  neval = (10*keyf + 1)*(2*neval + 1)
              else if (keyf .eq. 1) then
                  neval = 30*neval + 15
              end if

              return
           end if

           !!  Initialization  !!

           errmax = abserr
           maxerr = 1
           area = result
           errsum = abserr
           nrmax = 1
           iroff1 = 0
           iroff2 = 0

           !!  Main do loop  !!
    
           do last = 2,limit

               !!  Bisect the subinterval with the largest error estimate.  !!

               a1 = alist(maxerr)
               b1 = 0.5_mp*(alist(maxerr) + blist(maxerr))
               a2 = b1
               b2 = blist(maxerr)

               !!  Use 61-point rule only at this stage.  !!

               select case(keyf)
                      case(1)
                        call DQK15_4(f,a1,b1,area1,error1,resabs,defab1)
                        call DQK15_4(f,a2,b2,area2,error2,resabs,defab2)
                      case(6)
                        call DQK61_4(f,a1,b1,area1,error1,resabs,defab1)
                        call DQK61_4(f,a2,b2,area2,error2,resabs,defab2)
               end select 


               !!  Improve the previous approximations to the integral  !!
               !!  and error and test for accuracy.                     !!

               neval = neval + 1
               area12 = area1 + area2
               erro12 = error1 + error2
               errsum = errsum + erro12 - errmax
               area = area + area12 - rlist(maxerr)

               if (defab1 .ne. error1 .and. defab2 .ne. error2) then

                   if (abs(rlist(maxerr) - area12) .le. 0.1e-4_mp*abs(area12) &
                      .and. erro12 .ge. 0.99_mp*errmax) then
                      iroff1 = iroff1 + 1
                   end if

                   if (last .gt. 10 .and. erro12 .gt. errmax) then
                      iroff2 = iroff2 + 1
                   end if

               end if

               rlist(maxerr) = area1
               rlist(last)  = area2

               errbnd = max(epsabs,epsrel*abs(area))

               if (errsum .gt. errbnd) then 

                  !!  Test for roundoff error and eventually set  !!
                  !!  error flag.                                 !!
                  
                  if (iroff1 .ge. 6 .or. iroff2 .ge. 20) then 
                      ier = 2
                  end if

                  !!  Set error flag in the case that the number  !!
                  !!  of subintervals equals limit.               !!

                  if (last .eq. limit) then 
                     ier = 1
                  end if

                  !!  Set error flag in the case of bad integrand     !!   
                  !!  behaviour at a point of the integration range.  !!

                  if (max(abs(a1),abs(b2)) .le.  & 
                     (0.1e+1_mp + 0.1e+3_mp*epmach)*(abs(a2) + 0.1e+4_mp*uflow)) then 
                      ier = 3
                  end if
        
               end if 

               !!  Apend the newly-created intervals to the list.  !!

               if (error2 .gt. error1) then 
                   alist(maxerr) = a2
                   alist(last) = a1
                   blist(last) = b1
                   rlist(maxerr) = area2
                   rlist(last) = area1
                   elist(maxerr) = error2
                   elist(last) = error1
               else
                   alist(last) = a2
                   blist(maxerr) = b1
                   blist(last) = b2
                   elist(maxerr) = error1
                   elist(last) = error2
               end if

               !!  Call subroutine DQPSRT to maintain the descending ordering  !!
               !!  in the list of error estimates and select the subinterval   !!
               !!  with the largest error estimate (to be bisected next).      !!

               call DQPSRT_4(limit,last,maxerr,errmax,elist,iord,nrmax)

               !!  Exit do loop.  !!
 
               if (ier .ne. 0 .or. errsum .le. errbnd) then 
                  exit
               end if

           end do

           !!  Compare the final result. !!

           result = 0.0_mp 

           do k = 1,last
              result = result + rlist(k)
           end do
 
           abserr = errsum

           if (keyf .ne. 1) then
               neval = (10*keyf + 1)*(2*neval + 1)
           end if

           if (keyf .eq. 1) then 
               neval = 30*neval + 15
           end if

           return

       end Subroutine DQAGE4



       Subroutine DQK15_4(f,a,b,result,abserr,resabs,resasc)

           implicit none

           real(kind=mp), intent(in) :: a, b
           real(kind=mp), intent(out) :: abserr, resabs, resasc, result            

           real(kind=mp) :: absc, centr, dhlgth, fc, fsum, fval1,     &
                            fval2, hlgth, resg, resk, reskh

           real(kind=mp) :: fv1(7), fv2(7), wg(4), wgk(8), xgk(8)

           integer :: j, jtw, jtwm1

           real(kind=mp), external :: f

           !!-------------------------------------------------------------------------!!
           !!                                                                         !!
           !!  The abscissae and weights are given for the interval (-1,1).           !!
           !!  Because of symmetry only the positive abscissae and their              !!
           !!  corresponding weights are given.                                       !!
           !!                                                                         !!
           !!  xgk    - abscissae of the 15-point Kronrod rule                        !!
           !!         xgk(2), xgk(4), ...  abscissae of the 7-point                   !!
           !!         Gauss rule                                                      !!
           !!         xgk(1), xgk(3), ...  abscissae which are optimally              !!
           !!         added to the 7-point Gauss rule                                 !!
           !!                                                                         !!
           !!  wgk    - weights of the 15-point Kronrod rule                          !!
           !!                                                                         !!
           !!  wg     - weights of the 7-point Gauss rule                             !!
           !!                                                                         !!
           !!                                                                         !!
           !!  Gauss quadrature weights and Kronrod quadrature abscissae and weights  !!
           !!  as evaluated with 80 decimal digit arithmetic by L. W. Fullerton,      !!
           !!  Bell labs, Nov. 1981.                                                  !!
           !!                                                                         !!
           !!-------------------------------------------------------------------------!!

           data wg  (1) / 0.129484966168869693270611432679082_mp  /
           data wg  (2) / 0.279705391489276667901467771423780_mp  /
           data wg  (3) / 0.381830050505118944950369775488975_mp  /
           data wg  (4) / 0.417959183673469387755102040816327_mp  /

           data xgk (1) / 0.991455371120812639206854697526329_mp  /
           data xgk (2) / 0.949107912342758524526189684047851_mp  /
           data xgk (3) / 0.864864423359769072789712788640926_mp  /
           data xgk (4) / 0.741531185599394439863864773280788_mp  /
           data xgk (5) / 0.586087235467691130294144838258730_mp  /
           data xgk (6) / 0.405845151377397166906606412076961_mp  /
           data xgk (7) / 0.207784955007898467600689403773245_mp  /
           data xgk (8) / 0.000000000000000000000000000000000_mp  /

           data wgk (1) / 0.022935322010529224963732008058970_mp  /
           data wgk (2) / 0.063092092629978553290700663189204_mp  /
           data wgk (3) / 0.104790010322250183839876322541518_mp  /
           data wgk (4) / 0.140653259715525918745189590510238_mp  /
           data wgk (5) / 0.169004726639267902826583426598550_mp  /
           data wgk (6) / 0.190350578064785409913256402421014_mp  /
           data wgk (7) / 0.204432940075298892414161999234649_mp  /
           data wgk (8) / 0.209482141084727828012999174891714_mp  /

           !!  First executable statements of DQK15.  !!

           centr = 0.5_mp*(a + b)
           hlgth = 0.5_mp*(b - a)
           dhlgth = abs(hlgth)

           !!  Compute the 15-point Kronrod approximation to the  !!
           !!  integral, and estimate the absolute error.         !!

           fc = f(centr)
           resg = fc*wg(4)
           resk = fc*wgk(8)
           resabs = abs(resk)

           do j = 1,3
              jtw = j*2
              absc = hlgth*xgk(jtw)
              fval1 = f(centr - absc)
              fval2 = f(centr + absc)
              fv1(jtw) = fval1
              fv2(jtw) = fval2
              fsum = fval1 + fval2
              resg = resg + wg(j)*fsum
              resk  = resk + wgk(jtw)*fsum
              resabs = resabs + wgk(jtw)*(abs(fval1) + abs(fval2))
           end do

           do j = 1,4
              jtwm1 = j*2 - 1
              absc = hlgth*xgk(jtwm1)
              fval1 = f(centr - absc)
              fval2 = f(centr + absc)
              fv1(jtwm1) = fval1
              fv2(jtwm1) = fval2
              fsum = fval1 + fval2
              resk = resk + wgk(jtwm1)*fsum
              resabs = resabs + wgk(jtwm1)*(abs(fval1) + abs(fval2))
           end do

           reskh = resk*0.5_mp
           resasc = wgk(8)*abs(fc - reskh)

           do j = 1,7
              resasc = resasc + wgk(j)*(abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
           end do

           result = resk*hlgth
           resabs = resabs*dhlgth
           resasc = resasc*dhlgth
           abserr = abs((resk - resg)*hlgth)

           if (resasc .ne. 0.0_mp .and. abserr .ne. 0.0_mp) then
               abserr = resasc*min(0.1e+1_mp,(0.2e+3_mp*abserr/resasc)**1.5_mp)
           end if

           if (resabs .gt. uflow/(0.5e+2_mp*epmach)) then
               abserr = max((epmach*0.5e+2_mp)*resabs,abserr)
           end if

           return

       end Subroutine DQK15_4 



       Subroutine DQK61_4(f,a,b,result,abserr,resabs,resasc)
            
           implicit none

           real(kind=mp), intent(in) :: a, b
           real(kind=mp), intent(out) :: abserr, resabs, resasc, result

           real(kind=mp) :: dabsc, centr, dhlgth,            &        
                            fc, fsum, fval1, fval2, hlgth,   &
                            resg, resk, reskh

           real(kind=mp) :: fv1(30), fv2(30), wg(15), wgk(31), xgk(31)

           real(kind=mp), external :: f

           integer :: j, jtw, jtwm1
       
           !!---------------------------------------------------------!!
           !!                                                         !!
           !!  The abscissae and weights are given for the interval   !!
           !!  (-1,+1). Because of symmetry only the positive         !!  
           !!  abscissae and their corresponding weights are given.   !! 
           !!                                                         !!
           !!                                                         !!
           !!   xgk  -  abscissae of the 51-point Kronrod rule:       !!
           !!             xgk(2), xgk(4), ...                         !!
           !!                                                         !! 
           !!           abscissae of the 25-point Gauss rule:         !! 
           !!             xgk(1), xgk(3), ...                         !!
           !!                                                         !!
           !!           abscissae which are optimally added to        !!
           !!           the 25-point Gauss rule.                      !!
           !!                                                         !!
           !!   wgk  -  weights of the 51-point Kronrod rule.         !!
           !!                                                         !!
           !!   wg   -  weights of the 25-point Gauss rule.           !!
           !!                                                         !!
           !!   Gauss quadrature weights and Kronrod Quadrature       !!
           !!   abscissae and weights as evaluated with 80 decimal    !!
           !!   arithmetic by L.W. Fullerton, Bell Labs., Nov 1981.   !!
           !!                                                         !!
           !!---------------------------------------------------------!!
  
           data wg  (  1) / 0.007968192496166605615465883474674_mp  /
           data wg  (  2) / 0.018466468311090959142302131912047_mp  /
           data wg  (  3) / 0.028784707883323369349719179611292_mp  /
           data wg  (  4) / 0.038799192569627049596801936446348_mp  /
           data wg  (  5) / 0.048402672830594052902938140422808_mp  /
           data wg  (  6) / 0.057493156217619066481721689402056_mp  /
           data wg  (  7) / 0.065974229882180495128128515115962_mp  /
           data wg  (  8) / 0.073755974737705206268243850022191_mp  /
           data wg  (  9) / 0.080755895229420215354694938460530_mp  /
           data wg  ( 10) / 0.086899787201082979802387530715126_mp  /
           data wg  ( 11) / 0.092122522237786128717632707087619_mp  /
           data wg  ( 12) / 0.096368737174644259639468626351810_mp  /
           data wg  ( 13) / 0.099593420586795267062780282103569_mp  /
           data wg  ( 14) / 0.101762389748405504596428952168554_mp  /
           data wg  ( 15) / 0.102852652893558840341285636705415_mp  /

           data xgk (  1) / 0.999484410050490637571325895705811_mp  /
           data xgk (  2) / 0.996893484074649540271630050918695_mp  /
           data xgk (  3) / 0.991630996870404594858628366109486_mp  /
           data xgk (  4) / 0.983668123279747209970032581605663_mp  /
           data xgk (  5) / 0.973116322501126268374693868423707_mp  /
           data xgk (  6) / 0.960021864968307512216871025581798_mp  /
           data xgk (  7) / 0.944374444748559979415831324037439_mp  /
           data xgk (  8) / 0.926200047429274325879324277080474_mp  /
           data xgk (  9) / 0.905573307699907798546522558925958_mp  /
           data xgk ( 10) / 0.882560535792052681543116462530226_mp  /
           data xgk ( 11) / 0.857205233546061098958658510658944_mp  /
           data xgk ( 12) / 0.829565762382768397442898119732502_mp  /
           data xgk ( 13) / 0.799727835821839083013668942322683_mp  /
           data xgk ( 14) / 0.767777432104826194917977340974503_mp  /
           data xgk ( 15) / 0.733790062453226804726171131369528_mp  /
           data xgk ( 16) / 0.697850494793315796932292388026640_mp  /
           data xgk ( 17) / 0.660061064126626961370053668149271_mp  /
           data xgk ( 18) / 0.620526182989242861140477556431189_mp  /
           data xgk ( 19) / 0.579345235826361691756024932172540_mp  /
           data xgk ( 20) / 0.536624148142019899264169793311073_mp  /
           data xgk ( 21) / 0.492480467861778574993693061207709_mp  /
           data xgk ( 22) / 0.447033769538089176780609900322854_mp  /
           data xgk ( 23) / 0.400401254830394392535476211542661_mp  /
           data xgk ( 24) / 0.352704725530878113471037207089374_mp  /
           data xgk ( 25) / 0.304073202273625077372677107199257_mp  /
           data xgk ( 26) / 0.254636926167889846439805129817805_mp  /
           data xgk ( 27) / 0.204525116682309891438957671002025_mp  /
           data xgk ( 28) / 0.153869913608583546963794672743256_mp  /
           data xgk ( 29) / 0.102806937966737030147096751318001_mp  /
           data xgk ( 30) / 0.051471842555317695833025213166723_mp  /
           data xgk ( 31) / 0.000000000000000000000000000000000_mp  /

           data wgk (  1) / 0.001389013698677007624551591226760_mp  /
           data wgk (  2) / 0.003890461127099884051267201844516_mp  /
           data wgk (  3) / 0.006630703915931292173319826369750_mp  /
           data wgk (  4) / 0.009273279659517763428441146892024_mp  /
           data wgk (  5) / 0.011823015253496341742232898853251_mp  /
           data wgk (  6) / 0.014369729507045804812451432443580_mp  /
           data wgk (  7) / 0.016920889189053272627572289420322_mp  /
           data wgk (  8) / 0.019414141193942381173408951050128_mp  /
           data wgk (  9) / 0.021828035821609192297167485738339_mp  /
           data wgk ( 10) / 0.024191162078080601365686370725232_mp  /
           data wgk ( 11) / 0.026509954882333101610601709335075_mp  /
           data wgk ( 12) / 0.028754048765041292843978785354334_mp  /
           data wgk ( 13) / 0.030907257562387762472884252943092_mp  /
           data wgk ( 14) / 0.032981447057483726031814191016854_mp  /
           data wgk ( 15) / 0.034979338028060024137499670731468_mp  /
           data wgk ( 16) / 0.036882364651821229223911065617136_mp  /
           data wgk ( 17) / 0.038678945624727592950348651532281_mp  /
           data wgk ( 18) / 0.040374538951535959111995279752468_mp  /
           data wgk ( 19) / 0.041969810215164246147147541285970_mp  /
           data wgk ( 20) / 0.043452539701356069316831728117073_mp  /
           data wgk ( 21) / 0.044814800133162663192355551616723_mp  /
           data wgk ( 22) / 0.046059238271006988116271735559374_mp  /
           data wgk ( 23) / 0.047185546569299153945261478181099_mp  /
           data wgk ( 24) / 0.048185861757087129140779492298305_mp  /
           data wgk ( 25) / 0.049055434555029778887528165367238_mp  /
           data wgk ( 26) / 0.049795683427074206357811569379942_mp  /
           data wgk ( 27) / 0.050405921402782346840893085653585_mp  /
           data wgk ( 28) / 0.050881795898749606492297473049805_mp  /
           data wgk ( 29) / 0.051221547849258772170656282604944_mp  /
           data wgk ( 30) / 0.051426128537459025933862879215781_mp  /
           data wgk ( 31) / 0.051494729429451567558340433647099_mp  /
           
           !!-------------------------------------------------------!!
           !!                                                       !!
           !!     List of major variables                           !!          
           !!     -----------------------                           !!
           !!                                                       !!
           !!     centr  -  mid point of the interval               !! 
           !!     hlgth  -  half-length of the interval             !!
           !!     dabsc  -  abscissa                                !!
           !!     fval*  -  function value                          !!
           !!     resg   -  result of the 30-point gauss rule       !!
           !!     resk   -  result of the 61-point kronrod rule     !!
           !!     reskh  -  approximation to the mean value of f    !! 
           !!               over (a,b), i.e. to i/(b-a)             !! 
           !!                                                       !! 
           !!     Machine Dependent Constants                       !!
           !!     ---------------------------                       !!
           !!                                                       !!
           !!     epmach is the largest relative spacing.           !!
           !!     uflow is the smallest positive magnitude.         !! 
           !!                                                       !!
           !!-------------------------------------------------------!!
      
           centr = 0.5_mp*(b + a)
           hlgth = 0.5_mp*(b - a)
           dhlgth = abs(hlgth)

           !!  Compute the 61-point Kronrod approximation to the  !!
           !!  integral, and estimate the absolute error.         !!

           !!  First executable statement of DQK61.  !!

           resg = 0.0_mp 
           fc = f(centr)
           resk = wgk(31)*fc
           resabs = abs(resk)

           do j = 1,15
               jtw = j*2
               dabsc = hlgth*xgk(jtw)
               fval1 = f(centr - dabsc)
               fval2 = f(centr + dabsc)
               fv1(jtw) = fval1
               fv2(jtw) = fval2
               fsum = fval1 + fval2
               resg = resg + wg(j)*fsum
               resk = resk + wgk(jtw)*fsum
               resabs = resabs + wgk(jtw)*(abs(fval1) + abs(fval2))
           end do    

           do j = 1,15
               jtwm1 = j*2 - 1
               dabsc = hlgth*xgk(jtwm1)
               fval1 = f(centr - dabsc)
               fval2 = f(centr + dabsc)
               fv1(jtwm1) = fval1
               fv2(jtwm1) = fval2
               fsum = fval1 + fval2
               resk = resk + wgk(jtwm1)*fsum                    
               resabs = resabs + wgk(jtwm1)*(abs(fval1) + abs(fval2))
           end do 

           reskh = resk*0.5_mp
           resasc = wgk(31)*abs(fc - reskh)

           do j = 1,30,2
               resasc = resasc + wgk(j)*(abs(fv1(j) - reskh) + &
                        abs(fv2(j) - reskh))
               resasc = resasc + wgk(j+1)*(abs(fv1(j+1) - reskh) + &
                        abs(fv2(j+1) - reskh))
           end do

           result = resk*hlgth
           resabs = resabs*dhlgth
           resasc = resasc*dhlgth
           abserr = abs((resk - resg)*hlgth)

           if (resasc .ne. 0.0_mp .and. abserr .ne. 0.0_mp) then 
               abserr = resasc*min(0.1e+1_mp,(0.2e+3_mp*abserr/resasc)**1.5_mp)
           end if

           if (resabs .gt. uflow/(0.5e+2_mp*epmach)) then
              abserr = max((epmach*0.5e+2_mp)*resabs,abserr)
           end if

           return

       end Subroutine DQK61_4



       Subroutine DQPSRT_4(limit,last,maxerr,ermax,elist,iord,nrmax)
   
           !!----------------------------------------------------------------------!!
           !!                                                                      !!
           !!   Begin prologue DQPSRT                                              !!       
           !!                                                                      !!
           !!   Refer to  dqage,dqagie,dqagpe,dqawse                               !!
           !!   Routines called  (none)                                            !!
           !!   Revision date  810101   (yymmdd)                                   !!
           !!                                                                      !!
           !!   keywords  -   sequential sorting.                                  !!
           !!                                                                      !!
           !!   Author  Piessens,Robert,appl. math. & progr. div. - k.u.leuven     !!
           !!   de doncker,elise,appl. math. & progr. div. - k.u.leuven.           !!
           !!                                                                      !!
           !!                                                                      !!
           !!      Purpose                                                         !!
           !!     ---------                                                        !!
           !!                                                                      !!
           !!      This routine maintains the descending ordering in the           !!
           !!      list of the local error estimated resulting from the            !!
           !!      interval subdivision process. at each call two error            !!
           !!      estimates are inserted using the sequential search              !!
           !!      method, top-down for the largest error estimate and             !!
           !!      bottom-up for the smallest error estimate.                      !!
           !!                                                                      !!
           !!                                                                      !!
           !!      Description                                                     !!
           !!     -------------                                                    !!
           !!                                                                      !!
           !!      ordering routine                                                !!
           !!      standard fortran subroutine                                     !!
           !!                                                                      !!
           !!      Parameters:                                                     !!
           !!     -------------                                                    !!
           !!                                                                      !! 
           !!                  Meaning At Output                                   !!
           !!                 -------------------                                  !!
           !!                                                                      !!
           !!         limit  -   Integer                                           !! 
           !!                    Maximum number of error estimates the list        !!
           !!                    can contain                                       !!
           !!                                                                      !! 
           !!         last   -   Integer                                           !!
           !!                    Number of error estimates currently in the list.  !!
           !!                                                                      !!
           !!         maxerr -   Integer                                           !!
           !!                    maxerr points to the nrmax-th largest error       !! 
           !!                    estimate currently in the list                    !!
           !!                                                                      !!
           !!         ermax  -   Real                                              !! 
           !!                    nrmax-th largest error estimate.                  !!
           !!                    ermax = elist(maxerr)                             !!
           !!                                                                      !!
           !!         elist  -   Real                                              !!
           !!                    Vector of dimension last containing               !!
           !!                    the error estimates.                              !!
           !!                                                                      !!
           !!         iord   -   Integer                                           !!
           !!                    Vector of dimension last, the first k elements    !!
           !!                    of which contain pointers to the error            !!
           !!                    estimates, such that                              !!
           !!                    elist(iord(1)),...,  elist(iord(k))               !!
           !!                    form a decreasing sequence, with                  !!
           !!                    k = last if last.le.(limit/2+2), and              !!
           !!                    k = limit+1-last otherwise                        !!
           !!                                                                      !!
           !!         nrmax  -   Integer                                           !!
           !!                                                                      !!
           !!         maxerr =   iord(nrmax)                                       !!
           !!                                                                      !!
           !! end prologue  dqpsrt                                                 !!
           !!                                                                      !!
           !!----------------------------------------------------------------------!!

           implicit none

           real(kind=mp), intent(inout) :: ermax, elist(*)      
      
           real(kind=mp) :: errmax, errmin

           integer, intent(inout) :: iord(*), maxerr, nrmax
           integer, intent(in) :: limit

           integer :: i, ibeg, last, isucc, j, jbnd, jupbn, k, ido 

           !!  Check whether the list contains more than  !!
           !!  two error estimates.                       !!
            
           if (last .le. 2) then
               iord(1) = 1
               iord(2) = 2
               maxerr = iord(nrmax)
               ermax = elist(maxerr)
               return
           end if

           !!  This part of the routine is only executed if, due to a    !!
           !!  difficult integrand, subdivision increased the error      !!
           !!  estimate. In the normal case the insert procedure should  !!
           !!  start after the nrmax-th largest error estimate.          !!

           errmax = elist(maxerr)

           if (nrmax .ne. 1) then
               ido = nrmax - 1

               do i = 1, ido
                  isucc = iord(nrmax - 1)
                  
                  !!  Exit from do  loop.  !!
                  
                  if (errmax .le. elist(isucc)) then
                      exit
                  end if
            
                  iord(nrmax) = isucc
                  nrmax = nrmax - 1
               end do
           end if

           !!  Compute the number of elements in the list to be maintained  !!
           !!  in descending order. This number depends on the number of    !!
           !!  subdivisions still allowed.                                  !!

           jupbn = last

           if (last .gt. (limit/2 + 2)) then 
               jupbn = limit + 3 - last
           end if

           errmin = elist(last)

           !!  Insert errmax by traversing the list top-down, starting  !!
           !!  comparison from the element list elist(iord(nrmax + 1)). !!

           jbnd = jupbn - 1
           ibeg = nrmax + 1

           if (ibeg .le. jbnd) then 

               do i = ibeg, jbnd
                  isucc = iord(i)

                  !!  Exit from do loop.  !!

                  if (errmax .ge. elist(isucc)) then 
                      exit
                  end if

                  iord(i-1) = isucc
               end do

               if (errmax .ge. elist(isucc)) then 

                   !! Insert errmin by traversing the list bottom up.  !!

                   iord(i-1) = maxerr
                   k = jbnd

                   do j = i, jbnd
                      isucc = iord(k)

                      !!  Exit from do loop.  !!

                      if (errmin .lt. elist(isucc)) then 
                          exit
                      end if
           
                      iord(k+1) = isucc
                      k = k - 1
                   end do

                   if (errmin .lt. elist(isucc)) then
                       iord(k+1) = last
                   else
                       iord(i) = last
                   end if
          
               else
                   iord(jbnd) = maxerr
                   iord(jupbn) = last
               end if

           else
               iord(jbnd) = maxerr
               iord(jupbn) = last
           end if

           !!  Set maxerr and ermax.  !!

           maxerr = iord(nrmax)
           ermax = elist(maxerr)  
           
           return

       end Subroutine DQPSRT_4



        Subroutine Xerror(ier)

               !!-----------------------------------------------------------!!
               !!                                                           !!
               !!  If called, this subroutine will notify the user that an  !!
               !!  error has occurred during quadrature and will automat-   !!
               !!  ically end the program.                                  !!
               !!                                                           !!
               !!        If ier > 0, this subroutine will terminate. The    !!
               !!  estimates for the result and error are not reliable      !!
               !!  therfore it is assumed that the requested accuracy has   !!
               !!  not been achieved.                                       !!
               !!                                                           !!
               !!                                                           !!
               !!     ier = 1                                               !!
               !!            Maximum number of subdivisions allowed         !!
               !!            has been achieved. one can allow more          !!
               !!            subdivisions by increasing the value of        !!
               !!            limit (and taking the according dimension      !!
               !!            adjustments into account). however, if         !!
               !!            this yield no improvement it is advised        !!
               !!            to analyze the integrand in order to           !!
               !!            determine the integration difficulaties.       !!
               !!            if the position of a local difficulty can      !!
               !!            be determined (i.e.singularity,                !!
               !!            discontinuity within the interval) one         !!
               !!            will probably gain from splitting up the       !!
               !!            interval at this point and calling the         !!
               !!            integrator on the subranges. if possible,      !!
               !!            an appropriate special-purpose integrator      !!
               !!            should be used which is designed for           !!
               !!            handling the type of difficulty involved.      !!
               !!                                                           !!
               !!      ier = 2                                              !!
               !!            The occurrence of roundoff error is            !!
               !!            detected, which prevents the requested         !!
               !!            tolerance from being achieved.                 !!
               !!                                                           !!
               !!     ier = 3                                               !!
               !!            Extremely bad integrand behaviour occurs       !! 
               !!            at some points of the integration              !! 
               !!            interval.                                      !!
               !!                                                           !!
               !!     ier = 6                                               !!
               !!            The input is invalid, because                  !!
               !!            (epsabs.le.0 and                               !!
               !!            epsrel.lt.max(50*rel.mach.acc.,0.5d-28))       !!
               !!            or limit.lt.1 or lenw.lt.limit*4.              !!
               !!            result, abserr, neval, last are set            !!
               !!            to zero.                                       !!
               !!            except when lenw is invalid, iwork(1),         !!
               !!            work(limit*2+1) and work(limit*3+1) are        !!
               !!            set to zero, work(1) is set to a and           !!
               !!            work(limit+1) to b.                            !!
               !!                                                           !!
               !!-----------------------------------------------------------!!

               implicit none

               integer, intent(in) :: ier

               write(iout,10)

               select case(ier)
                  case(1)
                    write(iout,20)
                  case(2)
                    write(iout,30)
                  case(3)
                    write(iout,40)
                  case(6)
                    write(iout,50)
               end select

               stop  !!  Integration has failed therefore the code must  !!
                     !!  terminate here.                                 !!

               return
           
               !!-------------------------------------------------!!
               !!                                                 !!
               !!                FORMAT  STATMENTS                !!
               !!                                                 !!
               !!-------------------------------------------------!!

               10 format(/10x,'Sorry quadrature has failed - program terminated')

               20 format(/10x,'ier = 1'                                                       &
                         /10x,'Maximum number of subdivisions allowed has been achieved.'     &
                         /10x,'One can allow more subdivisions by increasing the value of'    &
                         /10x,'limit. However if this yields no improvement it is advised'    &
                         /10x,'to analyse the integrand in order to determine the '           &
                         /10x,'integration difficulties.')

               30 format(/10x,'ier = 2'                                                       &
                         /10x,'The occurrence of roundoff error is detected which prevents'   &
                         /10x,'the required tolerance from being achieved.')

              40 format(/10x,'ier = 3'                                                        &
                        /10x,'Extremely bad integrand behaviour occurs at some points'        &
                        /10x,'of the integration interval.')

               50 format(/10x,'ier = 6'                                                       &
                         /10x,'Input data is vaild.')

           end Subroutine Xerror


end Module Cross_Section


Program Ion_Atom

        !!----------------------------------------------------------------!!
        !!                                                                !!
        !!     Description                                                !!
        !!    -------------                                               !!
        !!                                                                !!
        !!  This is the main programming unit.                            !!
        !!                                                                !!
        !!                                                                !!
        !!     Details                                                    !!
        !!    ---------                                                   !!
        !!                                                                !!
        !!  The main program calls the subroutine Get_Data which reads    !!
        !!  the input data; this data then becomes accessable to all      !!
        !!  modules. The main prgram then determines the type of cross    !!
        !!  section to be evaluated via the variable 'my_section' and     !!
        !!  calls the subroutine Total_Differential which calculates the  !!
        !!  cross sections.                                               !!
        !!                                                                !!
        !!                                                                !!
        !!     Note                                                       !!
        !!    ------                                                      !!
        !!                                                                !!
        !!  High Performance Fortran directives are used in the main      !!
        !!  program to allow cross sections to be evaluated concurrently. !!
        !!  This can greatly speed-up execution time if the user has      !!
        !!  access to a parallel machine.                                 !!
        !!      If a HPF compiler is used then it should give a warning   !!
        !!  that the independent do-loops used in calling the subroutine  !!
        !!  Total_Differential are not pure. These warnings can be        !!
        !!  ignored since each cross section is calculated completely     !!
        !!  independent of all others thus avoiding data corruption.      !!
        !!                                                                !!
        !!----------------------------------------------------------------!!
        
        use Precision, only : mp
        use Global_Data
        use Cross_Section
        implicit none

        real(kind=mp) :: Differential_Cross_Section
        real(kind=mp), allocatable, dimension(:) :: ratio, Cross, Energy

        integer, parameter :: number_of_processors = 4  ! Need I say any more?

!HPF$   PROCESSORS procs(number_of_processors)
!HPF$   DISTRIBUTE (BLOCK) ONTO procs :: ratio, Cross, Energy 

        integer :: loop, points, stat

        !!  Create an output file for the raw data. This is very useful  !!
        !!  for plotting graphs.                                         !!
        
        open(unit=iout,file='he_data.dat',status='unknown',iostat=stat)
        if (stat .ne. 0) then
            write(iout,*)
            write(iout,*) "Error: could not create file he_data.dat"
            write(iout,*) "Check the read/write permissions on your directory."
            write(iout,*) 'Code has terminated.'
            stop
        end if

        open(unit=4,file='total_output.dat',status='unknown')
      
        call Machine_Accuracy  ! Determine machine dependent parameters.
        call Get_Data  ! Read the input data from file.

        if (my_section .eq. 0) then 
            points = int((Projectile_energy_finish - Projectile_energy_start)/    &    
                     Projectile_energy_increment) + 1 
        else
            points = int((Electron_Energy_Finish - Electron_Energy_Start)/     &
                     Electron_Energy_Step) + 1
        end if

        allocate (ratio(points), Cross(points), Energy(points))

        if (my_section .eq. 0) then
        
!HPF$   INDEPENDENT
            do loop = 1,points
                      Energy(loop) = Projectile_Energy_Start + (loop - 1)*  &
                      Projectile_Energy_Increment
            end do
            
        else
        
!HPF$   INDEPENDENT
            do loop = 1,points
                      Energy(loop) = Electron_Energy_Start + (loop - 1)*    &
                      Electron_Energy_Step
            end do
            
        end if
        
        if (my_section .eq. 0) then  ! Calculate total differential cross sections.
        
!HPF$   INDEPENDENT
            do loop = 1,points
               call Total_Differential(Energy(loop),Differential_Cross_Section)

               if (target .eq. 0 .or. target .eq. 2) then  !  For a hydrogen-like target we must
                                                           !  divide the cross section by two since
                                                           !  we have adopted the independent electron model.
                  Cross(loop) = Differential_Cross_Section*0.5_mp
               else
                  Cross(loop) = Differential_Cross_Section
               end if

            end do
            
         else  ! Calculate single or double differential cross sections.
         
!HPF$   INDEPENDENT
            do loop = 1,points
               call Total_Differential(Energy(loop),Differential_Cross_Section)
               ratio(loop) = (global_k/velocity)

               if (target .eq. 0 .or. target .eq. 2) then  !  For a hydrogen-like target we must 
                                                           !  divide the cross section by two since
                                                           !  we have adopted the independent electron model.
                  
                  Cross(loop) = pi*Differential_Cross_Section          
               else
                  Cross(loop) = 2.0_mp*pi*Differential_Cross_Section 
               end if

            end do
            
         end if

         !!  Write raw data to output file for plotting graphs.  !!
         
         if (my_section .eq. 0) then 
             do loop = 1,points
                write(4,*) Energy(loop), Cross(loop)
             end do
         else
             do loop = 1,points
                write(4,*) ratio(loop), Cross(loop)
             end do
         end if

        write(iout,100)
        write(iout,110)
        write(iout,120)
        write(iout,130)
        write(iout,150)
        write(iout,160)
        write(iout,165)

        select case(target)
            case(0)
              write(iout,166)
            case(1)
              write(iout,167)
            case(2)
              write(iout,168)
        end select

        if (no .eq. 1) then
            write(iout,170)
        else
            write(iout,180)
        end if

        select case(my_section)
            case(0)  ! Total cross section
               write(iout,190)
            case(1)  ! single differential cross section
               write(iout,200)
            case(2)  ! Total differential cross section
               write(iout,210)
        end select

        select case(my_section)
            case(0)  ! Total cross section
               write(iout,220) Projectile_Energy_Start
               write(iout,230) Projectile_Energy_Finish
               write(iout,240) Projectile_Energy_Increment
            case(1)  ! single differential cross sections
               write(iout,250) Electron_Energy_Start
               write(iout,260) Electron_Energy_Finish
               write(iout,270) Electron_Energy_Start
               write(iout,280) Projectile_Energy_Start
            case(2)  ! double differential cross section
               write(iout,250) Electron_Energy_Start
               write(iout,260) Electron_Energy_Finish
               write(iout,270) Electron_Energy_Step
               write(iout,280) Projectile_Energy_Start
               write(iout,290) angle
        end select

        write(iout,300)
        write(iout,310)
        write(iout,320) tolerance
        write(iout,330) acc

        write(iout,340)
        write(iout,350)

        if (my_section .eq. 0) then 
            write(iout,351)
        else
            write(iout,360)
        end if

        write(iout,370)

        if (my_section .eq. 0) then
            do loop = 1,points
               write(iout,*) Energy(loop), Cross(loop),'*1e-20sqr.metres/eV/sr' 
            end do
        else
            do loop = 1,points
               write(iout,*) Energy(loop), Cross(loop),'*1e-20sqr.metres/eV/sr'
            end do
        end if
     
     
100  format(//20x,'Ion-Atom-Wave')
110  format(19x,'---------------'/)
120  format(//10x,'Ion-Atom-Wave calculates cross sections for the single ionization')
130  format(10x,'of atomic and molecular hydrogen-like targets and helium-like targets')
150  format(10x,'by both light and heavy ion impact.'//)
160  format(20x,'Input data used for this run :')
165  format(19x,'------------------------------')
166  format(/7x,"The target is a hydrogen atom")
167  format(/7x,"The target is a helium atom")
168  format(/7x,"The target is a hydrogen molecule")
170  format(/7x,'The CDW model was used.')
180  format(7x,'The CDW-EIS model was used.')
190  format(/7x,'Total cross section evaluated for the following values :')
200  format(/7x,'single differential cross section evaluated for the following values :')
210  format(/7x,'double differential cross section evaluated for the following values :')

220  format(/10x,'Program started with the following projectile energy (keV) : ',e16.8)
230  format(10x,'Program ended with the following projectile energy (keV)   : ',e16.8)
240  format(10x,'The projectile energy was incremented by                   : ',e16.8)

250  format(/10x,'Program started with the following electron energy (eV) : ',e16.8)
260  format(10x,'Program ended with the following electron energy (eV)   : ',e16.8)
270  format(10x,'The electron energy was incremented by                  : ',e16.8)
280  format(10x,'Energy of the incoming ion (keV)                        : ',e16.8)
290  format(10x,'Angle of electron emission (in radians)                 : ',e16.8)

300  format(//20x,'Numerical details')
310  format(19x,'-------------------')
320  format(/10x,'Tolerance used throughout the code :',e16.8)
330  format(10x,'Tolerance used for the quadrature   :',e16.8)

340  format(//20x,'PROGRAM OUTPUT')
350  format(19x,'----------------')
351  format(/7x,'The first column contains the projectile energy (keV) and the second')
360  format(/7x,'The first column contains the electron energy (eV) and the second')
370  format(7x,'column contains the cross sections'//)

         stop

end Program Ion_Atom 
