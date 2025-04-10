!!!*************************************************************
! 文件/File: ion_atom_driver.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:17:02
!*************************************************************

!!   ion_atom_driver.f90
!!-------------------------------------------------------------!!
!!                                                             !!
!!     Description                                             !!
!!    -------------                                            !!
!!                                                             !!
!!  Get_Info is the driver program for the Ion-Atom programs   !!
!!  which requests the input data in a user-friendly fashion.  !!
!!  The output is written to unit 'iwrite' which is called     !!
!!  'info_total.dat'.                                          !!
!!                                                             !!
!!-------------------------------------------------------------!!



Module Precision

       implicit none
       integer, parameter :: mp = selected_real_kind(10,30)

end Module Precision


Program Get_Info

        use Precision
        implicit none
	
	real(kind=mp) :: ZP, Projectile_Energy_Start, Projectile_Energy_Finish,        &
       			 Projectile_Energy_Increment, Electron_Energy_Start,           &
                         Electron_Energy_Step, angle, tolerance, acc, dummy_variable,  &
                         tol_lower, Electron_Energy_Finish 

	integer :: target, no, my_section, error, check_file
	integer, parameter :: iwrite = 9  ! Output channel

        !! Evaluate the accuracy of the machine being used.  !!
        !! This helps us put a lower bound on the tolerance  !!
        !! that will be used for the code.                   !!

        tol_lower = epsilon(dummy_variable)  ! Lower bound on tolerance.

	!!  Read input values and check validity.  !!

        error = 0  ! Initialise error flag.

        do
           print*, ""
           print*, ""
           print*, " Please choose a target: "
           print*, "   0 = hydrogen"
           print*, "   1 = helium"
           print*, "   2 = molecular hydrogen"
           read*, target

           if (target .lt. 0 .or. target .gt. 2) then 
               error = error + 1
               call Failure(error)
           else
               exit
           end if
        end do

	error = 0  ! Initialise error flag.

	do
           print*, ""
           print*, ""
           print*, " Please enter which model you would like to use :"
           print*, ""
           print*, "  1  =  CDW"
           print*, "  2  =  CDW-EIS"
           read*, no

           if (no .lt. 1 .or. no .gt. 2) then 
               error = error + 1
               call Failure(error)
      	   else
	       exit
           end if
        
        end do

	error = 0
	do
           print*, ""
           print*, ""
           print*, " What type of cross section do you want to calculate :"
           print*, ""
           print*, "  0 = total cross section."
           print*, "  1 = single differential cross section."
           print*, "  2 = double differential cross section."
           read*, my_section

           if (my_section .lt. 0 .or. my_section .gt. 2) then 
	       error = error + 1
               call Failure(error)
           else
               exit
           end if

        end do

        error = 0
        do
           print*, ""
           print*, ""
           print*, " Enter projectile charge :"
           read*, ZP

           if (ZP .lt. 0.0_mp) then  ! Positive ions only.
               error = error + 1
               call Failure(error)
           else
               exit
           end if

        end do

        if (my_section .eq. 0) then   ! Need data for total cross section.
            
            error = 0
            do
               print*, ""
               print*, ""
               print*, " Enter starting value for projectile kinetic energy in keV/amu :"
               read*, Projectile_Energy_Start

               if (Projectile_Energy_Start .lt. 0.0_mp) then 
                   error = error + 1
                   call failure(error)
               else
                   exit
               end if

            end do

            error = 0
            do
               print*, ""
               print*, ""
               print*, " Enter upper limit for projectile kinetic energy in keV/amu :"
               read*, Projectile_Energy_Finish

               if (Projectile_Energy_Finish .lt. Projectile_Energy_Start &
                   .or. Projectile_Energy_Finish .lt. 0.0_mp) then
                   error = error + 1
                   call failure(error)
               else
                   exit
               end if

            end do

            error = 0
            do
               print*, ""
               print*, ""
               print*, " Enter increment value for projectile kinetic energy in keV/amu :"
               read*, Projectile_Energy_Increment
          
               if (Projectile_Energy_Increment .lt. 0.0_mp) then 
                   error = error + 1
                   call failure(error)
               else
                   exit
               end if

            end do

            Electron_Energy_Start = 1.0_mp       ! For total cross sections these values 
            Electron_Energy_Finish = 1.0_mp      ! are not needed so we supply dummy values.
            Electron_Energy_Step = 1.0_mp

	else   ! Need data for single and double differential cross sections.

            error = 0
            do
               print*, ""
               print*, ""
               print*, " Enter starting value for ejected electron energy in eV :"
               read*, Electron_Energy_Start

               if (Electron_Energy_Start .lt. 0.0_mp) then 
                   error = error + 1
                   call failure(error)
               else
                   exit
               end if

            end do

            error = 0
            do
               print*, ""
               print*, ""
               print*, " Enter upper limit for ejected electron energy in eV :"
               read*, Electron_Energy_Finish
               
               if (Electron_Energy_Finish .lt. Electron_Energy_Start .or.   &
                   Electron_Energy_Finish .lt. 0.0_mp) then 
                   error = error + 1
                   call failure(error)
               else
                   exit
               end if

            end do

            error = 0
	    do
               print*, ""
               print*, ""
               print*, " Enter increment value for ejected electron energy in eV :"
               read*, Electron_Energy_Step

               if (Electron_Energy_Step .lt. 0.0_mp) then 
                   error = error + 1
                   call failure(error)
               else
                   exit
               end if

            end do 

            error = 0
            do
               print*, ""
               print*, ""
               print*, " Enter energy for the projectile in keV/amu:"
               read*, Projectile_Energy_Start

               if (Projectile_Energy_Start .lt. 0.0_mp) then ! Positive ions only.
                   error = error + 1
                   call failure(error)
               else
                   exit
               end if

            end do

            Projectile_Energy_Finish = 1.0_mp       ! For single and double differential cross
            Projectile_Energy_Increment = 1.0_mp    ! sections these values are not needed so 
                                                    ! supply dummy values.
	end if

        if (my_section .eq. 2) then  ! For double differential cross sections only.

            print*, ""
            print*, ""
            print*, " Enter a value for the polar angle of emission for the"
            print*, " ejected electron (in radians please) :"
            read*, angle
           
        else   ! A value for 'angle' is not required for TDSC or 
               ! SDCS so we supply a dummy value.
            angle = 1.0_mp
        end if

        !!  Ask user to enter tolerances and check them with  !!
        !!  the lower bound set by tol_lower.                 !!
        
        error = 0
        do
           print*, ""
           print*, ""
           print*, " Enter a tolreance to be used throughout the calculation :"
           read*, tolerance

           if (tolerance .lt. tol_lower .or. tolerance .gt. 1.0_mp) then 
               error = error + 1
               call failure(error)
           else
               exit
           end if

        end do 
 
        error = 0
        do
           print*, ""
           print*, ""
           print*, " Enter a tolerance to used for the integration, preferably "
           print*, " greater than the previos tolerance :"
           read*, acc

           if (acc .lt. tol_lower .or. acc .gt. 1.0_mp) then 
               error = error + 1
               call failure(error)
           else
               exit
           end if

        end do
        
        !!  Finished collecting the information. Writing it to file now.  !!

        open(unit=iwrite, file='info_total.dat',status='unknown', iostat=check_file)
        if (check_file .ne. 0) then 
            print*, "Error could not write to file info_total.dat "
            print*, ""
            print*, "Check the read/write permissions on your directory."
            print*, "This program is terminating."
            stop
        end if

        write(iwrite,*) target
        write(iwrite,*) no
        write(iwrite,*) my_section
        write(iwrite,*) ZP
        write(iwrite,*) Projectile_Energy_Start
        write(iwrite,*) Projectile_Energy_Finish
        write(iwrite,*) Projectile_Energy_Increment
        write(iwrite,*) Electron_Energy_Start
        write(iwrite,*) Electron_Energy_Finish
        write(iwrite,*) Electron_Energy_Step
        write(iwrite,*) angle
        write(iwrite,*) tolerance
        write(iwrite,*) acc

Contains

       Subroutine Failure(count)

            !!  This subroutine is called when invalid data is entered.  !!
            !!  If the data is entered incorrectly more than five times  !!
            !!  then the code will terminate automatically.              !!
       
            integer, intent(in) :: count

            if (count .gt. 5) then 
               print*, ""
               print*, ""
               print*, " Five failures in a row, this program is terminating."
               stop
            else
               print*, ""
               print*, ""
               print*, " Sorry this value is invalid...  try again."
            end if

        end Subroutine Failure

end Program Get_Info
