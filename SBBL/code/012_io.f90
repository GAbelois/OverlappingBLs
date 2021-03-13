!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------     
module io
!-----------------------------------------------------------------------
! Modified to fit namelist feature, use allocatable attribution and
! >> remove the parameter and save attribution (Bicheng Chen 06/15/2016)
!-----------------------------------------------------------------------     
use types, only:rprec
use param
use intermediate
implicit none
save
! ---
integer,parameter :: num_hour_out = 1
! TOMAS CHOR included base in param.nml
integer :: nwrite
!integer,parameter::base=200,nwrite=base

!!!!  io_spec=.true. output plan-averaged spectrum
logical,parameter :: io_spec=.true., output_fields_3d_flag=.false.
integer,parameter :: spec_write_freqz=3000 !fields_3d_write_freqz=p_count*6
integer,parameter :: spec_write_start=1,spec_write_end=1000000
!integer,parameter::spec_write_start=1,spec_write_end=24*base
!! --------------------------------------------------------------------
!! The following block defines parameters for instantaneous slice output
!! inst_slice_freqz controls the output frequency
!! The 5 slices outputted every inst_slice_freqz (i.e. u,v,w,T,Cs in this order) ...
!! ... are saved in the 3rd dimension of inst_slice_array for inst_array_lim times ...
!! ... following which this array is outputted as a binary file and the process
!! ... starts over
!! Therefore, each file will contain inst_array_lim x-z slices of 5 variables
!! This information is important for post-processing of these binary files
!! --------------------------------------------------------------------

logical, parameter :: inst_slice_flag = .false.
integer, parameter :: num_vars = 4                  ! works currently only for u,v,w,T due to the size difference in Cs
!integer, parameter :: slice_inst = (nzt-1)/2    ! sets the value of the y node for the chosen x-z inst slice
!integer, parameter :: inst_slice_freqz = 5         ! 
!integer, parameter :: inst_array_lim=200
!real(kind=rprec), dimension(ldx,nz+1,num_vars*inst_array_lim) :: inst_slice_array
!integer:: inst_slice_counter

logical, parameter :: cumulative_time = .true.
character(*), parameter :: fcumulative_time = '../readin/total_time.dat'

integer, parameter :: n_avg_stats = 100             ! interval for updates in avg_stats
character(*), parameter :: end_hdr_avg = '# end header'

!! --------------------------------------------------------------------
!! The following block defines parameters for use in avgslice and theta_slice
!! --------------------------------------------------------------------
!integer,parameter :: average_dim_select_flag=1-(average_dim_num/2)
!! The variable average_dim_select_flag generates the following values based
!! on the value of average_dim_num in param.f90 :-
!! a) average_dim_num = 2 :
!! 	average_dim_select_flag = 0 ==> average the 3d array over x and y and output the z profile
!! b) average_dim_num = 1 :
!! 	average_dim_select_flag = 1 ==> average the 3d array over y and output the (x,z) profile
!integer, parameter :: dim1_size=average_dim_select_flag*(nxt-nz+1)+nz-1
!integer, parameter :: dim2_size=average_dim_select_flag*(nz-2)+1
!integer, parameter :: dim1_global=average_dim_select_flag*(nxt-nzt+1)+nzt-1
!integer, parameter :: dim2_global=average_dim_select_flag*(nzt-2)+1
!! --------------------------------------------------------------------
!! --------------------------------------------------------------------

! ---  time_spec>0 output time series spectrum (need additional calcu.)
integer, parameter :: time_spec = 0
integer :: n_obs, jt_total_init
integer, allocatable :: obs_pt(:,:)

! ---  io_mean=.true. output small domain time-averaged velocity
logical, parameter :: io_mean = .false.
integer :: jx_pls, jx_ple, width
integer :: jy_pls, jy_ple
real(rprec), dimension(:, :, :), allocatable :: mean_u, mean_v, mean_w,     &
                                                mean_u2, mean_v2, mean_w2

! --- mpi io
integer, dimension(3) :: gsize, lsize, start   ! dimensions of global and local variables 

!!!!  io_lambda2
!logical,parameter::io_lambda2=.false.
!real(rprec), dimension(:, :, :), allocatable :: lam2

! ---
contains
! --------------------------------------------------------------------
! file number 1-10 are used for temporary use
! 11-19 are basic output
! 20-40 are avgslice
! use >50 for debugging
! --------------------------------------------------------------------
    subroutine screen_display()
!-----------------------------------------------------------------------
!   echo the simulation parameters to make sure they are 
!   inputed correctly
!-----------------------------------------------------------------------
    use param
    implicit none
    ! --- Print computational environment on screen
    if ( rank == 0 ) then
        open(1,file='../output/screen.out')
        write(1,*) "++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	    write(1,*) "++--------------------------------------------------++"
	    write(1,*) "++          OCEANIC CANOPY FLOW WITH LES            ++"
	    write(1,*) "++--------------------------------------------------++"
	    write(1,*) "++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	    write(1,*) ""
	    write(1,*) ""
        
        write(1,*) "++--------------------------------------------------++"
        write(1,*) "++                 LES PARAMETERS                   ++"
        write(1,*) "++--------------------------------------------------++"
        if ( model .eq. 1 ) then  
            write(1,*) "    Smagrinsky model   " 
            write(1,*) "    Cs                =", cs
        else if ( model .eq. 2 ) then
            write(1,*) "    Dynamic smagrinsky model   "    
        else if ( model .eq. 3 ) then
            write(1,*) "    Scale-dependent dynamic smagrinsky model   " 
        else if ( model .eq. 4 ) then
            write(1,*) "    Lagrangian scale similar model   " 
        else if ( model .eq. 5 ) then
            write(1,*) "    Lagrangian scale-depandent model   "  
        end if
        write(1,*) "++--------------------------------------------------++"
        write(1,*) ""
        
        write(1,*) "++--------------------------------------------------++"
        write(1,*) "++                BACIC PARAMETERS                  ++"
        write(1,*) "++--------------------------------------------------++"
        write(1,*) "      nxt = ", nxt, "  nyt = ", nyt, "  NZ = ", nzt
        write(1,*) "++--------------------------------------------------++"
        write(1,*) "molecular viscosity (dimensional) = ", nu_molec
        write(1,*) "Number of timesteps               = ", nsteps
        write(1,*) "Time step size                    = ", dt
        write(1,*) "Number of processors              = ", npy*npz
        write(1,*) 'sampling stats every ', c_count, ' timesteps'
        write(1,*) 'writing stats every ', p_count, ' timesteps'
        write(1,*) "++--------------------------------------------------++"
        write(1,*) ""
        ! --- 
        close(1)
    end if
    ! ---
    end subroutine screen_display
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------          
    subroutine init_avg (file_avg, n_ccol, a_avg, n_avg, leaveopn)
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------   
    implicit none
    ! ---
    character(*), intent(in) :: file_avg
    integer, intent(in) :: n_ccol  !--num. coordz columns: x, y, etc.
    real(rprec), intent(out) :: a_avg(:, :)
    integer, intent(out) :: n_avg
    character (*), optional, intent (in) :: leaveopn
    ! ---
    character(128) :: buff
    logical :: exst, opn
    integer :: j
    real(rprec) :: z(n_ccol)
! ---
    inquire (file=file_avg, exist=exst, opened=opn)
    if (exst) then
        if (.not. opn) then
            open(1, file=file_avg)
            read(1, '(a)') buff

            if (buff(1:1) == '#') then
                read (buff(2:), *) n_avg
            else
                write (*, *) 'avg_stats: error'
                write (*, *) trim (file_avg), ' does not have expected format on line 1'
                stop  !--need to replace all stops with nice mpi exits
            end if
        end if

    !--skip data header lines here
        do
            read (1, '(a)') buff
            if (trim (buff) == trim (end_hdr_avg)) exit
        end do

        do j = 1, size (a_avg, 2)
            read (1, *) z, a_avg(:, j)  !--z is just placeholder here
        end do

        if (present (leaveopn)) then
            if (leaveopn /= 'yes') close (1)  !--case sensitive here
        else
            close (1)
        end if
    else
        n_avg = 0
        a_avg = 0._rprec
    end if
    ! ---
    end subroutine init_avg
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------     
    subroutine write_avg (file_avg, n_avg, x, a_avg, hdr, position)
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
    implicit none
    ! ---
    character(*), intent(in) :: file_avg
    integer, intent(in) :: n_avg
    real(rprec), intent(in) :: x(:, :)  !--coordz columns for x, y, etc
    real(rprec), intent(in) :: a_avg(:, :)

    character(*), optional, intent(in) :: hdr
    character(*), optional, intent(in) :: position

    character(64) :: r_fmt, fmt
    character(32) :: posn

    integer :: j

! --- check sizes compatible
    if (size (x, 2) /= size (a_avg, 2)) then
        write (*, *) 'write_avg: error with sizes of x, a_avg'
        stop
    end if

    if (present (position)) then
        posn = position
    else
        posn = 'rewind'
    end if

    open (1, file=file_avg, position=posn)
    if (trim (posn) /= 'append') then  !--case sensitive
        write (1, '(a,i0)') '# ', n_avg  !--not needed when appending
    end if

    if (present (hdr)) then
        !--write data header, if present
        write(1, '(a)') trim (hdr)
    end if

    !--write something to indicate end of header, always do this
    write(1, '(a)') end_hdr_avg

    !--define output format
    write(r_fmt, '(2(a,i0))') 'es', precision (1._rprec) + 7,  &
                               '.', precision (1._rprec)
    write(fmt, '(a,i0,3a)') '(', size (x, 1) + size (a_avg, 1),  &
                            '(1x,', trim (r_fmt), '))'
!--write to file
    do j = 1, size (a_avg, 2)
        write (1, fmt) x(:, j), a_avg(:, j)
    end do
    close (1)
    ! ---
    end subroutine write_avg
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
    subroutine post_spec(jt_local)
    
    use sim_param,only:u,v,w,theta,pcon
    use param
  use fft

  implicit none
  real(kind=rprec),dimension(4,nxt/2,nz-1)::spectra_uvwT
  real(kind=rprec),dimension(4,nxt/2,nzt-1)::spectra_uvwT_tot
  integer,intent(in)::jt_local
  integer::k,jz!,z

  real(kind=rprec),dimension(nz-1)::z   ! Ying changed z to array (09/30/2010)
  ! Chamecki changed z from integer to real (08/04/2006)
  real(kind=rprec),dimension(nzt-1)::z_tot

  character(len=64)::fname1
  !$if ($MPI)
  !$define $lbz 0
  integer :: recvcounts(npz)
  integer :: displs(npz)
  !$else
  !$define $lbz 1
  !$endif

  !$if ($MPI)
  recvcounts = size (z)
  displs = coord_of_rank * recvcounts
  do jz=1,nz-1
    z(jz) = (coordz*(nz-1)+jz-0.5_rprec)*dz*z_i
  enddo
  call mpi_gatherv (z(1), size (z), MPI_RPREC,&
     z_tot(1), recvcounts, displs,       &
     MPI_RPREC, 0, comm, ierr)
  !$else
  !do jz=1,nz-1
  !  z(jz)=(jz-0.5_rprec)*dz*z_i
  !  z_tot(jz)=z(jz)
  !enddo
  !$endif

  if((.not. USE_MPI) .or. (USE_MPI .and. coordz == 0)) then
    write(fname1,'(a,a)') path//'output/spec_x','.dat'
    open(82,file=fname1,form='formatted')
    do jz = 1,nzt-1
      write(82,*) (real(kx_2d(k,1)/z_i*z_tot(jz)),k=1,nxt/2)
    enddo
    close(82)
  endif

!write(fname1,'(a,a)') path//'output/spec_x','.dat'
!open(82,file=fname1,form='formatted')
  do jz=1,nz-1
!do jz=$lbz,nz-1
!   z=(jz-0.5_rprec)*dz*z_i
!   write(82,*) (real(2*pi/lx_tot*k/z_i*z),k=1,nxt/2-1)
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%COMMENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!! kx_2d=(2pi/lx_tot)*(0:nxt/2) => lx_tot is already non-dim by z_i
!! => k_x=[(2pi/lx_tot)*(0:nxt/2)/z_i] and
!! kz is given by = k_x*z => kz=[(2pi/lx_tot)*(0:nxt/2)*([1:nz]-0.5)*dz]
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%COMMENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!   write(82,*) (real(kx_2d(k,1)/z_i*z),k=1,nxt/2)
! Calculate spectrum now !! for each slice in MPI case !
!   call spectrum(u(:, :, jz), spectra_u(:,jz))
!   call spectrum(v(:, :, jz), spectra_v(:,jz))
!   call spectrum(w(:, :, jz), spectra_w(:,jz))
!   call spectrum(theta(:,:,jz),spectra_theta(:,jz))
    call spectrum(u(:, :, jz), spectra_uvwT(1,:,jz))
    call spectrum(v(:, :, jz), spectra_uvwT(2,:,jz))
    call spectrum(w(:, :, jz), spectra_uvwT(3,:,jz))
    call spectrum(theta(:, :, jz), spectra_uvwT(4,:,jz))

    ! Replace temperature spectra by Pollen concentration
    ! Chamecki - 08/10/2006
    IF (PCon_FLAG) call spectrum(PCon(:, :, jz,npcon), spectra_uvwT(4,:,jz))

  enddo
!   close(82)
!   print *,'spectra_sample_U',spectra_u(:,4)
!   print *,'spectra_sample_U2',spectra_uvwT(1,:,4)
!   print *,'spectra_sample_V',spectra_v(:,4)
!   print *,'spectra_sample_V2',spectra_uvwT(2,:,4)
!   print *,'spectra_sample_W',spectra_w(:,4)
!   print *,'spectra_sample_W2',spectra_uvwT(3,:,4)
!   print *,'spectra_sample_T',spectra_theta(:,4)
!   print *,'spectra_sample_T2',spectra_uvwT(4,:,4)
  !$if ($MPI)
  recvcounts = size (spectra_uvwT)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (spectra_uvwT(1, 1,1), size (spectra_uvwT), MPI_RPREC,&
     spectra_uvwT_tot(1, 1, 1), recvcounts, displs,       &
     MPI_RPREC, 0, comm, ierr)
  !$else
  !spectra_uvwT_tot=spectra_uvwT
  !$endif

  if((.not. USE_MPI) .or. (USE_MPI .and. coordz == 0)) then
    write(fname1,'(A,i6.6,A)')path//'output/spec_uvwT_',jt_local,'.bin'
    open(83,file=fname1,form='unformatted')

!write(fname1,'(a,i6.6,a)')path//'output/spec_u',jt_local,'.dat'
!open(1,file=fname1,form='formatted')
!write(fname2,'(A,i6.6,A)')path//'output/spec_v',jt_local,'.dat'
!open(2,file=fname2,form='formatted')
!write(fname3,'(A,i6.6,A)')path//'output/spec_w',jt_local,'.dat'
!open(3,file=fname3,form='formatted')
!write(fname4,'(A,i6.6,A)')path//'output/spec_t',jt_local,'.dat'
!open(4,file=fname4,form='formatted')

    write(83) real(spectra_uvwT_tot(:,1:nxt/2,:))
    close(83)
!do jz=1,nz
!   write(1,*)real(spectra_u(2:nxt/2,jz))
!   write(2,*)real(spectra_v(2:nxt/2,jz))
!   write(3,*)real(spectra_w(2:nxt/2,jz))
!   write(4,*)real(spectra_theta(2:nxt/2,jz))
!enddo
!close(1);close(2);close(3);close(4)
  end if
    ! ---
    end subroutine post_spec
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------    
    subroutine spectrum(u, spec)
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
    use fft
    implicit none
    ! ---
    real(kind=rprec), dimension(ldx, nynpy), intent(in) :: u
    real(kind=rprec), dimension(nxt/2), intent(out) :: spec  !--assumes Nyquist is 0
    ! ---
    integer :: jy, k 
    real(kind=rprec), dimension(nxt) :: vel_r, vel_c
    real(kind=rprec), dimension(nxt/2) :: tmp 

    integer*8, save :: plan
    logical, save :: init = .false.

    if (.not. init) then
        call dfftw_plan_r2r_1d(plan, nxt, vel_r, vel_c, FFTW_R2HC, FFTW_MEASURE)
        init = .true.
    end if

    ! initialize
    tmp(:) = 0._rprec
    do jy = 1, nynpy
        vel_r(:) = u(1:nxt,jy)/real(nxt,kind=rprec)
        ! check this normaliztion-part of forward
        ! call the fft
        call dfftw_execute_r2r(plan, vel_r, vel_c)
        ! compute magnitudes
        ! the 0.5 is the 1/2, all others are taken care of! (except maybe Nyquist)
        tmp(1) = tmp(1) + 0.5*vel_c(1)*vel_c(1)
    
        do k = 2, nxt/2
            tmp(k) = tmp(k) + vel_c(k)*vel_c(k) + vel_c(nxt+2-k)*vel_c(nxt+2-k)
            ! print *,'k,vel,spec',k,vel_c(k),tmp(k)
        end do
        !--assume Nyquist is 0
        !tmp(nxt/2+1)=tmp(nxt/2+1)+vel_c(nxt/2+1)*vel_c(nxt/2+1)
    end do
    tmp(:) = tmp(:)/real(nyt,kind=rprec) ! for average over nyt
    call mpi_allreduce(tmp, spec, nxt/2, mpi_rprec, mpi_sum, comm_ver, ierr)
    ! ---
    end subroutine spectrum  
!-----------------------------------------------------------------------
!   
!----------------------------------------------------------------------- 
!    subroutine timeseries_spec
!!-----------------------------------------------------------------------
!!   
!!-----------------------------------------------------------------------  
!    use sim_param,only:u,v,w,theta
!    implicit none
!    ! ---
!    integer :: jx, jy, jz, i
!    
!    if (mod(jt_total,time_spec)==0 .and. jt_total.gt.2000) then
!        jx = nxt/8
!        jy = nyt/2+1
!        jz = NZ/2
!!TSwrite(15)real(u(jx+nxt/24*2,jy,jz:jz+3)),real(u(jx+nxt/24*4,jy,jz:jz+3)),&
!!TS     real(u(jx+nxt/24*6,jy,jz:jz+3)),&
!!TS     real(v(jx+nxt/24*2,jy,jz:jz+3)),real(v(jx+nxt/24*4,jy,jz:jz+3)),&
!!TS     real(v(jx+nxt/24*6,jy,jz:jz+3)),&
!!TS     real(w(jx+nxt/24*2,jy,jz:jz+3)),real(w(jx+nxt/24*4,jy,jz:jz+3)),&
!!TS     real(w(jx+nxt/24*6,jy,jz:jz+3))
!    end if
!    ! ---
!    end subroutine timeseries_spec
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------   
    subroutine output_final(lun_opt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--lun_opt gives option for writing to different unit, and is used by
!  inflow_write
!--assumes the unit to write to (lun_default or lun_opt is already
!  open for sequential unformatted writing
!--this routine also closes the unit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    ! ---
    integer, intent(in), optional :: lun_opt  !--if present, write to unit lun
    ! ---
    integer, parameter :: lun_default = 17
    integer :: lun
    ! ---
    if (present (lun_opt)) then
        lun = lun_opt
    else
        lun = lun_default
    end if

    call checkpoint_final (lun)

    if (cumulative_time) then
    !--only do this for true final output, not intermediate recording
        open (1, file=fcumulative_time)
        write (1, *) nums
        close (1)
    end if
    ! ---
    end subroutine output_final
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
    subroutine checkpoint_final (fid_out)
    !use param, only : ldx,nyt,nz,nzt,theta_flag,PCon_FLAG
    use param
    use sim_param, only : u, v, w, RHSx, RHSy, RHSz, theta, pcon
    use sgsmodule, only : Cs_opt2, F_LM, F_MM, F_QN, F_NN, F_KX, F_XX,F_KX2,F_XX2
    use scalars_module, only : RHS_T,RHS_PCon,sgs_t3,deposition,Real_dep,Kc_t
    use bottombc, only : psi_m

    implicit none
    ! ---
    integer, intent (in) :: fid_out
    ! ---
    integer:: ipcon
    logical :: flag_dir
    integer(MPI_OFFSET_KIND) :: disp0, disp3d, disp2d
    character(100) :: cmd
    character(80) :: fname
    ! --- create directory for restart file
    inquire(directory=path_restart, exist=flag_dir)
    if (.not. flag_dir) then
        write(cmd,'(a)') 'mkdir '//path_restart
        call system(trim(cmd))
    end if
    ! ---
    disp3d = np  * sizeof(u(1:nxt, 1:nynpy, 1:nz-1))
    disp2d = npy * sizeof(u(1:nxt, 1:nynpy, 1))

    ! --- output velocity field
    write(fname, '(a,i8.8,a)') path_restart//'vel_tt', nums, '.out'
    disp0 = 0
    call write_file_mpi_3d(u(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
    disp0 = disp0+disp3d
    call write_file_mpi_3d(v(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
    disp0 = disp0+disp3d
    call write_file_mpi_3d(w(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
    disp0 = disp0+disp3d
    call write_file_mpi_3d(RHSx(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
    disp0 = disp0+disp3d
    call write_file_mpi_3d(RHSy(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
    disp0 = disp0+disp3d
    call write_file_mpi_3d(RHSz(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
    disp0 = disp0+disp3d
    call write_file_mpi_3d(Cs_opt2(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
    disp0 = disp0+disp3d
    call write_file_mpi_3d(F_LM(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
    disp0 = disp0+disp3d
    call write_file_mpi_3d(F_MM(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
    disp0 = disp0+disp3d
    call write_file_mpi_3d(F_QN(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
    disp0 = disp0+disp3d
    call write_file_mpi_3d(F_NN(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)

    ! --- output the temperature field
    if (theta_flag) then
        write(fname, '(a,i8.8,a)') path_restart//'temp_tt', nums, '.out'
        disp0 = 0
        call write_file_mpi_3d(theta(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
        disp0 = disp0+disp3d
        call write_file_mpi_3d(RHS_T(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
        disp0 = disp0+disp3d
        call write_file_mpi_xy(sgs_t3(1:nxt, 1:nynpy, 1), fname, disp0)
        disp0 = disp0+disp2d
        call write_file_mpi_xy(psi_m, fname, disp0)
    end if

    if (PCon_flag) then
        ! --- output the concentration
        write(fname, '(a,i8.8,a)') path_restart//'con_tt', nums, '.out'
        disp0 = 0
        do ipcon = 1, npcon
            call write_file_mpi_3d(Pcon(1:nxt, 1:nynpy, 1:nz-1, ipcon), fname, disp0)
            call write_file_mpi_3d(RHS_Pcon(1:nxt, 1:nynpy, 1:nz-1, ipcon), fname, disp0+npcon*disp3d)
            disp0 = disp0 + disp3d
        end do

        ! --- output related variables for concentration
        write(fname, '(a,i8.8,a)') path_restart//'con_diff_tt', nums, '.out'
        disp0 = 0
        call write_file_mpi_xy(Kc_t(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
        disp0 = disp0+disp3d
        call write_file_mpi_3d(F_KX(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
        disp0 = disp0+disp3d
        call write_file_mpi_3d(F_XX(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
        disp0 = disp0+disp3d
        call write_file_mpi_3d(F_KX2(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
        disp0 = disp0+disp3d
        call write_file_mpi_3d(F_XX2(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
        disp0 = disp0+disp3d
        call write_file_mpi_xy(deposition, fname, disp0)
        disp0 = disp0+disp2d
        call write_file_mpi_xy(Real_dep, fname, disp0)
    end if
    ! ---
    end subroutine checkpoint_final
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
! ---
end module io