!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!    
module param
!--------------------------------------------------------------------!
!  Purpose:
!       define the parameters contronl the simulation
!
!  Record of revisions:
!		date		programmer			description of change
!		====		==========			=====================
!
!  The control flag is flag_dynStress and flag_dynStokes
!--------------------------------------------------------------------!  
use types
use mpi
implicit none
save  
public
! --- Local File Setup
character(*), parameter :: path = '../'
character(*), parameter :: path_restart = '../restart/'

integer, parameter :: fid_param = 111
character(80), parameter :: fn_param = '../readin/param.nml'

! --- Debug
logical, parameter :: verbose = .false.     ! prints small stuff to screen
logical, parameter :: debug_bc = .true.     ! write lots of data/files

! --- Constants
real(rprec), parameter :: pi = 3.1415926535897932384626433_rprec
real(rprec), parameter :: vonk = .4_rprec
real(rprec), parameter :: g = 9.81_rprec
real(rprec), parameter :: na = 6.022_rprec*10._rprec**23._rprec
real(rprec), parameter :: r_gas = 8.3144621_rprec           ! ideal gas constant
real(rprec), parameter :: cd_ocean = 1.3_rprec*10._rprec**(-3._rprec)   ! drag coefficient of wind stress for wind velocity at 10m height
real(rprec), parameter :: rho_ocean = 1031._rprec           ! kg/m3
real(rprec), parameter :: rho_air = 1.225_rprec             ! kg/m3 and T=285.15K
real(rprec), parameter :: BOGUS = -1234567890._rprec
real(rprec), parameter :: TINYS = 1.e-30_rprec

! --- MPI Parameters
! --- Namelist variable list
namelist/mpi_param/ use_mpi, npy, npz
logical :: use_mpi = .false.            ! use mpi or not, always true
integer :: npy = int_missing          ! the number of processors in y direction of mpi
integer :: npz = int_missing          ! the number of processors in z direction of mpi
! --- Non-namelist variable list
integer :: status(mpi_status_size)
integer :: comm, comm_row, comm_col, comm_ver
integer :: np, global_rank, ierr
integer :: up, down, west, east, north, south
integer :: mpi_rprec, mpi_cprec
integer :: rank = int_missing
integer, dimension(2) :: dims
logical, dimension(2) :: periodic
integer :: coordy, coordz
! --- Di Yang added rank_of_coord, coord_of_rank
integer, allocatable, dimension(:,:) :: rank_of_coord
integer, allocatable, dimension(:) :: coord_of_rank
integer, allocatable, dimension(:) :: idy, idz, idy_ver, idz_col
integer, allocatable, dimension(:) :: rank_ver, rank_col
integer :: coord_ver, coord_col, ver_rank
! --- Output Setup
! --- Namelist variable list
namelist/output_control/ output, use_avgslice, base_time, average_dim_num, &
                         output_video, video_start, video_end, video_freq, &
                         flag_opt_3d, opt_3d_start, flag_opt_pre, flag_opt_eps,      & 
                         flag_opt_sgs, flag_opt_shear, flag_opt_frac,    & 
                         jx_opt_start, jy_opt_start, jz_opt_start, nxout, nyout, nzout
logical :: output = .true.
logical :: use_avgslice = .true.
integer :: base_time = 500          ! TOMAS CHOR included base_time in param.f90
integer :: average_dim_num = 2      ! if==2 then averages over x,y,t, if==1 then averages over y,t
logical :: output_video = .false.
integer :: video_start = int_missing
integer :: video_end = int_missing
integer :: video_freq = int_missing
logical :: flag_opt_3d = .true.     ! flag to output 3d field
integer :: opt_3d_start = 0         ! starting time step to output 3d field 
logical :: flag_opt_pre = .true.    ! flag to output pressure	
logical :: flag_opt_eps = .true.    ! flag to output dissipation
logical :: flag_opt_sgs = .false.   ! flag to output sgs stresses
logical :: flag_opt_shear = .false. ! flag to output horizontal shear gradient
logical :: flag_opt_frac = .false.  ! flag to output only a fraction of the 3d field
integer :: jx_opt_start = int_missing
integer :: jy_opt_start = int_missing
integer :: jz_opt_start = int_missing
integer :: nxout = int_missing
integer :: nyout = int_missing
integer :: nzout = int_missing

! --- used to output part of the domain data
integer, dimension(3) :: gsize_opt, lsize_opt, start_opt
integer :: bufsize

! --- Domain Parameters
! --- Namelist variable list
namelist/domain_param/ nxt, nyt, nzt, lx_tot, ly_tot, lz_tot, z_i
integer :: nxt = int_missing         ! the number of grid point of the whole domain in x-direction
integer :: nyt = int_missing         ! the number of grid point of the whole domain in y-direction
integer :: nzt = int_missing     ! the number of grid point of the whole domain in z-direction
real(rprec) :: lx_tot = float_missing   ! the length of the whole domain in x-direction
real(rprec) :: ly_tot = float_missing   ! the length of the whole domain in y-direction
real(rprec) :: lz_tot = float_missing   ! the length of the whole domain in z-direction, which is set to z_i conventionally
real(rprec) :: z_i = float_missing  ! the length scale for the simulation
! --- Non-namelist variable list
type(type_infoDm) :: info_dm        ! variable to save dimensional domain information
integer :: nz = int_missing         ! the grid point in sub-domain of each process in z-direction (our LES slice the whole domain only in z direction)
integer :: lhx                      ! the number of grid point for complex variable in x-direction
integer :: lhy                      ! the number of grid point for complex variable in y-direction 
integer :: ldx                      ! the number of grid point in x-direction after padding
integer :: ldy                      ! the number of grid point in y-direction after padding
integer :: nxt2                      ! the number of grid point for dealiasing in x-direction
integer :: nyt2                      ! the number of grid point for dealiasing in y-direction
integer :: lhx2                     ! the number of grid point for dealiasing complex variable in x-direction
integer :: ldx2                     ! the number of grid point dealiasing variable in x-direction after padding
integer :: lhy2                     ! the number of grid point for dealiasing complex variable in x-direction
integer :: ldy2                     ! the number of grid point dealiasing variable in x-direction after padding

integer :: nxnpy = int_missing         ! the grid point in sub-domain of each process in x-direction (our LES slice the whole domain only in x direction)
integer :: nynpy = int_missing         ! the grid point in sub-domain of each process in y-direction (our LES slice the whole domain only in x direction)
integer :: nxhnpy = int_missing
integer :: nxh2npy= int_missing 
integer :: nx2npy = int_missing
integer :: ny2npy = int_missing

real(rprec) :: ly = float_missing   ! then length of sub-domain of each process in y-direction
real(rprec) :: lz = float_missing   ! then length of sub-domain of each process in z-direction
real(rprec) :: dx, dy, dz           ! the spatial step in x,y,z-direction

real(rprec), allocatable, dimension(:) :: x, y, zuvp, zw       ! meshgrid variable
! --- Time Parameters
! --- Namelist variable list
namelist/time_param/ nsteps, dt, c_count, p_count, cs_count, &
                     flag_restart, nt_restart
integer :: nsteps = int_missing     ! the total number of time step of time integration
real(rprec) :: dt = float_missing   ! the time step
integer :: c_count = int_missing
integer :: p_count = int_missing
integer :: cs_count = int_missing   ! 
logical :: flag_restart = .false.   ! whether or not save the restart file (just in case unexpect termination of the job)
integer :: nt_restart = int_missing ! save the restart file for each nt_restart step

! --- Non-namelist variable list
type(type_infoTime) :: info_time    ! variable to save dimensional time information
real(rprec) :: tt                   ! physical time
integer :: jt                       ! the current time step
integer :: jt_total                 ! the total time step
integer :: nums
real(rprec), parameter :: tadv1 = 1.5_rprec         ! time advance parameters (AB2)
real(rprec), parameter :: tadv2 = 1._rprec - tadv1  ! 

! --- Flow Parameters
! --- Namelist variable list
namelist/flow_param/ u_star, u_scale, Pr, restart_vel, inilag, interp, DYN_init,        &
                     coriolis_forcing, freq_coriolis, ug_dim, vg_dim,                   &
                     nu_molec, walltype, molec, sgs, dns_bc,                            &
                     model, calc_lag_Ssim, Co, cs, ifilter,                             &
                     inflow, l_fringe,                                                  &
                     face_avg, read_inflow_file, write_inflow_file, jt_start_write,     &
                     ubc, damping_method, lbc_mom, ubc_mom,                             &
                     force_top_bot, use_mean_p_force, mean_p_force, use_force_angle,    &
                     force_angle, ocean_flag, flag_relaxForce,                          &
                     read_ic_profile, intp_vel

real(rprec) :: u_star  = float_missing   ! wind stress friction velocity
real(rprec) :: u_scale = float_missing   ! velocity scale
real(rprec) :: Pr = 0.4_rprec           ! turbulentg Prandtl number
logical :: restart_vel = .false.              ! flag to initialize velocity field
logical :: inilag = .true.              ! flag to initialize eddy viscosity
logical :: interp = .false.             ! 
integer :: DYN_init = int_missing       ! the time step to start using dynamical SGS model
logical :: coriolis_forcing = .false.   ! flag to use Coriolis forcing (always false for ocean_flag)
real(rprec) :: freq_coriolis = float_missing        ! Coriolis frequency
real(rprec) :: ug_dim = float_missing, vg_dim = float_missing   ! dimensional geostrophic wind in x,y-direction
real(rprec) :: nu_molec = 1.14e-5_rprec ! molecular viscosity
character(20) :: walltype = 'rough'     ! smooth or rough wall
logical(rprec) :: molec = .false.       ! 
logical(rprec) :: sgs = .true.
logical(rprec) :: dns_bc = .false.
integer :: model = 5                    ! subgrid-scale model type: 1->Smagorinsky; 2->Dynamic; 3->Scale depent; 4->Lagrangian scale-sim; 5->Lagrangian scale-dep
logical :: calc_lag_Ssim = .false.
real(rprec) :: Co = 0.16_rprec
real(rprec) :: cs = 0.2_rprec           ! Smagorinsky's constant
integer :: ifilter = 1                  ! test filter type: 1->cut off; 2->Gaussian; 3->Top-hat
logical :: inflow = .false.             ! flag to use prescribed inflow
real(rprec) :: l_fringe = float_missing ! length of blending zone (recycling)
real(rprec) :: face_avg = float_missing
logical :: read_inflow_file = .false.
logical :: write_inflow_file = .false.
integer :: jt_start_write = int_missing
integer :: ubc = 1
integer :: damping_method = 1
character(20) :: lbc_mom = 'wall'               ! 'wall' or 'stress free' (not used in ocean simulation)
character(20) :: ubc_mom = 'wall'               ! 'wall' or 'stress free' (for shallow or deep ocean simulation)
logical :: force_top_bot = .false.
logical :: use_mean_p_force = .false.
real(rprec) :: mean_p_force = float_missing     ! the mean pressure gradient acceleration: -1/rho*dp/dx
logical :: use_force_angle = .false.            ! Constant angle for the pressure forcing
real(rprec) :: force_angle = 0._rprec
logical :: ocean_flag = .false.                 ! flag to use ocean setup
logical :: flag_relaxForce = .false.            ! flag to use relaxation forcing (cellular flow forcing)
logical :: read_ic_profile = .false.            ! flag to initialize field with input profile as to minimize inertial oscillation
logical :: intp_vel = .false.                   ! interpolate the velocity field from coarse simulation

! --- Fringe method
integer :: jx_fringe, jx_relax
integer, dimension(3) :: gsize_fringe, lsize_fringe, start_fringe
real(rprec) :: lx_s, lx_pl
real(rprec), dimension(:), allocatable :: wgt 

! --- 
real(rprec), dimension(:, :, :), allocatable :: coef_ap, coef_bp, coef_cp
integer :: izs, ize
! --- canopy parameters
namelist/canopy_param/ flag_canopy, canopy_conf, use_field_scale, use_plant_scale,  &
                       lwest, least, lsouth, lnorth, LBB, LGL, &
                       canopy_dir, Cd_const, vel_exp 
logical :: flag_canopy = .true.
character(20) :: canopy_conf = 'block'  ! canopy configuration, "block" and "stripe"
logical :: use_field_scale = .true.
logical :: use_plant_scale = .false.
real(rprec) :: lwest  = float_missing        ! the upstream domain size of the canopy
real(rprec) :: least  = float_missing        ! the downstream domain size of the canopy
real(rprec) :: lsouth = float_missing        ! the south domain size of the canopy
real(rprec) :: lnorth = float_missing        ! the north domain size of the canopy
real(rprec) :: LBB = float_missing            ! the spacing between adjacent backbones
real(rprec) :: LGL = float_missing            ! the length of the growth line
real(rprec) :: canopy_dir = float_missing    ! orientation (degree) of the canopy relative to x-direction
real(rprec) :: Cd_const = float_missing
real(rprec) :: vel_exp = float_missing              ! velocity exponent

! --- Leaf Parameters
integer, parameter :: fid_leaf = 263
character(*), parameter :: fn_leaf = '../readin/leaf.nml'
namelist /leaf_param/ nc, hgt_leaf, lad_leaf
integer :: nc
real(rprec), dimension(100) :: hgt_leaf, lad_leaf

! --- Non-namelist variable list
real(rprec), dimension(:), allocatable :: ha
real(rprec), dimension(:), allocatable :: aleaf                         ! aleaf - leaf area index density
real(rprec), dimension(:, :, :), allocatable :: aleaf_x, aleaf_y
real(rprec), dimension(:,:,:), allocatable :: a_leafz_uv, a_leafz_w
real(rprec), dimension(:, :, :, :), allocatable :: SI, SS
real(rprec), dimension(:, :, :), allocatable :: Cd_uv, Cd_w
real(rprec), dimension(:, :, :), allocatable :: Cd_vex, Cd_cave
real(rprec), dimension(:, :, :), allocatable :: St, EI                  ! St - Stokes number, EI - Impaction coefficient
real(rprec) :: lai_data

real(rprec), parameter :: Pxx = 0.5_rprec	
real(rprec), parameter :: Pyy = 0.5_rprec	
real(rprec), parameter :: Pzz = 0.5_rprec

! --- Non-namelist variable list
real(rprec) :: ug = float_missing, vg = float_missing
real(rprec) :: coriol = float_missing
integer :: nnn = 2
integer :: BETA_FLAG = 1

! --- Temperature Parameters
! --- Namelist variable list
namelist/temperature_param/ theta_flag, restart_theta, spatial_flux_flag, OB_homog_flag, &
                            WALL_HOMOG_FLAG, model_sc, theta_init_time, inilag_sc, lbc, patch_flag, &
                            remote_flag, remote_homog_flag, remote_flux_homog_flag, &
                            remote_to_patch_flag, &
                            inv_strength, dTdz_top, theta_top, T_scale, wt_s, T_init, &
                            flag_dynWT, ttshift_dynWT, intp_theta
logical :: theta_flag = .false.     ! flag to use temperature field of not
logical :: restart_theta = .false.     ! flag to initialize concentration field
logical :: spatial_flux_flag = .false.  ! 
logical :: OB_homog_flag = .false.      !
logical :: WALL_HOMOG_FLAG = .false.    ! 
integer :: model_sc = 1                 ! subgrid-scale model type (scalar): 1->static prandtl; 2->Dynamic;
integer :: theta_init_time = int_missing      ! 
logical :: inilag_sc = .false.
integer :: lbc = 1
integer :: patch_flag = 1
integer :: remote_flag = 0
integer :: remote_homog_flag = 0
integer :: remote_flux_homog_flag = 0
logical :: remote_to_patch_flag = .false.
real(rprec) :: inv_strength = float_missing     ! the temperature gradient in inverse layer
real(rprec) :: dTdz_top = float_missing
real(rprec) :: theta_top = float_missing
real(rprec) :: T_scale = float_missing
real(rprec) :: wt_s = float_missing
real(rprec) :: T_init = float_missing
logical :: flag_dynWT = .false.                 ! dynamical kinematic heat flux
integer :: ttshift_dynWT = 0                    ! the time-step shift with respect to jt_total (ttshift_dynWT+jt_total)

logical :: intp_theta = .false.              ! interpolate the theta field from coarse simulation

! --- Non-namelist variable list
integer, parameter :: fid_wt = 260              ! the file id of fn_wt
character(*), parameter :: fn_wt = '../readin/wt.in'      ! the file name of kinematic heat flux when flag_dynWT is true
integer :: len_wt                               ! the length of wt_s

! --- Ocean Parameters (Added by Di Yang, Bicheng Chen add agl_Ust), Coriolis is always true
! --- Add rad_Ust to make agl_Ust more clear (Bicheng Chen 10/20/2016)
! --- Namelist variable list
namelist/ocean_param/ nu_ek, nu_sgs, alpha_w, prop_mixed,                                   &
                      flag_dynStress, agl_stress, stokes_flag, flag_dynStokes, flag_PMspec, &
                      amp_w, lambda_w, stokes_action, agl_Ust, Lat_dyn,                     &
                      flag_crossVel, udim_cross, vdim_cross
real(rprec) :: nu_ek = 1.16e-2_rprec                ! the eddy viscosity in Ekman layer (used only in initialization)
real(rprec) :: nu_sgs = 1.0e-3_rprec                ! 
real(rprec) :: alpha_w = 2.0e-4_rprec               ! the thermal expansion rate of water
real(rprec) :: prop_mixed = 1._rprec/3._rprec       ! the proportion of the mixed layer over the whole depth
logical     :: flag_dynStress = .true.              ! flag to use dynamics wind stress
integer     :: ttshift_dynStress = 0                ! the time-step shift with respect to jt_total (ttshift_dynStress+jt_total)
real(rprec) :: agl_stress = 0._rprec                ! the inclination angle of wind stress (w.r.t. x-direction)
logical     :: stokes_flag = .true.                 ! the flag to use Stokes drift velocity (due to wave)
logical     :: flag_dynStokes = .true.              ! flag to use dynamics Stokes
logical     :: flag_PMspec = .true.                 ! flag to use Pierson-Moskowitz Spectrum (Pierson and Moskowith, 1964) (only works while flag_PMspec is true)
real(rprec) :: amp_w = float_missing                ! the amplitude of the wave
real(rprec) :: lambda_w = float_missing             ! the wave length of the wave
character(50) :: stokes_action = 'donelan.pierson1987'
real(rprec) :: agl_Ust = float_missing              ! the inclination angle of wave propagation direction (w.r.t.x-direction)
real(rprec) :: Lat_dyn = 0.3_rprec
logical     :: flag_crossVel = .false.              ! cross velocity for the all domain
real(rprec) :: udim_cross = float_missing           ! dimensional cross velocity in x direction
real(rprec) :: vdim_cross = float_missing           ! dimensional cross velocity in y direction

! --- Non-namelist variable list
real(rprec) :: rad_stress = float_missing           ! the radians of inclination of the wind stress              ! 
real(rprec) :: u_cross = float_missing              ! the nondimensional cross velocity in x-direction
real(rprec) :: v_cross = float_missing              ! the nondimensional cross velocity in y-direction

integer, parameter :: fid_ustar = 258
character(*), parameter :: fn_ustar = '../readin/ustar.in'    ! the file name of friction velocity time series when flag_dynStress is true
real(rprec) :: ustar_dyn
integer :: len_ustar

integer, parameter :: fid_gwave = 259
character(*), parameter :: fn_gwave = '../readin/gwave.in'
integer :: len_gwave

! --- Concentration Parameters
! --- Namelist variable list
namelist/con_param/ pcon_flag, restart_pcon, n_con, pcon_scale, &
                    unit_pconScale, pcon_scheme, periodicbcx, periodicbcy, pcon_init, lbcp, &
                    pcon_sfc, models, sc, sc_boun, settling, pcon_acc, inflow_pcon, &
                    source_type, ini_src, end_src, active_pcon, pconsgs_acc_flag, model_psgsacc, cs_pa,  &
                    flag_endless, flag_uptake
logical :: pcon_flag = .false.              ! the flag to evolve concentration field(s)
logical :: restart_pcon = .false.               ! the flag to initialize the concentration field(s)
integer, target :: n_con                    ! the number of the concentration fields
integer :: pcon_init = 1                    ! the time step to initialize the concentration field
real(rprec) :: pcon_scale = 1._rprec        ! the scale of the concentration
character(20) :: unit_pconScale = 'NONE'    ! the unit of the pcon_scale (just a remind for myself mostly)
integer :: pcon_scheme = 3                  ! the scheme for concentration field
logical :: periodicbcx = .false.            ! the flag to use periodic boundary condition for concentration field in x-direction
logical :: periodicbcy = .false.            
integer :: lbcp = 1
real(rprec) :: pcon_sfc = 1._rprec
real(rprec) :: sfc_flux = 0._rprec
integer :: models = 1                       ! subgrid-scale model type (scalar): 1->static prandtl; 2->Dynamic;
real(rprec) :: sc = 0.8_rprec, sc_boun = 0.95_rprec
logical :: settling = .false.
logical :: pcon_acc = .false.
logical :: inflow_pcon = .false.
character(50) :: source_type = 'point'      ! source type: 'point', 'line', 'planar', and 'volume' 
integer :: ini_src = 1                      ! the initial time step to set up the source
integer :: end_src = 5000000                ! the final time step for source
logical :: active_pcon = .false.            ! flag to treat concentration field as active scalar
logical :: pconsgs_acc_flag = .false.
integer :: model_psgsacc = 1
real(rprec) :: cs_pa = 1._rprec*10._rprec**4._rprec
logical :: flag_endless = .false.
logical :: flag_uptake  = .false.           ! flag to determine if canopy uptake nutrients

type(type_infoCon) :: info_con
integer :: npcon = -1

! --- Particle Parameters
! --- Namelist variable list
namelist/particle_param/ vel_settling, ratio_dens, vol_spec, pcon_sfc_flux
real(rprec), dimension(:), allocatable, target :: vel_settling  ! settling(rise/slip) velocity of particles
real(rprec), dimension(:), allocatable, target :: ratio_dens    ! the density ratio between particl and fluid (water/air)
real(rprec), dimension(:), allocatable, target :: vol_spec      ! the specific volume (the reciprocal of the density)
real(rprec), dimension(:), allocatable, target :: pcon_sfc_flux
! --- Non-namelist variable list
real(rprec), dimension(:), allocatable :: settling_vel          ! nondimenional settling(rise/slip) velocity used in simulation
real(rprec), dimension(:), allocatable :: densratio_pcon        ! same thing for densratio_pcon and v_pcon
real(rprec), dimension(:), allocatable :: v_pcon
real(rprec), dimension(:), allocatable :: pcon_sfc_flux_nd
! --- Source Parameters
! --- fixed size for ix_src ..., need to change in the future
! --- Namelist variable list
namelist/number_source/ n_src
namelist/source_alignment/ src_align
namelist/source/ ix_src, iy_src, iz_src, x_src, y_src, z_src, row_src, col_src, rls_src
integer :: n_src
character(20) :: src_align                                      ! align in streamwise or lateral directions                                                ! the number of the source for each compound
integer, dimension(500) :: ix_src, iy_src, iz_src               ! the grid point of the source in x, y, z-direction
integer, dimension(500) :: x_src, y_src, z_src                  ! the coordinates of the source in x, y, z-direction
integer, dimension(500) :: row_src, col_src                     ! the row (y-direction) and the column (x-direction) of the source if ENDLESS is applied
real(rprec), dimension(500) :: rls_src                          ! the source strength for point source

! --- Non-namelist variable list
integer, parameter :: fid_src = 121                             ! the file id for source file
character(*), parameter :: fnPre_src = '../readin/source_'      ! the prefix of the file name for the source file (each compound has each own file)
character(*), parameter :: fnSuf_src = '.nml'                   ! the suffix  of the file name
type(src_rls), dimension(:), allocatable :: con_src             ! the variable to save source information

 !### ENDLESS Parameters ###
 !## Namelist variable list
 ! dm_endless - the variable to save the domain information for ENDLESS
 ! info_endless - the variable to save the infomation of ENDLESS
 ! flag_cell - the flag to use cellular flow for advection
namelist /endless_param/ dm_endless, info_endless, flag_cell
type(rp_dm) :: dm_endless
type(type_infoEndless) :: info_endless
logical :: flag_cell = .false.
 !## Non-namelist variable list
 ! map_conLoc - the varialbe to save the mapping relation (boundary relation)
 ! >>between different scalar fields
 ! usage_con - the usage information of the scalar fields when dynamical domain
 ! >>scheme is applied
 ! conMax - maximal concentration of each scalar field in sub-domain in each process
 ! conMax_all - maximal concentration of each scalar field over all process
type(type_map) :: map_conLoc
type(type_usageCon) :: usage_con
type(type_conMax) :: conMax, conMax_all
 !### END ENDLESS Parameters ###

 !### Cellular FLow Parameters ###
 !## Namelist variable list
 ! nmod_cell - the number of the mode of the cellular flow
 ! u0_cell - the velocity scale for cellular flow
 ! advU_cell - the phase velocity of the cellular flow in x-direction
 ! advV_cell - the phase velocity of the cellular flow in y-direction
 ! tt0_cell - the displacement time step of the phase of the cellular flow
 ! dph_cell - the depth of the cellular (not applied yet)
 ! d_cell - the wavelength of the cellular flow
 ! piPhix_cell - the phase diplacement of cellular flow in x-direction
 ! >>(normalized by pi)
 ! piPhiy_cell - the phase diplacement of cellular flow in y-direction
 ! >>(normalized by pi)
namelist /mode_cell/ nmod_cell
namelist /cellular_flow/ u0_cell, advU_cell, advV_cell, tt0_cell,&
   dph_cell, d_cell, piPhix_cell, piPhiy_cell
integer :: nmod_cell = int_missing
real(rprec), dimension(:), allocatable :: u0_cell
real(rprec), dimension(:), allocatable :: advU_cell, advV_cell
real(rprec) :: tt0_cell = float_missing
real(rprec), dimension(:), allocatable :: dph_cell
real(rprec), dimension(:), allocatable :: d_cell
real(rprec), dimension(:), allocatable :: piPhix_cell, piPhiy_cell
 !## Non-namelist variable list
type(type_paramCell) :: param_cell
 !### END Cellular FLow Parameters ###

 !### Relaxing Force Parameters ###
 !### >>Use cellular flow forcing
 !## Namelist variable list
 ! wn_cell - wave number of cellular flow
 ! type_filter - filter type, can be 'square', 'diagonal'
 ! t_relax - the relaxation time
namelist /mode_relaxForce/ nmod_cell
namelist /relax_force/ u0_cell, advU_cell, advV_cell, tt0_cell, dph_cell,&
   wn_cell, piPhix_cell, piPhiy_cell, wn_filter, type_filter, t_relax
integer, dimension(:), allocatable :: wn_cell
integer :: wn_filter
character(10) :: type_filter = 'square'
real(rprec) :: t_relax
 !## Non-namelist variable list
 ! the relaxation force in x-direction
 ! the relaxation force in y-direction
real(rprec), dimension(:, :, :), allocatable :: relaxForce_x, relaxForce_y
 !### END Relaxing Force Parameters ###

! --- Chemistry Parameters
! --- Namelist variable list
namelist/chem_param/ flag_chem, flag_bgCon, flag_mutReact, n_bgCon, &
                     fn_bgConProf, fn_bgChemRate, theta_bg, p_bgAtmos
logical :: flag_chem = .false.                  ! 
logical :: flag_bgCon = .false.
logical :: flag_mutReact = .false.
integer :: n_bgCon = int_missing
character(50) :: fn_bgConProf = 'conprof.bg'
character(50) :: fn_bgChemRate = 'chemrate.bg'
real(rprec) :: theta_bg = 298._rprec
real(rprec) :: p_bgAtmos = 1.01325_rprec*10._rprec**5._rprec
 
! --- Reorganize in the Furture
!!!!! --- leaf area index
!!!!real(rprec), parameter :: LAIp = 6.0
!!!!real(rprec), parameter :: LAIw = 2.9
!!!!! --- leaf area density
!!!!real(rprec), parameter :: a_leaf1 = 0.99
!!!!real(rprec), parameter :: a_leaf2 = 2.17
!!!!real(rprec), parameter :: a_leaf3 = 3.93
!!!!real(rprec), parameter :: a_leaf4 = 4.32
!!!!real(rprec), parameter :: a_leaf5 = 4.68
!!!!real(rprec), parameter :: a_leaf6 = 4.35
!!!!real(rprec), parameter :: a_leaf7 = 0.79
!!!!! --- height corresponding to leaf area density
!!!!real(rprec), parameter :: h_leaf1 = float_missing
!!!!real(rprec), parameter :: h_leaf2 = float_missing
!!!!real(rprec), parameter :: h_leaf3 = float_missing
!!!!real(rprec), parameter :: h_leaf4 = float_missing
!!!!real(rprec), parameter :: h_leaf5 = float_missing
!!!!real(rprec), parameter :: h_leaf6 = float_missing
!!!!real(rprec), parameter :: h_leaf7 = float_missing
integer, parameter :: diurnal_forcing_flag = 0, no_days = 1
logical, parameter :: jan_diurnal_run = .false., ug_diurnal_variation = .false.
logical, parameter :: GABLS_diurnal_test = .FALSE.
logical, parameter :: passive_scalar = .false., GABLS_test = .false.
logical, parameter :: test_phase = .FALSE., vec_map = .FALSE., smag_sc = .FALSE.
logical, parameter :: check_dt = .TRUE.
integer, parameter :: stencil_pts = 4
logical, parameter :: coarse_grain_flag = .FALSE.
real(rprec), parameter :: q_s1 = 11._rprec, q_s2 = 13._rprec
real(rprec), parameter :: q_mix = 12._rprec, q_top = 12._rprec, wq_s = 0.06
real(rprec), parameter::cap_thick = 80._rprec, z_decay = 1._rprec
! --- For the field experiment rag06
logical, parameter :: rag06 = .false.
logical, parameter :: fieldsize = .false.
logical, parameter :: fieldsize_stripe = .false.
! --- Starting and ending grid points of the field
! --- for fieldsize_stripe=.TRUE., only field_xs and field_xe are used
integer, parameter :: field_xs = 9 ! Starting x grid point
integer, parameter :: field_xe = 11 ! Ending x grid point
integer, parameter :: field_ys = 49 ! Starting y grid point
integer, parameter :: field_ye = 51 ! Ending y grid point
real(rprec), parameter :: V_pcon0 = 1._rprec/859.8708_rprec
! ---
end module param
!--------------------------------------------------------------------!
!                                                
!--------------------------------------------------------------------!


    