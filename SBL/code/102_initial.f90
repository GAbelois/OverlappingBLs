!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine init()
!--------------------------------------------------------------------!
!   initialize everything to make main program clean                                           
!--------------------------------------------------------------------!    
use param, only: theta_flag, pcon_flag, model, ubc, sgs, comm, ierr
use sim_param, only: u, v, w
use io, only: screen_display
use intermediate
use bottombc, only: patches
use fft
use test_filtermodule
use topbc, only: sponge, setsponge
implicit none
! ---
call init_namelist()            ! read all pararmeters from namelist 
! ---
call init_parameter()           ! domain variables initialization
! ---
call init_nondimensional()      ! normalize the variable
! ---
call allocate_flow_variable()   ! flow variables initialization
call allocate_output_variable() ! output variables initialization
! ---
call init_vel_field()
call init_temperature_field()
call init_concentration_field()
call data_exchange()            ! exchange data between neighboring zones
! ---
call patches()                  ! Initialize surface physical conditions
! --- formulate the fft plans
call init_fft_plan()            ! create fft plans and initialize the kx, ky arrays
! --- 
!call openfiles()                ! open output files
! --- initialize test filter
if ( sgs ) then
    call init_test_filter(2._rprec * filter_size, G_test)

    if (model == 3 .or. model == 5) then  !--scale dependent dynamic
        call init_test_filter(4._rprec * filter_size, G_test_test)
    end if
end if
! --- define the upper boundary condition (sponge damping layer)
if (ubc == 1) then
    call setsponge()
else
    sponge = 0._rprec
end if
! ---
call init_lad()
call init_coef_press()
! --- 
call filter_data(u)
call filter_data(v)
call filter_data(w)
! ---
! call screen_display()           ! echo the simulation parameters
! ---
end subroutine init
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine init_namelist()
!--------------------------------------------------------------------!
!     Read all the pararmeters in namelist                                        
!--------------------------------------------------------------------!
use param
implicit none
! ---
open(fid_param, file=fn_param, form='formatted', status='old')
read (fid_param, nml=domain_param)
read(fid_param, nml=output_control)
read (fid_param, nml=time_param)
read (fid_param, nml=flow_param)
read (fid_param, nml=canopy_param)
read (fid_param, nml=temperature_param)

read (fid_param, nml=ocean_param)  
read (fid_param, nml=con_param)

if (PCon_flag) then
    info_con%n_con => n_con
    allocate (vel_settling(n_con))
    allocate (ratio_dens(n_con))
    allocate (vol_spec(n_con))
    allocate (pcon_sfc_flux(n_con))
    read(fid_param, nml=particle_param)
    info_con%vel_settling => vel_settling
    info_con%ratio_dens => ratio_dens
    info_con%vol_spec => vol_spec
    info_con%pcon_sfc_flux => pcon_sfc_flux
end if
close (fid_param)

if ( flag_canopy ) then
    open(fid_leaf, file=fn_leaf, form='formatted', status='old')
    read(fid_leaf, nml=leaf_param)
    close(fid_leaf)
end if
! ---
end subroutine init_namelist
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine init_parameter()
!-----------------------------------------------------------------------
!   domain variables initialization (Bicheng Chen 07/02/2015)  
!-----------------------------------------------------------------------
use param
use scalars_module, only:hetero_count_out
use io, only:gsize, lsize, start
use intermediate, only:fr
implicit none
! ---
integer :: ind
! --- 
nz   = (nzt - 1)/npz + 1

nxt2  = 3*nxt/2
nyt2  = 3*nyt/2
lhx  = nxt/2 + 1
lhy  = nyt/2 + 1
ldx  = 2*lhx
ldy  = 2*lhy
lhx2 = nxt2/2 + 1
ldx2 = 2*lhx2
lhy2 = nyt2/2 + 1
ldy2 = 2*lhy2

nxnpy  = nxt / npy
nynpy  = nyt / npy
nx2npy = nxt2/ npy
ny2npy = nyt2/ npy
nxhnpy = nxnpy/2 + 1
nxh2npy= nx2npy/2 + 1

dx   = lx_tot / nxt
dy   = ly_tot / nyt
dz   = lz_tot / (nzt - 1)

ly   = ly_tot / npy 
lz   = lz_tot / npz

! --- initialize other time variables
hetero_count_out = p_count
! ---
! if (flag_restart) then
!     info_time%flag_restart = flag_restart
!     info_time%nt_restart = nt_restart
! end if
! ---
gsize(1) = nxt
gsize(2) = nyt
gsize(3) = nzt-1
lsize(1) = nxt
lsize(2) = nynpy
lsize(3) = nz-1
start(1) = 0
start(2) = nynpy * coordy
start(3) = (nz-1) * coordz
! ---
if ( mod(nxt2, npy) /= 0 .and. mod(nyt2, npy) /= 0 ) then
    write(*, *) 'Grid numbers and number of processors in y are not compatible'
    call mpi_finalize(ierr)
    stop
end if

fr=(1._rprec/real(p_count,kind=rprec))*real(c_count,kind=rprec)
! ---
end subroutine init_parameter
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine init_nondimensional()
!-----------------------------------------------------------------------
!   normalize the variables
!-----------------------------------------------------------------------
use types, only:rprec
use param
use stokes_drift
use intermediate
implicit none
! --- 
ly = ly / z_i
lz = lz / z_i

lx_tot = lx_tot/z_i
ly_tot = ly_tot/z_i
lz_tot = lz_tot/z_i

dx = dx/z_i
dy = dy/z_i
dz = dz/z_i

call init_meshgrid()
! --- 
dt = dt * u_scale / z_i
u_star = u_star / u_scale
! --- canopy normalization
if ( flag_canopy ) then
    lwest  = lwest / z_i
    least  = least / z_i
    lsouth = lsouth / z_i
    lnorth = lnorth / z_i
    LBB = LBB / z_i
    LGL = LGL / z_i
end if
! --- Coriolis effect
if (coriolis_forcing) then
    coriol = freq_coriolis*z_i/u_scale
    ug = ug_dim/u_scale
    vg = vg_dim/u_scale
end if
! --- Mean pressure gradient force
if (use_mean_p_force) then
    if (mean_p_force == float_missing) then
        mean_p_force = 1._rprec / lz_tot
    else
        mean_p_force = mean_p_force / u_scale**2 * z_i
    end if
else 
    mean_p_force = 0._rprec
end if

! ---
if (ocean_flag .and. theta_flag) then
    alpha_w = alpha_w * T_scale     ! normalize thermal expansion rate
end if

! --- Ocean flag
call init_stokes_drift

if (ocean_flag) then
    coriol = freq_coriolis*z_i/u_scale
    ug = ug_dim/u_scale
    vg = vg_dim/u_scale
    !- if use dynamical stress, read friction velocity and its angle
    if (flag_dynStress) then
        inquire (iolength=len_ustar) ustar_dyn, agl_stress
        open (unit=fid_ustar, file=fn_ustar, access='direct', recl=len_ustar)
        read (unit=fid_ustar, rec=1) ustar_dyn, agl_stress
        close (fid_ustar)
        ustar_dyn = ustar_dyn / u_scale
    end if 
    rad_stress = agl_stress/180._rprec*pi
    
    ! >>Cross flow over all domain
    if (flag_crossVel) then
        u_cross = udim_cross / u_scale
        v_cross = vdim_cross / u_scale
    end if
end if
! ---
end subroutine init_nondimensional
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine init_meshgrid()
!-----------------------------------------------------------------------
!    
!-----------------------------------------------------------------------
use param
implicit none
! ---
integer :: i, jx, jy, jz
integer :: ipy1, ipy2, ipz1, ipz2
integer :: jy_opt1, jy_opt2, jz_opt1, jz_opt2
integer :: flagy, flagz
real(rprec), dimension(:, :, :), allocatable :: temp
! ---
allocate(x(ldx), y(nynpy), zuvp(0:nz), zw(0:nz))
! --- 
do jz = 0, nz
    zw(jz)   = (coordz*(nz-1) + jz - 1) * dz
    zuvp(jz) = (coordz*(nz-1) + jz - 0.5_rprec) * dz
end do

do jx = 1, ldx 
    x(jx) = (jx - 1) * dx
end do

do jy = 1, nynpy
    y(jy) = (coordy*nynpy + jy - 1) * dy
end do
! --- for fringe method
if ( inflow ) then
    l_fringe = l_fringe / z_i
    ! do jx = 2, nxt      ! errors caused by round-up could occur, made some special treatment here
    !     if ( x(jx) .lt. lx_tot - l_fringe .and. x(jx+1) .ge. lx_tot - l_fringe ) then
    !         jx_fringe = jx - 1
    !         jx_relax  = nxt - jx_fringe
    !     end if
    ! end do
    jx_relax  = int(l_fringe/dx)
    jx_fringe = nxt - jx_relax      ! jx_relax is readin from the param.nml

    gsize_fringe(1) = jx_relax
    gsize_fringe(2) = nyt
    gsize_fringe(3) = nzt-1
    lsize_fringe(1) = jx_relax
    lsize_fringe(2) = nynpy
    lsize_fringe(3) = nz-1
    start_fringe(1) = 0
    start_fringe(2) = nynpy * coordy
    start_fringe(3) = (nz-1) * coordz

    if ( read_inflow_file ) then
        allocate(wgt(jx_relax))
        lx_s  = lx_tot - l_fringe       ! Stevens et al. (2014) Renewable Energy
        lx_pl = lx_tot - 0.25_rprec*l_fringe
    
        do i = 1, jx_relax
            jx = jx_fringe + i
            if ( x(jx) .lt. lx_pl ) then
                wgt(i) = 0.5_rprec*( 1._rprec - cos(pi*(x(jx) - lx_s) / (lx_pl - lx_s)) )
            else
                wgt(i) = 1._rprec    
            end if
        end do
    end if

    ! if ( rank == 0 ) then
    !     write(*,*) "jx_relax and jx_fringe are:", jx_relax, jx_fringe
    ! end if
end if
! ---
if ( flag_opt_3d .and. flag_opt_frac ) then
    ipy1 = (jy_opt_start-1)/nynpy;  ipy2 = (jy_opt_start-1+nyout-1)/nynpy
    ipz1 = (jz_opt_start-1)/(nz-1); ipz2 = (jz_opt_start-1+nzout-1)/(nz-1)
    jy_opt1 = jy_opt_start - ipy1*nynpy;  jy_opt2 = (jy_opt_start+nyout-1) - ipy2*nynpy; 
    jz_opt1 = jz_opt_start - ipz1*(nz-1); jz_opt2 = (jz_opt_start+nzout-1) - ipz2*(nz-1)
    
    if (ipy1 == ipy2) then
        flagy = 1
    else
        flagy = 0
    end if

    if (ipz1 == ipz2) then
        flagz = 1
    else
        flagz = 0
    end if
    ! ---
    if (coordz == ipz1) then
        lsize_opt(3) = flagz*jz_opt2 + (1-flagz)*(nz-1) - jz_opt1 + 1
        start_opt(3) = jz_opt1 - 1 + (nz-1)*ipz1
    else if ( coordz > ipz1 .and. coordz < ipz2 ) then
        lsize_opt(3) = nz - 1
        start_opt(3) = (nz - 1) * coordz
    else if (coordz == ipz2) then
        lsize_opt(3) = jz_opt2 - (flagz*jz_opt1 + (1-flagz)*1) + 1
        start_opt(3) = flagz*(jz_opt1-1) + (nz-1)*ipz2
    else 
        lsize_opt(3) = nzt-1
        start_opt(3) = 0
    end if

    if (coordy == ipy1) then
        lsize_opt(2) = flagy*jy_opt2 + (1-flagy)*nynpy - jy_opt1 + 1
        start_opt(2) = jy_opt1 - 1 + nynpy*ipy1
    else if (coordy > ipy1 .and. coordy < ipy2) then
        lsize_opt(2) = nynpy
        start_opt(2) = nynpy*coordy
    else if (coordy == ipy2) then
        lsize_opt(2) = jy_opt2 - (flagy*jy_opt1 + (1-flagy)*1) + 1
        start_opt(2) = flagy*(jy_opt1-1) + nynpy*ipy2
    else
        lsize_opt(2) = nyt
        start_opt(2) = 0 
    end if
            
    gsize_opt(1) = nxout
    gsize_opt(2) = nyt
    gsize_opt(3) = nzt-1
    lsize_opt(1) = nxout
    start_opt(1) = 0

    ! --- 
    if ( lsize_opt(2) == nyt .or. lsize_opt(3) == nzt-1 ) then
        bufsize = 0
        if (lsize_opt(2) == nyt) then
            lsize_opt(2) = nynpy
            start_opt(2) = nynpy*coordy
        end if

        if (lsize_opt(3) == nzt-1) then
            lsize_opt(3) = nz - 1
            start_opt(3) = (nz - 1) * coordz
        end if
    else 
        allocate(temp(1:lsize_opt(1), 1:lsize_opt(2), 1:lsize_opt(3)))
        bufsize = size(temp(1:lsize_opt(1), 1:lsize_opt(2), 1:lsize_opt(3)))
        deallocate(temp)
    end if
end if
! ---
end subroutine init_meshgrid
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine allocate_output_variable()
!-----------------------------------------------------------------------
!   output variables initialization (Bicheng Chen 06/13/2016)
!-----------------------------------------------------------------------
use io
use param !, only:nxt, nyt, nz, nzt, use_avgslice, average_dim_num, USE_MPI
use intermediate !, only:average_dim_select_flag, dim1_size, dim2_size, &
                 !     dim1_global, dim2_global, avg_out
implicit none
! ---
if (use_avgslice) then
    allocate (avg_u(nxt, nz-1), avg_v(nxt, nz-1), avg_w(nxt, nz-1), avg_p(nxt, nz-1),   &
              avg_dudt(nxt, nz-1), avg_dvdt(nxt, nz-1), avg_dwdt(nxt, nz-1),            &
              avg_dudz(nxt, nz-1), avg_dvdz(nxt, nz-1))
    allocate (avg_u2(nxt, nz-1), avg_v2(nxt, nz-1), avg_w2(nxt, nz-1), avg_p2(nxt, nz-1), &
              avg_uw(nxt, nz-1), avg_vw(nxt, nz - 1), avg_uv(nxt, nz - 1), avg_wp(nxt, nz-1),   &
              avg_txx(nxt, nz-1), avg_txy(nxt, nz-1), avg_txz(nxt, nz-1),               &
              avg_tyy(nxt, nz-1), avg_tyz(nxt, nz-1), avg_tzz(nxt, nz-1))
    allocate (avg_u3(nxt, nz-1), avg_v3(nxt, nz-1), avg_w3(nxt, nz-1))

    allocate (avg_u2_fluc(nxt, nz-1), avg_v2_fluc(nxt, nz-1), avg_uv_fluc(nxt, nz-1),   &
              avg_u2w_fluc(nxt, nz-1), avg_v2w_fluc(nxt, nz-1), avg_dissip(nxt, nz-1),  &
              avg_utxz_fluc(nxt, nz-1), avg_vtyz_fluc(nxt, nz-1), avg_wtzz(nxt, nz-1), avg_Cs(nxt, nz-1))

    avg_u = 0._rprec; avg_v = 0._rprec; avg_w = 0._rprec; avg_p = 0._rprec
    avg_dudt = 0._rprec; avg_dvdt = 0._rprec; avg_dwdt = 0._rprec
    avg_dudz = 0._rprec; avg_dvdz = 0._rprec
    avg_u2 = 0._rprec; avg_v2 = 0._rprec; avg_w2 = 0._rprec; avg_p2 = 0._rprec
    avg_uw = 0._rprec; avg_vw = 0._rprec; avg_uv = 0._rprec; avg_wp = 0._rprec
    avg_txx = 0._rprec; avg_txy = 0._rprec; avg_txz = 0._rprec
    avg_tyy = 0._rprec; avg_tyz = 0._rprec; avg_tzz = 0._rprec
    avg_u3 = 0._rprec; avg_v3 = 0._rprec; avg_w3 = 0._rprec
    avg_u2_fluc = 0._rprec; avg_v2_fluc = 0._rprec; avg_uv_fluc = 0._rprec
    avg_u2w_fluc = 0._rprec; avg_v2w_fluc = 0._rprec; avg_dissip = 0._rprec
    avg_utxz_fluc = 0._rprec; avg_vtyz_fluc = 0._rprec; avg_wtzz = 0._rprec; avg_Cs = 0._rprec
                
    if (theta_flag) then 
        allocate (avg_theta(nxt, nz-1), avg_theta2(nxt, nz-1), avg_dTdz(nxt, nz-1),     &
                  avg_sgsTheta(nxt, nz-1), avg_wTheta(nxt, nz-1), avg_Nut(nxt, nz-1),   &
                  avg_T2_fluc(nxt, nz-1))
        
        avg_theta = 0._rprec; avg_theta2 = 0._rprec; avg_dTdz = 0._rprec
        avg_sgsTheta = 0._rprec; avg_wTheta = 0._rprec; avg_Nut = 0._rprec
        avg_T2_fluc = 0._rprec
    end if
    
    if (PCon_FLAG) then 
        allocate (avg_PCon(nxt, nz-1, npcon), avg_PCon2(nxt, nz-1, npcon), avg_dPCondz(nxt, nz-1, npcon),   &
                  avg_sgsPCon1(nxt, nz-1, npcon), avg_sgsPCon3(nxt, nz-1, npcon),   &
                  avg_uPCon(nxt, nz-1, npcon), avg_wPCon(nxt, nz-1, npcon),   &
                  avg_Kct(nxt, nz-1), avg_Cs2Sc(nxt, nz-1))

        avg_PCon = 0._rprec; avg_PCon2 = 0._rprec; avg_dPCondz = 0._rprec; 
        avg_sgsPCon1 = 0._rprec; avg_sgsPCon3 = 0._rprec; 
        avg_uPCon = 0._rprec; avg_wPCon = 0._rprec;
        avg_Kct = 0._rprec; avg_Cs2Sc = 0._rprec
    end if
end if
! ---
end subroutine allocate_output_variable
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine allocate_flow_variable()
!-----------------------------------------------------------------------
!   flow variables initialization
!! >>Add dynamical wind stress and dynamical Stokes drift
!-----------------------------------------------------------------------
use param
use sim_param
use fft
use test_filtermodule, only:G_test, G_test_test
use sgsmodule
use bottombc, only:ustar_avg, zo, z_os, patch, d0,          &
                   T_s, q_s, phi_m, psi_m, phi_h, psi_h,    &
                   zo_PCon, PCon_s
use topbc, only:sponge
use intermediate
use scalars_module
implicit none
! ---
integer :: jx, jy, ind

integer, dimension(npz) :: sendcounts, recvcounts, displs

real(rprec), dimension(:, :, :, :), allocatable :: k_ttlCR_tot
character(20), dimension(:), allocatable :: unt_bg      ! initial unit of background concentration
character(20) :: unt_to  
real(rprec), dimension(:, :), allocatable :: con_bg
real(rprec), dimension(:, :, :, :), allocatable :: con_bg_3d
real(rprec), dimension(:, :), allocatable :: k_bg       ! chemical reaction rate of specified species
integer, parameter :: fid_bg = 822
namelist/conprof/  unt_bg, unt_to, con_bg
namelist/chemrate/ k_bg
! ---
allocate (u(ldx, nynpy, 0:nz), v(ldx, nynpy, 0:nz), w(ldx, nynpy, 0:nz))
allocate (ke(ldx, nynpy, 0:nz), ke_temp(ldx, nynpy, 0:nz), &
          dragx(ldx, nynpy, 0:nz), dragy(ldx, nynpy, 0:nz), dragz(ldx, nynpy, 0:nz))
allocate (dudx(ldx, nynpy, 0:nz), dudy(ldx, nynpy, 0:nz), dudz(ldx, nynpy, 0:nz),        &
          dvdx(ldx, nynpy, 0:nz), dvdy(ldx, nynpy, 0:nz), dvdz(ldx, nynpy, 0:nz),        &
          dwdx(ldx, nynpy, 0:nz), dwdy(ldx, nynpy, 0:nz), dwdz(ldx, nynpy, 0:nz),        &
          RHSx(ldx, nynpy, 0:nz), RHSy(ldx, nynpy, 0:nz), RHSz(ldx, nynpy, 0:nz),        &
          RHSx_f(ldx, nynpy, 0:nz), RHSy_f(ldx, nynpy, 0:nz), RHSz_f(ldx, nynpy, 0:nz))
allocate (dudt(ldx, nynpy, 0:nz), dvdt(ldx, nynpy, 0:nz), dwdt(ldx, nynpy, 0:nz))
allocate (dkedx(ldx, nynpy, 0:nz), dkedy(ldx, nynpy, 0:nz), dkedz(ldx, nynpy, 0:nz))
allocate (txx(ldx, nynpy, 0:nz), txy(ldx, nynpy, 0:nz), txz(ldx, nynpy, 0:nz),           &
          tyy(ldx, nynpy, 0:nz), tyz(ldx, nynpy, 0:nz), tzz(ldx, nynpy, 0:nz))
allocate (divtx(ldx, nynpy, 0:nz), divty(ldx, nynpy, 0:nz), divtz(ldx, nynpy, 0:nz))
allocate (p(ldx, nynpy, 0:nz))
allocate (dpdx(ldx, nynpy, nz), dpdy(ldx, nynpy, nz), dpdz(ldx, nynpy, nz))

allocate (kx(lhx), ky(lhx))
allocate (kx_2d(lhx, nyt), ky_2d(lhx, nyt), k2(lhx, nyt), notNyquist(lhx, nyt))
allocate (kx_2d_mpi(nxhnpy, nyt), ky_2d_mpi(nxhnpy, nyt), k2_mpi(nxhnpy, nyt))
allocate (dealias(nxhnpy, nyt), notNyquist_mpi(nxhnpy, nyt))
allocate (G_test(nxhnpy, nyt), G_test_test(nxhnpy, nyt))

allocate (Nu_t(ldx, nynpy, 0:nz), dissip(ldx, nynpy, 0:nz), magS(ldx, nynpy, 0:nz),     &
          u_lag(ldx, nynpy, 0:nz),v_lag(ldx, nynpy, 0:nz), w_lag(ldx, nynpy, 0:nz))
allocate (F_LM(ldx, nynpy, nz), F_MM(ldx, nynpy, nz), F_QN(ldx, nynpy, nz),             &
          F_NN(ldx, nynpy, nz), Beta(ldx, nynpy, nz), Betaclip(ldx, nynpy, nz))
allocate (Beta_avg(nz), Betaclip_avg(nz))
allocate (Cs_opt2(ldx, nynpy, nz), Cs_opt2_avg(ldx, nynpy, nz), Cs_Ssim(ldx, nynpy, nz), &
          epsilon_lag(ldx, nynpy, nz), epsilon_lag2(ldx, nynpy, nz),                  &
          xlag(ldx, nynpy, nz), ylag(ldx, nynpy, nz), zlag(ldx, nynpy, nz))
allocate (F_KX(nxt, nynpy, nz), F_XX(nxt, nynpy, nz), F_KX2(nxt, nynpy, nz), F_XX2(nxt, nynpy, nz))

allocate (ustar_avg(nxt, nynpy), zo(nxt, nynpy), z_os(nxt, nynpy), patch(nxt, nynpy), d0(nxt, nynpy))
allocate (sponge(0:nz))

allocate (cross_x(ldx, nynpy, nz), cross_y(ldx, nynpy, nz), cross_z(ldx, nynpy, nz))
allocate (vort1(ldx, nynpy, 0:nz), vort2(ldx, nynpy, 0:nz), vort3(ldx, nynpy, 0:nz))
allocate (beta_scal(ldx, nynpy, 0:nz))

allocate (cross_x_big(nxt2, ny2npy, nz), cross_y_big(nxt2, ny2npy, nz), cross_z_big(nxt2, ny2npy, nz))
allocate (u_big(nxt2, ny2npy, 0:nz), v_big(nxt2, ny2npy, 0:nz), w_big(nxt2, ny2npy, 0:nz))
allocate (vort1_big(nxt2, ny2npy, 0:nz), vort2_big(nxt2, ny2npy, 0:nz), vort3_big(nxt2, ny2npy, 0:nz))

! --- press_stag_array
allocate (rH_x(ldx, nynpy, 0:nz), rH_y(ldx, nynpy, 0:nz), rH_z(ldx, nynpy, 0:nz))
allocate (H_x(nxhnpy, nyt, 0:nz), H_y(nxhnpy, nyt, 0:nz), H_z(nxhnpy, nyt, 0:nz))
allocate (rtopw(ldx, nynpy), rbottomw(ldx, nynpy))
allocate (topw(nxhnpy, nyt), bottomw(nxhnpy, nyt))

! --- SGS model
if ( sgs ) then
    if (model == 3) then
        allocate (L11(ldx, nynpy), L12(ldx, nynpy), L13(ldx, nynpy),                        &
                  L22(ldx, nynpy), L23(ldx, nynpy), L33(ldx, nynpy),                        &
                  Q11(ldx, nynpy), Q12(ldx, nynpy), Q13(ldx, nynpy),                        &
                  Q22(ldx, nynpy), Q23(ldx, nynpy), Q33(ldx, nynpy), S_bar(ldx, nynpy),     &
                  S11_bar(ldx, nynpy), S12_bar(ldx, nynpy), S13_bar(ldx, nynpy),            &
                  S22_bar(ldx, nynpy), S23_bar(ldx, nynpy), S33_bar(ldx, nynpy),            &
                  S_S11_bar(ldx, nynpy), S_S12_bar(ldx, nynpy), S_S13_bar(ldx, nynpy),      &
                  S_S22_bar(ldx, nynpy), S_S23_bar(ldx, nynpy), S_S33_bar(ldx, nynpy),      &
                  S_hat(ldx, nynpy),       &
                  S11_hat(ldx, nynpy), S12_hat(ldx, nynpy), S13_hat(ldx, nynpy),            &
                  S22_hat(ldx, nynpy), S23_hat(ldx, nynpy), S33_hat(ldx, nynpy),            &
                  S_S11_hat(ldx, nynpy), S_S12_hat(ldx, nynpy), S_S13_hat(ldx, nynpy),      &
                  S_S22_hat(ldx, nynpy), S_S23_hat(ldx, nynpy), S_S33_hat(ldx, nynpy),      &
                  u_bar(ldx, nynpy), v_bar(ldx, nynpy), w_bar(ldx, nynpy),                  &
                  u_hat(ldx, nynpy), v_hat(ldx, nynpy), w_hat(ldx, nynpy), S(ldx, nynpy))
        allocate (beta_sd(nz))
        M11 => Q11; M12 => Q12; M13 => Q13; M22 => Q22; M23 => Q23; M33 => Q33
    else if (model == 4) then
        allocate (u_pr(nxt, nynpy), v_pr(nxt, nynpy), w_pr(nxt, nynpy))
        allocate (w_nod(ldx, nynpy), S(ldx, nynpy), tempos(ldx, nynpy),                 &    
                  L11(ldx, nynpy), L12(ldx, nynpy), L13(ldx, nynpy),                    &
                  L22(ldx, nynpy), L23(ldx, nynpy), L33(ldx, nynpy),                    &
                  Q11(ldx, nynpy), Q12(ldx, nynpy), Q13(ldx, nynpy),                    &
                  Q22(ldx, nynpy), Q23(ldx, nynpy), Q33(ldx, nynpy),                    &
                  M11(ldx, nynpy), M12(ldx, nynpy), M13(ldx, nynpy),                    &
                  M22(ldx, nynpy), M23(ldx, nynpy), M33(ldx, nynpy),                    &
                  N11(ldx, nynpy), N12(ldx, nynpy), N13(ldx, nynpy),                    &
                  N22(ldx, nynpy), N23(ldx, nynpy), N33(ldx, nynpy),                    & 
                  LM(ldx, nynpy), MM(ldx, nynpy), QN(ldx, nynpy), NN(ldx, nynpy),       &
                  Tn(ldx, nynpy), epsi(ldx, nynpy), dumfac(ldx, nynpy), S_bar(ldx, nynpy),  &
                  S11_bar(ldx, nynpy), S12_bar(ldx, nynpy), S13_bar(ldx, nynpy),        &
                  S22_bar(ldx, nynpy), S23_bar(ldx, nynpy), S33_bar(ldx, nynpy),        &
                  S_S11_bar(ldx, nynpy), S_S12_bar(ldx, nynpy), S_S13_bar(ldx, nynpy),  &
                  S_S22_bar(ldx, nynpy), S_S23_bar(ldx, nynpy), S_S33_bar(ldx, nynpy),  &
                  S_hat(ldx, nynpy),    & 
                  S11_hat(ldx, nynpy), S12_hat(ldx, nynpy), S13_hat(ldx, nynpy),        &
                  S22_hat(ldx, nynpy), S23_hat(ldx, nynpy), S33_hat(ldx, nynpy),        &
                  S_S11_hat(ldx, nynpy), S_S12_hat(ldx, nynpy), S_S13_hat(ldx, nynpy),  &
                  S_S22_hat(ldx, nynpy), S_S23_hat(ldx, nynpy), S_S33_hat(ldx, nynpy),  &
                  u_bar(ldx, nynpy), v_bar(ldx, nynpy), w_bar(ldx, nynpy),              &
                  u_hat(ldx, nynpy), v_hat(ldx, nynpy), w_hat(ldx, nynpy), fourbeta(ldx, nynpy))
        allocate (xp(ldx, nynpy, 0:nz), yp(ldx, nynpy, 0:nz), zp(ldx, nynpy, 0:nz),     &
                  u_temp(ldx, nynpy, 0:nz), v_temp(ldx, nynpy, 0:nz))
        allocate (FF_LM(nxt + 2, nynpy + 2, nz + 2), FF_MM(nxt + 2, nynpy + 2, nz + 2))
  
    else if (model == 5) then
        allocate (visc(ldx, nynpy, nz), Cs_opt2_2d(ldx, nynpy, nz), Cs_opt2_4d(ldx, nynpy, nz))
        allocate (S(ldx, nynpy), tempos(ldx, nynpy),                                        & 
                  L11(ldx, nynpy), L12(ldx, nynpy), L13(ldx, nynpy),                        &
                  L22(ldx, nynpy), L23(ldx, nynpy), L33(ldx, nynpy),                        &
                  Q11(ldx, nynpy), Q12(ldx, nynpy), Q13(ldx, nynpy),                        &
                  Q22(ldx, nynpy), Q23(ldx, nynpy), Q33(ldx, nynpy),                        &
                  M11(ldx, nynpy), M12(ldx, nynpy), M13(ldx, nynpy),                        &
                  M22(ldx, nynpy), M23(ldx, nynpy), M33(ldx, nynpy),                        &
                  N11(ldx, nynpy), N12(ldx, nynpy), N13(ldx, nynpy),                        &
                  N22(ldx, nynpy), N23(ldx, nynpy), N33(ldx, nynpy),                        & 
                  LM(ldx, nynpy), MM(ldx, nynpy), QN(ldx, nynpy), NN(ldx, nynpy),           &
                  Tn(ldx, nynpy), epsi(ldx, nynpy), dumfac(ldx, nynpy), S_bar(ldx, nynpy),  &
                  S11_bar(ldx, nynpy), S12_bar(ldx, nynpy), S13_bar(ldx, nynpy),            &
                  S22_bar(ldx, nynpy), S23_bar(ldx, nynpy), S33_bar(ldx, nynpy),            &
                  S_S11_bar(ldx, nynpy), S_S12_bar(ldx, nynpy), S_S13_bar(ldx, nynpy),      & 
                  S_S22_bar(ldx, nynpy), S_S23_bar(ldx, nynpy), S_S33_bar(ldx, nynpy),      &
                  S_hat(ldx, nynpy),        & 
                  S11_hat(ldx, nynpy), S12_hat(ldx, nynpy), S13_hat(ldx, nynpy),            &
                  S22_hat(ldx, nynpy), S23_hat(ldx, nynpy), S33_hat(ldx, nynpy),            &
                  S_S11_hat(ldx, nynpy), S_S12_hat(ldx, nynpy), S_S13_hat(ldx, nynpy),      &
                  S_S22_hat(ldx, nynpy), S_S23_hat(ldx, nynpy), S_S33_hat(ldx, nynpy),      &  
                  u_bar(ldx, nynpy), v_bar(ldx, nynpy), w_bar(ldx, nynpy),                  &
                  u_hat(ldx, nynpy), v_hat(ldx, nynpy), w_hat(ldx, nynpy))
        allocate (LMvert(nz), MMvert(nz), QNvert(nz), NNvert(nz))
        allocate (xp(ldx, nynpy, 0:nz), yp(ldx, nynpy, 0:nz), zp(ldx, nynpy, 0:nz),         &
                  u_temp(ldx, nynpy, 0:nz), v_temp(ldx, nynpy, 0:nz))
        allocate (FF_LM(nxt + 2, nynpy + 2, nz + 2), FF_MM(nxt + 2, nynpy + 2, nz + 2),     &
                  FF_QN(nxt + 2, nynpy + 2, nz + 2), FF_NN(nxt + 2, nynpy + 2, nz + 2),     &
                  Beta_t(nxt + 2, nynpy + 2, nz + 2)) 
    end if
end if

! --- Temperature field
allocate (phi_m(nxt, nynpy), psi_m(nxt, nynpy), phi_h(nxt, nynpy), psi_h(nxt, nynpy))
allocate (L(nxt, nynpy), wstar(nxt, nynpy))

if (theta_flag) then
    allocate (T_s(nxt, nynpy), q_s(nxt, nynpy))
    allocate (theta(ldx, nynpy, 0:nz), q(ldx, nynpy, 0:nz))
    allocate (Pr_(ldx, nynpy, 0:nz), dTdz(ldx, nynpy, 0:nz), dqdz(ldx, nynpy, 0:nz),         &
              RHS_Tf(ldx, nynpy, 0:nz), RHS_T(ldx, nynpy, 0:nz), RHS_qf(ldx, nynpy, 0:nz),   &
              RHS_q(ldx, nynpy, 0:nz), sgs_t3(ldx, nynpy, 0:nz), sgs_q3(ldx, nynpy, 0:nz))    
    allocate (T_s_filtered(ldx, nynpy))
    
    ! - SGS model
    if (model_sc == 4) then
        allocate (FF_KX(nxt + 2, nynpy + 2, nz + 2), FF_XX(nxt + 2, nynpy + 2, nz + 2))
    else if (model_sc == 5) then
        allocate (FF_KX(nxt + 2, nynpy + 2, nz + 2), FF_XX(nxt + 2, nynpy + 2, nz + 2), &
                  FF_KX2(nxt + 2, nynpy + 2, nz + 2), FF_XX2(nxt + 2, nynpy + 2, nz + 2))
    end if
    !- use dynamical kinematic heat flux (Bicheng Chen)
    if (flag_dynWT) then
        inquire (iolength=len_wt) wt_s
        open (unit=fid_wt, file=fn_wt, access='direct', recl=len_wt)
        read (unit=fid_wt, rec=1) wt_s
        close (fid_wt)
    end if
end if

! --- concentration variables initialization
allocate (beta_pcon(ldx, 0:nynpy+1, 0:nz))  ! allocate the variable always needed
if (pcon_flag) then
    npcon = info_con%n_con
    allocate (zo_PCon(nxt, nynpy), PCon_s(nxt, nynpy))
    allocate (P_surf_flux(nxt, nynpy), deposition(nxt, nynpy),  &
              P_surf_flux_dep(nxt, nynpy), Real_dep(nxt, nynpy))
    allocate (matrix_x(nxt, nxt), matrix_y(nyt, nyt))
    allocate (dvector_x(nxt), dvector_y(nyt))

    allocate (settling_vel(npcon), densratio_pcon(npcon), v_pcon(npcon), pcon_sfc_flux_nd(npcon))
    if (settling) then
        settling_vel = info_con%vel_settling(1:npcon)/u_scale
        densratio_pcon = info_con%ratio_dens(1:npcon)
        v_pcon = info_con%vol_spec(1:npcon)*pcon_scale
    else 
        settling_vel = 0._rprec; densratio_pcon = 0._rprec
        v_pcon = 0._rprec
    end if
    pcon_sfc_flux_nd = info_con%pcon_sfc_flux(1:npcon)/u_scale/pcon_scale

    allocate (Kc_t(ldx, 0:nynpy+1, 0:nz), dPCondz(ldx, 0:nynpy+1, 0:nz, npcon),          &
              sgs_PCon3(ldx, 0:nynpy+1, 0:nz, npcon), res_PCon3(ldx, 0:nynpy+1, 0:nz, npcon),   &
              sgs_PCon1(ldx, 0:nynpy+1, 0:nz, npcon), Cs2Sc(ldx, 0:nynpy+1, 0:nz))
    
    !- allocate concentration
    allocate (PCon(ldx, 0:nynpy+1, 0:nz, npcon))
    allocate (RHS_PConf(ldx, 0:nynpy+1, 0:nz, npcon), RHS_PCon(ldx, 0:nynpy+1, 0:nz, npcon))

    call read_source_info()
end if

! ---
Cs_opt2_avg = 0._rprec
! ---
end subroutine allocate_flow_variable
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine init_vel_field()
!-----------------------------------------------------------------------
!   initialize velocity field
!-----------------------------------------------------------------------
use types, only:rprec
use param
use sim_param
implicit none
! ---
if ( restart_vel ) then
    if ( inflow .and. read_inflow_file .and. interp ) then
        call read_precursor_field()
    else 
        call read_flow_field()
    end if 
else
    if ( read_ic_profile ) then
        call ic_read()
    else
        if (theta_flag) then
            if ( intp_vel ) then
                call create_fft_plan_intp()
                call read_coarse_field()
                call intp_flow_field()
            else
                call ic_scal()
            end if
        else
            call ic()
        end if
    end if
end if

!- add cross velocity to the all velocity field (Bicheng Chen 09/07/2016)
if (flag_crossVel .and. ocean_flag) then
    u(1:nxt, 1:nynpy, 0:nz) = u(1:nxt, 1:nynpy, 0:nz) + u_cross
    v(1:nxt, 1:nynpy, 0:nz) = v(1:nxt, 1:nynpy, 0:nz) + v_cross
end if

! ---
end subroutine init_vel_field
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine init_temperature_field()
!-----------------------------------------------------------------------
!   initialize temperature field
!-----------------------------------------------------------------------
use types, only:rprec
use param
use sim_param
use scalars_module, only:RHS_T, sgs_t3
use bottombc, only:psi_m
implicit none
! ---
integer(MPI_OFFSET_KIND) :: disp0, disp2d, disp3d
real(rprec), dimension(nxt/2, nynpy, nz-1) :: theta_tmp, RHS_T_tmp
real(rprec), dimension(nxt/2, nynpy) :: sgs_t3_tmp, psi_m_tmp
! ---
character(80) :: fname
logical :: exst
! ---
if (theta_flag .and. restart_theta .and. restart_vel) then

    if ( inflow .and. read_inflow_file .and. interp ) then    
        write (fname, '(a,i8.8,a)') path_restart//'temp_tt', nums, '.out'
        
        disp0 = 0
        disp2d = npy * sizeof(sgs_t3_tmp(1:nxt/2, 1:nynpy))
        disp3d = np  * sizeof(theta_tmp(1:nxt/2, 1:nynpy, 1:nz-1))

        call read_file_mpi_3d_inflow(theta_tmp(1:nxt/2, 1:nynpy, 1:nz-1), fname, disp0)
        disp0 = disp0 + disp3d
        call read_file_mpi_3d_inflow(RHS_T_tmp(1:nxt/2, 1:nynpy, 1:nz-1), fname, disp0)
        disp0 = disp0 + disp3d
        call read_file_mpi_xy_inflow(sgs_t3_tmp(1:nxt/2, 1:nynpy), fname, disp0)
        disp0 = disp0 + disp2d
        call read_file_mpi_xy_inflow(psi_m_tmp, fname, disp0)

        theta(1:nxt/2, 1:nynpy, 1:nz-1) = theta_tmp(1:nxt/2, 1:nynpy, 1:nz-1)
        RHS_T(1:nxt/2, 1:nynpy, 1:nz-1) = RHS_T_tmp(1:nxt/2, 1:nynpy, 1:nz-1)
        sgs_t3(1:nxt/2, 1:nynpy, 1) = sgs_t3_tmp(1:nxt/2, 1:nynpy)
        psi_m(1:nxt/2, 1:nynpy) = psi_m_tmp(1:nxt/2, 1:nynpy)

        theta(nxt/2+1:nxt, 1:nynpy, 1:nz-1) = theta_tmp(1:nxt/2, 1:nynpy, 1:nz-1)
        RHS_T(nxt/2+1:nxt, 1:nynpy, 1:nz-1) = RHS_T_tmp(1:nxt/2, 1:nynpy, 1:nz-1)
        sgs_t3(nxt/2+1:nxt, 1:nynpy, 1) = sgs_t3_tmp(1:nxt/2, 1:nynpy)
        psi_m(nxt/2+1:nxt, 1:nynpy) = psi_m_tmp(1:nxt/2, 1:nynpy)
    else
        write (fname, '(a,i8.8,a)') path_restart//'temp_tt', nums, '.out'
        disp0 = 0
        disp2d = npy * sizeof(sgs_t3(1:nxt, 1:nynpy, 1))
        disp3d = np  * sizeof(theta(1:nxt, 1:nynpy, 1:nz-1))

        call read_file_mpi_3d(theta(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
        disp0 = disp0 + disp3d
        call read_file_mpi_3d(RHS_T(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
        disp0 = disp0 + disp3d
        call read_file_mpi_xy(sgs_t3(1:nxt, 1:nynpy, 1), fname, disp0)
        disp0 = disp0 + disp2d
        call read_file_mpi_xy(psi_m, fname, disp0)
    end if

else if (theta_flag .and. restart_theta .and. (.not. restart_vel)) then !LV1
    print *, "Cannot initialize temperature field with initializing velocity."
    print *, "Stop"
    stop
else if ((.not. theta_flag) .and. restart_theta .and. restart_vel) then !LV1
    print *, "Cannot initialize temperature field with theta_flag."
    print *, "Stop"
    stop
else
    !- the initialization of temperature without file is in ic_scal
end if
! ---
end subroutine init_temperature_field
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine init_concentration_field()
!-----------------------------------------------------------------------
!   concentration variables initialization
!-----------------------------------------------------------------------
use types, only:rprec
use param
use sim_param
use scalars_module, only:deposition, Real_dep, RHS_PCon, Kc_t
use sgsmodule, only:F_KX, F_XX, F_KX2, F_XX2
implicit none
! ---
integer(MPI_OFFSET_KIND) :: disp0, disp2d, disp3d
! ---
character(80) :: fname
logical :: exst
integer :: ipcon
! ---
if (PCon_flag .and. restart_pcon) then
    write (fname, '(a,i8.8,a)') path_restart//'con_tt', nums, '.out'
    disp0 = 0
    disp2d = npy * sizeof(deposition)
    disp3d = np  * sizeof(PCon(1:nxt, 1:nynpy, 1:nz-1, 1))
    do ipcon = 1, npcon
        call read_file_mpi_3d(PCon(1:nxt, 1:nynpy, 1:nz-1, ipcon), fname, disp0)
        call read_file_mpi_3d(RHS_Pcon(1:nxt, 1:nynpy, 1:nz-1, ipcon), fname, disp0+npcon*disp3d)
        disp0 = disp0 + disp3d
    end do

    write (fname, '(a,i8.8,a)') path_restart//'con_diff_tt', nums, '.out'
    disp0 = 0
    call read_file_mpi_3d(Kc_t(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
    disp0 = disp0+disp3d
    call read_file_mpi_3d(F_KX(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
    disp0 = disp0+disp3d
    call read_file_mpi_3d(F_XX(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
    disp0 = disp0+disp3d
    call read_file_mpi_3d(F_KX2(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
    disp0 = disp0+disp3d
    call read_file_mpi_3d(F_XX2(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
    disp0 = disp0+disp3d
    call read_file_mpi_xy(deposition, fname, disp0)
    disp0 = disp0+disp2d
    call read_file_mpi_xy(Real_dep, fname, disp0)
else if ((.not. pcon_flag) .and. restart_pcon) then
    write(*,*) "Cannot initialize concentration field without pcon_flag. Stop."
    stop
else if (PCon_flag .and. .not. restart_pcon) then
    !- Always initialize pollen concentration with zeros
    PCon(:, :, :, :) = 0._rprec
    if ( coordz == 0 ) then
        PCon(:, :, 0, :) = BOGUS
    end if
end if 
    
if (PCon_flag .and. .not. restart_pcon .and. inflow .and. write_inflow_file) then
    call nutrient_init_read()
else if (PCon_flag .and. .not. restart_pcon .and. .not. inflow .and. read_inflow_file) then
    call nutrient_inflow_read()
end if
! ---
end subroutine init_concentration_field
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine read_source_info
!-----------------------------------------------------------------------
!   initialize point sources and single sources
!-----------------------------------------------------------------------
use param
use scalars_module
implicit none
! ---
character(80) :: fname
integer :: ipcon
logical :: exst
! ---
select case (source_type)
case ("point")
    !-- allocate variable
    allocate (con_src(npcon))
    !-- read source release data
    do ipcon = 1, info_con%n_con
        write(fname, "(a,i0,a)") trim(fnPre_src), ipcon, trim(fnSuf_src)
        inquire(file=fname, exist=exst)
        if (exst) then
            !--- readin temp variable
            open(fid_src, file=fname, form='formatted', status='old')
            read(fid_src, nml=number_source)
            read(fid_src, nml=source)
            close(fid_src)
            !--- copy to real source variable
            con_src(ipcon)%n = n_src
            con_src(ipcon)%ix(1:n_src) = ix_src(1:n_src)
            con_src(ipcon)%iy(1:n_src) = iy_src(1:n_src)
            con_src(ipcon)%iz(1:n_src) = iz_src(1:n_src)
            con_src(ipcon)%rls(1:n_src) = rls_src(1:n_src)
        end if
    end do
  
    !-- calculate the specific cpu coord for point sources and normalize the source release
    do ipcon = 1, npcon
        con_src(ipcon)%icpu_y = int(con_src(ipcon)%iy/nynpy)
        con_src(ipcon)%icpu_z = int((con_src(ipcon)%iz - 1)/(nz - 1))
        con_src(ipcon)%iy_cpu = con_src(ipcon)%iy - con_src(ipcon)%icpu_y*nynpy
        con_src(ipcon)%iz_cpu = con_src(ipcon)%iz - con_src(ipcon)%icpu_z*(nz - 1)
        con_src(ipcon)%rls = con_src(ipcon)%rls/(dx*z_i*dy*z_i*dz*z_i)/(PCon_scale/(z_i/u_scale))
    end do
case ("line")
    !-- allocate variable
    allocate (con_src(npcon))
    !-- read source release data
    do ipcon = 1, info_con%n_con
        write(fname, "(a,i0,a)") trim(fnPre_src), ipcon, trim(fnSuf_src)
        inquire(file=fname, exist=exst)
        if (exst) then
            !--- read in temp variable
            open(fid_src, file=fname, form='formatted', status='old')
            read(fid_src, nml=number_source)
            read(fid_src, nml=source_alignment)
            read(fid_src, nml=source)
            close(fid_src)
            !--- copy to real source variable
            select case (src_align)
            case ('lateral')
                con_src(ipcon)%n = n_src
                con_src(ipcon)%ix(1:n_src) = ix_src(1:n_src)
                ! con_src(ipcon)%iy(1:n_src) = iy_src(1:n_src)
                con_src(ipcon)%iz(1:n_src) = iz_src(1:n_src)
                con_src(ipcon)%rls(1:n_src) = rls_src(1:n_src)
                ! con_src(ipcon)%xs(1:n_src) = x_src(1:n_src)
                ! ! con_src(ipcon)%ys(1:n_src) = y_src(1:n_src)
                ! con_src(ipcon)%zs(1:n_src) = z_src(1:n_src)
                ! con_src(ipcon)%rls(1:n_src) = rls_src(1:n_src)
            case ('streamwise')
                con_src(ipcon)%n = n_src
                ! con_src(ipcon)%ix(1:n_src) = ix_src(1:n_src)
                con_src(ipcon)%iy(1:n_src) = iy_src(1:n_src)
                con_src(ipcon)%iz(1:n_src) = iz_src(1:n_src)
                con_src(ipcon)%rls(1:n_src) = rls_src(1:n_src)
                ! ! con_src(ipcon)%xs(1:n_src) = x_src(1:n_src)
                ! con_src(ipcon)%ys(1:n_src) = y_src(1:n_src)
                ! con_src(ipcon)%zs(1:n_src) = z_src(1:n_src)
                ! con_src(ipcon)%rls(1:n_src) = rls_src(1:n_src)
            case default
                write (*, *) 'invalid source alignment type'
                stop
            end select
        end if
    end do

    !-- calculate the specific cpu coord for point sources and normalize the source release
    do ipcon = 1, npcon
        con_src(ipcon)%icpu_y = int(con_src(ipcon)%iy/nynpy)
        con_src(ipcon)%icpu_z = int((con_src(ipcon)%iz - 1)/(nz - 1))
        con_src(ipcon)%iy_cpu = con_src(ipcon)%iy - con_src(ipcon)%icpu_y*nynpy
        con_src(ipcon)%iz_cpu = con_src(ipcon)%iz - con_src(ipcon)%icpu_z*(nz - 1)
        con_src(ipcon)%rls = con_src(ipcon)%rls/(dx*z_i*dy*z_i*dz*z_i)/(PCon_scale/(z_i/u_scale))
    end do
case ("planar")
    !-- allocate variable
    allocate (con_src(npcon))
    !-- read source release data
    do ipcon = 1, info_con%n_con
        write (fname, "(a,i0,a)") trim(fnPre_src), ipcon, trim(fnSuf_src)
        inquire(file=fname, exist=exst)
        if (exst) then
            !--- readin temp variable
            open(fid_src, file=fname, form='formatted', status='old')
            read(fid_src, nml=number_source)
            read(fid_src, nml=source)
            close(fid_src)
            !--- copy to real source variable
            con_src(ipcon)%n = n_src
            con_src(ipcon)%iz(1:n_src) = iz_src(1:n_src)
            con_src(ipcon)%rls(1:n_src) = rls_src(1:n_src)
        end if
    end do

    !-- calculate the specific cpu coord for point sources and normalize the source release
    do ipcon = 1, npcon
        con_src(ipcon)%icpu_z = int((con_src(ipcon)%iz - 1)/(nz - 1))
        con_src(ipcon)%iz_cpu = con_src(ipcon)%iz - con_src(ipcon)%icpu_z*(nz - 1)
        con_src(ipcon)%rls = con_src(ipcon)%rls/(dx*z_i*dy*z_i*dz*z_i)/(PCon_scale/(z_i/u_scale))
    end do
case default
    write(*,*) "Source type NOT supprted"
    stop
end select
! ---
end subroutine read_source_info
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------