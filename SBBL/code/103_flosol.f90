!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
subroutine flosol
!-----------------------------------------------------------------------
!   solve all the governing equations
!----------------------------------------------------------------------- 
use param
use sim_param
use scalars_module, only:beta_scal, obukhov, theta_all_in_one,    &
                         pollen_all_in_one2, beta_pcon
use scalars_module2
use derivatives
use intermediate
use topbc, only:sponge
use io
use stokes_drift, only:ust, vst
use canopy, only:nutrient_inflow 
implicit none
! ---
integer :: jz, ipcon
! ---
call newstep

! --- Call obukhov to calculate the MO functions
call obukhov()

! --- Calculate the spatial derivatives

call ddx(dudx, u)
call ddy(dudy, u)
call ddx(dvdx, v)
call ddy(dvdy, v)
call ddx(dwdx, w)
call ddy(dwdy, w)
call ddz_uv(dudz, u)    ! on exit of ddz_uv, have dudz, dvdz at 1:nz, except bottom process has 2:nz, top processor has 1:nz-1
call ddz_uv(dvdz, v)
call ddz_w(dwdz, w)     ! on exit of ddz_w, have dwdz at 0:nz-1, and bottom process has 1:nz-1

! --- calculate wall stress and calculate derivatives at bottom and top boundary
if (dns_bc) then
    if ( coordz == 0 ) then
        call wallstress_dns_bottom()
    else if ( coordz == npz-1 ) then
        call wallstress_dns_top()
    end if
else
    if ( coordz == 0 ) then
        call wallstress_bottom() !--provides txz, tyz, dudz, dvdz at jz=1 and nz
    else if ( coordz == npz-1 ) then
        call wallstress_top()
    end if
end if

! --- compute turbulent viscosity (const.)
if (dns_bc .and. molec) then
    call dns_stress(txx, txy, txz, tyy, tyz, tzz)
else
    call sgs_stag()         !--MPI: txx, txy, tyy, tzz at 1:nz-1; txz, tyz at 1:nz
end if

!--exchange ghost-node information for tij
!--send stuff up to ghost nodes
!--move this into sgs_stag?
call mpi_sendrecv(txz(:, :, nz - 1), ldx*nynpy, MPI_RPREC, up,   4, &
                  txz(:, :, 0),      ldx*nynpy, MPI_RPREC, down, 4, &
                  comm, status, ierr)
call mpi_sendrecv(tyz(:, :, nz - 1), ldx*nynpy, MPI_RPREC, up,   5, &
                  tyz(:, :, 0),      ldx*nynpy, MPI_RPREC, down, 5, &
                  comm, status, ierr)
call mpi_sendrecv(tzz(:, :, nz - 1), ldx*nynpy, MPI_RPREC, up,   6, &
                  tzz(:, :, 0),      ldx*nynpy, MPI_RPREC, down, 6, &
                  comm, status, ierr)
call mpi_sendrecv(tzz(:, :, 1),      ldx*nynpy, MPI_RPREC, down, 7, &
                  tzz(:, :, nz),     ldx*nynpy, MPI_RPREC, up,   7, &
                  comm, status, ierr)
! compute divergence of SGS shear stresses
! note: the divt's and the diagonal elements of t are equivalenced!
!--actually, they are not equivalenced in this version

! --- provides divtz 1:nz-1
call divstress_uv(divtx, txx, txy, txz)
call divstress_uv(divty, txy, tyy, tyz)
call divstress_w(divtz, txz, tyz, tzz)  ! --- provides divtz 1:nz-1, except 1:nz at top process 

! --- provides RHS{x,y,z} 1:nz-1
call convec(RHSx, RHSy, RHSz)

call mpi_sendrecv(divtx(:, :, nz - 1), ldx*nynpy, MPI_RPREC, up,   6, &
                  divtx(:, :, 0),      ldx*nynpy, MPI_RPREC, down, 6, &
                  comm, status, ierr)
call mpi_sendrecv(divtx(:, :, 1),      ldx*nynpy, MPI_RPREC, down, 7, &
                  divtx(:, :, nz),     ldx*nynpy, MPI_RPREC, up,   7, &
                  comm, status, ierr)
! ---
if (theta_flag .and. (jt .ge. theta_init_time)) then
    call theta_all_in_one
else
    beta_scal = 0._rprec
end if
! --- Time advance for pollen Chamecki - 08/01/2006
beta_pcon = 0._rprec    !DY Initialize beta_pcon to be zero

if (PCon_flag .and. (jt .ge. PCon_init)) then
    if (pcon_scheme == 1) then      ! Pseudo-spectral
        !!$ IF (PCon_scheme==1) CALL pollen_all_in_one
        print *, "Pseudo-spectral method for particle, stop!"
        stop
    end if
      ! QUICK and SMART
    if (pcon_scheme == 2 .or. pcon_scheme == 3) call pollen_all_in_one2
end if

! --- Compute preliminary RHS matrices for pressure calculation
RHSx(:, :, 1:nz - 1) = -RHSx(:, :, 1:nz - 1) - divtx(:, :, 1:nz - 1)
RHSy(:, :, 1:nz - 1) = -RHSy(:, :, 1:nz - 1) - divty(:, :, 1:nz - 1)
RHSz(:, :, 1:nz - 1) = -RHSz(:, :, 1:nz - 1) - divtz(:, :, 1:nz - 1)

dudt(:, :, 1:nz - 1) = -divtx(:, :, 1:nz - 1)
dvdt(:, :, 1:nz - 1) = -divty(:, :, 1:nz - 1)
dwdt(:, :, 1:nz - 1) = -divtz(:, :, 1:nz - 1)

! ---
if (theta_flag .and. (.not. passive_scalar)) then
    !--add buoyancy term...only valid for theta
    RHSz(:, :, 1:nz - 1) = RHSz(:, :, 1:nz - 1) + beta_scal(:, :, 1:nz - 1)
end if

! Di Yang Add buoyancy term due to droplets
if (active_pcon) then
    RHSz(:, :, 1:nz - 1) = RHSz(:, :, 1:nz - 1) + beta_pcon(:, :, 1:nz - 1)
end if
! ---
if (ocean_flag) then
    ! Enforce coriolis forcing for ocean simulation without geostrophic wind vector
    do jz = 1, nz - 1
        !-BC Changed by Bicheng Chen for direction Stokes drift
        RHSx(:, :, jz) = RHSx(:, :, jz) + coriol*(v(:, :, jz) + vst(jz)) - coriol*vg
        RHSy(:, :, jz) = RHSy(:, :, jz) - coriol*(u(:, :, jz) + ust(jz)) + coriol*ug
        dudt(:, :, jz) = dudt(:, :, jz) + coriol*(v(:, :, jz) + vst(jz)) - coriol*vg
        dvdt(:, :, jz) = dvdt(:, :, jz) - coriol*(u(:, :, jz) + ust(jz)) + coriol*ug
    end do
    !- Add the cross flow effect through pressure gradient force (09/07/2016)
    if (flag_crossVel) then
        RHSx(:, :, 1:nz - 1) = RHSx(:, :, 1:nz - 1) - coriol*v_cross
        RHSy(:, :, 1:nz - 1) = RHSy(:, :, 1:nz - 1) + coriol*u_cross
        dudt(:, :, 1:nz - 1) = dudt(:, :, 1:nz - 1) - coriol*v_cross
        dvdt(:, :, 1:nz - 1) = dvdt(:, :, 1:nz - 1) + coriol*u_cross
    end if
else
    if (coriolis_forcing) then
        RHSx(:, :, 1:nz - 1) = RHSx(:, :, 1:nz - 1) + coriol*(v(:, :, 1:nz - 1) - vg)
        RHSy(:, :, 1:nz - 1) = RHSy(:, :, 1:nz - 1) + coriol*(ug - u(:, :, 1:nz - 1))
        dudt(:, :, 1:nz - 1) = dudt(:, :, 1:nz - 1) + coriol*(v(:, :, 1:nz - 1) - vg)
        dvdt(:, :, 1:nz - 1) = dvdt(:, :, 1:nz - 1) + coriol*(ug - u(:, :, 1:nz - 1))
    end if
end if


call add_damping_layer()            ! add damping terms to the momentum RHS near the top
call canopy_model()                 ! calculate canopy drag
call Admas_Bashforth()              ! time advancement

! --- solve Poisson equation for pressure
!--do we ever need p itself, or only its gradient? -> probably
!  do not need to store p
!--provides p, dpdx, dpdy at 0:nz-1
call press_stag_array(p, dpdx, dpdy) 
call mpi_sendrecv(p(1, 1, nz - 1), ldx*nynpy, MPI_RPREC, up,   2, &
                  p(1, 1, 0),      ldx*nynpy, MPI_RPREC, down, 2, &
                  comm, status, ierr)
!--calculate dpdz here
!--careful, p is not dimensioned the same as the others
dpdz(1:nxt, 1:nynpy, 1:nz - 1) = (p(1:nxt, 1:nynpy, 1:nz - 1) - &
                                  p(1:nxt, 1:nynpy, 0:nz - 2))/dz
! ---
!--if really wanted to, could avoid storing pressure gradients
!  just add them directly to RHS in press_stag
RHSx(:, :, 1:nz - 1) = RHSx(:, :, 1:nz - 1) - dpdx(:, :, 1:nz - 1)
RHSy(:, :, 1:nz - 1) = RHSy(:, :, 1:nz - 1) - dpdy(:, :, 1:nz - 1)
RHSz(:, :, 1:nz - 1) = RHSz(:, :, 1:nz - 1) - dpdz(:, :, 1:nz - 1)

dudt(:, :, 1:nz - 1) = dudt(:, :, 1:nz - 1) - dpdx(:, :, 1:nz - 1)
dvdt(:, :, 1:nz - 1) = dvdt(:, :, 1:nz - 1) - dpdy(:, :, 1:nz - 1)
dwdt(:, :, 1:nz - 1) = dwdt(:, :, 1:nz - 1) - dpdz(:, :, 1:nz - 1)

! --- provides u, v, w at 1:nz
call project()
call data_exchange()
! ---
nums = nums + 1
call energy()
! ---
if ( coordz == 0 ) then
    ke(:,:,0) = 0._rprec
end if

call mpi_sendrecv(ke(1, 1, nz-1), ldx*nynpy, MPI_RPREC, up,   1,   &
                  ke(1, 1, 0)   , ldx*nynpy, MPI_RPREC, down, 1,   &
                  comm, status, ierr)
call mpi_sendrecv(ke(1, 1, 1) , ldx*nynpy, MPI_RPREC, down, 2,     &
                  ke(1, 1, nz), ldx*nynpy, MPI_RPREC, up,   2,     &
                  comm, status, ierr)
ke_temp = ke

call filter_data(ke_temp)
call ddx(dkedx, ke_temp)
call ddy(dkedy, ke_temp)
call ddz_uv(dkedz,ke)

dudt(:, :, 1:nz - 1) = dudt(:, :, 1:nz - 1) + dkedx(:, :, 1:nz - 1)
dvdt(:, :, 1:nz - 1) = dvdt(:, :, 1:nz - 1) + dkedy(:, :, 1:nz - 1)
dwdt(:, :, 1:nz - 1) = dwdt(:, :, 1:nz - 1) + dkedz(:, :, 1:nz - 1)

call mpi_sendrecv(dudt(1, 1, nz - 1), ldx*nynpy, MPI_RPREC, up,   1,   &
                  dudt(1, 1, 0),      ldx*nynpy, MPI_RPREC, down, 1,   &
                  comm, status, ierr)

call mpi_sendrecv(dvdt(1, 1, nz - 1), ldx*nynpy, MPI_RPREC, up,   2,   &
                  dvdt(1, 1, 0),      ldx*nynpy, MPI_RPREC, down, 2,   &
                  comm, status, ierr)

call mpi_sendrecv(dwdt(1, 1, nz - 1), ldx*nynpy, MPI_RPREC, up,   3,   &
                  dwdt(1, 1, 0),      ldx*nynpy, MPI_RPREC, down, 3,   &
                  comm, status, ierr)
!--------
call filter_data(u)
call filter_data(v)
call filter_data(w)
! call avg_stats() !--only does something once every n_avg_stats steps
if ( use_avgslice .and. mod(jt, c_count)==0 ) then
    call flow_slice()
    
    if ( theta_flag .and. (jt .ge. theta_init_time) ) then
        call theta_slice()
    end if
        
    if ( PCon_flag .and. (jt .ge. PCon_init) ) then
        call pcon_slice()
    end if
end if
! ---
!call screen_diagnose()
if ( flag_restart .and. mod(jt, nt_restart) == 0 ) then
    call mpi_sendrecv(RHSx(1, 1, 1),  ldx*nynpy, MPI_RPREC, down, 4, &
                      RHSx(1, 1, nz), ldx*nynpy, MPI_RPREC, up,   4, &
                      comm, status, ierr)
    call mpi_sendrecv(RHSy(1, 1, 1),  ldx*nynpy, MPI_RPREC, down, 5, &
                      RHSy(1, 1, nz), ldx*nynpy, MPI_RPREC, up,   5, &
                      comm, status, ierr)
    call mpi_sendrecv(RHSz(1, 1, 1),  ldx*nynpy, MPI_RPREC, down, 6, &
                      RHSz(1, 1, nz), ldx*nynpy, MPI_RPREC, up,   6, &
                      comm, status, ierr)
    call output_final()
end if

if (flag_opt_3d) then
    if ( nums .ge. opt_3d_start .and. mod(nums, base_time) == 0 ) then
        if (flag_opt_frac) then
            call output_field_fraction(nums)
        else
            call output_field(nums)
        end if
    end if
end if

if ( output_video ) then
    if ( nums .ge. video_start .and. nums .le. video_end .and. mod(nums, video_freq) == 0 ) then
        call check_field(nums)
    end if
end if

if ( inflow .and. write_inflow_file ) then
    if (nums .ge. jt_start_write) then
        call output_fringe(nums-jt_start_write)
    end if
end if

if (.not. inflow .and. read_inflow_file ) then
    ! --- can be used to input inflow tracer profile, but without fluctuations
    if (PCon_flag .and. .not. restart_pcon.and. (jt .ge. PCon_init)) then	
        do ipcon = 1, npcon	
            do jz = 1, nz-1		
                PCon(1, 1:nynpy, jz, ipcon) = nutrient_inflow(jz, ipcon)		
            end do	
        end do	
    end if
end if
! ---
end subroutine flosol
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
subroutine newstep
!-----------------------------------------------------------------------
!   prepare for next step
!-----------------------------------------------------------------------
use types, only: rprec
use param, only: tt, dt
use sim_param, only: RHSx, RHSy, RHSz, RHSx_f, RHSy_f, RHSz_f
implicit none
! ---
real(kind=rprec) :: CFL, visc_stab
! ---
call check_cfl(CFL, visc_stab)

tt = tt + dt

! --- save previous time's right-hand-sides for Adams-Bashforth Integration
!  (In subroutine "STEP" use first order time advancement on first time step).
RHSx_f = RHSx
RHSy_f = RHSy
RHSz_f = RHSz
! ---
end subroutine newstep
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
subroutine Admas_Bashforth
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
use types, only:rprec
use param
use sim_param
! use scalars_module, only: scale_us
implicit none
! ---
real(kind=rprec) :: force 
! ---
if ( jt == 1 .and. (.not. restart_vel) ) then
    ! if restart_vel, then this is read from the initialization file
    ! else for the first step put RHS_f=RHS
    !--i.e. at first step, take an Euler step
    RHSx_f = RHSx
    RHSy_f = RHSy
    RHSz_f = RHSz
end if

!--calculate u^(*) (intermediate vel field)
!  at this stage, p, dpdx_i are from previous time step
!  (assumes old dpdx has NOT been added to RHSx_f, etc)
!  we add force (mean press forcing) here so that u^(*) is as close
!  to the final velocity as possible

if (use_mean_p_force) then
    force = mean_p_force
else
    force = 0._rprec
end if
! ! --- Added for the Ragweed field
! if (rag06 .and. PCon_FLAG) force = force*(scale_us**2)

if (use_force_angle) then
    u(:, :, 1:nz - 1) = u(:, :, 1:nz - 1) + dt *    &
                    (tadv1*RHSx(:, :, 1:nz - 1) + tadv2*RHSx_f(:, :, 1:nz - 1) + &
                     dragx(:, :, 1:nz - 1) + force*cos(force_angle))
    v(:, :, 1:nz - 1) = v(:, :, 1:nz - 1) + dt *    &
                    (tadv1*RHSy(:, :, 1:nz - 1) + tadv2*RHSy_f(:, :, 1:nz - 1) + &
                     dragy(:, :, 1:nz - 1) + force*sin(force_angle))
    w(:, :, 1:nz - 1) = w(:, :, 1:nz - 1) + dt *    &
                    (tadv1*RHSz(:, :, 1:nz - 1) + tadv2*RHSz_f(:, :, 1:nz - 1) + &
                     dragz(:, :, 1:nz - 1))

    dudt(:, :, 1:nz - 1) = dudt(:, :, 1:nz - 1) + dragx(:, :, 1:nz - 1) + mean_p_force*cos(force_angle)
    dvdt(:, :, 1:nz - 1) = dvdt(:, :, 1:nz - 1) + dragy(:, :, 1:nz - 1) + mean_p_force*sin(force_angle)
    dwdt(:, :, 1:nz - 1) = dwdt(:, :, 1:nz - 1) + dragz(:, :, 1:nz - 1)
else
    u(:, :, 1:nz - 1) = u(:, :, 1:nz - 1) + dt *    &
                    (tadv1*RHSx(:, :, 1:nz - 1) + tadv2*RHSx_f(:, :, 1:nz - 1) + &
                     dragx(:, :, 1:nz - 1) + force)
    v(:, :, 1:nz - 1) = v(:, :, 1:nz - 1) + dt *    &
                    (tadv1*RHSy(:, :, 1:nz - 1) + tadv2*RHSy_f(:, :, 1:nz - 1) + &
                     dragy(:, :, 1:nz - 1))
    w(:, :, 1:nz - 1) = w(:, :, 1:nz - 1) + dt *    &
                    (tadv1*RHSz(:, :, 1:nz - 1) + tadv2*RHSz_f(:, :, 1:nz - 1) + &
                     dragz(:, :, 1:nz - 1))
                    
    dudt(:, :, 1:nz - 1) = dudt(:, :, 1:nz - 1) + dragx(:, :, 1:nz - 1) + force
    dvdt(:, :, 1:nz - 1) = dvdt(:, :, 1:nz - 1) + dragy(:, :, 1:nz - 1)
    dwdt(:, :, 1:nz - 1) = dwdt(:, :, 1:nz - 1) + dragz(:, :, 1:nz - 1)
end if

! --- after this point, u,v,w at jz = 0 are not useful, until updated
u(:, :, 0) = BOGUS
v(:, :, 0) = BOGUS
w(:, :, 0) = BOGUS
! ---
end subroutine Admas_Bashforth
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine energy
!--------------------------------------------------------------------!
!   calculate the total kinetic energy                                      
!--------------------------------------------------------------------!    
use types, only:rprec
use param
use sim_param, only:u, v, w, ke
implicit none
! ---
integer :: jx, jy, jz, jz_min
real(kind=rprec) :: kea, denom, temp_w, ke_global
! ---
jz_min = 1

kea = 0._rprec
denom = nxt * nynpy * (nz - jz_min)
    
do jz = jz_min, nz - 1
do jy = 1, nynpy
do jx = 1, nxt
    temp_w = .5_rprec*(w(jx, jy, jz) + w(jx, jy, jz + 1))
    !--MPI: assumes w(jz=0) is in sync here
    ke(jx, jy, jz) = .5_rprec*(u(jx, jy, jz)**2 + v(jx, jy, jz)**2 + temp_w**2)
    kea = kea + ke(jx, jy, jz)/denom
end do
end do
end do
! ---
call mpi_reduce(kea, ke_global, 1, MPI_RPREC, MPI_SUM, 0, comm, ierr)
if (rank == 0) then
    open(13,file=path//'output/check_ke.out',status="unknown",position="append")
    kea = ke_global/np
    write(13, *) nums, kea
    close(13)
end if
! ---
end subroutine energy
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------!  