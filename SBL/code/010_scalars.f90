!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
module scalars_module
!-----------------------------------------------------------------------
! HUMIDITY subroutines in place but not yet turned on !!
! Modified to fit namelist feature, use allocatable attribution and
! >> remove the parameter attribution (Bicheng Chen 06/12/2016)
!-----------------------------------------------------------------------
use types, only:rprec
use param
use stokes_drift, only: ust, vst, ust_x, vst_y
!!BC add k_ttlCR by Bicheng Chen for chemical reaction with background species
use sim_param,only:u,v,w,theta,q,dudt,dvdt,dwdt, pcon, k_ttlCR
use bottombc !Includes patches subroutine
use sgsmodule,only:Nu_t,magS
use derivatives
implicit none
save
!public
! ---
integer, parameter:: tag_counter = 200
logical, parameter:: SCALAR_DEBUG=.false.

real(rprec), dimension(:, :, :), allocatable :: beta_scal, Pr_
real(rprec), dimension(:, :, :), allocatable :: dTdz, dqdz ! Only ones needed for output
!! Might need to add x and y derivatives here in case they need to be outputted
!! Right now they are in the "scalar"_all_in_one routines below !!
real(rprec), dimension(:, :, :), allocatable :: RHS_Tf, RHS_T, RHS_qf, RHS_q
real(rprec), dimension(:, :, :), allocatable :: sgs_t3, sgs_q3 !defines the surface sgs flux

integer :: len_dcon, rec_dcon
integer, parameter :: fid_dcon = 257
character(80) :: fn_dcon

real(rprec), dimension(:, :), allocatable :: L, wstar !defines obukhov length and convective vel scale, w_star
real(rprec), dimension(:, :), allocatable :: T_s_filtered !filtered T_s for calc of wT_s

! Now define local u_scale in bottombc.f90 and calculate in wallstress.f90 and use that value
! everywhere else
integer, parameter:: obukhov_output=0 !Controls whether obukhov variables are outputted by theta_slice
integer, parameter:: wt_s_vector_dim1=no_days*86400/300+1

!real(kind=rprec),dimension(floor(no_days*86400._rprec/300._rprec)+1,1) :: wt_s_vector
real(kind=rprec),dimension(wt_s_vector_dim1,1) :: wt_s_vector
! Variables for heterogeneity analysis
! hetero_array_freqz = number of time steps equivalent to 20 seconds
!integer,parameter:: hetero_array_freqz=int(20/dt_dim),hetero_count_out=p_count
integer,parameter :: hetero_array_freqz=100
integer :: hetero_count_out
integer,save::time_ind

! Variables added for pollen
! Chamecki - 08/01/2006
real(rprec), dimension(:, :, :), allocatable  :: Kc_t    ! 3D matrix of SGS diffusivity for pollen
real(rprec), dimension(:, :, :, :), allocatable  :: dPCondz    ! Vertical derivatives (needed also for output)

real(rprec), dimension(:, :), allocatable   :: P_surf_flux   ! Surface pollen flux
real(rprec), dimension(:, :), allocatable   :: deposition    ! Surface pollen deposition (this is actually the net flux=deposition-source)
real(rprec), dimension(:, :), allocatable   :: P_surf_flux_dep  ! Surface pollen flux for Cr=0 everywhere
real(rprec), dimension(:, :), allocatable   :: Real_dep      ! This is the real deposition (using Cr=0 everywhere)
REAL(rprec) :: flux_out, flux_out_prev

! Vars for interpolation of velocity field
real(rprec), dimension(:, :), allocatable :: matrix_x      ! Matrix for x direction
real(rprec), dimension(:, :), allocatable :: matrix_y      ! Matrix for y direction
real(rprec), dimension(:), allocatable :: dvector_x     ! Derivative vector for x direction
real(rprec), dimension(:), allocatable :: dvector_y     ! Derivative vector for y direction

!BC modified to use new parameter file (Bicheng Chen 05/29/2015)
!REAL(kind=rprec),DIMENSION(ldx,nyt,$lbz:nz,npcon)  :: RHS_PConf,RHS_PCon ! RHS for PCon equation
real(rprec), dimension(:, :, :, :), allocatable :: RHS_PConf, RHS_PCon ! RHS for PCon equation
real(rprec), dimension(:, :, :, :), allocatable  :: sgs_PCon3     ! Defines the sgs vertical flux
real(rprec), dimension(:, :, :, :), allocatable  :: res_PCon3     ! Defines the resolved vertical flux
real(rprec), dimension(:, :, :, :), allocatable  :: sgs_PCon1     ! Defines the sgs x-direction flux
real(rprec), dimension(:, :, :), allocatable  :: Cs2Sc         ! Dynamic coefficient for scalar equation
!
! Variables for pollen balance
real(kind=rprec):: released,airborne,deposited,gone

! Scaling factors for time-varying u* and Co
real(kind=rprec):: scale_us, scale_Co
real(kind=rprec), dimension(3) :: scale_Co3
real(rprec), dimension(:, :, :), allocatable :: beta_pcon

real(rprec), dimension(:, :, :), allocatable :: scalar_x,scalar_y,scalar_z      ! Interpolated scalar in x, y and z directions
real(rprec), dimension(:, :, :), allocatable :: u_int,v_int,w_int 	            ! Interpolated velocity field for convective term
real(rprec), dimension(:, :), allocatable :: ghost_x0,ghost_xLx                 ! Ghost nodes at x=-dx/2 and x=lx_tot+dx/2
real(rprec), dimension(:, :), allocatable :: ghost_y0,ghost_yLy                 ! Ghost nodes at y=-dy/2 and y=ly_tot+dy/2
real(rprec), dimension(:, :), allocatable :: ghost_z0    

real(rprec), dimension(:, :, :), allocatable :: magS_int                        ! Interpolated |S|   
real(rprec) :: delta                        ! Filter size
    
!!!!!!!--------------------------------------------
!! Part I of the scalar files - contains the basic subroutines
!! Also look at scalars_module2.f90 for other subroutines !!
!! CONTAINS subroutines :
!! scalar_RHS_calc - computes the RHS side of the scalar evolution equation
!! calc_beta - computes the buoyancy term for temperature
!! step_scalar - time evolves the scalar
!! obukhov - computes the obukhov similarity terms for use in scalars,wallstress and derivwall
!! Authored by Vijayant Kumar
!! Last modified - April 24, 2004
!!!!!!!--------------------------------------------
!
contains
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine theta_all_in_one
!-----------------------------------------------------------------------
!   Performs the various steps for theta calculation
!-----------------------------------------------------------------------   
    use topbc, only:sponge
    use test_filtermodule
    implicit none
    ! ---  
    real(kind=rprec) :: wt_s_current, theta_mn
    integer :: jz
    real(kind=rprec), dimension(ldx,nynpy,0:nz) :: dTdx, dTdy

    if ( coordz == 0 ) then
        !- Add dynamical heat flux
        if (flag_dynWT) then
            open(unit=fid_wt, file=fn_wt, access='direct',recl=len_wt)
            read(unit=fid_wt, rec=nums+ttshift_dynWT+1) wt_s
            close(fid_wt)
        end if
        wt_s_current = wt_s
    end if

    ! Right now set the Prandtl num matrix equal to a constant Prandtl
    Pr_ = Pr

    call filter_data(theta)
    call ddx(dTdx, theta)
    call ddy(dTdy, theta)
    call ddz_uv(dTdz, theta)

    if ( coordz == npz-1 ) then
        dTdz(:,:,nz) = dTdz_top/T_scale*z_i ! Valid for temperature
    end if

    ! Need to synchronize w and dTdz across the processors for uvp node
    ! based computation in scalar_RHS_calc (calculation of RHS_m)
    call mpi_sendrecv(w(1, 1, 1),  ldx*nynpy, MPI_RPREC, down, tag_counter+1,       &
                      w(1, 1, nz), ldx*nynpy, MPI_RPREC, up,   tag_counter+1,       &
                      comm, status, ierr)
    call mpi_sendrecv(dTdz(1, 1, 1),  ldx*nynpy, MPI_RPREC, down, tag_counter+3,    &
                      dTdz(1, 1, nz), ldx*nynpy, MPI_RPREC, up,   tag_counter+3,    &
                      comm, status, ierr)

    ! Also need to synchronize Nu_t across the processors for computation in scalar_RHS_calc (calculation of RHS_m)
    call mpi_sendrecv(Nu_t(1, 1, 1),  ldx*nynpy, MPI_RPREC, down, tag_counter+3,    &
                      Nu_t(1, 1, nz), ldx*nynpy, MPI_RPREC, up,   tag_counter+3,    &
                      comm, status, ierr)
    call mpi_sendrecv(Nu_t(1, 1, nz-1), ldx*nynpy, MPI_RPREC, up,   tag_counter+4,  &
                      Nu_t(1, 1, 0),    ldx*nynpy, MPI_RPREC, down, tag_counter+4,  &
                      comm, status, ierr)

    RHS_Tf = RHS_T

    ! Perform test filtering of T_s for calculation of surf fluxes
    if ( (nums .eq. theta_init_time) .and. (lbc .eq. 0) ) then
        !print *,'T_s b4 filtering',sqrt(sum((T_s-sum(T_s)/float(nxt*nyt))**2))/float(nxt*nyt)
        T_s_filtered(1:nxt,1:nynpy) = T_s
        call test_filter(T_s_filtered, G_test)
        T_s = T_s_filtered(1:nxt,1:nynpy)
        !print *,'T_s after filtering',sqrt(sum((T_s-sum(T_s)/float(nxt*nyt))**2))/float(nxt*nyt)
    end if

    call scalar_RHS_calc(theta,dTdx,dTdy,dTdz,T_s,z_os,RHS_T,sgs_t3,wt_s_current)

    if (ubc==1 .and. damping_method/=1) then !add the damping term to the scalar equation
        do jz = 1, nz-1
            call calc_hor_avg3(theta(1:nxt,1:nynpy,jz), theta_mn)
            RHS_T(1:nxt,1:nynpy,jz) = RHS_T(1:nxt,1:nynpy,jz)-0.5_rprec*(sponge(jz)+sponge(jz+1))*&
                            (theta(1:nxt,1:nynpy,jz)-theta_mn)
        end do
    end if

    ! call check_RHST(jt)
    ! Calculates the buoyancy term which gets added to the vertical momentum equation
    call calc_beta(theta)

    if (jt .eq. theta_init_time .and. (.not. restart_theta)) then
        RHS_Tf = RHS_T
    end if

    call step_scalar(theta,RHS_T,RHS_Tf)

    call mpi_sendrecv(theta(1, 1, 1),  ldx*nynpy, MPI_RPREC, down, tag_counter+7,   &
                      theta(1, 1, nz), ldx*nynpy, MPI_RPREC, up,   tag_counter+7,   &
                      comm, status, ierr)
    call mpi_sendrecv(theta(1, 1, nz-1), ldx*nynpy, MPI_RPREC, up,   tag_counter+8,     &
                      theta(1, 1, 0),    ldx*nynpy, MPI_RPREC, down, tag_counter+8,     &
                      comm, status, ierr)
    ! ---
    end subroutine theta_all_in_one
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
    subroutine calc_beta (scalar)
    ! This calculates the buoyancy term (beta_scal) to be added to the vertical
    ! momentum equation for temperature
    ! Authored by Vijayant Kumar
    ! Last updated April 14, 2004
    implicit none
    ! ---
    real(kind=rprec), dimension(ldx,nynpy,0:nz), intent(in) :: scalar
    real(kind=rprec), dimension(0:nz-1) :: scalar_bar
    real(kind=rprec) :: g_hat, above, below
    integer :: i, j, k, jz_min

    !..Non-dimensionalize gravity
    g_hat = g*(z_i/(u_scale**2))
    beta_scal = 0._rprec

    !..Note Beta is stored on W nodes, but Theta is on UVP nodes
    !....We do not time-advance the ground nodes, so start at k=2
    ! VK: Inserted the averaging code inside this file itself
    ! rather than doing it in prof

    scalar_bar = 0._rprec
    call calc_hor_avg(scalar, scalar_bar)
    
    !....We do not time-advance the ground nodes, so start at k=2
    !.. For the MPI case, this means that we start from jz=2 for
    !.. coordz=0 and jz=1 otherwise... enable by an if statment

    if ( coordz == 0 ) then
        jz_min = 2
    else
        jz_min = 1
    end if

    !BC added by Bicheng Chen for thermal expansion coef for ocean
    if (ocean_flag) then
        do k = jz_min, nz-1
        do j = 1, nynpy
        do i = 1, nxt
            above = alpha_w*(scalar(i,j,k)-scalar_bar(k))
            below = alpha_w*(scalar(i,j,k-1)-scalar_bar(k-1))
            beta_scal(i,j,k) = g_hat*(above + below)/2._rprec
        end do
        end do
        end do
    else
        do k = jz_min, nz-1
        do j = 1, nynpy
        do i = 1, nxt
            above = (scalar(i,j,k)-scalar_bar(k))/scalar_bar(k)
            below = (scalar(i,j,k-1)-scalar_bar(k-1))/scalar_bar(k-1)
            beta_scal(i,j,k) = g_hat*(above + below)/2._rprec
        end do
        end do
        end do
    end if
    ! ---
    end subroutine calc_beta
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
    subroutine scalar_RHS_calc(scalar,dsdx,dsdy,dsdz,S_Surf,z_os,RHS,   &
                               sgs_vert,surf_flux_current)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------     
    use test_filtermodule
    use intermediate, only:u_big, v_big, w_big
    implicit none
    ! ---
    real(kind=rprec), dimension(ldx,nynpy,0:nz) :: scalar
    real(kind=rprec), dimension(ldx,nynpy,0:nz) :: dsdx, dsdy, dsdz
    real(kind=rprec), dimension(nxt,nynpy) :: S_Surf, z_os
    real(kind=rprec), dimension(ldx,nynpy,0:nz) :: RHS
    real(kind=rprec), dimension(ldx,nynpy,0:nz) :: sgs_vert
    real(kind=rprec) :: surf_flux_current

    integer :: i, j, k, jz, jz_min, jz_max, ubc_jz
    real(kind=rprec) :: const, crap2
    real(kind=rprec), dimension(ldx,nynpy,0:nz)   :: temp, dtemp
    real(kind=rprec), dimension(nxt2, ny2npy, 0:nz) :: dsdx_m,dsdy_m,dsdz_m
    real(kind=rprec), dimension(nxt2, ny2npy, 0:nz) :: RHS_m
    !cVK - changed the dimensions for RHS_m,u_m etc. to ldx2
    !cVK as otherwise it causes segmentation errors
    real(kind=rprec),dimension(nxt,nynpy) :: ustar_local, surf_flux
    real(kind=rprec),dimension(nxt,nyt) :: surf_flux_temp
    real(kind=rprec),dimension(ldx,nynpy) :: scalar_node_1 ! test filtered and used for computing surface flux
!   real(kind=rprec),dimension(nxt,nyt):: psi_h,phi_h
!   real(kind=rprec),dimension(nxt,nyt):: psi_h,phi_h,ustar2
    real(kind=rprec),dimension (ptypes):: ustar, wt_s2
    character (64) :: fname_hetero

    character(*), parameter :: fmt_5168="(1400(E14.5))"
    character(*), parameter :: fmt_5169="(1400(E17.10))"
    ! ---
    if (coordz == 0) then
        ustar = 0._rprec

        if (patch_flag==1) then
            if (spatial_flux_flag) then
                if (jt .eq. theta_init_time) then
                    if ( coord_ver == 0 ) then
                        open(unit=77,file='../readin/spatial_flux.dat',status='unknown')
                        do j = 1, nyt
                            read(77,fmt_5168) (surf_flux_temp(i,j), i=1,nxt)
                        end do
                        !The readin surf_flux is dimensional - so non-dimensionalize
                        surf_flux_temp = surf_flux_temp/u_scale/T_scale
                    end if
                    
                    do i = 1, nxt
                        call scatter_y(surf_flux_temp(i,:), surf_flux(i,:))
                    end do
                else if (jt .gt. theta_init_time) then
                    surf_flux = sgs_t3(1:nxt,1:nynpy,1)
                end if

                ustar_local = ustar_avg !set ustar as value computed in obukhov

                do j = 1, nynpy
                do i = 1, nxt
                    dsdz(i,j,1) = -phi_h(i,j)*surf_flux(i,j)/   &
                                (ustar_local(i,j)*vonk*dz/2._rprec) !set the gradient at the first point
                end do
                end do
            else    ! SPATIAL FLUX IF BLOCK CONTINUES !block 102 conts.
                do k = 1, ptypes
                    wt_s2(k) = (-1.)**(k+1)*surf_flux_current
                    !VK Creates patches of type ptypes with alternating signs of heat flux
                end do

!c This computes the average value of the scalar
!c S_surf refers to the value assigned to the surface according to
!c routine patches.f based on number of patches, and parameters in
!c dimen.h
! Also calculated is a local averaged value of u_scale from the subgrid
! stress at the wall

                do j = 1, nynpy
                do i = 1, nxt
                do k = 1, ptypes
                    if ( patch(i,j) == k ) then
!                       ustar(patch(i,j))=ustar(patch(i,j))+(txz(i,j,1)**2+&
!                                           tyz(i,j,1)**2)**.25/patchnum(patch(i,j))
                        ustar(patch(i,j)) = ustar(patch(i,j)) +     &
                                            ustar_avg(i,j)/patchnum(patch(i,j))
                    end if
                end do
                end do
                end do

                ! distribute ustar value to patches
                do j = 1, nynpy
                do i = 1, nxt
                do k = 1, ptypes
                    if ( patch(i,j) == k ) then
!                       ustar2(i,j) = ustar(k)
                        ustar_local(i,j) = ustar(k)
                    end if
                end do
                end do
                end do

                ! Compute surface flux and dsdz at z=DZ/2
                do j = 1, nynpy
                do i = 1, nxt
                    ! lbc=1 is used for prescribing the surface flux while
                    ! lbc=0 has been used to prescribe the temperature
                    if ( lbc==1 .and. scalar(1,1,1)<2 ) then
                        do k = 1, ptypes
                            if ( patch(i,j) == k ) then
                                ! The u_scale is coming from dimen.h = Ug for coriolis and is not
                                ! the local u_scale computed from stress at the surface.
                                surf_flux(i,j) = wt_s2(k)/T_scale/u_scale
                            end if
                        end do
                    else if ( lbc==0 .and. scalar(1,1,1)<2 ) then
                        ustar_local = ustar_avg
                        surf_flux(i,j)=(S_Surf(i,j)-scalar(i,j,1))*vonk*ustar_local(i,j)&
                                        /(dlog(dz/(2._rprec*z_os(i,j)))-psi_h(i,j))

                        !  equ. 4.28. in brutsaert
                        ! Impose the gradient value at the first node.
                        ! Note that ustar2 is a variable coming from sgs stresses at the wall,
                        ! ustar2 has a single value for all the nodes in a single patch.
                        ! phi_h is obtained from the routine obukhov.f, surf_flux is as computed above
                        ! and vonk and dz are constants. Everything is in dimensionless form
                    end if

                    ! Now we have the lowest dsdz on the UVP nodes all others on w nodes
                    dsdz(i,j,1) = -phi_h(i,j)*surf_flux(i,j)/(ustar_local(i,j)*vonk*dz/2._rprec)
                end do
                end do
            end if
        end if

        ! Do the heterogeneity analysis
        if ((lbc .eq. 0) .and. (nums .gt. theta_init_time)) then !ONLY FOR SURF BC = surface temperature !IF BLOCK 200
!           if (jt .eq. hetero_array_freqz) then !If Block 200.1 STARTS
            if (nums .eq. theta_init_time+1) then !If Block 200.1 STARTS
                print *,'Hetero vars: hetero_array_freqz,hetero_count_out = ',hetero_array_freqz,hetero_count_out
                if (mod(hetero_count_out,hetero_array_freqz) .NE. 0) then !BLOCK 200.1.1 STARTS
                    print *,'Hetero count out not exactly divisible by hetero_array_freqz .. QUITTING'
                    print *,'Please change the values in scalars_module.f90'
                    STOP
                    time_ind = 1
                end if  ! BLOCK 200.1.1 ENDS
            end if ! BLOCK 200.1 ENDS

            if (mod(jt,hetero_array_freqz) .eq. 0) then !BLOCK 200.2 STARTS
                print *, 'Writing hetero fields out at jt = ',jt
                write (fname_hetero, '(a,i6.6,a)') path//'result/fields_3d/hetero_data',time_ind,'.bin'
                open(77,file=fname_hetero,form='unformatted',position='append')
                write(77) nums,real(scalar_node_1),real(ustar_local),&
                          real(psi_h), real(surf_flux), real(phi_h)
                close(77)
                if (mod(jt,hetero_count_out) .eq. 0) then !BLOCK 200.2.1 STARTS
                    print *, 'Changing file name string counter at jt = ',jt
                    time_ind = time_ind + 1
                end if !BLOCK 200.2.1 ENDS
            end if !BLOCK 200.2 ENDS
        end if !IF BLOCK 200 ENDS
    end if

    do jz = 0, nz
        call zero_padding(dsdx(:,:,jz), dsdx_m(:,:,jz))
        call zero_padding(dsdy(:,:,jz), dsdy_m(:,:,jz))
        call zero_padding(dsdz(:,:,jz), dsdz_m(:,:,jz))
    end do
    
    if ( coordz == 0 ) then
        jz_min = 2
    else
        jz_min = 1
    end if
    
    if ( coordz == npz-1 ) then
        jz_max = nz - 2
    else
        jz_max = nz - 1
    end if

    ubc_jz = nz-1

! Now compute the RHS term of the filtered scalar equation.
! Note that this is the advection term with the scalar as
! the diffusion term has been thrown away. This is done step
! by step for each of the expanded arrays from dealias1 separately
! for the first node & last node AND the rest of the nodes.

! xxxxx ------Comments valid for MPI case only ---------XXXX
! For MPI case, all the nodes have fluid nodes (1:nz-1) except for
! near-the wall processes (2:nz-1 for w nodes and 1:nz-1 for uvp nodes)
! and the top nodes (1:nz)
! The following loop executes from 2:nz-1 for the near-wall process
! and 1:nz-1 for the other processes. Also, the subsequent 2 loops
! take care of the first node for the near-wall process (coordz = 0)
! and the topmost node for the top process (coordz = npz-1).
! Note also that we need to make ghost data available for dTdz and
! w for the topmost node (jz=n) within each process and therefore,
! this synchronization (MPI_SENDRECV()) has been done in the subroutine
! theta_all_in_one ()
! xxxxx --------- MPI Comment block ends ------------------XXXX

    do k = jz_min, nz-1
    do j = 1, ny2npy
    do i = 1, nxt2
        RHS_m(i,j,k) = (u_big(i,j,k)+ust(k))*dsdx_m(i,j,k)+   &
                       (v_big(i,j,k)+vst(k))*dsdy_m(i,j,k)+   &
                       (w_big(i,j,k)*dsdz_m(i,j,k)+w_big(i,j,k+1)*dsdz_m(i,j,k+1))/2._rprec
    end do
    end do
    end do

    if ( coordz == 0 ) then
        do j = 1, ny2npy
        do i = 1, nxt2
            RHS_m(i,j,1) = (u_big(i,j,1)+ust(1))*dsdx_m(i,j,1)&
                         + (v_big(i,j,1)+vst(1))*dsdy_m(i,j,1)&
                         + (0.5_rprec*w_big(i,j,2))*dsdz_m(i,j,2)
        end do
        end do
    end if

    ! if ( coordz == npz-1 ) then
    !     do j = 1, ny2npy
    !     do i = 1, nxt2
    !         RHS_m(i,j,nz-1) = (u_big(i,j,nz-1)+ust(nz-1))*dsdx_m(i,j,nz-1)    &
    !                         + (v_big(i,j,nz-1)+vst(nz-1))*dsdy_m(i,j,nz-1)    &
    !                         + (0.5_rprec*w_big(i,j,nz-1))*dsdz_m(i,j,nz-1)
    !     end do
    !     end do
    ! end if

    do jz = 1, nz - 1
        call zero_unpadding(RHS_m(:,:,jz), RHS(:,:,jz))
    end do

!c...Now buildxing the SGS part of the RHS.
! Here the sgs_term for scalars is built up using Nu_t from sgs_stag_W.f
! and dividing it by the turbulent Prandtl # specified in dimen.h
!c....Note: Since we bring the Convective term to RHS its sign changes.
!c....Below "Temp" is used for SGS flux; its divergence is added to RHS
!VK.. Nu_t is on w nodes everywhere except at z=dz/2.. while
!VK dsdx is on uvp nodes.. so, interpolate Nu_t as we want temp to
!VK be on uvp nodes
! All this valid only till Pr_ is a constant..
! This will need a makeover once Pr_ becomes dynamic as well...

    if ( coordz == 0 ) then
    ! at jz=1, both Nu_t and dsdx are on uvp nodes .. no need for interp
        do j = 1, nynpy
        do i = 1, nxt
            temp(i,j,1)=(1./Pr_(i,j,1))*Nu_t(i,j,1)*dsdx(i,j,1)
        end do
        end do
    end if
    
    if ( coordz == npz-1 ) then
    ! at jz=1, both Nu_t and dsdx are on uvp nodes .. no need for interp
        do j = 1, nynpy
        do i = 1, nxt
            temp(i,j,nz-1)=(1./Pr_(i,j,nz-1))*Nu_t(i,j,nz)*dsdx(i,j,nz-1)
        end do
        end do
    end if

    do k = jz_min, jz_max
    do j = 1, nynpy
    do i = 1, nxt
!       temp(i,j,k)=(1./Pr_(i,j,k))*Nu_t(i,j,k)*dsdx(i,j,k)
        temp(i,j,k)=(1./Pr_(i,j,k))*0.5_rprec*(Nu_t(i,j,k)+Nu_t(i,j,k+1))*dsdx(i,j,k)
    end do
    end do
    end do

    call ddx(dtemp, temp)

    if ( coordz == 0 ) then
    ! at jz=1, both Nu_t and dsdy are on uvp nodes .. no need for interp
        do j = 1, nynpy
        do i = 1, nxt
            RHS(i,j,1) = (-1.*RHS(i,j,1) + dtemp(i,j,1))
            temp(i,j,1)= (1./Pr_(i,j,1))*Nu_t(i,j,1)*dsdy(i,j,1)
        end do
        end do
    end if
    
    if ( coordz == npz-1 ) then
        do j = 1, nynpy
        do i = 1, nxt
            RHS(i,j,nz-1) = (-1.*RHS(i,j,nz-1) + dtemp(i,j,nz-1))
            temp(i,j,nz-1)= (1./Pr_(i,j,nz-1))*Nu_t(i,j,nz)*dsdy(i,j,nz-1)
        end do
        end do
    end if

    do k = jz_min, jz_max
    ! Nu_t is on w nodes and dsdy is on uvp nodes for jz=2 to nz
        do j = 1, nynpy
        do i = 1, nxt
            RHS(i,j,k) = (-1.*RHS(i,j,k) + dtemp(i,j,k))
            temp(i,j,k)=(1./Pr_(i,j,k))*0.5_rprec*(Nu_t(i,j,k)+Nu_t(i,j,k+1))*dsdy(i,j,k)
!           temp(i,j,k)=(1./Pr_(i,j,k))*Nu_t(i,j,k)*dsdy(i,j,k)
        end do
        end do
    end do

    call ddy(dtemp, temp)

  !...Use MO flux at wall for the scalar sgs term !
  ! Note that the total contribution to the scalar sgs term at
  ! the first node comes from the surface flux computed above from
  ! the specified heat flux, wt_s

    if ( coordz == 0 ) then
        do j = 1, nynpy
        do i = 1, nxt
            RHS(i,j,1) = RHS(i,j,1) + dtemp(i,j,1)
            temp(i,j,1) = -1.*surf_flux(i,j)
            sgs_vert(i,j,1) = surf_flux(i,j)
        end do
        end do
    end if

  ! Note sgs_vert is -1*temp because sgs_vert is modeled as -Nu_t*dsdz/Pr
  ! while temp is the same term but w/o the minus sign due to the additional
  ! minus outside the scalar fluctuation flux term in RHS
  ! need to run this loop nz due to the nature of the differenetiation in ddz_w

    do k = jz_min, nz
    do j = 1, nynpy
    do i = 1, nxt
        RHS(i,j,k) = RHS(i,j,k) + dtemp(i,j,k)
!       temp(i,j,k)=(1./Pr_(i,j,k))*0.5*(Nu_t(i,j,k)+Nu_t(i,j,k-1))*dsdz(i,j,k)
        temp(i,j,k)=(1./Pr_(i,j,k))*Nu_t(i,j,k)*dsdz(i,j,k)
        sgs_vert(i,j,k)=-1.*temp(i,j,k)
    end do
    end do
    end do

  ! The SGS_z flux is on the W nodes, but DDZ_W will put it back on UVP nodes!
  ! Also note that sgs_vert(i,j,k) influences the computations in
  ! OBUKHOV.f and is not involved in any computations in this routine.
  ! sgs_t3(i,j,1) (<w'theta'> is used for computing wt at the surface in OBUKHOV)

    call ddz_w(dtemp, temp)

    do k = 1, ubc_jz
    do j = 1, nynpy
    do i = 1, nxt
        RHS(i,j,k) = RHS(i,j,k) + dtemp(i,j,k)
    end do
    end do
    end do
    ! ---
    end subroutine scalar_RHS_calc
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
    subroutine step_scalar(scalar,RHS_pre,RHS_post)
!-----------------------------------------------------------------------
!cVK - This routine moves the scalar field (scalar in this case)
!cVK - forward in time using the scalar from previous time step
!cVK - and the RHS terms from the previous two time steps
!cVK - using second order Adams-Bashforth scheme   
!-----------------------------------------------------------------------  
    implicit none
    ! ---
    real(kind=rprec), dimension(ldx,nynpy,0:nz) :: scalar, RHS_pre, RHS_post
    integer :: i,j,k

! Note that in the last staments within this file, we set the value
! of scalar at the topmost node based on prescribed bc (inv_strength)
! and so, we couldx save some computation by only performing
! the scalar computation till Nz-1 global node...

    do k = 1, nz-1
    do j = 1, nynpy
    do i = 1, nxt
        scalar(i,j,k) = scalar(i,j,k)+dt*(1.5_rprec*RHS_pre(i,j,k)-0.5_rprec*RHS_post(i,j,k))
    end do
    end do
    end do

!VK Note that this subroutine was designed to be working with a set of scalars (incl.
!VK temperature and humidity and so, these boundary conditions as given below shouldx
!VK be interpreted in the right context and not just for temperature
!VK For example, the first if block refers to a condition with humidity while the
!VK second and third statements are focussed to temperature

    if ( coordz == npz-1 ) then
! if MPI - then clicks and is used for the process dealing wih the top nodes
! else in general is used for temp bc
! Note that L_z*npz is the total vertical extent of the domain for the MPI and
! non-MPI cases ..(with MPI in place, we can not use L_z>z_i anymore)
        if ( (lz*npz) > 1._rprec ) then ! for temperature and non-neutral case
!     if (scalar(2,2,2)<2._rprecand.L_z>z_i) then ! for temperature and non-neutral case
!         print *,'setting scalar value at topmost global node for coordz = ',coordz
!         scalar(:,:,Nz)=scalar(:,:,Nz-1)+inv_strength/T_scale*z_i*dz !ubc
            scalar(:,:,nz) = scalar(:,:,nz-1) + dTdz_top/T_scale*z_i*dz !ubc
! inv_strength - refers to the slope of the inversion ~ 0.003K/Km for temperature
        else ! for everything else - neutral and passive scalars (may be modified depending on need)
            scalar(:,:,nz) = scalar(:,:,nz-1)
        end if
    end if
    ! ---
    end subroutine step_scalar
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
    subroutine obukhov
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
    use types, only: rprec
    use param, only: path, walltype
    use sim_param, only: u, v, theta, dudz
    use bottombc !Includes patches subroutine, phi_m,psi_m,phi_h,psi_h,ustar_avg
    use test_filtermodule
    implicit none
    ! ---
    real(kind=rprec), dimension(ldx,nynpy) :: wt_avg, theta_avg, u1, v1
    real(kind=rprec), dimension(nxt,nynpy) :: x, zeta, u_avg ! wstar, L already defined above
    real(kind=rprec) :: g_,wt_,ustar_,theta_,L_,zo_,fr,wstar_avg
    real(kind=rprec), save :: obuk_L,obuk_ustar,obuk_phi_m,obuk_phi_h,obuk_psi_m
    real(kind=rprec), save :: obuk_wt_sfc,obuk_psi_h,obuk_zo,obuk_wstar
    
    integer :: jx, jy
    real(kind=rprec) :: U_tmp, psi_m_tmp, phi_m_tmp, phi_h_tmp, psi_h_tmp, T_s_tmp
    character(*), parameter :: fmt_5168 = "(E14.5,9(1x,E14.5))"
    ! ---
    if ( ocean_flag ) then
        if (coordz == npz-1 .and. ubc_mom == 'wall') then
            psi_m = 0._rprec !psi_m is being read in initial.f90
            phi_m = 1._rprec
            psi_h = 0._rprec
            phi_h = 1._rprec
            
            select case (walltype)
            case('smooth')
                ustar_avg = sum( sqrt(u(1:nxt,1:nynpy,nz-1)**2+v(1:nxt,1:nynpy,nz-1)**2) ) / float(nxt*nynpy)   &
                        / ( dlog(0.5_rprec*dz/(nu_molec/u_scale/z_i))/vonK + B )
            case('rough')
                ustar_avg = sum( sqrt(u(1:nxt,1:nynpy,nz-1)**2+v(1:nxt,1:nynpy,nz-1)**2) ) / float(nxt*nynpy)   &
                        * vonK / dlog(0.5_rprec*dz/zo)
            case default
                write(*,*) 'invalid wall type'
                stop
            end select
            L = 0._rprec
            wstar = 0._rprec
        end if
    else
        if ( coordz == 0 ) then
            if ( (.not. theta_flag) .or. jt .lt. theta_init_time ) then
                if (.not. restart_theta) psi_m = 0._rprec !psi_m is being read in initial.f90
                phi_m = 1._rprec
                psi_h = 0._rprec
                phi_h = 1._rprec
                
                select case (walltype)
                case('smooth')
                    ustar_avg = sum( sqrt(u(1:nxt,1:nynpy,1)**2+v(1:nxt,1:nynpy,1)**2) ) / float(nxt*nynpy)   &
                            / ( dlog(0.5_rprec*dz/(nu_molec/u_scale/z_i))/vonK + B )
                case('rough')
                    ustar_avg = sum( sqrt(u(1:nxt,1:nynpy,1)**2+v(1:nxt,1:nynpy,1)**2) ) / float(nxt*nynpy)   &
                            * vonK / ( dlog(0.5_rprec*dz/zo) - sum(psi_m(1:nxt,1:nynpy))/float(nxt*nynpy) )
                case default
                    write(*,*) 'invalid wall type'
                    stop
                end select
                L = 0._rprec
                wstar = 0._rprec
                return
            end if

            ! Do the following only for jt .eq. theta_init_time
            ! This becomes curcial for the supercomputing runs as we need to break the
            ! simulation into smaller chunks and thereby need an accurate way to continue
            ! the simulation from the vel_sc.out file..
            ! therefore, added sgs_t3(:,:,1) i.e. the surface flux to the list of output variables
            ! in io.f90
            if (jt .eq. theta_init_time) then
                obuk_L=0._rprec;obuk_ustar=0._rprec;obuk_wstar=0._rprec;obuk_phi_m=0._rprec
                obuk_phi_h=0._rprec;obuk_psi_m=0._rprec;obuk_psi_h=0._rprec;obuk_zo=0._rprec
                if (.not. restart_theta) then
                    sgs_t3(:,:,1) = wt_s / u_scale / T_scale
                end if
            else if (flag_dynWT .and. (jt > theta_init_time)) then !LV1
                !- Add dynamical heat flux (Bicheng Chen 10/25/2016)
                open(unit=fid_wt, file=fn_wt, access='direct',recl=len_wt)
                read(unit=fid_wt, rec=nums+ttshift_dynWT+1) wt_s
                close(fid_wt)
                sgs_t3(:,:,1) = wt_s / u_scale / T_scale
            end if

            !  nondimensionalize g
            g_ = g/(u_scale**2/z_i)

            theta_avg = theta(:,:,1)
            wt_avg = sgs_t3(:,:,1) ! We need only the surface flux - defined by sgs
            call calc_hor_avg3(dlog(zo(1:nxt,1:nynpy)), zo_)
            zo_ = exp(zo_)
            ! averages over x-y plane @ z = 1
            call calc_hor_avg3(sgs_t3(1:nxt,1:nynpy,1), wt_)   
            call calc_hor_avg3(sqrt(u(1:nxt,1:nynpy,1)**2+v(1:nxt,1:nynpy,1)**2), U_tmp)
            call calc_hor_avg3(psi_m(1:nxt,1:nynpy), psi_m_tmp)

            
            ustar_ = U_tmp*vonK/(dlog(0.5_rprec*dz/zo_)-psi_m_tmp)
            call calc_hor_avg3(theta_avg(1:nxt,1:nynpy), theta_)

            if ((patch_flag==1 .and. num_patch==1) .OR. (OB_homog_flag)) then
                do jx = 1, nxt
                do jy = 1, nynpy
                    wt_avg(jx,jy) = wt_
                    ustar_avg(jx,jy) = ustar_
                    theta_avg(jx,jy) = theta_
                end do
                end do
            else
                u1 = u(:,:,1)
                v1 = v(:,:,1)
                call test_filter(u1,G_test)
                call test_filter(v1,G_test)
                do jx = 1, nxt
                do jy = 1, nynpy
                    u_avg(jx,jy) = sqrt(u1(jx,jy)**2+v1(jx,jy)**2)
                end do
                end do
                ustar_avg(1:nxt,:)=u_avg(:,:)*vonK/(dlog(0.5_rprec*dz/zo(:,:))-psi_m(:,:))
            end if

            ! Compute Obukhov Length
            do jx = 1, nxt
            do jy = 1, nynpy
                L(jx,jy)=-ustar_avg(jx,jy)**3/(vonk*g_/theta_avg(jx,jy)*wt_avg(jx,jy))

                ! w_star is defined as [(g/<T_0>)*<w'T'>*z_i]**(1/3) (refer to
                ! Nieuwstadt et al., Turbulent Shear flows, 1991)
                ! Therefore, for our case, where we are computing the non-dimensional w_star,
                ! the formula transforms to [(g_nd/<T_0_nd>)*<w'T'>_nd]**(1/3) where the suffix
                ! _nd refers to being non-dimensionalized using Ug (coriolis velocity), T_scale (300K)
                ! and z_i
                wstar(jx,jy)=sign((g_/theta_avg(jx,jy)*abs(wt_avg(jx,jy)))**(1./3.),wt_avg(jx,jy))

                ! wstar(jx,jy)=sign((g_/theta_avg(jx,jy)*abs(wt_avg(jx,jy))*z_i)**(1./3.),wt_avg(jx,jy))
                ! The above is the earlier wrong formula where there has been this additional z_i which for
                ! the post-processing means division by 10 as the usual value of z_i=1000 which with cube root on
                ! w_star just becomes a multiplication by 10. So, in case you feel that the value of w_star is quite
                ! high and the data is dated before Dec. 7th, 2004, please make sure to divide by 10
                
                ! --- for unstable conditions
                if ((L(jx,jy) < 0._rprec) .and. (wt_avg(jx,jy) .ne. 0._rprec)) then
                    x(jx,jy)=(1._rprec-16._rprec*dz/2._rprec/L(jx,jy))**.25_rprec
                    psi_m(jx,jy)=2._rprec*dlog((1.+x(jx,jy))/2._rprec)+&
                        dlog((1._rprec+x(jx,jy)**2)/2._rprec)-2._rprec*datan(x(jx,jy))+pi/2._rprec
                    psi_h(jx,jy) = 2._rprec*dlog((1._rprec+x(jx,jy)**2)/2._rprec)
                    phi_m(jx,jy) = x(jx,jy)**(-1)
                    phi_h(jx,jy) = x(jx,jy)**(-2)
                else if ((L(jx,jy) > 0._rprec) .and. (wt_avg(jx,jy) .ne. 0._rprec)) then
                    ! Implementing new formulations for phi and psi for stable case
                    ! using Cheng & Brutsaert (2004): source - Brutsaert's book from
                    ! Marc's Hydrology course
                    ! the new relations are from the GABLS study
                    zeta(jx,jy) = 0.5_rprec*dz/L(jx,jy)

                    ! %%%%%%%%%%%%%%%%%%%%%%%% Formulations used in Brutsaert (1982) %%%%%%%%%%%%%%%%%%%%
                    !  phi_m(jx,jy)=1._rprec+5.0_rprec*zeta(jx,jy)
                    !  phi_h(jx,jy)=1._rprec+5.0_rprec*zeta(jx,jy)
                    !  psi_m(jx,jy)=-1._rprec*5.0_rprec*zeta(jx,jy)
                    !  psi_h(jx,jy)=-1._rprec*5.0_rprec*zeta(jx,jy)

                    phi_m(jx,jy) = 1._rprec+4.8_rprec*zeta(jx,jy)
                    phi_h(jx,jy) = 1._rprec+7.8_rprec*zeta(jx,jy)
                    psi_m(jx,jy) = -1._rprec*4.8_rprec*zeta(jx,jy)
                    psi_h(jx,jy) = -1._rprec*7.8_rprec*zeta(jx,jy)
                else
                    psi_m(jx,jy) = 0._rprec
                    psi_h(jx,jy) = 0._rprec
                    phi_m(jx,jy) = 1._rprec
                    phi_h(jx,jy) = 1._rprec
                end if
            end do
            end do

            L_ = -(ustar_**3)/(vonk*(g_/theta_)*wt_)
            wstar_avg = sign((g_/theta_*abs(wt_))**(1./3.),wt_)
        end if
    end if    
    ! ---
    end subroutine obukhov
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
!! Here is the finite volumes discretization
!! For now, all the routines are specific
!! Afer everything is working, I shouldx clean up all the pollen
!! routines and put them together in a nice and organized final version.
!
!! Main routine for pollen time advance
!! Chamecki - 03/29/2007
!
    subroutine pollen_all_in_one2
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
    use test_filtermodule
    implicit none
    ! ---
    integer :: j, ipcon, ind
    real(kind=rprec) :: sim_time, data_time, data_time_prev
    real(kind=rprec) :: data_ustar, data_ustar_prev
    real(kind=rprec), dimension(3)     :: data_co, data_co_prev
    real(kind=rprec), dimension(nxt,nynpy) :: PCon_sfc_temp

    ! Calculate vertical derivative just as diagnostic variable
    ! Note: in the finite volume formualtion this is never used to update the concentration field
    do ipcon = 1, npcon
        call ddz_uv_con(dPCondz(:,:,:,ipcon),PCon(:,:,:,ipcon))
    end do
    ! Need to synchronize w and dPCondz across the processors for uvp node
    ! based computation in scalar_RHS_calc (calculation of RHS_m)
    call mpi_sendrecv(w(1, 1, 1),  ldx*nynpy, MPI_RPREC, down, tag_counter+1,       &
                      w(1, 1, nz), ldx*nynpy, MPI_RPREC, up,   tag_counter+1,       &
                      comm, status, ierr)
    call mpi_sendrecv(dudt(1, 1, 1),  ldx*nynpy, MPI_RPREC, down, tag_counter+2,    &
                      dudt(1, 1, nz), ldx*nynpy, MPI_RPREC, up,   tag_counter+2,    &
                      comm, status, ierr)
    call mpi_sendrecv(dvdt(1, 1, 1),  ldx*nynpy, MPI_RPREC, down, tag_counter+3,    &
                      dvdt(1, 1, nz), ldx*nynpy, MPI_RPREC, up,   tag_counter+3,    &
                      comm, status, ierr)
    call mpi_sendrecv(dwdt(1, 1, 1),  ldx*nynpy, MPI_RPREC, down, tag_counter+4,    &
                      dwdt(1, 1, nz), ldx*nynpy, MPI_RPREC, up,   tag_counter+4,    &
                      comm, status, ierr)
    ! Also need to synchronize Nu_t across the processors for computation in scalar_RHS_calc (calculation of RHS_m)
    call mpi_sendrecv(Nu_t(1, 1, 1),  ldx*nynpy, MPI_RPREC, down, tag_counter+5,    &
                      Nu_t(1, 1, nz), ldx*nynpy, MPI_RPREC, up,   tag_counter+5,    &
                      comm, status, ierr)
    call mpi_sendrecv(Nu_t(1, 1, nz-1), ldx*nynpy, MPI_RPREC, up,   tag_counter+6,  &
                      Nu_t(1, 1, 0),    ldx*nynpy, MPI_RPREC, down, tag_counter+6,  &
                      comm, status, ierr)
    call mpi_sendrecv(magS(1, 1, 1),  ldx*nynpy, MPI_RPREC, down, tag_counter+7,    &
                      magS(1, 1, nz), ldx*nynpy, MPI_RPREC, up,   tag_counter+7,    &
                      comm, status, ierr)
    call mpi_sendrecv(magS(1, 1, nz-1), ldx*nynpy, MPI_RPREC, up,   tag_counter+8,  &
                      magS(1, 1, 0),    ldx*nynpy, MPI_RPREC, down, tag_counter+8,  &
                      comm, status, ierr)

    if (nums >= ini_src .and. nums < end_src) then
        PCon_sfc_temp = PCon_s
    else
        PCon_sfc_temp = 0._rprec
    end if

    RHS_PConf = RHS_PCon        ! Store previous time step RHS for AB2 time integration

    ! Calculate RHS of equation
    ! The total flux out of the domain is also calculated inside the routine
    call PCon_RHS_calc(PCon,dPCondz,PCon_sfc_temp,zo_PCon,RHS_PCon,sgs_PCon3,res_PCon3,pcon_sfc_flux_nd,sgs_PCon1)

    ! For first time step, RHS(n-1)=RHS(n) for time integration
    if ((jt .eq. PCon_init) .or. (jt .eq. ini_src) .and. (.not. restart_pcon)) then
        RHS_PConf = RHS_PCon
    end if

    if (active_pcon) call calcbeta_pcon(PCon)       ! buoyancy term added to vertical momentum equation
    call step_pollen2(PCon,RHS_PCon,RHS_PConf)      ! Time integration using AB2

    do ipcon = 1, npcon
        call data_exchange_y(PCon(:,:,1:nz,ipcon), periodicbcy)
        call mpi_sendrecv(PCon(:, :, 1,  ipcon), ldx*(nynpy+2), MPI_RPREC, down, tag_counter+7,     &
                          PCon(:, :, nz, ipcon), ldx*(nynpy+2), MPI_RPREC, up,   tag_counter+7,     &
                          comm, status, ierr)
        call mpi_sendrecv(PCon(:, :, nz-1, ipcon), ldx*(nynpy+2), MPI_RPREC, up,   tag_counter+8,   &
                          PCon(:, :, 0, ipcon),    ldx*(nynpy+2), MPI_RPREC, down, tag_counter+8,   &
                          comm, status, ierr)
    end do
    
    ! call check_mass_balance()
    ! ---
    end subroutine pollen_all_in_one2
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
    subroutine check_mass_balance()
!-----------------------------------------------------------------------
!   Check mass balance
!-----------------------------------------------------------------------
    implicit none
    ! ---
    real(kind=rprec) :: airborne_global, deposited_global, gone_global
    integer :: ipcon
    ! --- Pollen Grains released (point source) - uses AB2
    if (nums >= ini_src) then 
        select case (source_type)
        case ("point")
            released = 0.d0
            do ipcon = 1, npcon
                released = released + 1.5_rprec*((nums-ini_src+1)*dt&
                                    *sum(con_src(ipcon)%rls)*dx*dy*dz)  &
                         - 0.5_rprec*((nums-ini_src)*dt&
                                    *sum(con_src(ipcon)%rls)*dx*dy*dz)
            end do
        case ("planar")
            released = 1.5_rprec*((nums-ini_src+1)*dt&
                            *sum(con_src(ipcon)%rls)*nxt*nyt*dx*dy*dz)  &
                     - 0.5_rprec*((nums-ini_src)*dt&   
                            *sum(con_src(ipcon)%rls)*nxt*nyt*dx*dy*dz)

        end select
    end if

    if (nums >= end_src) then
        select case (source_type)
        case ("point")
            released = 0.d0
            do ipcon = 1, npcon
                released = released + ((end_src-ini_src)*dt*sum(con_src(ipcon)%rls)*dx*dy*dz)
            end do
        case ("planar")
            released = ((end_src-ini_src)*dt*sum(con_src(ipcon)%rls)*nxt*nyt*dx*dy*dz)
        end select
    end if
    ! --- Pollen Grains airborne and deposited in the valid domain
    airborne = (sum(PCon(1:nxt,1:nynpy,1:nz-1,1))*dx*dy*dz)
    deposited= (sum(deposition(1:nxt,1:nynpy))*dx*dy)
    
    ! --- Pollen grains out of domain (#)
    if (periodicbcx .and. periodicbcy) then
        gone = 0._rprec
    else
        gone = gone+1.5_rprec*flux_out-0.5_rprec*flux_out_prev
    end if

    call mpi_reduce (airborne,airborne_global,1,MPI_RPREC,MPI_SUM,0,comm,ierr)
    call mpi_reduce (deposited,deposited_global,1,MPI_RPREC,MPI_SUM,0,comm,ierr)
    ! call mpi_reduce (gone, gone_global, 1, MPI_RPREC,MPI_SUM,0,comm,ierr)
    
    if ( rank == 0 ) then
        open(1,file=path//'output/pollen_balance.out',status="unknown",position="append")
        write(1,'(I8,6G18.8)') nums, (released), (airborne_global), &                          
                               (deposited_global), (gone), &                                     
                               ((released-airborne_global-deposited_global-gone))
        close(1)
    end if
    ! ---
    end subroutine check_mass_balance
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
    subroutine calcbeta_pcon (PCon)
    ! This calculates the buoyancy term due to oil concentration (beta_pcon)
    ! to be added to the vertical momentum equation
    ! Authored by Di Yang
    implicit none
    real(kind=rprec),dimension(ldx,0:nynpy+1,0:nz,npcon),intent(in)::PCon
    integer::i, j, k, ipcon   
    real(kind=rprec)::g_hat
    !real(kind=rprec),dimension(npcon):: V_pcon

    !..Non-dimensionalize gravity
    g_hat=g*(z_i/(u_scale**2))
    beta_pcon=0._rprec

    !  do ipcon=1,npcon
    !!!     V_pcon(ipcon) = V_pcon0/2._rprec**(ipcon-1)
    !     V_pcon(ipcon) = V_pcon0
    !  enddo

    do k = 0, nz-1
    do j = 1, nynpy
    do i = 1, nxt
        beta_pcon(i,j,k) = 0._rprec
        do ipcon = 1, npcon
            beta_pcon(i,j,k) = beta_pcon(i,j,k)+g_hat*    &
                    (densratio_pcon(ipcon)-1._rprec)*V_pcon(ipcon)*PCon(i,j,k,ipcon)
        end do
    end do
    end do
    end do
    ! ---
    end subroutine calcbeta_pcon
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
    subroutine PCon_RHS_calc(scalar,dsdz,S_Surf,z_os,RHS,sgs_vert,  &
                              res_vert,sfc_flux,sgs_x)
!-----------------------------------------------------------------------
!   This is the finite volumes version of RHS calculation
!   This routine include both QUICK and SMART schemes with periodic and inflow/outflow bc's.
!   The main reference for the discretization of the advection term is the paper by
!     Waterson and Deconinck (1995) cited by Xie et al JoT (2004).
!-----------------------------------------------------------------------
!DY Changed by Di Yang for pconsgs_acc_flag=.true.
    use sim_param, only: dudx,dvdy,dwdz,dudy,dudz,dvdx,dvdz,  &
                         dwdx,dwdy,dpdx,dpdy,dpdz
    implicit none
    ! --- BC added by Bicheng Chen for index ind, jnd and n_end
    real(kind=rprec), dimension(ldx,0:nynpy+1,0:nz,npcon)  :: scalar
    real(kind=rprec), dimension(ldx,0:nynpy+1,0:nz,npcon)  :: dsdz
    real(kind=rprec), dimension(nxt,nynpy)             :: S_Surf,z_os
    real(kind=rprec), dimension(ldx,0:nynpy+1,0:nz,npcon)  :: RHS
    real(kind=rprec), dimension(ldx,0:nynpy+1,0:nz,npcon)  :: sgs_vert                  ! Store SGS vertical pollen flux
    real(kind=rprec), dimension(ldx,0:nynpy+1,0:nz,npcon)  :: res_vert                  ! Store resolved vertical pollen flux
    real(kind=rprec), dimension(npcon), intent(in) :: sfc_flux  
    real(kind=rprec), dimension(ldx,0:nynpy+1,0:nz,npcon)  :: sgs_x                     ! Store SGS x-direction pollen flux

    real(kind=rprec), dimension(ldx,nynpy+1,nz)        :: u_int0,v_int0,w_int0 	      ! Interpolated velocity field without settling velocity
    real(kind=rprec), dimension(nxt,nynpy)  :: surf_flux_m                  ! Previous surface flux for correct calculation of deposition
    real(kind=rprec), dimension(nxt,nynpy)  :: surf_flux_src_m              ! Previous surface flux for total flux
    integer :: i, j, k, jz, ipcon, jx, jy, jz_min
    ! ---
    allocate (scalar_x(nxt+1,nynpy+1,nz), scalar_y(nxt+1,nynpy+1,nz), scalar_z(nxt+1,nynpy+1,nz))
    allocate (u_int(ldx,nynpy+1,nz), v_int(ldx,nynpy+1,nz), w_int(ldx,nynpy+1,nz))
    allocate (ghost_x0(nynpy,nz-1), ghost_xLx(nynpy,nz-1))
    allocate (ghost_y0(nxt,nz-1), ghost_yLy(nxt,nz-1))
    allocate (ghost_z0(ldx,nynpy))
    allocate (magS_int(ldx,nynpy,0:nz) )
    ! --- Step 1 - Interpolations for velocity
    call vel_intp(u_int0, v_int0, w_int0)
    
    ! This is the new interpolation for |S| - center of finite volumes (interpolation only in z)
    if ( coordz == npz-1 ) then
        ! To avoid problems in the upper boundary
        magS(1:nxt,1:nynpy,nz) = magS(1:nxt,1:nynpy,nz-1)
    end if

    magS_int(1:nxt,1:nynpy,2:nz-1)=(magS(1:nxt,1:nynpy,2:nz-1)+magS(1:nxt,1:nynpy,3:nz))/2._rprec

    if ( coordz == 0 ) then
        ! The first level does not need any interpolation
        magS_int(1:nxt,1:nynpy,1) = magS(1:nxt,1:nynpy,1)
    else
        magS_int(1:nxt,1:nynpy,1) = (magS(1:nxt,1:nynpy,1)+magS(1:nxt,1:nynpy,0))/2._rprec
    end if

    call mpi_sendrecv(magS_int(1, 1, 1),  ldx*nynpy, MPI_RPREC, down, tag_counter+5,    &
                      magS_int(1, 1, nz), ldx*nynpy, MPI_RPREC, up,   tag_counter+5,    &
                      comm, status, ierr)
    call mpi_sendrecv(magS_int(1, 1, nz-1), ldx*nynpy, MPI_RPREC, up,   tag_counter+6,  &
                      magS_int(1, 1, 0),    ldx*nynpy, MPI_RPREC, down, tag_counter+6,  &
                      comm, status, ierr)

    do ipcon = 1, npcon
        ! settling_vel=settling_vel0*2._rprec**(ipcon-1)

        ! Add settling velocity to interpolated vertical velocity
        ! Note: w_int @ k=1 and k=nz is zero (imposed B.C.)
        !       there is no need to add the settling for these two locations since
        !       on the boundaries the settling is accounted for together with the
        !       vertical turbulent diffusion (they add to zero at the top and to the
        !       theoretical equilibrium profile at the surface).

        u_int = u_int0
        v_int = v_int0
        w_int = w_int0  

        if (settling) call vel_settling_correction(ipcon)
        ! 
        ! Step 2 - Prepare bottom boundary condition
        !          
        ! Store previous surface flux to integrate deposition
        surf_flux_m = P_surf_flux
        surf_flux_src_m = P_surf_flux_dep
        call calc_srf_dep(ipcon, scalar(:,1:nynpy,:,:), S_surf, sfc_flux(ipcon), dsdz(:,1:nynpy,:,:))      

        ! For first time step, RHS(n-1)=RHS(n) for time integration
        if ((jt .eq. PCon_init) .or. (jt .eq. ini_src)) then
            surf_flux_m = P_surf_flux
            surf_flux_src_m = P_surf_flux_dep
        end if

        ! Calculate pollen deposition (time integration using AB2 scheme)
        deposition = deposition&
                   + (-dt*(1.5_rprec*P_surf_flux-0.5_rprec*surf_flux_m))

        ! Calculate pollen source
        Real_dep = Real_dep&
                 + (-dt*(1.5_rprec*P_surf_flux_dep-0.5_rprec*surf_flux_src_m))

        ! --- Step 3 - Prepare boundary conditions (face values or ghost nodes)
        call enforce_periodic_bndc(ipcon, scalar)
        
        ! --- Step 4 - Interpolation for particle concentration
        ! Calculate the interpolated concentrations
        ! This is actually where QUICK, SMART, etc are different
        if ( PCon_scheme == 2 ) then            ! QUICK scheme
            call QUICK_scheme(ipcon, scalar)
        else if ( PCon_scheme == 3 ) then       ! SMART scheme
            call SMART_scheme(ipcon, scalar)
        end if

        call sgs_sc_stag(ipcon, scalar)

        ! --- Step 5 - Assemble RHS
        ! Now it is just standard centered finite differences
        ! Add the terms one by one...
        call advection(ipcon, RHS)

        ! Store SGS vertical flux
        res_vert(1:nxt,1:nynpy,1:nz,ipcon)=w_int(1:nxt,1:nynpy,1:nz)*scalar_z(1:nxt,1:nynpy,1:nz)
        call sgs_diffusion(ipcon, scalar, RHS, sgs_vert, sgs_x)

        if(pconsgs_acc_flag) call sgs_acceleration(ipcon, scalar, RHS)
        call calc_source_term(ipcon, scalar, RHS)

        ! Calculate total flux at x=lx_tot for mass balance
        ! Save previous value for time integration
        flux_out_prev = flux_out
        flux_out = 0
    
        do j = 1, nynpy
            flux_out = flux_out +       &
                (sum((u_int(nxt+1,j,1:nz-1)+ust(1:nz-1))*scalar_x(nxt+1,j,1:nz-1))*dt*dy*dz)
        end do

        do i = 1, nxt
            flux_out = flux_out + &
                (sum((v_int(i,1,1:nz-1)+vst(1:nz-1))*scalar_x(i,1,1:nz-1))*dt*dx*dz)
        end do
    end do

    if (flag_chem .and. (jt>=ini_src .and. jt<end_src)) then
        if (flag_bgCon) then
            RHS = RHS - k_ttlCR * PCon
        end if
    end if

    ! ---
    deallocate (scalar_x, scalar_y, scalar_z)
    deallocate (u_int, v_int, w_int)
    deallocate (ghost_x0, ghost_xLx)
    deallocate (ghost_y0, ghost_yLy)
    deallocate (ghost_z0)
    deallocate (magS_int)
    ! ---
    end subroutine PCon_RHS_calc
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
    subroutine step_pollen2(scalar,RHS_pre,RHS_post)
!-----------------------------------------------------------------------
!   This is a simplified version of step_scalar()
!   Chamecki 10/04/2006
!-----------------------------------------------------------------------        
    implicit none
    ! ---
    integer :: i, j, k, ipcon  
    real(kind=rprec), dimension(ldx,0:nynpy+1,0:nz,npcon) :: scalar, RHS_pre, RHS_post
    ! ---
    do ipcon = 1, npcon
        do k = 1, nz-1
        do j = 1, nynpy
        do i = 1, nxt
            scalar(i,j,k,ipcon) = scalar(i,j,k,ipcon) + dt *    &
            (1.5_rprec*RHS_pre(i,j,k,ipcon)-0.5_rprec*RHS_post(i,j,k,ipcon))
        end do
        end do
        end do

        if ( coordz == npz-1 ) then
            scalar(:,:,nz,ipcon) = scalar(:,:,nz-1,ipcon)
        end if
    end do
    ! ---
    end subroutine step_pollen2
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
    subroutine vel_intp(u_int0, v_int0, w_int0)
!-----------------------------------------------------------------------
!   Conservative interpolation
!-----------------------------------------------------------------------
    use sim_param, only: dudx, dvdy
    implicit none
    ! ---
    real(kind=rprec), dimension(ldx,nynpy+1,nz), intent(out) :: u_int0, v_int0, w_int0 
    integer :: i, j, k
    ! real(kind=rprec), dimension(ldx,nz) :: temp1, temp0m
    real(kind=rprec), dimension(nxnpy, nyt) :: dvdy_tmp, v_tmp
    real(kind=rprec), dimension(nxnpy, nyt+1) :: v_int0_tmp 
    ! ---3a - u-component
    ! Assemble inverted matrix for x-direction
    matrix_x = 0._rprec
    do i = 1, (nxt-1)
        do j = 1, i
            matrix_x(i,j) = -(nxt-i)*(1._rprec/nxt)
        end do
    
        do j = (i+1), nxt
            matrix_x(i,j) = i*(1._rprec/nxt)
        end do
    end do
    matrix_x(nxt,:) = 1._rprec/nxt

    ! Calculate u-component of velocity (Loop over all yz-locations)
    do j = 1, nynpy
    do k = 1, nz
        ! Assemble derivative vector
        dvector_x(1:nxt-1) = dx*dudx(1:nxt-1,j,k)
        dvector_x(nxt) = sum(u(1:nxt,j,k))

        ! Calculate velocity field
        do i = 1, nxt
            u_int0(i,j,k) = sum(matrix_x(:,i)*dvector_x(:))
        end do
    end do
    end do

    ! Use periodicity
    u_int0(nxt+1,1:nynpy,1:nz) = u_int0(1,1:nynpy,1:nz)

    ! ---3b - v-component
    ! Assemble inverted matrix for y-direction
    matrix_y = 0._rprec
    do i = 1, (nyt-1)
        do j = 1, i
            matrix_y(i,j) = -(nyt-i)*(1._rprec/nyt)
        end do
    
        do j = (i+1), nyt
            matrix_y(i,j) = i*(1._rprec/nyt)
        end do
    end do
    matrix_y(nyt,:) = 1._rprec/nyt

    ! Calculate v-component of velocity (Loop over all xz-locations)
    do k = 1, nz
        call transpose_phy_y_to_x(dvdy(1:nxt,:,k), dvdy_tmp)
        call transpose_phy_y_to_x(v(1:nxt,:,k), v_tmp)
        do i = 1, nxnpy
            ! Assemble derivative vector
            dvector_y(1:nyt-1) = dy*dvdy_tmp(i,1:nyt-1)
            dvector_y(nyt) = sum(v_tmp(i,1:nyt))

            ! Calculate velocity field
            do j = 1, nyt
                v_int0_tmp(i,j) = sum(matrix_y(:,j)*dvector_y(:))
            end do
        end do
        ! Use periodicity
        v_int0_tmp(1:nxnpy,nyt+1) = v_int0_tmp(1:nxnpy,1)
        call transpose_phy2_x_to_y(v_int0_tmp, v_int0(1:nxt,:,k))
    end do

    ! ---3c - w-component
    ! w is already on the right position
    w_int0(1:nxt,1:nynpy,1:nz) = w(1:nxt,1:nynpy,1:nz)
    ! ---
    end subroutine vel_intp
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
    subroutine vel_settling_correction(ipcon)
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
    implicit none
    ! ---
    integer, intent(in) :: ipcon
    real(rprec) :: temp
    real(kind=rprec), dimension(ldx,nz) :: temp0, temp0m, temp1, temp1m
    ! ---
    w_int(1:nxt,1:nynpy,2:nz-1) = w_int(1:nxt,1:nynpy,2:nz-1)-settling_vel(ipcon)
    if (coordz > 0) then
        w_int(1:nxt,1:nynpy,1)  = w_int(1:nxt,1:nynpy,1)-settling_vel(ipcon)
    end if
    
    if (coordz < npz-1) then
        w_int(1:nxt,1:nynpy,nz) = w_int(1:nxt,1:nynpy,nz)-settling_vel(ipcon)
    end if

    if (PCon_acc) then
        ! Multiply by taup=ws/g (use dimensionless gravity)
        ! Add contribution to u_int (linear interpolation to faces)
        temp = settling_vel(ipcon)/(9.81D0/(u_scale**2/z_i))

        u_int(2:nxt,1:nynpy,1:nz)=u_int(2:nxt,1:nynpy,1:nz)&
            -0.5D0*(dudt(1:nxt-1,1:nynpy,1:nz)+dudt(2:nxt,1:nynpy,1:nz))*temp
        u_int(1,1:nynpy,1:nz)=u_int(1,1:nynpy,1:nz)&
            -0.5D0*(dudt(nxt,1:nynpy,1:nz)+dudt(1,1:nynpy,1:nz))*temp
        u_int(nxt+1,1:nynpy,1:nz)=u_int(nxt+1,1:nynpy,1:nz)&
            -0.5D0*(dudt(nxt,1:nynpy,1:nz)+dudt(1,1:nynpy,1:nz))*temp

        ! Add contribution to v_int (linear interpolation to faces)
        v_int(1:nxt,2:nynpy,1:nz)=v_int(1:nxt,2:nynpy,1:nz)&
            -0.5D0*(dvdt(1:nxt,1:nynpy-1,1:nz)+dvdt(1:nxt,2:nynpy,1:nz))*temp

        temp1(:,:)  = dvdt(:,1,1:nz)
        temp1m(:,:) = dvdt(:,nynpy,1:nz)

        call mpi_sendrecv(temp1m(:,:), ldx*nz, mpi_rprec, north, 1, &
                          temp0(:,:),  ldx*nz, mpi_rprec, south, 1, & 
                          comm, status, ierr)
        call mpi_sendrecv(temp1(:,:),  ldx*nz, mpi_rprec, south, 2, &
                          temp0m(:,:), ldx*nz, mpi_rprec, north, 2, & 
                          comm, status, ierr)
        
        v_int(1:nxt,1,1:nz) = v_int(1:nxt,1,1:nz) - 0.5D0*(temp0(1:nxt,1:nz)+temp1(1:nxt,1:nz))*temp
        v_int(1:nxt,nynpy+1,1:nz) = v_int(1:nxt,nynpy+1,1:nz)-0.5D0*(temp1m(1:nxt,1:nz)+temp0m(1:nxt,1:nz))*temp


        ! Add contribution to w_int (no interpolation needed - see note above about settling velocity)
        w_int(1:nxt,1:nynpy,2:nz-1)=w_int(1:nxt,1:nynpy,2:nz-1)&
            -dwdt(1:nxt,1:nynpy,2:nz-1)*temp

        if ( coordz > 0 ) then
            w_int(1:nxt,1:nynpy,1)=w_int(1:nxt,1:nynpy,1)-dwdt(1:nxt,1:nynpy,1)*temp
        end if

        if ( coordz < npz-1 ) then
            w_int(1:nxt,1:nynpy,nz)=w_int(1:nxt,1:nynpy,nz)-dwdt(1:nxt,1:nynpy,nz)*temp
        end if
    end if
    ! ---
    end subroutine vel_settling_correction
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
    subroutine calc_srf_dep(ipcon, scalar, S_surf, surf_flux_current, dsdz)
!-----------------------------------------------------------------------
! Compute surface flux and at z=dz/2 using similarity theory
!BC modified by BC to improve the computational efficiency 05/28/2015
!- eliminate do loop and use array operations   
!-----------------------------------------------------------------------
    implicit none
    ! ---
    integer, intent(in) :: ipcon
    real(kind=rprec), dimension(ldx,nynpy,0:nz,npcon), intent(in) :: scalar
    real(kind=rprec), dimension(nxt,nynpy), intent(in) :: S_surf
    real(kind=rprec), intent(in) :: surf_flux_current 
    real(kind=rprec), dimension(ldx,nynpy,0:nz,npcon), intent(out) :: dsdz
    ! ---
    real(kind=rprec), dimension(nxt, nynpy) :: term1,term_alpha,term_xi,&
                                               term_tao,term_y,term_r,term_omega     
    real(kind=rprec), dimension(nxt,nynpy)  :: ustar_local
    real(kind=rprec), dimension(ptypes)     :: sfc_flx_patch                ! Surface flux for each patch
    real(kind=rprec), dimension(ptypes)     :: ustar_patch                  ! u* averaged for each patch
    integer :: i, j, k
    ! --- For now the imposed surface flux is the same for all patches
    do k = 1, ptypes
        sfc_flx_patch(k) = surf_flux_current
    end do

    ! Initialize patch-averaged u*
    ustar_patch = 0._rprec

    ! Calculate patch-averaged u*
    do j = 1, nynpy
    do i = 1, nxt
    do k = 1, ptypes
        if (patch(i,j)==k) then
            ustar_patch(patch(i,j)) = ustar_patch(patch(i,j))&
                                    + ustar_avg(i,j)/patchnum(patch(i,j))
        end if
    end do
    end do
    end do

    ! Assign u* to each node based on patch type
    do j = 1, nynpy
    do i = 1, nxt
    do k = 1, ptypes
        if (patch(i,j)==k) then
            ustar_local(i,j) = ustar_patch(k)
        end if
    end do
    end do
    end do

    if (lbcp == 1) then
        ! lbc=1 is used for prescribing the surface flux while
        ! Assign surface flux based on patch type
        do k = 1, ptypes
            where (patch==k)
                P_surf_flux = sfc_flx_patch(k)
            end where
        end do
        ! If imposed flux, src and total are the same
        P_surf_flux_dep = P_surf_flux
    else if (lbcp == 0) then
        ! Use average u* instead of local value based on patch type
        ustar_local = ustar_avg
        ! Calculate surface flux for neutral atmospheric stability
        if (settling) then
            term1 = ((z_os-d0)/(dz/2._rprec-d0))&
                ** (Sc_boun*settling_vel(ipcon)/(vonk*ustar_local))
            P_surf_flux = ((S_Surf*term1-scalar(1:nxt,:,1,ipcon)) &
                * settling_vel(ipcon))/(1._rprec-term1)
            ! Calculate effects of atmospheric stability
            if (wt_s .ne. 0._rprec) then
                term_alpha = Sc_boun*settling_vel(ipcon)/vonk/ustar_local
                term_xi = (dz/2._rprec-d0)/L
                term_omega = 1._rprec
                where (L .lt. 0.0_rprec) ! For unstable atmosphere
                    term_tao = 16._rprec*term_xi*(0.5_rprec-term_alpha)&
                             - (1._rprec+term_alpha)
                    term_y = 2._rprec*term_alpha  &
                           / ((term_tao**2._rprec  &
                           - 4._rprec*term_alpha*16._rprec*term_xi  &
                           * (1._rprec+term_alpha-0.5_rprec))**0.5_rprec  &
                           - term_tao)
                    term_r = term_y**2._rprec/term_alpha  &
                           + (1._rprec-term_y)**2._rprec  &
                           - 0.5_rprec*(16._rprec*term_xi)**2._rprec  &
                           / (1._rprec-16._rprec*term_xi*term_y)**2._rprec  &
                           * term_y**2._rprec / term_alpha  &
                           * (1._rprec-term_y)**2._rprec
                    term_omega = (1._rprec+term_alpha) &
                               **(1._rprec+term_alpha-0.5_rprec)  &
                               * term_r**(-0.5_rprec)  &
                               * (term_y/term_alpha)**term_alpha  &
                               * (1._rprec-term_y)  &
                               * (1._rprec-16._rprec*term_xi*term_y)**(-0.5_rprec)
                else where (L .gt. 0._rprec) ! For stable atmosphere
                    term_omega = 1._rprec+5._rprec*(term_alpha/(term_alpha+1._rprec))*term_xi
                end where
                P_surf_flux = P_surf_flux / term_omega
            end if
            ! Calculate real deposition surface flux for neutral atmospheric
            !- stability
            !- This is the surface flux if there is Cr=0 everywhere in the domain
            !- and is used to characterize the deposition
            P_surf_flux_dep =&
                ((0._rprec-scalar(1:nxt,:,1,ipcon))*settling_vel(ipcon))&
                /(1._rprec-term1)
        else
            term1 = dlog((dz/2._rprec-d0)/(z_os-d0))
            P_surf_flux = ((S_Surf-scalar(1:nxt,:,1,ipcon))*vonk*ustar_local)&
                          /(Sc_boun*term1)
            ! Calculate real deposition surface flux for neutral atmospheric
            !- stability
            !- This is the surface flux if there is Cr=0 everywhere in the domain
            !- and is used to characterize the deposition
            P_surf_flux_dep = ((0._rprec-scalar(1:nxt,:,1,ipcon))&
                            * vonk*ustar_local)/(Sc_boun*term1)
        end if
    end if

    ! Calculate vertical derivative at first vertical node
    ! Note: in the finite volume formualtion this is never used to update the concentration field
    if ( coordz == 0 ) then
        if (settling) then
            dsdz(1:nxt,:,1,ipcon) = (-Sc_boun*P_surf_flux &
                            - settling_vel(ipcon)*scalar(1:nxt,:,1,ipcon))&
                            / (ustar_local*vonk*(dz/2._rprec-d0))
        else
            dsdz(1:nxt,:,1,ipcon) = -(Sc_boun*P_surf_flux) &
                            / (ustar_local*vonk*(dz/2._rprec-d0))
        end if
    end if

    ! Do not allow positive fluxes for point source case
    ! i.e. fluxes from the ground when PCon becomes negative due to Gibbs
    ! if (sourcetype == "point") then
    !     where(P_surf_flux > 0) P_surf_flux = 0._rprec
    ! end if
    ! ---
    end subroutine calc_srf_dep
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
    subroutine enforce_periodic_bndc(ipcon, scalar)
!-----------------------------------------------------------------------
!   For periodic boundary condition in x and y
!   The ghost nodes are just explicit enforcement of periodicity   
!-----------------------------------------------------------------------
    implicit none
    ! ---
    integer, intent(in) :: ipcon   
    real(kind=rprec), dimension(ldx,0:nynpy+1,0:nz,npcon), intent(in) :: scalar
    ! ---
    if (periodicbcx) then
        ghost_x0(1:nynpy,1:nz-1) = scalar(nxt,1:nynpy,1:nz-1,ipcon)
        ghost_xLx(1:nynpy,1:nz-1)= scalar(1,  1:nynpy,1:nz-1,ipcon)
        ! The face values have to be calculated in the interpolation scheme below
    else
        ! For inflow condition/outflow
        ! The ghost nodes enforce zero derivatives
        ! ghost_xLx(1:nyt,1:nz-1)=scalar(nxt,1:nyt,1:nz-1,ipcon)

        ! The face value enforces zero concentration at the inlet
        ! scalar_x(1,1:nyt,1:nz-1)=0._rprec
        ! BC modified by Bicheng Chen, not an inlet
        ! scalar_x(1,1:nyt,1:nz-1)=scalar(1,1:nyt,1:nz-1,ipcon)

        ! The ghost node is also set to be zero
        ! The alternative approach wouldx be to calculate its value based
        ! on linear extrapolation for the advective term and 2nd order
        ! differences for the diffusion term
        ! ghost_x0(1:nyt,1:nz-1)=0._rprec
        ! ghost_x0(1:nyt,1:nz-1)=scalar(1,1:nyt,1:nz-1,ipcon)
        ! END BC modified by Bicheng Chen, not an inlet

        ! The other face values are set based on a "rough" approximation
        ! I have to rethink this approach
        ! scalar_x(nxt+1,1:nyt,1:nz-1)=scalar(nxt,1:nyt,1:nz-1,ipcon)

        ! BC added by Bicheng Chen for determining boundary is a inlet or
        ! outlet; inlet scalar_x=0, outlet scalar_x=scalar_bdy
        ! x0 bdy
        ghost_x0(1:nynpy,1:nz-1) = scalar(1,1:nynpy,1:nz-1,ipcon)*&
            (1-sign(1._rprec, u_int(1,1:nynpy,1:nz-1)+ust_x(1:nynpy,1:nz-1)))/2
        scalar_x(1,1:nynpy,1:nz-1) = ghost_x0(1:nynpy,1:nz-1)
        ! xLx bdy
        ghost_xLx(1:nynpy,1:nz-1) = scalar(nxt,1:nynpy,1:nz-1,ipcon)*&
            (1-sign(1._rprec, -(u_int(1,1:nynpy,1:nz-1)+ust_x(1:nynpy,1:nz-1))))/2
        scalar_x(nxt+1,1:nynpy,1:nz-1) = ghost_xLx(1:nynpy,1:nz-1)
        ! END BC added by Bicheng Chen for determining boundary is a inlet or
        ! outlet; inlet scalar_x=0, outlet scalar_x=scalar_bdy
    end if
    !BC END modified by Bicheng Chen for replicated domain

    !BC modified by Bicheng Chen for replicated domain
    
    
    if ( periodicbcy ) then
        ghost_y0(1:nxt,1:nz-1) = scalar(1:nxt,0,1:nz-1,ipcon)
        ghost_yLy(1:nxt,1:nz-1)= scalar(1:nxt,nynpy+1,1:nz-1,ipcon)
    else
        ghost_y0(1:nxt,1:nz-1) = scalar(1:nxt,1,1:nz-1,ipcon)*&
                (1-sign(1._rprec, v_int(1:nxt,1,1:nz-1)+vst_y(1:nxt,1:nz-1)))/2
        scalar_y(1:nxt,1,1:nz-1) = ghost_y0(1:nxt,1:nz-1)

        ghost_yLy(1:nxt,1:nz-1) = scalar(1:nxt,nynpy,1:nz-1,ipcon)*&
                (1-sign(1._rprec, -(v_int(1:nxt,1,1:nz-1)+vst_y(1:nxt,1:nz-1))))/2
        scalar_y(1:nxt,nynpy+1,1:nz-1) = ghost_yLy(1:nxt,1:nz-1)
        
        ! END BC added by Bicheng Chen for determining boundary is a inlet or
        ! outlet; inlet scalar_y=0, outlet scalar_y=scalar_bdy
    end if

    ghost_z0(1:ldx,1:nynpy) = 0._rprec
    ! ---
    end subroutine enforce_periodic_bndc
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
    subroutine sgs_sc_stag(ipcon, scalar)
!-----------------------------------------------------------------------
!   In these 3 cases assign Kc_t=Nu_t/Sc (Sc from param.f90)
!   Case 1 - initialization (jt<cs_count)
!   Case 2 - No dynamic Sc (model_sc==1)
!   Case 3 - Initial part of the run (jt<DYN_init)
!-----------------------------------------------------------------------
    implicit none
    ! ---
    integer, intent(in) :: ipcon   
    real(kind=rprec), dimension(ldx,0:nynpy+1,0:nz,npcon), intent(in) :: scalar
    ! ---
    if (jt<cs_count .or. model_sc==1 .or. jt<DYN_init) then
        ! This is the new interpolation for nu_t - center of finite volumes (interpolation in z)
        ! The first level does not need any interpolation
        if ( coordz == 0 ) then
            Kc_t(1:nxt,1:nynpy,1) = Nu_t(1:nxt,1:nynpy,1)
        else
            Kc_t(1:nxt,1:nynpy,1) = (Nu_t(1:nxt,1:nynpy,1)+Nu_t(1:nxt,1:nynpy,2))/2._rprec
        end if

        ! From now on, interpolates in z
        Kc_t(1:nxt,1:nynpy,2:nz-2)=(Nu_t(1:nxt,1:nynpy,2:nz-2)+Nu_t(1:nxt,1:nynpy,3:nz-1))/2._rprec

        ! The top one is only copied
        if ( coordz == npz-1 ) then
            Kc_t(1:nxt,1:nynpy,nz-1)=Nu_t(1:nxt,1:nynpy,nz-1)
        else
            Kc_t(1:nxt,1:nynpy,nz-1)=(Nu_t(1:nxt,1:nynpy,nz-1)+Nu_t(1:nxt,1:nynpy,nz))/2._rprec
        end if
        ! Divide by constant Schmidt number
        Kc_t = Kc_t/Sc
    else
        ! Calculate dynamic coefficient Cs2/Sc
        if (mod(jt,cs_count)==0) then
            ! Plane Averaged Scale Invariant (PASI);
            if (model_sc==2) then
                call dynamic_sc_pasi(scalar(:,1:nynpy,1:nz,ipcon),scalar_x,scalar_y,scalar_z,u_int,v_int,w_int)
            ! Plane Averaged Scale Dependent (PASD);
            else if (model_sc==3) then
                call dynamic_sc_pasd(scalar(:,1:nynpy,1:nz,ipcon),scalar_x,scalar_y,scalar_z,u_int,v_int,w_int)
            ! Lagrangian Averaged Scale Invariant (LASI);
            else if (model_sc==4) then
                call dynamic_sc_lasi(scalar(:,1:nynpy,1:nz,ipcon),scalar_x,scalar_y,scalar_z,u_int,v_int,w_int)
            ! Lagragian Averaged Scale Dependent (LASD);
            else if (model_sc==5) THEN
                call dynamic_sc_lasd(scalar(:,1:nynpy,1:nz,ipcon),scalar_x,scalar_y,scalar_z,u_int,v_int,w_int)
            end if
        end if
        ! Filter size
        delta = (dx*dy*dz)**(1._rprec/3._rprec)
        ! Calculate eddy-diffusivity
        Kc_t = Cs2Sc*(delta**2)*magS_int
    end if

    call data_exchange_y(Kc_t(:,:,1:nz), .true.)
    
    call mpi_sendrecv(Kc_t(:, :, 1),  ldx*(nynpy+2), MPI_RPREC, down, tag_counter+9,  &
                      Kc_t(:, :, nz), ldx*(nynpy+2), MPI_RPREC, up,   tag_counter+9,  &
                      comm, status, ierr)
    call mpi_sendrecv(Kc_t(:, :, nz-1), ldx*(nynpy+2), MPI_RPREC, up,   tag_counter+10,   &
                      Kc_t(:, :, 0),    ldx*(nynpy+2), MPI_RPREC, down, tag_counter+10,   &
                      comm, status, ierr)
    ! ---
    end subroutine sgs_sc_stag
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
    subroutine advection(ipcon, RHS)
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
    implicit none
    ! ---
    integer, intent(in) :: ipcon  
    real(kind=rprec), dimension(ldx,0:nynpy+1,0:nz,npcon), intent(inout) :: RHS 
    integer :: i, j, k
    ! --- 1.1) Advection => -u*dC/dx
    do k = 1, nz-1
    do j = 1, nynpy
    do i = 1, nxt
        RHS(i,j,k,ipcon)=-(1._rprec/dx)*(u_int(i+1,j,k)*scalar_x(i+1,j,k)-u_int(i,j,k)*scalar_x(i,j,k)&
                        +ust(k)*(scalar_x(i+1,j,k)-scalar_x(i,j,k)))
    end do
    end do
    end do

    ! 1.2) Advection => -v*dC/dy
    do k = 1, nz-1
    do j = 1, nynpy
    do i = 1, nxt
        RHS(i,j,k,ipcon)=RHS(i,j,k,ipcon)-(1._rprec/dy)*&
            (v_int(i,j+1,k)*scalar_y(i,j+1,k)-v_int(i,j,k)*scalar_y(i,j,k)&
            +vst(k)*(scalar_y(i,j+1,k)-scalar_y(i,j,k)))
    end do
    end do
    end do

    ! 1.3) Advection => -w*dC/dz
    do k = 1, nz-1
    do j = 1, nynpy
    do i = 1, nxt
        RHS(i,j,k,ipcon)=RHS(i,j,k,ipcon)-(1._rprec/dz)*    &
            (w_int(i,j,k+1)*scalar_z(i,j,k+1)-w_int(i,j,k)*scalar_z(i,j,k))
    end do
    end do
    end do
    ! ---
    end subroutine advection
!-----------------------------------------------------------------------
!   
!----------------------------------------------------------------------- 
    subroutine sgs_diffusion(ipcon, scalar, RHS, sgs_vert, sgs_x)
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
    implicit none
    ! ---
    integer, intent(in) :: ipcon   
    real(kind=rprec), dimension(ldx,0:nynpy+1,0:nz,npcon), intent(in) :: scalar
    real(kind=rprec), dimension(ldx,0:nynpy+1,0:nz,npcon), intent(inout) :: RHS 
    real(kind=rprec), dimension(ldx,0:nynpy+1,0:nz,npcon), intent(out) :: sgs_vert                  ! Store SGS vertical pollen flux
    real(kind=rprec), dimension(ldx,0:nynpy+1,0:nz,npcon), intent(out) :: sgs_x 
    integer :: i, j, k
    ! --- 2.1) SGS diffusion => d((nu/Sc)*dc/dx)/dx
    do k = 1, nz-1
    do j = 1, nynpy
        do i = 2, nxt-1
            RHS(i,j,k,ipcon)=RHS(i,j,k,ipcon)+(1._rprec/dx**2)*(((Kc_t(i,j,k)+Kc_t(i+1,j,k))/2._rprec)*(scalar(i+1,j,k,ipcon)-scalar(i,j,k,ipcon)) &
            -((Kc_t(i,j,k)+Kc_t(i-1,j,k))/2._rprec)*(scalar(i,j,k,ipcon)-scalar(i-1,j,k,ipcon)))
        end do
        RHS(1,j,k,ipcon)=RHS(1,j,k,ipcon)+(1._rprec/dx**2)*(((Kc_t(1,j,k)+Kc_t(2,j,k))/2._rprec)*(scalar(2,j,k,ipcon)-scalar(1,j,k,ipcon)) &
            -((Kc_t(1,j,k)+Kc_t(nxt,j,k))/2._rprec)*(scalar(1,j,k,ipcon)-ghost_x0(j,k)))
        RHS(nxt,j,k,ipcon)=RHS(nxt,j,k,ipcon)+(1._rprec/dx**2)*(((Kc_t(nxt,j,k)+Kc_t(1,j,k))/2._rprec)*(ghost_xLx(j,k)-scalar(nxt,j,k,ipcon)) &
            -((Kc_t(nxt,j,k)+Kc_t(nxt-1,j,k))/2._rprec)*(scalar(nxt,j,k,ipcon)-scalar(nxt-1,j,k,ipcon)))
    end do
    end do

    ! 2.2) SGS diffusion => d((nu/Sc)*dc/dy)/dy
    do k = 1, nz-1
    do i = 1, nxt
        ! do j = 2, nynpy-1
        do j = 1, nynpy
            RHS(i,j,k,ipcon)=RHS(i,j,k,ipcon)+(1._rprec/dy**2)*(((Kc_t(i,j,k)+Kc_t(i,j+1,k))/2._rprec)*(scalar(i,j+1,k,ipcon)-scalar(i,j,k,ipcon)) &
            -((Kc_t(i,j,k)+Kc_t(i,j-1,k))/2._rprec)*(scalar(i,j,k,ipcon)-scalar(i,j-1,k,ipcon)))
        end do

        ! RHS(i,1,k,ipcon)=RHS(i,1,k,ipcon)+(1._rprec/dy**2)*(((Kc_t(i,1,k)+Kc_t(i,2,k))/2._rprec)*(scalar(i,2,k,ipcon)-scalar(i,1,k,ipcon)) &
        !     -((Kc_t(i,1,k)+Kc_t(i,nynpy,k))/2._rprec)*(scalar(i,1,k,ipcon)-ghost_y0(i,k)))
        ! RHS(i,nynpy,k,ipcon)=RHS(i,nynpy,k,ipcon)+(1._rprec/dy**2)*(((Kc_t(i,nynpy,k)+Kc_t(i,1,k))/2._rprec)*(ghost_yLy(i,k)-scalar(i,nynpy,k,ipcon)) &
        !     -((Kc_t(i,nynpy,k)+Kc_t(i,nynpy-1,k))/2._rprec)*(scalar(i,nynpy,k,ipcon)-scalar(i,nynpy-1,k,ipcon)))
    end do
    end do

    ! 2.3) SGS diffusion => d((nu/Sc)*dc/dz)/dz
    do j = 1, nynpy
    do i = 1, nxt
        do k = 2, nz-2
            RHS(i,j,k,ipcon)=RHS(i,j,k,ipcon)+(1._rprec/dz**2)*(((Kc_t(i,j,k+1)+Kc_t(i,j,k))/2._rprec)*(scalar(i,j,k+1,ipcon)-scalar(i,j,k,ipcon)) &
            -((Kc_t(i,j,k)+Kc_t(i,j,k-1))/2._rprec)*(scalar(i,j,k,ipcon)-scalar(i,j,k-1,ipcon)))

            ! Store SGS vertical flux
            sgs_vert(i,j,k,ipcon)=-(1._rprec/dz)*((Kc_t(i,j,k)+Kc_t(i,j,k-1))/2._rprec)*(scalar(i,j,k,ipcon)-scalar(i,j,k-1,ipcon))
        end do

        if ( coordz == 0 ) then
            ! Impose surface flux at bottom of element 1
            RHS(i,j,1,ipcon)=RHS(i,j,1,ipcon)+(1._rprec/dz**2)*(((Kc_t(i,j,2)+Kc_t(i,j,1))/2._rprec)*(scalar(i,j,2,ipcon)-scalar(i,j,1,ipcon))) &
            +(1._rprec/dz)*(P_surf_flux(i,j))
            ! Store surface flux
            sgs_vert(i,j,1,ipcon)=P_surf_flux(i,j)
        else
            RHS(i,j,1,ipcon)=RHS(i,j,1,ipcon)+(1._rprec/dz**2)*(((Kc_t(i,j,2)+Kc_t(i,j,1))/2._rprec)*(scalar(i,j,2,ipcon)-scalar(i,j,1,ipcon)) &
            -((Kc_t(i,j,1)+Kc_t(i,j,0))/2._rprec)*(scalar(i,j,1,ipcon)-scalar(i,j,0,ipcon)))
            sgs_vert(i,j,1,ipcon)=-(1._rprec/dz)*((Kc_t(i,j,1)+Kc_t(i,j,0))/2._rprec)*(scalar(i,j,1,ipcon)-scalar(i,j,0,ipcon))
        end if

        if ( coordz == npz-1 ) then
        ! Impose zero flux at top of element nz-1
            RHS(i,j,nz-1,ipcon)=RHS(i,j,nz-1,ipcon)+(1._rprec/dz**2)*(0._rprec &
            -((Kc_t(i,j,nz-1)+Kc_t(i,j,nz-2))/2._rprec)*(scalar(i,j,nz-1,ipcon)-scalar(i,j,nz-2,ipcon)))
        else
            RHS(i,j,nz-1,ipcon)=RHS(i,j,nz-1,ipcon)+(1._rprec/dz**2)*(((Kc_t(i,j,nz)+Kc_t(i,j,nz-1))/2._rprec)*(scalar(i,j,nz,ipcon)-scalar(i,j,nz-1,ipcon)) &
            -((Kc_t(i,j,nz-1)+Kc_t(i,j,nz-2))/2._rprec)*(scalar(i,j,nz-1,ipcon)-scalar(i,j,nz-2,ipcon)))
        end if

        sgs_vert(i,j,nz-1,ipcon)=-(1._rprec/dz)*((Kc_t(i,j,nz-1)+Kc_t(i,j,nz-2))/2._rprec)*(scalar(i,j,nz-1,ipcon)-scalar(i,j,nz-2,ipcon))

        do k = 1, nz
            ! Store SGS x-direction flux
            sgs_x(i,j,k,ipcon)=-(1._rprec/dx)*((Kc_t(i+1,j,k)+Kc_t(i,j,k)/2._rprec)*(scalar(i+1,j,k,ipcon)-scalar(i,j,k,ipcon)))
        end do
    end do
    end do
    ! ---
    end subroutine sgs_diffusion
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
    subroutine sgs_acceleration(ipcon, scalar, RHS)
!-----------------------------------------------------------------------
!   Added by Di Yang for sgs acceleration
!-----------------------------------------------------------------------
    use sim_param, only: dudx,dvdy,dwdz,dudy,dudz,dvdx,dvdz,  &
                         dwdx,dwdy,dpdx,dpdy,dpdz
    implicit none
    ! ---
    integer, intent(in) :: ipcon
    real(kind=rprec), dimension(ldx,0:nynpy+1,0:nz,npcon), intent(in) :: scalar
    real(kind=rprec), dimension(ldx,0:nynpy+1,0:nz,npcon), intent(inout) :: RHS 
    ! ---
    real(kind=rprec), dimension(0:nxt+1,0:nynpy+1,0:nz) :: acc_mag         ! acceleration magnitude
    real(kind=rprec), dimension(0:nxt+1,0:nynpy+1,  nz) :: S11, S12, S22, S33, S13, S23
    real(kind=rprec) :: ux, uy, uz, vx, vy, vz, wx, wy, wz
    integer :: i, j, k, jx, jy, jz, jz_min

    real(kind=rprec), dimension(0:nxt+1, 0:nz-1) :: temp1, temp2, temp3, temp4  
    ! ---
    
    delta=(dx*dy*dz)**(1._rprec/3._rprec)   ! Filter size

    select case(model_psgsacc)
    case (1)    
        !DY model 1: based on gradient of acceleration magnitude
        acc_mag(1:nxt,1:nynpy,0:nz-1) = sqrt(dudt(1:nxt,1:nynpy,0:nz-1)**2+dvdt(1:nxt,1:nynpy,0:nz-1)**2 &
            +0.25D0*(dwdt(1:nxt,1:nynpy,0:nz-1)+dwdt(1:nxt,1:nynpy,0+1:nz))**2)
        !--this is only required b/c of the unnatural placement of all strains
        !  onto w-nodes, be careful not to overwrite nz on top process with garbage
        !--dwdz(jz=0) is already known, except at bottom process (OK)
        !call mpi_sendrecv (dwdz(1, 1, nz-1), ldx*nynpy, MPI_RPREC, up, 1,  &
        !                   dwdz(1, 1, 0), ldx*nynpy, MPI_RPREC, down, 1,   &
        !                   comm, status, ierr)
        call mpi_sendrecv(acc_mag(1, 1, 1),  ldx*nynpy, MPI_RPREC, down, 2,   &
                          acc_mag(1, 1, nz), ldx*nynpy, MPI_RPREC, up,   2,   &
                          comm, status, ierr)
        if ( coordz == npz-1 ) then
            acc_mag(1:nxt,1:nynpy,nz) = sqrt(dudt(1:nxt,1:nynpy,nz)**2+dvdt(1:nxt,1:nynpy,nz)**2 &
                +dwdt(1:nxt,1:nynpy,nz)**2)
        end if
        acc_mag(0,1:nynpy,0:nz-1) = acc_mag(nxt,1:nynpy,0:nz-1)
        acc_mag(nxt+1,1:nynpy,0:nz-1) = acc_mag(1,1:nynpy,0:nz-1)
        call data_exchange_y(acc_mag(0:nxt+1,0:nynpy+1,0:nz-1), .true.)

        do k=1,nz-1
        do j=1,nynpy
        do i=1,nxt
            RHS(i,j,k,ipcon)=RHS(i,j,k,ipcon)-((1._rprec/dx**2)*((acc_mag(i+1,j,k)-acc_mag(i,j,k))*scalar_x(i+1,j,k) &
                -(acc_mag(i,j,k)-acc_mag(i-1,j,k))*scalar_x(i,j,k)) &
                +(1._rprec/dy**2)*((acc_mag(i,j+1,k)-acc_mag(i,j,k))*scalar_y(i,j+1,k) &
                -(acc_mag(i,j,k)-acc_mag(i,j-1,k))*scalar_y(i,j,k))) &
                *settling_vel(ipcon)/(9.81D0/(u_scale**2/z_i))*Cs_pa*delta
        end do
        end do
        end do
        !DY Add z direction
        do j=1,nynpy
        do i=1,nxt
            do k=2,nz-1
                RHS(i,j,k,ipcon)=RHS(i,j,k,ipcon)-(1._rprec/dz**2)*((acc_mag(i,j,k+1)-acc_mag(i,j,k))*scalar_z(i,j,k+1) &
                -(acc_mag(i,j,k)-acc_mag(i,j,k-1))*scalar_z(i,j,k)) &
                *settling_vel(ipcon)/(9.81D0/(u_scale**2/z_i))*Cs_pa*delta
            end do
            if ( coordz == 0 ) then
                RHS(i,j,1,ipcon)=RHS(i,j,1,ipcon)-(1._rprec/dz**2)*((acc_mag(i,j,2)-acc_mag(i,j,1))*scalar_z(i,j,2) &
                -(acc_mag(i,j,1)-0.5D0*(3._rprec*acc_mag(i,j,1)-acc_mag(i,j,2)))*scalar_z(i,j,1)) &
                *settling_vel(ipcon)/(9.81D0/(u_scale**2/z_i))*Cs_pa*delta
            else
                RHS(i,j,1,ipcon)=RHS(i,j,1,ipcon)-(1._rprec/dz**2)*((acc_mag(i,j,2)-acc_mag(i,j,1))*scalar_z(i,j,2) &
                -(acc_mag(i,j,1)-acc_mag(i,j,0))*scalar_z(i,j,1)) &
                *settling_vel(ipcon)/(9.81D0/(u_scale**2/z_i))*Cs_pa*delta
            end if
        end do
        end do
    
    case(2)
        !DY model 2: based on gradient of Q (=|S|^2-|omega|^2)/4.
        !DY To save memory, save Q into acc_mag
        if ( coordz == 0 ) then
        ! calculate Sij on w-nodes
        ! calculate |S| on w-nodes
            ! at level z=dz/2.  (Where dudz and dvdz are on UVP nodes)
            do jy=1,nynpy
            do jx=1,nxt
                ux=dudx(jx,jy,1)  ! uvp-node
                uy=dudy(jx,jy,1)  ! uvp-node
                uz=dudz(jx,jy,1)  ! uvp-node
                vx=dvdx(jx,jy,1)  ! uvp-node
                vy=dvdy(jx,jy,1)  ! uvp-node
                vz=dvdz(jx,jy,1)  ! uvp-node
                ! special case
                wx=0.5_rprec*(dwdx(jx,jy,1)+dwdx(jx,jy,2))  ! uvp-node
                wy=0.5_rprec*(dwdy(jx,jy,1)+dwdy(jx,jy,2))  ! uvp-node
                wz=dwdz(jx,jy,1)  ! uvp-node
                S11(jx,jy,1)=ux          ! uvp-node
                S12(jx,jy,1)=0.5_rprec*(uy+vx) ! uvp-node
                ! taken care of with wall stress routine
                S13(jx,jy,1)=0.5_rprec*(uz+wx) ! uvp
                S22(jx,jy,1)=vy          ! uvp-node
                ! taken care of with wall stress routine
                S23(jx,jy,1)=0.5_rprec*(vz+wy) ! uvp
                S33(jx,jy,1)=wz          ! uvp-node
            end do
            end do
            jz_min = 2
        else
            jz_min = 1
        end if

        !--this is only required b/c of the unnatural placement of all strains
        !  onto w-nodes, be careful not to overwrite nz on top process with garbage
        !--dwdz(jz=0) is already known, except at bottom process (OK)
        !call mpi_sendrecv (dwdz(1, 1, nz-1), ldx*nynpy, MPI_RPREC, up, 1,  &
        !                   dwdz(1, 1, 0), ldx*nynpy, MPI_RPREC, down, 1,   &
        !                   comm, status, ierr)
        call mpi_sendrecv(dwdz(1, 1, 1),  ldx*nynpy, MPI_RPREC, down, 2,  &
                          dwdz(1, 1, nz), ldx*nynpy, MPI_RPREC, up,   2,  &
                          comm, status, ierr)

        ! calculate derivatives/strain on w-nodes
        !--in MPI version, calculating up to nz saves some interprocess exchange
        !  later but note that dwdz is not provided w/o some communication
        !  (unless its the top process)
        do jz=jz_min, nz
        do jy=1,nynpy
        do jx=1,nxt
            ux=0.5_rprec*(dudx(jx,jy,jz) + dudx(jx,jy,jz-1))  ! w-node
            uy=0.5_rprec*(dudy(jx,jy,jz) + dudy(jx,jy,jz-1))  ! w-node
            uz=dudz(jx,jy,jz)  ! w-node
            vx=0.5_rprec*(dvdx(jx,jy,jz) + dvdx(jx,jy,jz-1))  ! w-node
            vy=0.5_rprec*(dvdy(jx,jy,jz) + dvdy(jx,jy,jz-1))  ! w-node
            vz=dvdz(jx,jy,jz)  ! w-node
            wx=dwdx(jx,jy,jz)  ! w-node
            wy=dwdy(jx,jy,jz)  ! w-node
            wz=0.5_rprec*(dwdz(jx,jy,jz) + dwdz(jx,jy,jz-1))  ! w-node
            S11(jx,jy,jz)=ux          ! w-node
            S12(jx,jy,jz)=0.5_rprec*(uy+vx) ! w-node
            S13(jx,jy,jz)=0.5_rprec*(uz+wx) ! w-node
            S22(jx,jy,jz)=vy          ! w-node
            S23(jx,jy,jz)=0.5_rprec*(vz+wy) ! w-node
            S33(jx,jy,jz)=wz          ! w-node
        end do
        end do
        end do

        do k=1,nz
        do j=1,nynpy
        do i=1,nxt
            if ( coordz == 0 .and. (k == 1) ) then
            !--dwdy(jz=1) shouldx be 0, so we couldx use this
            !--dwdx(jz=1) shouldx be 0, so we couldx use this
                acc_mag(i,j,1) = (0.5_rprec*(dwdy(i,j,1)+dwdy(i,j,2))-dvdz(i,j,1))**2 &
                +(dudz(i,j,1)-0.5_rprec*(dwdx(i,j,1)+dwdx(i,j,2)))**2 &
                +(dvdx(i,j,1)-dudy(i,j,1))**2
            else
                acc_mag(i,j,k) = (dwdy(i,j,k)-dvdz(i,j,k))**2+(dudz(i,j,k)-dwdx(i,j,k))**2 &
                    +(dvdx(i,j,k)-dudy(i,j,k))**2
            end if
        end do
        end do
        end do

        do k=0,nz
        do j=1,nynpy
        do i=1,nxt
            acc_mag(i,j,k) = (S11(i,j,k)**2 + S22(i,j,k)**2 + S33(i,j,k)**2 &
                + 2._rprec*(S12(i,j,k)**2 + S13(i,j,k)**2 + S23(i,j,k)**2) &
                - acc_mag(i,j,k))/4._rprec
        end do
        end do
        end do
        acc_mag(0,1:nynpy,0:nz) = acc_mag(nxt,1:nynpy,0:nz)
        acc_mag(nxt+1,1:nynpy,0:nz) = acc_mag(1,1:nynpy,0:nz)
        call data_exchange_y(acc_mag(0:nxt+1,0:nynpy+1,0:nz-1), .true.)

        !DY Calculate x and y directions first
        
        do k=1,nz-1
        do j=1,nynpy
        do i=1,nxt
            RHS(i,j,k,ipcon)=RHS(i,j,k,ipcon)+  &
            ((1._rprec/dx**2)*((acc_mag(i+1,j,k)-acc_mag(i,j,k))*scalar_x(i+1,j,k) &
                -(acc_mag(i,j,k)-acc_mag(i-1,j,k))*scalar_x(i,j,k)) &
                +(1._rprec/dy**2)*((acc_mag(i,j+1,k)-acc_mag(i,j,k))*scalar_y(i,j+1,k) &
                -(acc_mag(i,j,k)-acc_mag(i,j-1,k))*scalar_y(i,j,k))) &
                *settling_vel(ipcon)/(9.81D0/(u_scale**2/z_i))*Cs_pa*delta**2
        end do
        end do
        end do
        !DY Add z direction
        do j=1,nynpy
        do i=1,nxt
            do k=2,nz-1
                RHS(i,j,k,ipcon)=RHS(i,j,k,ipcon)+  &
                (1._rprec/dz**2)*((acc_mag(i,j,k+1)-acc_mag(i,j,k))*scalar_z(i,j,k+1) &
                -(acc_mag(i,j,k)-acc_mag(i,j,k-1))*scalar_z(i,j,k)) &
                *settling_vel(ipcon)/(9.81D0/(u_scale**2/z_i))*Cs_pa*delta**2
            end do
            
            if ( coordz == 0 ) then
                RHS(i,j,1,ipcon)=RHS(i,j,1,ipcon)+  &
                (1._rprec/dz**2)*((acc_mag(i,j,2)-acc_mag(i,j,1))*scalar_z(i,j,2) &
                -(acc_mag(i,j,1)-0.5D0*(3._rprec*acc_mag(i,j,1)-acc_mag(i,j,2)))*scalar_z(i,j,1)) &
                *settling_vel(ipcon)/(9.81D0/(u_scale**2/z_i))*Cs_pa*delta**2
            else
                RHS(i,j,1,ipcon)=RHS(i,j,1,ipcon)+  &
                (1._rprec/dz**2)*((acc_mag(i,j,2)-acc_mag(i,j,1))*scalar_z(i,j,2) &
                -(acc_mag(i,j,1)-acc_mag(i,j,0))*scalar_z(i,j,1)) &
                *settling_vel(ipcon)/(9.81D0/(u_scale**2/z_i))*Cs_pa*delta**2
            end if
        end do
        end do

    case(3)
        !DY model 3: based on pressure Hessian
        !DY To save memory, save pressure Hessian into S_ij

        !DY Calculate x direction
        do k=1,nz
        do j=1,nynpy
        do i=1,nxt
            if(i.eq.1) then
                if(k.eq.nz) then
                    S11(1,j,k) = (1._rprec/dx)*   &
                    (0.5D0*(3._rprec*dpdx(1,j,k-1)-dpdx(1,j,k-2)) &
                    -0.5D0*(3._rprec*dpdx(nxt,j,k-1)-dpdx(nxt,j,k-2)))
                    S12(1,j,k) = (1._rprec/dx)*   &
                    (0.5D0*(3._rprec*dpdy(1,j,k-1)-dpdy(1,j,k-2)) &
                    -0.5D0*(3._rprec*dpdy(nxt,j,k-1)-dpdy(nxt,j,k-2)))
                else
                    S11(1,j,k) = (1._rprec/dx)*(dpdx(1,j,k)-dpdx(nxt,j,k))
                    S12(1,j,k) = (1._rprec/dx)*(dpdy(1,j,k)-dpdy(nxt,j,k))
                end if
            else
                if(k.eq.nz) then
                    S11(i,j,k) = (1._rprec/dx)*   &
                    (0.5D0*(3._rprec*dpdx(i,j,k-1)-dpdx(i,j,k-2)) &
                    -0.5D0*(3._rprec*dpdx(i-1,j,k-1)-dpdx(i-1,j,k-2)))
                    S12(i,j,k) = (1._rprec/dx)*   &
                    (0.5D0*(3._rprec*dpdy(i,j,k-1)-dpdy(i,j,k-2)) &
                    -0.5D0*(3._rprec*dpdy(i-1,j,k-1)-dpdy(i-1,j,k-2)))
                else
                    S11(i,j,k) = (1._rprec/dx)*(dpdx(i,j,k)-dpdx(i-1,j,k))
                    S12(i,j,k) = (1._rprec/dx)*(dpdy(i,j,k)-dpdy(i-1,j,k))
                end if
            end if
            
            if(k.eq.nz) then
                if(i.eq.nxt) then
                    S13(i,j,k) = 1._rprec/dx*(dpdz(1,j,k)-dpdz(i,j,k))
                else
                    S13(i,j,k) = 1._rprec/dx*(dpdz(i+1,j,k)-dpdz(i,j,k))
                end if
            else
                if(i.eq.nxt) then
                    S13(i,j,k) = 0.5_rprec/dx*((dpdz(1,j,k+1)+dpdz(1,j,k)) &
                    -(dpdz(i,j,k+1)+dpdz(i,j,k)))
                else
                    S13(i,j,k) = 0.5_rprec/dx*((dpdz(i+1,j,k+1)+dpdz(i+1,j,k)) &
                    -(dpdz(i,j,k+1)+dpdz(i,j,k)))
                end if
            end if
        end do
        end do
        end do
        
        S11(0,1:nynpy,1:nz) = S11(nxt,1:nynpy,1:nz)
        S11(nxt+1,1:nynpy,1:nz) = S11(1,1:nynpy,1:nz)
        call data_exchange_y(S11(0:nxt+1,0:nynpy+1,1:nz), .true.)

        S12(0,1:nynpy,1:nz) = S12(nxt,1:nynpy,1:nz)
        S12(nxt+1,1:nynpy,1:nz) = S12(1,1:nynpy,1:nz)
        call data_exchange_y(S12(0:nxt+1,0:nynpy+1,1:nz), .true.)

        S13(0,1:nynpy,1:nz) = S13(nxt,1:nynpy,1:nz)
        S13(nxt+1,1:nynpy,1:nz) = S13(1,1:nynpy,1:nz)
        call data_exchange_y(S13(0:nxt+1,0:nynpy+1,1:nz), .true.)

        !DY Calculate x direction
        !DY Save x component into acc_mag to save memory
        do k=1,nz-1
        do j=1,nynpy
        do i=1,nxt
            if(i.eq.1) then
                acc_mag(1,j,k) = (1._rprec/dx)* &
                ((1._rprec/dx)*(scalar(2,j,k,ipcon)-scalar(1,j,k,ipcon))*S11(2,j,k) &
                -(1._rprec/dx)*(scalar_x(2,j,k)-scalar_x(1,j,k))*S11(1,j,k))
            else if(i.eq.nxt) then
                acc_mag(nxt,j,k) = (1._rprec/dx) &
                *((1._rprec/dx)*(scalar_x(nxt+1,j,k)-scalar_x(nxt,j,k))*S11(nxt,j,k) &
                -(1._rprec/dx)*(scalar(nxt,j,k,ipcon)-scalar(nxt-1,j,k,ipcon))*S11(nxt-1,j,k))
            else
                acc_mag(i,j,k) = (1._rprec/dx)* &
                ((1._rprec/dx)*(scalar(i+1,j,k,ipcon)-scalar(i,j,k,ipcon))*S11(i+1,j,k) &
                -(1._rprec/dx)*(scalar(i,j,k,ipcon)-scalar(i-1,j,k,ipcon))*S11(i,j,k))
            end if
            
            if(j.eq.1) then
                acc_mag(i,j,k) = acc_mag(i,j,k)+(1._rprec/dx/dy) &
                *((scalar_x(i+1,j+1,k)-scalar_x(i+1,j,k))*S12(i+1,j,k) &
                -(scalar_x(i,j+1,k)-scalar_x(i,j,k))*S12(i,j,k))
            else if(j.eq.nynpy) then
                acc_mag(i,j,k) = acc_mag(i,j,k)+(1._rprec/dx/dy) &
                *((scalar_x(i+1,j,k)-scalar_x(i+1,j-1,k))*S12(i+1,j,k) &
                -(scalar_x(i,j,k)-scalar_x(i,j-1,k))*S12(i,j,k))
            else
                acc_mag(i,j,k) = acc_mag(i,j,k)+(0.5_rprec/dx/dy) &
                *((scalar_x(i+1,j+1,k)-scalar_x(i+1,j-1,k))*S12(i+1,j,k) &
                -(scalar_x(i,j+1,k)-scalar_x(i,j-1,k))*S12(i,j,k))
            end if
            
            if(k.eq.1) then
                acc_mag(i,j,k) = acc_mag(i,j,k)+(1._rprec/dx/dz) &
                *((scalar_x(i+1,j,k+1)-scalar_x(i+1,j,k))*S13(i+1,j,k) &
                -(scalar_x(i,j,k+1)-scalar_x(i,j,k))*S13(i,j,k))
            else if(k.eq.nz) then
                acc_mag(i,j,k) = acc_mag(i,j,k)+(1._rprec/dx/dz) &
                *((scalar_x(i+1,j,k)-scalar_x(i+1,j,k-1))*S13(i+1,j,k) &
                -(scalar_x(i,j,k)-scalar_x(i,j,k-1))*S13(i,j,k))
            else
                acc_mag(i,j,k) = acc_mag(i,j,k)+(0.5_rprec/dx/dz) &
                *((scalar_x(i+1,j,k+1)-scalar_x(i+1,j,k-1))*S13(i+1,j,k) &
                -(scalar_x(i,j,k+1)-scalar_x(i,j,k-1))*S13(i,j,k))
            end if
            RHS(i,j,k,ipcon) = RHS(i,j,k,ipcon)-  &
            acc_mag(i,j,k)*settling_vel(ipcon)/(9.81D0/(u_scale**2/z_i))*Cs_pa*delta**2
        end do
        end do
        end do

        !DY Calculate y direction
        do k=1,nz
        do j=1,nynpy
        do i=1,nxt
            if(j.eq.1) then
                if(k.eq.nz) then
                    S22(i,1,k) = (1._rprec/dy)*   &
                    (0.5D0*(3._rprec*dpdy(i,1,k-1)-dpdy(i,1,k-2)) &
                    -0.5D0*(3._rprec*dpdy(i,nynpy,k-1)-dpdy(i,nynpy,k-2)))
                    S12(i,1,k) = (1._rprec/dy)*   &
                    (0.5D0*(3._rprec*dpdx(i,1,k-1)-dpdx(i,1,k-2)) &
                    -0.5D0*(3._rprec*dpdx(i,nynpy,k-1)-dpdx(i,nynpy,k-2)))
                else
                    S22(i,1,k) = (1._rprec/dy)*(dpdy(i,1,k)-dpdy(i,nynpy,k))
                    S12(i,1,k) = (1._rprec/dy)*(dpdx(i,1,k)-dpdx(i,nynpy,k))
                end if
            else
                if(k.eq.nz) then
                    S22(i,j,k) = (1._rprec/dy)* &
                    (0.5D0*(3._rprec*dpdy(i,j,k-1)-dpdy(i,j,k-2)) &
                    -0.5D0*(3._rprec*dpdy(i,j-1,k-1)-dpdy(i,j-1,k-2)))
                    S12(i,j,k) = (1._rprec/dy)* &
                    (0.5D0*(3._rprec*dpdx(i,j,k-1)-dpdx(i,j,k-2)) &
                    -0.5D0*(3._rprec*dpdx(i,j-1,k-1)-dpdx(i,j-1,k-2)))
                else
                    S22(i,j,k) = (1._rprec/dy)*(dpdy(i,j,k)-dpdy(i,j-1,k))
                    S12(i,j,k) = (1._rprec/dy)*(dpdx(i,j,k)-dpdx(i,j-1,k))
                end if
            end if
            
            if(k.eq.nz) then
                if(j.eq.nynpy) then
                    S23(i,j,k) = 1._rprec/dy*(dpdz(i,1,k)-dpdz(i,j,k))
                else
                    S23(i,j,k) = 1._rprec/dy*(dpdz(i,j+1,k)-dpdz(i,j,k))
                end if
            else
                if(j.eq.nynpy) then
                    S23(i,j,k) = 0.5_rprec/dy*((dpdz(i,1,k+1)+dpdz(i,1,k)) &
                    -(dpdz(i,j,k+1)+dpdz(i,j,k)))
                else
                    S23(i,j,k) = 0.5_rprec/dy*((dpdz(i,j+1,k+1)+dpdz(i,j+1,k)) &
                    -(dpdz(i,j,k+1)+dpdz(i,j,k)))
                end if
            end if
        end do
        end do
        end do
        
        S12(0,1:nynpy,1:nz)       = S12(nxt,1:nynpy,1:nz)
        S12(nxt+1,1:nynpy,1:nz)   = S12(1,1:nynpy,1:nz)
        call data_exchange_y(S12(0:nxt+1,0:nynpy+1,1:nz), .true.)

        S22(0,1:nynpy,1:nz)       = S22(nxt,1:nynpy,1:nz)
        S22(nxt+1,1:nynpy,1:nz)   = S22(1,1:nynpy,1:nz)
        call data_exchange_y(S22(0:nxt+1,0:nynpy+1,1:nz), .true.)

        S23(0,1:nynpy,1:nz)       = S23(nxt,1:nynpy,1:nz)
        S23(nxt+1,1:nynpy,1:nz)   = S23(1,1:nynpy,1:nz)
        call data_exchange_y(S23(0:nxt+1,0:nynpy+1,1:nz), .true.)

        !DY Calculate y direction
        !DY Save y component into acc_mag to save memory
        do k=1,nz-1
        do j=1,nynpy
        do i=1,nxt
            if(j.eq.1) then
                acc_mag(i,1,k) = (1._rprec/dy)*     &
                ((1._rprec/dy)*(scalar(i,2,k,ipcon)-scalar(i,1,k,ipcon))*S22(i,2,k) &
                -(1._rprec/dy)*(scalar_y(i,2,k)-scalar_y(i,1,k))*S22(i,1,k))
            else if(j.eq.nynpy) then
                acc_mag(i,nynpy,k) = (1._rprec/dy) &
                *((1._rprec/dy)*(scalar_y(i,nynpy+1,k)-scalar_y(i,nynpy,k))*S22(i,nynpy,k) &
                -(1._rprec/dy)*(scalar(i,nynpy,k,ipcon)-scalar(i,nynpy-1,k,ipcon))*S22(i,nynpy-1,k))
            else
                acc_mag(i,j,k) = (1._rprec/dy)*     &
                ((1._rprec/dy)*(scalar(i,j+1,k,ipcon)-scalar(i,j,k,ipcon))*S22(i,j+1,k) &
                -(1._rprec/dy)*(scalar(i,j,k,ipcon)-scalar(i-1,j,k,ipcon))*S22(i,j,k))
            end if
            
            if(i.eq.1) then
                acc_mag(i,j,k) = acc_mag(i,j,k)+(1._rprec/dx/dy) &
                *((scalar_y(i+1,j+1,k)-scalar_y(i,j+1,k))*S12(i,j+1,k) &
                -(scalar_y(i+1,j,k)-scalar_y(i,j,k))*S12(i,j,k))
            else if(i.eq.nxt) then
                acc_mag(i,j,k) = acc_mag(i,j,k)+(1._rprec/dx/dy) &
                *((scalar_y(i,j+1,k)-scalar_y(i-1,j+1,k))*S12(i,j+1,k) &
                -(scalar_y(i,j,k)-scalar_y(i-1,j,k))*S12(i,j,k))
            else
                acc_mag(i,j,k) = acc_mag(i,j,k)+(0.5_rprec/dx/dy) &
                *((scalar_y(i+1,j+1,k)-scalar_y(i-1,j+1,k))*S12(i,j+1,k) &
                -(scalar_y(i+1,j,k)-scalar_y(i-1,j,k))*S12(i,j,k))
            end if
            
            if(k.eq.1) then
                acc_mag(i,j,k) = acc_mag(i,j,k)+(1._rprec/dy/dz) &
                *((scalar_y(i,j+1,k+1)-scalar_y(i,j+1,k))*S23(i,j+1,k) &
                -(scalar_y(i,j,k+1)-scalar_y(i,j,k))*S23(i,j,k))
            else if(k.eq.nz) then
                acc_mag(i,j,k) = acc_mag(i,j,k)+(1._rprec/dy/dz) &
                *((scalar_y(i,j+1,k)-scalar_y(i,j+1,k-1))*S23(i,j+1,k) &
                -(scalar_y(i,j,k)-scalar_y(i,j,k-1))*S23(i,j,k))
            else
                acc_mag(i,j,k) = acc_mag(i,j,k)+(0.5_rprec/dy/dz) &
                *((scalar_y(i,j+1,k+1)-scalar_y(i,j+1,k-1))*S23(i,j+1,k) &
                -(scalar_y(i,j,k+1)-scalar_y(i,j,k-1))*S23(i,j,k))
            end if
            RHS(i,j,k,ipcon) = RHS(i,j,k,ipcon)-acc_mag(i,j,k)* &
                settling_vel(ipcon)/(9.81D0/(u_scale**2/z_i))*Cs_pa*delta**2
        end do
        end do
        end do

        !DY Calculate z direction
        do k=1,nz
        do j=1,nynpy
        do i=1,nxt
            if(k.eq.1) then
                S13(i,j,k) = (1._rprec/dz)*(dpdx(i,j,2)-dpdx(i,j,1))
                S23(i,j,k) = (1._rprec/dz)*(dpdy(i,j,2)-dpdy(i,j,1))
                S33(i,j,k) = (1._rprec/dz)*(dpdz(i,j,2)-dpdz(i,j,1))
            else if(k.eq.nz) then
                S13(i,j,k) = (1._rprec/dz)*(dpdx(i,j,nz-1)-dpdx(i,j,nz-2))
                S23(i,j,k) = (1._rprec/dz)*(dpdy(i,j,nz-1)-dpdy(i,j,nz-2))
                S33(i,j,k) = (1._rprec/dz)*(dpdz(i,j,nz-1)-dpdz(i,j,nz-2))
            else
                S13(i,j,k) = (1._rprec/dz)*(dpdx(i,j,k)-dpdx(i,j,k-1))
                S23(i,j,k) = (1._rprec/dz)*(dpdy(i,j,k)-dpdy(i,j,k-1))
                S33(i,j,k) = (1._rprec/dz)*(dpdz(i,j,k)-dpdz(i,j,k-1))
            end if
        end do
        end do
        end do
        
        S13(0,1:nynpy,1:nz)       = S13(nxt,1:nynpy,1:nz)
        S13(nxt+1,1:nynpy,1:nz)   = S13(1,1:nynpy,1:nz)
        call data_exchange_y(S13(0:nxt+1,0:nynpy+1,1:nz), .true.)

        S23(0,1:nynpy,1:nz)       = S23(nxt,1:nynpy,1:nz)
        S23(nxt+1,1:nynpy,1:nz)   = S23(1,1:nynpy,1:nz)
        call data_exchange_y(S23(0:nxt+1,0:nynpy+1,1:nz), .true.)

        S33(0,1:nynpy,1:nz)       = S33(nxt,1:nynpy,1:nz)
        S33(nxt+1,1:nynpy,1:nz)   = S33(1,1:nynpy,1:nz)
        call data_exchange_y(S33(0:nxt+1,0:nynpy+1,1:nz), .true.)

        !DY Calculate z direction
        !DY Save z component into acc_mag to save memory
        do k = 1, nz-1
        do j = 1, nynpy
        do i = 1, nxt
            if(i.eq.1) then
                acc_mag(1,j,k) = (1._rprec/dz) &
                *((1._rprec/dx)*(scalar_z(i+1,j,k+1)-scalar_z(i,j,k+1))*S13(i,j,k+1) &
                - (1._rprec/dx)*(scalar_z(i+1,j,k)  -scalar_z(i,j,k))  *S13(i,j,k))
            else if(i.eq.nxt) then
                acc_mag(nxt,j,k) = (1._rprec/dz) &
                *((1._rprec/dx)*(scalar_z(i,j,k+1)-scalar_z(i-1,j,k+1))*S13(i,j,k+1) &
                -(1._rprec/dx)*(scalar_z(i,j,k)-scalar_z(i-1,j,k))*S13(i,j,k))
            else
                acc_mag(i,j,k) = (1._rprec/dz) &
                *((0.5_rprec/dx)*(scalar_z(i+1,j,k+1)-scalar_z(i-1,j,k+1))*S13(i,j,k+1) &
                -(0.5_rprec/dx)*(scalar_z(i+1,j,k)-scalar_z(i-1,j,k))*S13(i,j,k))
            end if
            
            if(j.eq.1) then
                acc_mag(i,j,k) = acc_mag(i,j,k)+(1._rprec/dy/dz) &
                *((scalar_z(i,j+1,k+1)-scalar_z(i,j,k+1))*S23(i,j,k+1) &
                -(scalar_z(i,j+1,k)-scalar_z(i,j,k))*S23(i,j,k))
            else if(j.eq.nynpy) then
                acc_mag(i,j,k) = acc_mag(i,j,k)+(1._rprec/dy/dz) &
                *((scalar_z(i,j,k+1)-scalar_z(i,j-1,k+1))*S23(i,j,k+1) &
                -(scalar_z(i,j,k)-scalar_z(i,j-1,k))*S23(i,j,k))
            else
                acc_mag(i,j,k) = acc_mag(i,j,k)+(0.5_rprec/dy/dz) &
                *((scalar_z(i,j+1,k+1)-scalar_z(i,j-1,k+1))*S23(i,j,k+1) &
                -(scalar_z(i,j+1,k)-scalar_z(i,j-1,k))*S23(i,j,k))
            end if
            
            if(k.eq.1) then
                acc_mag(i,j,k) = acc_mag(i,j,k)+(1._rprec/dz) &
                *((1._rprec/dz)*(scalar(i,j,k+1,ipcon)-scalar(i,j,k,ipcon))*S33(i,j,k+1) &
                -(1._rprec/dz)*(scalar(i,j,k+1,ipcon)-scalar(i,j,k,ipcon))*S33(i,j,k))
            else if(k.eq.nz-1) then
                acc_mag(i,j,k) = acc_mag(i,j,k)+(1._rprec/dz) &
                *((1._rprec/dz)*(scalar(i,j,k,ipcon)-scalar(i,j,k-1,ipcon))*S33(i,j,k+1) &
                -(1._rprec/dz)*(scalar(i,j,k,ipcon)-scalar(i,j,k-1,ipcon))*S33(i,j,k))
            else
                acc_mag(i,j,k) = acc_mag(i,j,k)+(1._rprec/dz**2) &
                *((scalar(i,j,k+1,ipcon)-scalar(i,j,k,ipcon))*S33(i,j,k+1) &
                -(scalar(i,j,k,ipcon)-scalar(i,j,k-1,ipcon))*S33(i,j,k))
            end if
            RHS(i,j,k,ipcon) = RHS(i,j,k,ipcon)-    &
                acc_mag(i,j,k)*settling_vel(ipcon)/(9.81D0/(u_scale**2/z_i))*Cs_pa*delta**2
        end do
        end do
        end do

    case default
        write (*,*) 'scalar_module: invalid model_psgsacc number'
        stop
    end select
    ! ---
    end subroutine sgs_acceleration
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
    subroutine calc_source_term(ipcon, scalar, RHS)
!-----------------------------------------------------------------------
!   Add source and sink terms if needed
!-----------------------------------------------------------------------
    use canopy, only: Vmax, Km
    implicit none
    ! ---
    integer, intent(in) :: ipcon
    real(kind=rprec), dimension(ldx,0:nynpy+1,0:nz,npcon), intent(in) :: scalar  
    real(kind=rprec), dimension(ldx,0:nynpy+1,0:nz,npcon), intent(inout) :: RHS 
    ! ---
    integer :: ind
    ! ---
    if (jt >= ini_src .and. jt < end_src) then 
        select case (source_type)
        case ("point")
            do ind = 1, con_src(ipcon)%n
                if (coordy == con_src(ipcon)%icpu_y(ind) .and. coordz == con_src(ipcon)%icpu_z(ind)) then
                    RHS(con_src(ipcon)%ix(ind),con_src(ipcon)%iy_cpu(ind),&
                        con_src(ipcon)%iz_cpu(ind),ipcon)= &
                        RHS(con_src(ipcon)%ix(ind),con_src(ipcon)%iy_cpu(ind),&
                            con_src(ipcon)%iz_cpu(ind),ipcon) + &
                            con_src(ipcon)%rls(ind)
                end if
            end do
        case ("line")
            do ind = 1, con_src(ipcon)%n
                select case (src_align)
                case ('lateral')
                    if (coordz == con_src(ipcon)%icpu_z(ind)) then
                        RHS(con_src(ipcon)%ix(ind),1:nynpy,&
                            con_src(ipcon)%iz_cpu(ind),ipcon)= &
                            RHS(con_src(ipcon)%ix(ind),1:nynpy,&
                                con_src(ipcon)%iz_cpu(ind),ipcon) + &
                                con_src(ipcon)%rls(ind)
                    end if
                case ('streamwise')
                    if (coordy == con_src(ipcon)%icpu_y(ind) .and.  &
                        coordz == con_src(ipcon)%icpu_z(ind)) then
                        RHS(1:nxt,con_src(ipcon)%iy_cpu(ind),&
                            con_src(ipcon)%iz_cpu(ind),ipcon)= &
                            RHS(1:nxt,con_src(ipcon)%iy_cpu(ind),&
                                con_src(ipcon)%iz_cpu(ind),ipcon)+ &
                                con_src(ipcon)%rls(ind)
                    end if
                case default
                    write (*, *) 'invalid source alignment type'
                    stop
                end select
            end do
        case("planar")
            do ind = 1, con_src(ipcon)%n
                if (coordz == con_src(ipcon)%icpu_z(ind)) then
                    RHS(1:nxt,1:nynpy,con_src(ipcon)%iz_cpu(ind),ipcon) = &
                        RHS(1:nxt,1:nynpy,con_src(ipcon)%iz_cpu(ind),ipcon) + &
                        con_src(ipcon)%rls(ind)
                end if
            end do
        case default
            write(*,*) "invalid source type"
        end select
    end if 

    ! --- sink term
    if ( flag_canopy .and. flag_uptake ) then
        if (use_field_scale) then
            RHS(1:nxt,1:nynpy,1:nz-1,ipcon) = RHS(1:nxt,1:nynpy,1:nz-1,ipcon) -     &
                        z_i/u_scale*a_leafz_uv(1:nxt,1:nynpy,1:nz-1)/z_i*Vmax *     &
                max(scalar(1:nxt,1:nynpy,1:nz-1,ipcon),0._rprec)/(Km+max(scalar(1:nxt,1:nynpy,1:nz-1,ipcon),0._rprec))
        end if
    end if
    ! ---
    end subroutine calc_source_term
!-----------------------------------------------------------------------
!   
!----------------------------------------------------------------------- 
    subroutine QUICK_scheme(ipcon, scalar)
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
    implicit none
    ! ---
    integer, intent(in) :: ipcon   
    real(kind=rprec), dimension(ldx,0:nynpy+1,0:nz,npcon), intent(in) :: scalar
    ! ---
    integer :: i, j, k
    real(kind=rprec), dimension(ldx,0:nz) :: temp1, temp2, temp3, temp4
    ! ---
    ! Interpolated concentrations in x
    do k = 1, nz-1
    do j = 1, nynpy
        ! interior nodes
        do i = 3, nxt-1
            ! Choose the appropriate upwind direction
            if ((u_int(i,j,k)+ust(k)) >= 0._rprec) then
                scalar_x(i,j,k)=(3._rprec*scalar(i,j,k,ipcon)+  &
                    6._rprec*scalar(i-1,j,k,ipcon)-scalar(i-2,j,k,ipcon))/8._rprec
            else
                scalar_x(i,j,k)=(3._rprec*scalar(i-1,j,k,ipcon)+    &
                    6._rprec*scalar(i,j,k,ipcon)-scalar(i+1,j,k,ipcon))/8._rprec
            end if
        end do

        ! nodes using BC

        ! Second node (i=2)
        ! Choose the appropriate upwind direction
        if ((u_int(2,j,k)+ust(k)) >= 0._rprec) then
            scalar_x(2,j,k)=(3._rprec*scalar(2,j,k,ipcon)+  &
                    6._rprec*scalar(1,j,k,ipcon)-ghost_x0(j,k))/8._rprec
        else
            scalar_x(2,j,k)=(3._rprec*scalar(1,j,k,ipcon)+  &
                    6._rprec*scalar(2,j,k,ipcon)-scalar(3,j,k,ipcon))/8._rprec
        end if

        ! Before last node (i=nxt)
        ! Choose the appropriate upwind direction
        if ((u_int(nxt,j,k)+ust(k)) >= 0._rprec) then
            scalar_x(nxt,j,k)=(3._rprec*scalar(nxt,j,k,ipcon)+      &
                    6._rprec*scalar(nxt-1,j,k,ipcon)-scalar(nxt-2,j,k,ipcon))/8._rprec
        else
            scalar_x(nxt,j,k)=(3._rprec*scalar(nxt-1,j,k,ipcon)+    &
                    6._rprec*scalar(nxt,j,k,ipcon)-ghost_xLx(j,k))/8._rprec
        end if
    end do
    end do

    ! Interpolated concentrations in y
    !BC added by Bicheng Chen for direction of Stokes drift
    do k = 1, nz-1
    do i = 1, nxt
        ! interior nodes
        do j = 3, nynpy-1
        ! Choose the appropriate upwind direction
            if (v_int(i,j,k)+vst(k) >= 0._rprec) then
                scalar_y(i,j,k)=(3._rprec*scalar(i,j,k,ipcon)&
                    +6._rprec*scalar(i,j-1,k,ipcon)-scalar(i,j-2,k,ipcon))/8._rprec
            else
                scalar_y(i,j,k)=(3._rprec*scalar(i,j-1,k,ipcon)&
                    +6._rprec*scalar(i,j,k,ipcon)-scalar(i,j+1,k,ipcon))/8._rprec
            end if
        end do

        ! nodes using BC

        ! Second node (j=2)
        ! Choose the appropriate upwind direction
        if (v_int(i,2,k)+vst(k) >= 0._rprec) then
            scalar_y(i,2,k)=(3._rprec*scalar(i,2,k,ipcon)&
                +6._rprec*scalar(i,1,k,ipcon)-ghost_y0(i,k))/8._rprec
        else
            scalar_y(i,2,k)=(3._rprec*scalar(i,1,k,ipcon)&
                +6._rprec*scalar(i,2,k,ipcon)-scalar(i,3,k,ipcon))/8._rprec
        end if

        ! Before last node (j=nyt)
        ! Choose the appropriate upwind direction
        if (v_int(i,nynpy,k)+vst(k) >= 0._rprec) then
            scalar_y(i,nynpy,k)=(3._rprec*scalar(i,nynpy,k,ipcon)&
                +6._rprec*scalar(i,nynpy-1,k,ipcon)-scalar(i,nynpy-2,k,ipcon))/8._rprec
        else
            scalar_y(i,npy,k)=(3._rprec*scalar(i,npy-1,k,ipcon)&
                +6._rprec*scalar(i,nynpy,k,ipcon)-ghost_yLy(i,k))/8._rprec
        end if
    end do
    end do

    ! Interpolated concentrations in z
    ! Note: interpolated values at k=1 and k=nz are not needed, since the fluxes are
    ! imposed through BC and w'=0 anyway
    do k = 3, nz-2
    do j = 1, nynpy
    do i = 1, nxt
        ! Choose the appropriate upwind direction
        if (w_int(i,j,k) >= 0._rprec) then
            scalar_z(i,j,k)=(3._rprec*scalar(i,j,k,ipcon)+  &
                    6._rprec*scalar(i,j,k-1,ipcon)-scalar(i,j,k-2,ipcon))/8._rprec
        else
            scalar_z(i,j,k)=(3._rprec*scalar(i,j,k-1,ipcon)+    &
                    6._rprec*scalar(i,j,k,ipcon)-scalar(i,j,k+1,ipcon))/8._rprec
        end if
    end do
    end do
    end do

    ! For the nodes k=2 and k=nz-1
    do j = 1, nynpy
    do i = 1, nxt
        ! k=2
        ! Choose the appropriate upwind direction
        if (w_int(i,j,2) >= 0._rprec) then
            ! scalar_z(i,j,2)=(3._rprec*scalar(i,j,2)+6._rprec*scalar(i,j,1)-scalar(i,j,0))/8._rprec
            ! For now use a simple interpolation in this case
            if ( coordz == 0 ) then
                scalar_z(i,j,2)=(scalar(i,j,2,ipcon)+scalar(i,j,1,ipcon))/2._rprec
            else
                scalar_z(i,j,2)=(3._rprec*scalar(i,j,2,ipcon)+  &
                        6._rprec*scalar(i,j,1,ipcon)-scalar(i,j,0,ipcon))/8._rprec
            end if
        else
            scalar_z(i,j,2)=(3._rprec*scalar(i,j,1,ipcon)+  &
                    6._rprec*scalar(i,j,2,ipcon)-scalar(i,j,3,ipcon))/8._rprec
        end if

        ! k=nz-1
        ! Choose the appropriate upwind direction
        if (w_int(i,j,nz-1) >= 0._rprec) then
            scalar_z(i,j,nz-1)=(3._rprec*scalar(i,j,nz-1,ipcon)+    &
                    6._rprec*scalar(i,j,nz-2,ipcon)-scalar(i,j,nz-3,ipcon))/8._rprec
        else
            if ( coordz == npz-1 ) then
                ! scalar_z(i,j,nz-1)=(3._rprec*scalar(i,j,nz-2)+6._rprec*scalar(i,j,nz-1)-scalar(i,j,nz))/8._rprec
                scalar_z(i,j,nz-1)=(3._rprec*scalar(i,j,nz-2,ipcon)+    &
                        6._rprec*scalar(i,j,nz-1,ipcon)-0._rprec)/8._rprec
            else
                scalar_z(i,j,nz-1)=(3._rprec*scalar(i,j,nz-2,ipcon)+    &
                        6._rprec*scalar(i,j,nz-1,ipcon)-scalar(i,j,nz,ipcon))/8._rprec
            end if
        end if
    end do
    end do

    ! Compute face nodes for periodic boundary conditions
    if (periodicbcx) then
        do k = 1, nz-1
        do j = 1, nynpy
            ! First node (i=1)
            ! Choose the appropriate upwind direction
            if ((u_int(1,j,k)+ust(k)) >= 0._rprec) then
                scalar_x(1,j,k)=(3._rprec*scalar(1,j,k,ipcon)+  &
                    6._rprec*scalar(nxt,j,k,ipcon)-scalar(nxt-1,j,k,ipcon))/8._rprec
            else
                scalar_x(1,j,k)=(3._rprec*scalar(nxt,j,k,ipcon)+    &
                    6._rprec*scalar(1,j,k,ipcon)-scalar(2,j,k,ipcon))/8._rprec
            end if

            ! Last node (i=nxt+1)
            ! Choose the appropriate upwind direction
            if ((u_int(nxt+1,j,k)+ust(k)) >= 0._rprec) then
                scalar_x(nxt+1,j,k)=(3._rprec*scalar(1,j,k,ipcon)+  &
                    6._rprec*scalar(nxt,j,k,ipcon)-scalar(nxt-1,j,k,ipcon))/8._rprec
            else
                scalar_x(nxt+1,j,k)=(3._rprec*scalar(nxt,j,k,ipcon)+    &
                    6._rprec*scalar(1,j,k,ipcon)-scalar(2,j,k,ipcon))/8._rprec
            end if
        end do
        end do
    end if

    if (periodicbcy) then
        temp3(:,:) = scalar(:,nynpy,:,ipcon)
        temp4(:,:) = scalar(:,nynpy-1,:,ipcon)

        call mpi_sendrecv(temp3(:,:), ldx*(nz+1), mpi_rprec, north, 7,  &
                          temp1(:,:), ldx*(nz+1), mpi_rprec, south, 7,  &
                          comm, status, ierr)
        call mpi_sendrecv(temp4(:,:), ldx*(nz+1), mpi_rprec, north, 8,  &
                          temp2(:,:), ldx*(nz+1), mpi_rprec, south, 8,  &
                          comm, status, ierr)

        !-BC added by Bicheng Chen for direction of Stokes drift
        if ( coordy == 0 ) then
            do k = 1, nz-1
            do i = 1, nxt
                ! First node (j=1)
                ! Choose the appropriate upwind direction
                if (v_int(i,1,k)+vst(k) >= 0._rprec) then
                    scalar_y(i,1,k)=(3._rprec*scalar(i,1,k,ipcon)+  &
                        6._rprec*temp1(i,k) - temp2(i,k))/8._rprec
                else
                    scalar_y(i,1,k)=(3._rprec*temp1(i,k)+   &
                        6._rprec*scalar(i,1,k,ipcon)-scalar(i,2,k,ipcon))/8._rprec
                end if
            end do
            end do
        end if


        temp3(:,:) = scalar(:,1,:,ipcon)
        temp4(:,:) = scalar(:,2,:,ipcon)

        call mpi_sendrecv(temp3(:,:), ldx*(nz+1), mpi_rprec, north, 7,  &
                          temp1(:,:), ldx*(nz+1), mpi_rprec, south, 7,  &
                          comm, status, ierr)
        call mpi_sendrecv(temp4(:,:), ldx*(nz+1), mpi_rprec, north, 8,  &
                          temp2(:,:), ldx*(nz+1), mpi_rprec, south, 8,  &
                          comm, status, ierr)
        
        if ( coordy == npy-1 ) then 
            do k = 1, nz-1
            do i = 1, nxt
                ! Last node (j=nyt+1)
                ! Choose the appropriate upwind direction
                if (v_int(i,nynpy+1,k)+vst(k) >= 0._rprec) then
                    scalar_y(i,nynpy+1,k)=(3._rprec*temp1(i,k)+    &
                        6._rprec*scalar(i,nynpy,k,ipcon)-scalar(i,nynpy-1,k,ipcon))/8._rprec
                else     
                    scalar_y(i,nynpy+1,k)=(3._rprec*scalar(i,nynpy,k,ipcon)+  &
                        6._rprec*temp1(i,k)-temp2(i,k))/8._rprec
                end if
            end do
            end do
        end if    
    end if
    ! ---
    end subroutine QUICK_scheme
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
    subroutine SMART_scheme(ipcon, scalar)
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
    implicit none
    ! ---
    integer, intent(in) :: ipcon   
    real(kind=rprec), dimension(ldx,0:nynpy+1,0:nz,npcon), intent(in) :: scalar
    ! ---
    integer :: i, j, k
    real(kind=rprec) :: rx,ry,rz                     ! Gradient ratios in x, y and z
    real(kind=rprec), dimension(nxt+1,nz) :: temp1, temp3
    real(kind=rprec), dimension(ldx,0:nz) :: temp2, temp4
    ! --- Interpolated concentrations in x

    do k = 1, nz-1
    do j = 1, nynpy
        ! interior nodes
        do i = 3, nxt-1
            ! Choose the appropriate upwind direction
            if ((u_int(i,j,k)+ust(k)) >= 0._rprec) then
                ! To avoid problems of divisions by zero
                if ((scalar(i-1,j,k,ipcon)-scalar(i-2,j,k,ipcon))==0._rprec) then
                    scalar_x(i,j,k)=scalar(i-1,j,k,ipcon)
                else
                ! Gradient ratio
                    rx=(scalar(i,j,k,ipcon)-scalar(i-1,j,k,ipcon))&
                        /(scalar(i-1,j,k,ipcon)-scalar(i-2,j,k,ipcon))
                ! Bound the gradient
                    rx=MAX(0._rprec,MIN(MIN(2._rprec*rx,0.75_rprec*rx&
                            +0.25_rprec),4._rprec))
                ! Intrepolate the value
                    scalar_x(i,j,k)=scalar(i-1,j,k,ipcon)&
                        +rx*(scalar(i-1,j,k,ipcon)-scalar(i-2,j,k,ipcon))/2._rprec
                end if
            else
                ! To avoid problems of divisions by zero
                if ((scalar(i,j,k,ipcon)-scalar(i+1,j,k,ipcon))==0._rprec) then
                    scalar_x(i,j,k)=scalar(i,j,k,ipcon)
                else
                    ! Gradient ratio
                    rx=(scalar(i-1,j,k,ipcon)-scalar(i,j,k,ipcon))&
                        /(scalar(i,j,k,ipcon)-scalar(i+1,j,k,ipcon))
                    ! Bound the gradient
                    rx=MAX(0._rprec,MIN(MIN(2._rprec*rx,0.75_rprec*rx&
                        +0.25_rprec),4._rprec))
                    ! Intrepolate the value
                    scalar_x(i,j,k)=scalar(i,j,k,ipcon)&
                        +rx*(scalar(i,j,k,ipcon)-scalar(i+1,j,k,ipcon))/2._rprec
                end if
            end if
        end do

        ! nodes using BC

        ! Second node (i=2)
        ! Choose the appropriate upwind direction
        if ((u_int(2,j,k)+ust(k)) >= 0._rprec) then
            ! To avoid problems of divisions by zero
            if ((scalar(1,j,k,ipcon)-ghost_x0(j,k))==0._rprec) then
                scalar_x(2,j,k)=scalar(1,j,k,ipcon)
            else
                ! Gradient ratio
                rx=(scalar(2,j,k,ipcon)-scalar(1,j,k,ipcon))&
                    /(scalar(1,j,k,ipcon)-ghost_x0(j,k))
                ! Bound the gradient
                rx=MAX(0._rprec,MIN(MIN(2._rprec*rx,0.75_rprec*rx+0.25_rprec),4._rprec))
                ! Intrepolate the value
                scalar_x(2,j,k)=scalar(1,j,k,ipcon)&
                    +rx*(scalar(1,j,k,ipcon)-ghost_x0(j,k))/2._rprec
            end if
        else
            ! To avoid problems of divisions by zero
            if ((scalar(2,j,k,ipcon)-scalar(3,j,k,ipcon))==0._rprec) then
                scalar_x(2,j,k)=scalar(2,j,k,ipcon)
            else
            ! Gradient ratio
                rx=(scalar(1,j,k,ipcon)-scalar(2,j,k,ipcon))&
                    /(scalar(2,j,k,ipcon)-scalar(3,j,k,ipcon))
            ! Bound the gradient
                rx=MAX(0._rprec,MIN(MIN(2._rprec*rx,0.75_rprec*rx+0.25_rprec),4._rprec))
            ! Intrepolate the value
                scalar_x(2,j,k)=scalar(2,j,k,ipcon)&
                    +rx*(scalar(2,j,k,ipcon)-scalar(3,j,k,ipcon))/2._rprec
            end if
        end if

        ! Before last node (i=nxt)
        ! Choose the appropriate upwind direction
        if ((u_int(nxt,j,k)+ust(k)) >= 0._rprec) then
            ! To avoid problems of divisions by zero
            if ((scalar(nxt-1,j,k,ipcon)&
                -scalar(nxt-2,j,k,ipcon))==0._rprec) then
                scalar_x(nxt,j,k)=scalar(nxt-1,j,k,ipcon)
            else
                ! Gradient ratio
                rx=(scalar(nxt,j,k,ipcon)-scalar(nxt-1,j,k,ipcon))&
                    /(scalar(nxt-1,j,k,ipcon)-scalar(nxt-2,j,k,ipcon))
                ! Bound the gradient
                rx=MAX(0._rprec,MIN(MIN(2._rprec*rx,0.75_rprec*rx+0.25_rprec),4._rprec))
                ! Intrepolate the value
                scalar_x(nxt,j,k)=scalar(nxt-1,j,k,ipcon)&
                    +rx*(scalar(nxt-1,j,k,ipcon)-scalar(nxt-2,j,k,ipcon))/2._rprec
            end if
        else
            ! To avoid problems of divisions by zero
            if ((scalar(nxt,j,k,ipcon)-ghost_xLx(j,k))==0._rprec) then
                scalar_x(nxt,j,k)=scalar(nxt,j,k,ipcon)
            else
            ! Gradient ratio
                rx=(scalar(nxt-1,j,k,ipcon)-scalar(nxt,j,k,ipcon))&
                    /(scalar(nxt,j,k,ipcon)-ghost_xLx(j,k))
            ! Bound the gradient
                rx=MAX(0._rprec,MIN(MIN(2._rprec*rx,0.75_rprec*rx+0.25_rprec),4._rprec))
            ! Intrepolate the value
                scalar_x(nxt,j,k)=scalar(nxt,j,k,ipcon)&
                    +rx*(scalar(nxt,j,k,ipcon)-ghost_xLx(j,k))/2._rprec
            end if
        end if
    end do
    end do   
    
    ! --- Interpolated concentrations in y
    do k = 1, nz-1
    do i = 1, nxt
        ! interior nodes
        do j = 3, nynpy-1
        ! choose the appropriate upwind direction
            if (v_int(i,j,k)+vst(k) >= 0._rprec) then
            ! to avoid problems of divisions by zero
                if ((scalar(i,j-1,k,ipcon)-scalar(i,j-2,k,ipcon))==0._rprec) then
                    scalar_y(i,j,k)=scalar(i,j-1,k,ipcon)
                else
                ! gradient ratio
                    ry=(scalar(i,j,k,ipcon)-scalar(i,j-1,k,ipcon))/&
                        (scalar(i,j-1,k,ipcon)-scalar(i,j-2,k,ipcon))
                ! bound the gradient
                    ry=max(0._rprec,min(min(2._rprec*ry,0.75_rprec*ry&
                        +0.25_rprec),4._rprec))
                ! intrepolate the value
                    scalar_y(i,j,k)=scalar(i,j-1,k,ipcon)&
                        +ry*(scalar(i,j-1,k,ipcon)-scalar(i,j-2,k,ipcon))/2._rprec
                end if
            else
            ! to avoid problems of divisions by zero
                if ((scalar(i,j,k,ipcon)-scalar(i,j+1,k,ipcon))==0._rprec) then
                    scalar_y(i,j,k)=scalar(i,j,k,ipcon)
                else
                ! gradient ratio
                    ry=(scalar(i,j-1,k,ipcon)-scalar(i,j,k,ipcon))&
                        /(scalar(i,j,k,ipcon)-scalar(i,j+1,k,ipcon))
                ! bound the gradient
                    ry=max(0._rprec,min(min(2._rprec*ry,0.75_rprec*ry&
                        +0.25_rprec),4._rprec))
                ! intrepolate the value
                    scalar_y(i,j,k)=scalar(i,j,k,ipcon)&
                        +ry*(scalar(i,j,k,ipcon)-scalar(i,j+1,k,ipcon))/2._rprec
                end if
            end if
        end do
    end do
    end do
        
    ! nodes using BC
    do k = 1, nz-1
    do i = 1, nxt
        ! Second node (j=2)
        ! Choose the appropriate upwind direction
        if (v_int(i,2,k)+vst(k) >= 0._rprec) then
        ! to avoid problems of divisions by zero
            ! if ( coordy == 0 ) then
            if ( scalar(i,1,k,ipcon) - ghost_y0(i,k) == 0._rprec ) then
                scalar_y(i,2,k) = scalar(i,1,k,ipcon)
            else
                ! if ((scalar(i,1,k,ipcon)-scalar(i,0,k,ipcon))==0._rprec) then
                !     scalar_y(i,2,k) = scalar(i,1,k,ipcon)
                ! else
                ! gradient ratio
                    ry=(scalar(i,2,k,ipcon)-scalar(i,1,k,ipcon))&
                        /(scalar(i,1,k,ipcon)-ghost_y0(i,k))
                ! bound the gradient
                    ry=max(0._rprec,min(min(2._rprec*ry,0.75_rprec*ry&
                        +0.25_rprec),4._rprec))
                ! intrepolate the value
                    scalar_y(i,2,k)=scalar(i,1,k,ipcon)&
                        +ry*(scalar(i,1,k,ipcon)-ghost_y0(i,k))/2._rprec
                ! end if
            end if
        else
        ! to avoid problems of divisions by zero
            if ((scalar(i,2,k,ipcon)-scalar(i,3,k,ipcon))==0._rprec) then
                scalar_y(i,2,k) = scalar(i,2,k,ipcon)
            else
            ! gradient ratio
                ry=(scalar(i,1,k,ipcon)-scalar(i,2,k,ipcon))&
                    /(scalar(i,2,k,ipcon)-scalar(i,3,k,ipcon))
            ! bound the gradient
                ry=max(0._rprec,min(min(2._rprec*ry,0.75_rprec*ry&
                    +0.25_rprec),4._rprec))
            ! intrepolate the value
                scalar_y(i,2,k)=scalar(i,2,k,ipcon)&
                    +ry*(scalar(i,2,k,ipcon)-scalar(i,3,k,ipcon))/2._rprec
            end if
        end if

        ! Before last node (j=nyt)
        ! Choose the appropriate upwind direction
        if (v_int(i,nynpy,k)+vst(k) >= 0._rprec) then
            ! to avoid problems of divisions by zero
            if ((scalar(i,nynpy-1,k,ipcon)-scalar(i,nynpy-2,k,ipcon))==0._rprec) then
                scalar_y(i,nynpy,k)=scalar(i,nynpy-1,k,ipcon)
            else
            ! gradient ratio
                ry=(scalar(i,nynpy,k,ipcon)-scalar(i,nynpy-1,k,ipcon))&
                    /(scalar(i,nynpy-1,k,ipcon)-scalar(i,nynpy-2,k,ipcon))
            ! bound the gradient
                ry=max(0._rprec,min(min(2._rprec*ry,0.75_rprec*ry&
                    +0.25_rprec),4._rprec))
            ! interpolate the value
                scalar_y(i,nynpy,k)=scalar(i,nynpy-1,k,ipcon)&
                    +ry*(scalar(i,nynpy-1,k,ipcon)-scalar(i,nynpy-2,k,ipcon))/2._rprec
            end if
        else
            ! to avoid problems of divisions by zero
            if ((scalar(i,nynpy,k,ipcon)-ghost_yLy(i,k))==0._rprec) then
                scalar_y(i,nynpy,k) = scalar(i,nynpy,k,ipcon)
            else
            ! gradient ratio
                ry=(scalar(i,nynpy-1,k,ipcon)-scalar(i,nynpy,k,ipcon))&
                    /(scalar(i,nynpy,k,ipcon)-ghost_yLy(i,k))
            ! bound the gradient
                ry=max(0._rprec,min(min(2._rprec*ry,0.75_rprec*ry&
                    +0.25_rprec),4._rprec))
            ! intrepolate the value
                scalar_y(i,nynpy,k)=scalar(i,nynpy,k,ipcon)&
                    +ry*(scalar(i,nynpy,k,ipcon)-ghost_yLy(i,k))/2._rprec
            end if
        end if
    end do
    end do
    
    ! --- Interpolated concentrations in z
    ! Note: interpolated values at k=1 and k=nz are not needed, since the fluxes are
    ! imposed through BC and w'=0 anyway
    DO k = 3, nz-2
    DO j = 1, nynpy
    DO i = 1, nxt
        ! Choose the appropriate upwind direction
        IF (w_int(i,j,k) >= 0._rprec) THEN
        ! To avoid problems of divisions by zero
            IF ((scalar(i,j,k-1,ipcon)-scalar(i,j,k-2,ipcon))==0._rprec) THEN
                scalar_z(i,j,k)=scalar(i,j,k-1,ipcon)
            ELSE
            ! Gradient ratio
                rz=(scalar(i,j,k,ipcon)-scalar(i,j,k-1,ipcon))&
                    /(scalar(i,j,k-1,ipcon)-scalar(i,j,k-2,ipcon))
            ! Bound the gradient
                rz=MAX(0._rprec,MIN(MIN(2._rprec*rz,0.75_rprec*rz&
                    +0.25_rprec),4._rprec))
            ! Intrepolate the value
                scalar_z(i,j,k)=scalar(i,j,k-1,ipcon)&
                    +rz*(scalar(i,j,k-1,ipcon)-scalar(i,j,k-2,ipcon))/2._rprec
            END IF
        ELSE
        ! To avoid problems of divisions by zero
            IF ((scalar(i,j,k,ipcon)-scalar(i,j,k+1,ipcon))==0._rprec) THEN
                scalar_z(i,j,k)=scalar(i,j,k,ipcon)
            ELSE
            ! Gradient ratio
                rz=(scalar(i,j,k-1,ipcon)-scalar(i,j,k,ipcon))&
                    /(scalar(i,j,k,ipcon)-scalar(i,j,k+1,ipcon))
            ! Bound the gradient
                rz=MAX(0._rprec,MIN(MIN(2._rprec*rz,0.75_rprec*rz&
                    +0.25_rprec),4._rprec))
            ! Intrepolate the value
                scalar_z(i,j,k)=scalar(i,j,k,ipcon)&
                    +rz*(scalar(i,j,k,ipcon)-scalar(i,j,k+1,ipcon))/2._rprec
            END IF
        END IF
    END DO
    END DO
    END DO

    ! For the nodes k=2 and k=nz-1
    DO j=1,nynpy
    DO i=1,nxt
    ! k=2
    ! Choose the appropriate upwind direction
        IF (w_int(i,j,2) >= 0._rprec) THEN
            if ( coordz == 0 ) then
            ! For now use a simple interpolation in this case
                scalar_z(i,j,2)=(scalar(i,j,2,ipcon)+scalar(i,j,1,ipcon))/2._rprec
            else
            ! To avoid problems of divisions by zero
                IF ((scalar(i,j,1,ipcon)-scalar(i,j,0,ipcon)==0._rprec)) THEN
                    scalar_z(i,j,2)=scalar(i,j,1,ipcon)
                ELSE
                ! Gradient ratio
                    rz=(scalar(i,j,2,ipcon)-scalar(i,j,1,ipcon))&
                      /(scalar(i,j,1,ipcon)-scalar(i,j,0,ipcon))
                ! Bound the gradient
                    rz=MAX(0._rprec,MIN(MIN(2._rprec*rz,0.75_rprec*rz&
                      +0.25_rprec),4._rprec))
                ! Interpolate the value
                    scalar_z(i,j,2)=scalar(i,j,1,ipcon)&
                      +rz*(scalar(i,j,1,ipcon)-scalar(i,j,0,ipcon))/2._rprec
                ENDIF
            endif
        ELSE
        ! To avoid problems of divisions by zero
            IF ((scalar(i,j,2,ipcon)-scalar(i,j,3,ipcon))==0._rprec) THEN
                scalar_z(i,j,2)=scalar(i,j,2,ipcon)
            ELSE
            ! Gradient ratio
                rz=(scalar(i,j,1,ipcon)-scalar(i,j,2,ipcon))&
                  /(scalar(i,j,2,ipcon)-scalar(i,j,3,ipcon))
            ! Bound the gradient
                rz=MAX(0._rprec,MIN(MIN(2._rprec*rz,0.75_rprec*rz&
                  +0.25_rprec),4._rprec))
            ! Intrepolate the value
                scalar_z(i,j,2)=scalar(i,j,2,ipcon)&
                    +rz*(scalar(i,j,2,ipcon)-scalar(i,j,3,ipcon))/2._rprec
            END IF
        END IF

        ! k=nz-1
        ! Choose the appropriate upwind direction
        IF (w_int(i,j,nz-1) >= 0._rprec) THEN
        ! To avoid problems of divisions by zero
            IF ((scalar(i,j,nz-2,ipcon)-scalar(i,j,nz-3,ipcon))==0._rprec) THEN
                scalar_z(i,j,nz-1)=scalar(i,j,nz-2,ipcon)
            ELSE
            ! Gradient ratio
                rz=(scalar(i,j,nz-1,ipcon)-scalar(i,j,nz-2,ipcon))&
                  /(scalar(i,j,nz-2,ipcon)-scalar(i,j,nz-3,ipcon))
            ! Bound the gradient
                rz=MAX(0._rprec,MIN(MIN(2._rprec*rz,0.75_rprec*rz&
                  +0.25_rprec),4._rprec))
            ! Interpolate the value
                scalar_z(i,j,nz-1)=scalar(i,j,nz-2,ipcon)&
                  +rz*(scalar(i,j,nz-2,ipcon)-scalar(i,j,nz-3,ipcon))/2._rprec
            END IF
        ELSE
        ! To avoid problems of divisions by zero
            IF ((scalar(i,j,nz-1,ipcon)-scalar(i,j,nz,ipcon))==0._rprec) THEN
                scalar_z(i,j,nz-1)=scalar(i,j,nz-1,ipcon)
            ELSE
            ! Gradient ratio
                rz=(scalar(i,j,nz-2,ipcon)-scalar(i,j,nz-1,ipcon))&
                  /(scalar(i,j,nz-1,ipcon)-scalar(i,j,nz,ipcon))
            ! Bound the gradient
                rz=MAX(0._rprec,MIN(MIN(2._rprec*rz,0.75_rprec*rz&
                    +0.25_rprec),4._rprec))
            ! Intrepolate the value
                scalar_z(i,j,nz-1)=scalar(i,j,nz-1,ipcon)&
                    +rz*(scalar(i,j,nz-1,ipcon)-scalar(i,j,nz,ipcon))/2._rprec
            END IF
        END IF
    END DO
    END DO

    ! ---
    call mpi_sendrecv(scalar(:,1:nynpy,nz-2,ipcon),ldx*nynpy,MPI_RPREC,up,&
                tag_counter+4,ghost_z0(1,1),ldx*nynpy,MPI_RPREC,down,&
                tag_counter+4,comm,status,ierr)
    if ( coordz == 0 ) then
        scalar_z(1:nxt,1:nynpy,1)=0._rprec
        !scalar_z(1:nxt,1:nyt,2)=(scalar(1:nxt,1:nyt,2)+scalar(1:nxt,1:nyt,1))/2._rprec
    else
        do j = 1, nynpy
        do i = 1, nxt
        ! k=1
        ! Choose the appropriate upwind direction
            if (w_int(i,j,1) >= 0._rprec) then
                if ((scalar(i,j,0,ipcon)-ghost_z0(i,j)==0._rprec)) then
                    scalar_z(i,j,1) = scalar(i,j,0,ipcon)
                else
                ! Gradient ratio
                    rz=(scalar(i,j,1,ipcon)-scalar(i,j,0,ipcon))/(scalar(i,j,0,ipcon)-ghost_z0(i,j))
                ! Bound the gradient
                    rz=MAX(0._rprec,MIN(MIN(2._rprec*rz,0.75_rprec*rz+0.25_rprec),4._rprec))
                ! Interpolate the value
                    scalar_z(i,j,1)=scalar(i,j,0,ipcon)+rz*(scalar(i,j,0,ipcon)-ghost_z0(i,j))/2._rprec
                end if
            else
                ! To avoid problems of divisions by zero
                if ((scalar(i,j,1,ipcon)-scalar(i,j,2,ipcon))==0._rprec) then
                    scalar_z(i,j,1)=scalar(i,j,1,ipcon)
                else
                ! Gradient ratio
                    rz=(scalar(i,j,0,ipcon)-scalar(i,j,1,ipcon))/(scalar(i,j,1,ipcon)-scalar(i,j,2,ipcon))
                ! Bound the gradient
                    rz=MAX(0._rprec,MIN(MIN(2._rprec*rz,0.75_rprec*rz+0.25_rprec),4._rprec))
                ! Intrepolate the value
                    scalar_z(i,j,1)=scalar(i,j,1,ipcon)+rz*(scalar(i,j,1,ipcon)-scalar(i,j,2,ipcon))/2._rprec
                end if
            end if
        end do
        end do
    end if
        
    call mpi_sendrecv(scalar_z(1,1,1), (nxt+1)*(nynpy+1),MPI_RPREC,down,tag_counter+1, &
                      scalar_z(1,1,nz),(nxt+1)*(nynpy+1),MPI_RPREC,up,  tag_counter+1, &
                      comm,status,ierr)

    if ( coordz == npz-1 ) then
        scalar_z(1:nxt,1:nynpy,nz) = scalar_z(1:nxt,1:nynpy,nz-1)
    end if
    
    ! Compute face nodes for periodic boundary conditions
    if (periodicbcx) then
        ! Interpolated concentrations in x
        DO k=1,nz-1
        DO j=1,nynpy
        ! First node (i=1)
        ! Choose the appropriate upwind direction
            IF ((u_int(1,j,k)+ust(k)) >= 0._rprec) THEN
            ! To avoid problems of divisions by zero
                IF ((scalar(nxt,j,k,ipcon)&
                    -scalar(nxt-1,j,k,ipcon))==0._rprec) THEN
                    scalar_x(1,j,k)=scalar(nxt,j,k,ipcon)
                ELSE
                ! Gradient ratio
                    rx=(scalar(1,j,k,ipcon)-scalar(nxt,j,k,ipcon))&
                      /(scalar(nxt,j,k,ipcon)-scalar(nxt-1,j,k,ipcon))
                ! Bound the gradient
                    rx=MAX(0._rprec,MIN(MIN(2._rprec*rx,0.75_rprec*rx&
                        +0.25_rprec),4._rprec))
                ! Intrepolate the value
                    scalar_x(1,j,k)=scalar(nxt,j,k,ipcon)&
                        +rx*(scalar(nxt,j,k,ipcon)-scalar(nxt-1,j,k,ipcon))/2._rprec
                END IF
            ELSE
            ! To avoid problems of divisions by zero
                IF ((scalar(1,j,k,ipcon)-scalar(2,j,k,ipcon))==0._rprec) THEN
                    scalar_x(1,j,k)=scalar(1,j,k,ipcon)
                ELSE
                ! Gradient ratio
                    rx=(scalar(nxt,j,k,ipcon)-scalar(1,j,k,ipcon))&
                      /(scalar(1,j,k,ipcon)-scalar(2,j,k,ipcon))
                ! Bound the gradient
                    rx=MAX(0._rprec,MIN(MIN(2._rprec*rx,0.75_rprec*rx&
                        +0.25_rprec),4._rprec))
                ! Intrepolate the value
                    scalar_x(1,j,k)=scalar(1,j,k,ipcon)&
                        +rx*(scalar(1,j,k,ipcon)-scalar(2,j,k,ipcon))/2._rprec
                END IF
            END IF
            ! Last node (i=nxt+1) - this is the same as the first node!
            scalar_x(nxt+1,j,k)=scalar_x(1,j,k)
        END DO
        END DO
    end if

    if (periodicbcy) then 
        temp4(:,:) = scalar(:,nynpy-1,:,ipcon)
        call mpi_sendrecv(temp4(:,:), ldx*(nz+1), mpi_rprec, north, 8,  &
                          temp2(:,:), ldx*(nz+1), mpi_rprec, south, 8,  &
                          comm, status, ierr)
        ! if ( coordy == 0 ) then
        !     scalar_y(1:nxt,1,1:nz) = 0._rprec
        ! else
            do k = 1, nz-1
            do i = 1, nxt
                ! First node (j=1)
                ! Choose the appropriate upwind direction
                if (v_int(i,1,k)+vst(k) >= 0._rprec) then
                ! To avoid problems of divisions by zero    
                    if ( scalar(i,0,k,ipcon)-temp2(i,k) == 0._rprec ) then
                        scalar_y(i,1,k) = scalar(i,0,k,ipcon)
                    else
                        ! Gradient ratio
                        ry=(scalar(i,1,k,ipcon)-scalar(i,0,k,ipcon))/(scalar(i,0,k,ipcon)-temp2(i,k))
                        ! Bound the gradient
                        ry=MAX(0._rprec,MIN(MIN(2._rprec*ry,0.75_rprec*ry+0.25_rprec),4._rprec))
                        ! Intrepolate the value
                        scalar_y(i,1,k)=scalar(i,0,k,ipcon)+ry*(scalar(i,0,k,ipcon)-temp2(i,k))/2._rprec
                    end if
                else
                    if ((scalar(i,1,k,ipcon)-scalar(i,2,k,ipcon))==0._rprec) then
                        scalar_y(i,1,k) = scalar(i,1,k,ipcon)
                    else
                    ! Gradient ratio
                        ry= (scalar(i,0,k,ipcon)-scalar(i,1,k,ipcon))/   &
                        (scalar(i,1,k,ipcon)-scalar(i,2,k,ipcon))
                    ! Bound the gradient
                        ry=MAX(0._rprec,MIN(MIN(2._rprec*ry,0.75_rprec*ry+0.25_rprec),4._rprec))
                    ! Intrepolate the value
                        scalar_y(i,1,k)=scalar(i,1,k,ipcon)+ry*(scalar(i,1,k,ipcon)-scalar(i,2,k,ipcon))/2._rprec 
                    end if
                end if
            end do
            end do
        ! end if
        
        temp3(:,:) = scalar_y(:,1,:)
        call mpi_sendrecv(temp3(:,:), (nxt+1)*nz, mpi_rprec, south, 2,  &
                          temp1(:,:), (nxt+1)*nz, mpi_rprec, north, 2,  &
                          comm, status, ierr)
        scalar_y(:,nynpy+1,:) = temp1(:,:)
    end if
        
    ! if (periodicbcy) then   

    !     temp3(:,:) = scalar(:,nynpy,:,ipcon)
    !     temp4(:,:) = scalar(:,nynpy-1,:,ipcon)

    !     call mpi_sendrecv(temp3(:,:), ldx*(nz+1), mpi_rprec, north, 7,  &
    !                       temp1(:,:), ldx*(nz+1), mpi_rprec, south, 7,  &
    !                       comm, status, ierr)
    !     call mpi_sendrecv(temp4(:,:), ldx*(nz+1), mpi_rprec, north, 8,  &
    !                       temp2(:,:), ldx*(nz+1), mpi_rprec, south, 8,  &
    !                       comm, status, ierr)
        
    !     if ( coordy == 0 ) then           
    !         do k = 1, nz-1
    !         do i = 1, nxt
    !             ! First node (j=1)
    !             ! Choose the appropriate upwind direction
    !             if (v_int(i,1,k)+vst(k) >= 0._rprec) then
    !             ! To avoid problems of divisions by zero
    !                 if ((temp1(i,k)-temp2(i,k))==0._rprec) then
    !                     scalar_y(i,1,k) = temp1(i,k)
    !                 else
    !                     ! Gradient ratio
    !                     ry=(scalar(i,1,k,ipcon)-temp1(i,k))/(temp1(i,k)-temp2(i,k))
    !                     ! Bound the gradient
    !                     ry=MAX(0._rprec,MIN(MIN(2._rprec*ry,0.75_rprec*ry+0.25_rprec),4._rprec))
    !                     ! Intrepolate the value
    !                     scalar_y(i,1,k)=temp1(i,k)+ry*(temp1(i,k)-temp2(i,k))/2._rprec
    !                 end if
    !             else
    !                 ! To avoid problems of divisions by zero
    !                 if ((scalar(i,1,k,ipcon)-scalar(i,2,k,ipcon))==0._rprec) then
    !                     scalar_y(i,1,k) = scalar(i,1,k,ipcon)
    !                 else
    !                 ! Gradient ratio
    !                     ry= (temp1(i,k)-scalar(i,1,k,ipcon))/   &
    !                         (scalar(i,1,k,ipcon)-scalar(i,2,k,ipcon))
    !                 ! Bound the gradient
    !                     ry=MAX(0._rprec,MIN(MIN(2._rprec*ry,0.75_rprec*ry+0.25_rprec),4._rprec))
    !                 ! Intrepolate the value
    !                     scalar_y(i,1,k)=scalar(i,1,k,ipcon)+ry*(scalar(i,1,k,ipcon)-scalar(i,2,k,ipcon))/2._rprec
    !                 end if
    !             end if
    !         end do
    !         end do
    !     end if

    !     call mpi_sendrecv(scalar_y(:,1,:),       (nxt+1)*nz, mpi_rprec, south, 2,  &
    !                       scalar_y(:,nynpy+1,:), (nxt+1)*nz, mpi_rprec, north, 2,  &
    !                       comm, status, ierr)
    ! end if
    ! ---
    end subroutine SMART_scheme
! ---
end module scalars_module
