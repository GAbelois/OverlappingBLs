!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!    
subroutine ic()
!--------------------------------------------------------------------!
! Log profile that is modified to flatten at z=z_i
!--------------------------------------------------------------------!
use types, only:rprec
use param
use sim_param, only:u, v, w
use bottombc
implicit none
! ---
real(kind=rprec), parameter :: alpha = 0.5_rprec
integer :: jx, jy, jz, seed, jz_abs

real(kind=rprec), dimension(nz) :: ubar, vbar
real(kind=rprec) :: rms, noise, arg, arg2
real(kind=rprec) :: z, w_star, ran3
! ---
! if ((inflow) .and. (.not. read_inflow_file)) then !--no turbulence
!     u = face_avg
!     v = 0._rprec
!     w = 0._rprec
! else
    !BC uncorelate w_star from wt_s by Bicheng Chen
    w_star = (9.81_rprec/T_init*0.06_rprec*z_i)**(1._rprec/3._rprec)
    !w_star=(9.81_rprec/T_init*wt_s*z_i)**(1._rprec/3._rprec)
    !      T_star=wt_s/w_star
    !      q_star=T_star
    !write(*,*) 'Modified Log Profile for IC'
    
    if ( ocean_flag ) then
        call ekman_layer(ubar, vbar)
        ! do jz = 1, nz-1
        !     if ( ubc_mom == 'wall' ) then               ! shallow ocean flow
        !         select case (walltype)
        !         case('smooth')
        !             arg2 = (lz_tot - zuvp(jz)) / (nu_molec/u_scale/z_i)
        !             arg = (1._rprec/vonk)*log(arg2) + B !-1./(2.*vonk*z_i*z_i)*z*z
        !         case('rough')
        !             arg2 = (lz_tot - zuvp(jz)) /(zo1/z_i)
        !             arg = (1._rprec/vonk)*log(arg2)
        !             !arg = 2._rprec*((lz_tot - zuvp(jz))/lz_tot)**0.3
        !         case default
        !             write(*,*) 'invalid wall type'
        !             stop
        !         end select
        !     else if ( ubc_mom == 'stress free' ) then   ! deep ocean flow
        !         ! waits to be implemented
        !         arg2 = zuvp(jz) /(zo1/z_i)
        !         arg = 1._rprec
        !         !arg = (1._rprec/vonk)*log(arg2)
        !     end if
        
        !     if ( coriolis_forcing ) then
        !         ubar(jz) = ug
        !         vbar(jz) = vg
        !     else
        !         ubar(jz) = arg
        !         vbar(jz) = 0._rprec
        !     end if
        
        !     !if ((coriolis_forcing) .and. (zuvp(jz) .gt. (.5_rprec*z_i))) ubar(jz) = ug 
        ! end do
    else   
        do jz = 1, nz-1 
            if ( ubc_mom == 'wall' ) then   ! channel flow with bottom and top walls
                select case (walltype)      ! bottom boundary
                case('smooth')
                    if ( zuvp(jz) .le. lz_tot/2 ) then
                        arg2 = zuvp(jz) / (nu_molec/u_scale/z_i)
                    else
                        arg2 = (lz_tot-zuvp(jz)) / (nu_molec/u_scale/z_i)
                    end if
                    arg = (1._rprec/vonk)*log(arg2) + B !-1./(2.*vonk*z_i*z_i)*z*z
                case('rough')
                    ! IC in equilibrium with rough surface (rough dominates in effective zo)
                    ! arg2=z/(sum(zo)/float(nxt*nyt))
                    if ( zuvp(jz) .le. lz_tot/2 ) then
                        arg2 = zuvp(jz)/(zo1/z_i)
                    else
                        arg2 = (lz_tot-zuvp(jz))/(zo1/z_i)
                    end if
                    arg = (1._rprec/vonk)*log(arg2)
                case default
                    write(*,*) 'invalid wall type'
                    stop
                end select
            else if ( ubc_mom == 'stress free' ) then   ! atmospheric boundary layer flow
                select case (walltype)      ! bottom boundary
                case('smooth')
                    arg2 = zuvp(jz) / (nu_molec/u_scale/z_i)
                    arg = (1._rprec/vonk)*log(arg2) + B !-1./(2.*vonk*z_i*z_i)*z*z
                case('rough')
                    ! IC in equilibrium with rough surface (rough dominates in effective zo)
                    ! arg2=z/(sum(zo)/float(nxt*nyt))
                    arg2 = zuvp(jz)/(zo1/z_i)
                    arg = (1._rprec/vonk)*log(arg2)
                case default
                    write(*,*) 'invalid wall type'
                    stop
                end select
            end if
                
            if ( coriolis_forcing ) then
                ubar(jz) = ug
                vbar(jz) = vg
                ! Note that ug and vg have already been non-dimensionalized in param.f90
                !           ubar(jz)=arg/30._rprec
            else
                ubar(jz) = arg
                vbar(jz) = 0._rprec
            end if

            if ((coriolis_forcing) .and. (zuvp(jz) .gt. (.5_rprec*z_i))) ubar(jz) = ug 
        end do
    end if
    
    rms = 1._rprec
    do jz = 1, nz-1
        jz_abs = coordz*(nz - 1) + jz
        seed = -80 - jz_abs !--trying to make consistent init for MPI
        do jy = 1, nynpy
        do jx = 1, nxt
            !...Ran3 returns uniform RV between 0 and 1. (sigma_rv=0.289)
            !...Taking std dev of vel as 1 at all heights
            if (zuvp(jz) .le. lz_tot) then
                noise = rms/.289_rprec*(ran3(seed) - .5_rprec)
                u(jx, jy, jz) = alpha*noise*abs(0.5_rprec*lz_tot - zuvp(jz)) + ubar(jz)
                noise = rms/.289_rprec*(ran3(seed) - 0.5_rprec)
                !v(jx, jy, jz) = alpha*noise*abs(0.5_rprec*lz_tot - zuvp(jz)) + vbar(jz)
                v(jx, jy, jz) = vbar(jz)
                noise = rms/.289_rprec*(ran3(seed) - .5_rprec)
                !w(jx, jy, jz) = alpha*noise*abs(0.5_rprec*lz_tot - zuvp(jz))
                w(jx, jy, jz) = 0._rprec
            else
                noise = rms/.289_rprec*(ran3(seed) - .5_rprec)
                u(jx, jy, jz) = noise*w_star/u_scale*.01_rprec + ubar(jz)
                noise = rms/.289_rprec*(ran3(seed) - 0.5_rprec)
                v(jx, jy, jz) = noise*w_star/u_scale*.01_rprec + vbar(jz)
                noise = rms/.289_rprec*(ran3(seed) - 0.5_rprec)
                w(jx, jy, jz) = noise*w_star/u_scale*.01_rprec
            end if
        end do
        end do
    end do

    !...BC for W
    if ( coordz == 0 ) then
        w(1:nxt, 1:nynpy, 1) = 0._rprec
    end if
    if ( coordz == npz - 1 ) then
        w(1:nxt, 1:nynpy, nz) = 0._rprec
    endif

    !...BC for U, V
    if ( coordz == npz - 1 ) then
        if ( ubc_mom == 'wall' ) then
            u(1:nxt, 1:nynpy, nz) = BOGUS
            v(1:nxt, 1:nynpy, nz) = BOGUS
        else if ( ubc_mom == 'stress free' ) then
            u(1:nxt, 1:nynpy, nz) = u(1:nxt, 1:nynpy, nz - 1)
            v(1:nxt, 1:nynpy, nz) = v(1:nxt, 1:nynpy, nz - 1)
        end if
    end if
! end if
! ---
end subroutine ic
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!    
subroutine ic_dns()
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!    
use types, only:rprec
use param
use sim_param, only:u, v, w
implicit none
! ---  
real(kind=rprec), dimension(nz)::ubar
real(kind=rprec) :: rms, temp
integer :: jx, jy, jz, seed
real(kind=rprec) :: z

real(rprec), external :: ran3
! ---
if (inflow) then
    ! uniform flow case:
    ubar = face_avg
else
    seed = -112

    ! calculate height of first uvp point in wall units
    ! lets do a laminar case (?)
    ! *** cyan we shouldn't impose laminar flow profile using u_scale as velocity scale
    ! ***      because u_scale in laminar and turbulent scenarios are completely different
    do jz = 1, nz
        ubar(jz) = (u_scale*z_i/nu_molec)*zuvp(jz)*(1._rprec - .5_rprec*zuvp(jz)) ! non-dimensional
    end do
end if

! rms=0.0001 seems to work in some cases
! the "default" rms of a unif variable is 0.289
rms = 0.2_rprec
do jz = 1, nz
do jy = 1, nynpy
do jx = 1, nxt
    u(jx, jy, jz) = ubar(jz) + (rms/.289_rprec)*(ran3(seed) - .5_rprec)/u_scale
    v(jx, jy, jz) = 0._rprec + (rms/.289_rprec)*(ran3(seed) - .5_rprec)/u_scale
    w(jx, jy, jz) = 0._rprec + (rms/.289_rprec)*(ran3(seed) - .5_rprec)/u_scale
end do
end do
end do

! make sure w-mean is 0
temp = 0._rprec
do jz = 1, nz
do jy = 1, nynpy
do jx = 1, nxt
    temp = temp + w(jx, jy, jz)
end do
end do
end do
temp = temp/(nxt*nynpy*nz)

do jz = 1, nz
do jy = 1, nynpy
do jx = 1, nxt
    w(jx, jy, jz) = w(jx, jy, jz) - temp
end do
end do
end do

w(:, :, 1)  = 0._rprec
w(:, :, nz) = 0._rprec
u(:, :, nz) = u(:, :, nz - 1)
v(:, :, nz) = v(:, :, nz - 1)
! ---
end subroutine ic_dns
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine ic_scal()
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
use types, only:rprec
use param
use sim_param, only:u, v, w, theta, q
use bottombc
implicit none
! ---

real(kind=rprec), dimension(nz) :: ubar, vbar, wbar
real(kind=rprec) :: rms, noise, arg, arg2, theta_mean, ustar_ratio
real(kind=rprec) :: z, w_star, T_star, q_star, ran3, bv_freq
real(kind=rprec) :: z_turb_limit, perturb_height_factor, z_inv
integer :: jx, jy, jz, seed, jz_abs

character(*), parameter :: fmt_7781 = "('z, ubar, vbar, wbar, Tbar:',5(1x, F9.4))"
! ---
if ( lbc .eq. 0 ) then
    theta_mean = T_s_min - 0.001*T_s_min !T_s_min is dimensional while T_s is non-dimensional
else 
    theta_mean = T_init
end if

if (wt_s .lt. 0._rprec) then
    perturb_height_factor = 0.30_rprec
    z_inv = 0.30_rprec*z_i
else
    perturb_height_factor = 0.3_rprec
    z_inv = 0.57_rprec*z_i

    if (ocean_flag) then
        z_inv = prop_mixed*z_i
    end if
end if
z_turb_limit = perturb_height_factor*z_i

if (wt_s .eq. 0.0_rprec) then
    ! w_star is of O(1) with z_i=500 and wt_s=0.06 for atmos
    ! For ocean, ~5e-5 is better
    w_star = (g/theta_mean*5.3691d-5*z_i)**(1._rprec/3._rprec)
    T_star = 5.3691d-5/w_star
    q_star = T_star
else
    w_star = sign((g/theta_mean*abs(wt_s)*z_i)**(1._rprec/3._rprec), wt_s)
    T_star = wt_s/w_star
    q_star = T_star
end if

if (ocean_flag) then
    if ( ubc_mom == 'wall' ) then
        call init_profile(ubar, vbar, z_inv)
    else
        ! call ekman_layer(ubar, vbar)
        ubar = ug
        vbar = vg
    end if
    wbar = 0._rprec
else
    do jz = 1, nz
        ustar_ratio = ug*vonk/log(z_inv/zo1)
        arg2 = zuvp(jz)/(zo1/z_i)
        arg = ustar_ratio*(1._rprec/vonk)*log(arg2) !-1./(2.*vonk*z_i*z_i)*z*z
        if (coriolis_forcing) then
            ubar(jz) = ug
            vbar(jz) = vg
            wbar(jz) = 0._rprec
        else
            ubar(jz) = arg
            vbar(jz) = 0._rprec
            wbar(jz) = 0._rprec
        end if

        if (zuvp(jz) .gt. z_inv/z_i) then
            ubar(jz) = ug
        end if
    end do
end if 

rms = 1._rprec
bv_freq = (g/T_scale*inv_strength)**(1._rprec/2._rprec)
do jz = 1, nz-1
    jz_abs = coordz*(nz - 1) + jz
    z = (coordz*(nz-1) + jz - 0.5_rprec)*dz*z_i
    seed = -80 - jz_abs !--trying to make consistent init for MPI
    do jy = 1, nynpy
    do jx = 1, nxt
        !...Ran3 returns uniform RV between 0 and 1. (sigma_rv=0.289)
        !...Taking std dev of vel as 1 at all heights
        if ( z .le. z_inv ) then
            noise = rms/.289_rprec*(ran3(seed) - .5_rprec)
            u(jx, jy, jz) = noise*(1._rprec - z/z_i)*w_star/u_scale + ubar(jz)
            noise = rms/.289_rprec*(ran3(seed) - .5_rprec)
            v(jx, jy, jz) = noise*(1._rprec - z/z_i)*w_star/u_scale + vbar(jz)
            noise = rms/.289_rprec*(ran3(seed) - .5_rprec)
            w(jx, jy, jz) = noise*(1._rprec - z/z_i)*w_star/u_scale + wbar(jz)
            noise = rms/0.289_rprec*(ran3(seed) - .5_rprec)
            theta(jx, jy, jz) = (theta_mean + 10._rprec*noise*(1 - z/z_i)*T_star)/T_scale
            !!!noise = rms/0.289_rprec*(ran3(seed) - .5_rprec)
            !!!q(jx, jy, jz) = q_mix + 50._rprec*noise*(1 - z/z_i)*q_star
        else
            u(jx, jy, jz) = ubar(jz)
            v(jx, jy, jz) = vbar(jz)
            w(jx, jy, jz) = wbar(jz)
            theta(jx, jy, jz) = (theta_mean + (z - z_inv)*inv_strength)/T_scale
            !!q(jx, jy, jz) = q_mix
        end if

        if (coordz == 0 .and. (jz == 1 .or. jz == 2)) then
            noise = rms/.289_rprec*(ran3(seed) - .5_rprec)
            u(jx, jy, jz) = u(jx, jy, jz) + noise*w_star/u_scale
            v(jx, jy, jz) = v(jx, jy, jz) + noise*w_star/u_scale
            w(jx, jy, jz) = w(jx, jy, jz) + noise*w_star/u_scale
            theta(jx, jy, jz) = theta(jx, jy, jz) + 10._rprec*noise*T_star/T_scale
        end if
    end do
    end do
end do
! ---
if ( coordz == 0 ) then
    w(1:nxt, 1:nynpy, 1) = 0._rprec
end if
if ( coordz == npz - 1 ) then
    w(1:nxt, 1:nynpy, nz) = 0._rprec
endif

if ( coordz == npz - 1 ) then
    if ( ubc_mom == 'wall' ) then
        u(1:nxt, 1:nynpy, nz) = BOGUS
        v(1:nxt, 1:nynpy, nz) = BOGUS
        theta(1:nxt, 1:nynpy, nz) = BOGUS
    else if ( ubc_mom == 'stress free' ) then
        u(1:nxt, 1:nynpy, nz) = u(1:nxt, 1:nynpy, nz - 1)
        v(1:nxt, 1:nynpy, nz) = v(1:nxt, 1:nynpy, nz - 1)
        theta(1:nxt, 1:nynpy, nz) = theta(1:nxt, 1:nynpy, nz - 1) + dTdz_top/T_scale*z_i*dz
    end if
end if
! ---
end subroutine ic_scal
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine ic_read()
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
use types, only:rprec
use param
use sim_param, only:u, v, w, theta, q
use bottombc
implicit none
! ---
real, dimension(nzt-1) :: ubar_tot, vbar_tot, wbar_tot, tbar_tot
real(kind=rprec), dimension(nz-1) :: ubar, vbar, wbar, tbar
real(kind=rprec) :: rms, noise, arg, theta_mean
real(kind=rprec) :: z, w_star, T_star, q_star, ran3
integer :: jx, jy, jz, seed, jz_abs
logical :: exst

! --- read initial profile to minimize inertial oscillation
if ( coordz == 0 ) then
    inquire(file='../readin/ic.dat', exist=exst)
    if (exst) then
        open(9, file='../readin/ic.dat')
        do jz = 1, nzt-1
            read(9,*) ubar_tot(jz), vbar_tot(jz), tbar_tot(jz)
        end do
        close(9)
    else
        write(*, *) 'file ../readin/ic.dat NOT found'
        stop
    end if
end if
call scatter_z(ubar_tot, ubar)
call scatter_z(vbar_tot, vbar)
! call scatter_z(wbar_tot, wbar)
call scatter_z(tbar_tot, tbar)

! if ( rank == 0 ) then
!     open(9, file='../output/check_ic_profile.csv')
!     write(9,*) 'z,u,v,T'
!     do jz = 1, nzt-1
!         z = -(real(jz)-0.5_rprec)*dz*z_i
!         write(9,'(e15.7, a, e15.7, a, e15.7, a, e15.7)')  &
!                 z, ',', ubar_tot(jz), ',', vbar_tot(jz), ',',    &
!                 tbar_tot(jz)  
!     end do
!     close(9) 
! end if

if ( lbc .eq. 0 ) then
    theta_mean = T_s_min - 0.001*T_s_min !T_s_min is dimensional while T_s is non-dimensional
else 
    theta_mean = T_init
end if

if (wt_s .eq. 0.0_rprec) then
    ! w_star is of O(1) with z_i=500 and wt_s=0.06 for atmos
    ! For ocean, ~5e-5 is better
    w_star = (g/theta_mean*5.3691d-5*z_i)**(1._rprec/3._rprec)
    T_star = 5.3691d-5/w_star
    q_star = T_star
else
    w_star = sign((g/theta_mean*abs(wt_s)*z_i)**(1._rprec/3._rprec), wt_s)
    T_star = wt_s/w_star
    q_star = T_star
end if


rms = 1._rprec
do jz = 1, nz-1
    jz_abs = coordz*(nz - 1) + jz
    z = (coordz*(nz-1) + jz - 0.5_rprec)*dz*z_i
    seed = -80 - jz_abs !--trying to make consistent init for MPI
    do jy = 1, nynpy
    do jx = 1, nxt
        !...Ran3 returns uniform RV between 0 and 1. (sigma_rv=0.289)
        !...Taking std dev of vel as 1 at all heights
        if (coordz .eq. 0 .and. (jz == 1 .or. jz == 2)) then
            noise = rms/.289_rprec*(ran3(seed) - .5_rprec)
            u(jx, jy, jz) = noise*(1._rprec - z/z_i)*w_star/u_scale + ubar(jz)
            noise = rms/.289_rprec*(ran3(seed) - .5_rprec)
            v(jx, jy, jz) = noise*(1._rprec - z/z_i)*w_star/u_scale + vbar(jz)
            noise = rms/.289_rprec*(ran3(seed) - .5_rprec)
            w(jx, jy, jz) = noise*(1._rprec - z/z_i)*w_star/u_scale
            noise = rms/0.289_rprec*(ran3(seed) - .5_rprec)
            theta(jx, jy, jz) = noise*(1._rprec - z/z_i)*T_star/T_scale + tbar(jz) 
        else
            u(jx, jy, jz) = ubar(jz)
            v(jx, jy, jz) = vbar(jz)
            w(jx, jy, jz) = 0._rprec
            theta(jx, jy, jz) = tbar(jz)
        end if
    end do
    end do
end do
! ---
if ( coordz == 0 ) then
    w(1:nxt, 1:nynpy, 1) = 0._rprec
end if
if ( coordz == npz - 1 ) then
    w(1:nxt, 1:nynpy, nz) = 0._rprec
endif

if ( coordz == npz - 1 ) then
    if ( ubc_mom == 'wall' ) then
        u(1:nxt, 1:nynpy, nz) = BOGUS
        v(1:nxt, 1:nynpy, nz) = BOGUS
        theta(1:nxt, 1:nynpy, nz) = BOGUS
    else if ( ubc_mom == 'stress free' ) then
        u(1:nxt, 1:nynpy, nz) = u(1:nxt, 1:nynpy, nz - 1)
        v(1:nxt, 1:nynpy, nz) = v(1:nxt, 1:nynpy, nz - 1)
        theta(1:nxt, 1:nynpy, nz) = theta(1:nxt, 1:nynpy, nz - 1) + dTdz_top/T_scale*z_i*dz
    end if
end if
! ---
end subroutine ic_read
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine nutrient_init_read()
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
use types, only:rprec
use param
use sim_param, only:PCon
use bottombc
implicit none
! ---
real, dimension(:), allocatable :: PCon_inflow, hgt_PCon
real(kind=rprec), dimension(:), allocatable :: nutrient_init
integer :: nz_pcon, jz, jh
logical :: exst
real(kind=rprec) :: c0, c1
character(80) :: fname
integer :: ipcon
! --- read initial profile to minimize inertial oscillation
allocate(nutrient_init(nz-1))
do ipcon = 1, npcon
    nutrient_init = 0.0_rprec
    write(fname, '(a,i0,a)') '../readin/nutrient_', ipcon, '.dat'
    inquire(file=fname, exist=exst)
    if (exst) then
        open(9, file=fname)
        read(9,*) nz_pcon
        allocate(PCon_inflow(nz_pcon), hgt_PCon(nz_pcon))
        do jz = 1, nz_pcon
            read(9,*) hgt_PCon(jz), PCon_inflow(jz)
        end do
        close(9)

        do jz = 1, nz-1
        do jh = 1, nz_pcon-1
            if ( zuvp(jz)*z_i .ge. -hgt_PCon(jh) .and. zuvp(jz)*z_i .le. -hgt_PCon(jh+1) ) then
                c0 = (-hgt_PCon(jh+1)-zuvp(jz)*z_i) / (-hgt_PCon(jh+1)-(-hgt_PCon(jh)))
                c1 = 1._rprec - c0
                nutrient_init(jz) = (c0*PCon_inflow(jh) + c1*PCon_inflow(jh+1))*1000._rprec ! unit conversion from liter-1 to m-3 
            end if
        end do
        end do

        deallocate(PCon_inflow, hgt_PCon)
    else
        write(*, *) 'file ../readin/nutrient.dat NOT found'
        stop
    end if

    do jz = 1, nz-1
        PCon(:, 1:nynpy, jz, ipcon) = nutrient_init(jz)
    end do
end do
deallocate(nutrient_init)
! ---
end subroutine nutrient_init_read
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine nutrient_inflow_read()
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
use types, only:rprec
use param
use sim_param, only:PCon
use bottombc
use canopy, only:nutrient_inflow
implicit none
! ---
real, dimension(:), allocatable :: PCon_inflow, hgt_PCon
integer :: nz_pcon, jz, jh
logical :: exst
real(kind=rprec) :: c0, c1
character(80) :: fname
integer :: ipcon
! --- read initial profile to minimize inertial oscillation
allocate(nutrient_inflow(nz-1, npcon))
nutrient_inflow = 0.0_rprec
do ipcon = 1, npcon
    write(fname, '(a,i0,a)') '../readin/nutrient_', ipcon, '.dat'
    inquire(file=fname, exist=exst)
    if (exst) then
        open(9, file=fname)
        read(9,*) nz_pcon
        allocate(PCon_inflow(nz_pcon), hgt_PCon(nz_pcon))
        do jz = 1, nz_pcon
            read(9,*) hgt_PCon(jz), PCon_inflow(jz)
        end do
        close(9)

        do jz = 1, nz-1
        do jh = 1, nz_pcon-1
            if ( zuvp(jz)*z_i .ge. -hgt_PCon(jh) .and. zuvp(jz)*z_i .le. -hgt_PCon(jh+1) ) then
                c0 = (-hgt_PCon(jh+1)-zuvp(jz)*z_i) / (-hgt_PCon(jh+1)-(-hgt_PCon(jh)))
                c1 = 1._rprec - c0
                nutrient_inflow(jz,ipcon) = (c0*PCon_inflow(jh) + c1*PCon_inflow(jh+1))*1000._rprec ! unit conversion from liter-1 to m-3 
            end if
        end do
        end do

        deallocate(PCon_inflow, hgt_PCon)
    else
        write(*, *) 'file ../readin/nutrient.dat NOT found'
        stop
    end if

    do jz = 1, nz-1
        PCon(1, 1:nynpy, jz, ipcon) = nutrient_inflow(jz, ipcon)
    end do
end do
! deallocate(nutrient_inflow)
! ---
end subroutine nutrient_inflow_read
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
function ran3(idum)
!-----------------------------------------------------------------------
!   random number generator
!-----------------------------------------------------------------------
use types, only:rprec
implicit none
! ---
integer(kind=4) :: idum
real(kind=rprec)    :: ran3
! ---
integer, parameter :: mbig  = 1000000000
integer, parameter :: mseed = 161803398
integer, parameter :: mz    = 0
real(kind=rprec), parameter :: fac = 1.d0/mbig
!parameter (mbig=4000000.,mseed=1618033.,mz=0.,fac=1./mbig)
! ---
integer :: i, iff, ii, inext, inextp, k
integer :: mj, mk, ma(55)
save iff, inext, inextp, ma
data iff /0/
! ---  
if( idum .lt. 0 .or. iff .eq. 0 ) then
    iff = 1
    mj = mseed - iabs(idum)
    mj = mod(mj, mbig)
    ma(55) = mj
    mk = 1
    do i = 1, 54
        ii = mod(21*i, 55)
        ma(ii) = mk
        mk = mj - mk
        if(mk .lt. mz) mk = mk + mbig
        mj = ma(ii)
    end do
    
    do k = 1, 4
        do i = 1, 55
            ma(i) = ma(i) - ma(1+mod(i+30,55))
            if(ma(i) .lt. mz) ma(i) = ma(i) + mbig
        end do
    end do
    inext = 0
    inextp = 31
    idum = 1
end if
  
inext = inext + 1
if(inext .eq. 56) inext = 1
inextp = inextp + 1
if(inextp .eq. 56) inextp = 1
mj = ma(inext) - ma(inextp)
if(mj .lt. mz) mj = mj + mbig
ma(inext) = mj
ran3 = mj * fac
! ---
return
end function ran3
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine read_flow_field
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
use types, only:rprec
use param
use sim_param
use sgsmodule
use io, only:fcumulative_time
implicit none
! ---
character(80) :: fname
logical :: exst
integer(MPI_OFFSET_KIND) :: disp0, disp3d 
! ---
inquire (file=fcumulative_time, exist=exst)
if (exst) then
    open(1, file=fcumulative_time)
    read(1, *) nums
    close(1)
else
    write(*, *) "Restart file NOT found"
    nums = 0
end if

disp3d = np  * sizeof(u(1:nxt, 1:nynpy, 1:nz-1))

write (fname, '(a,i8.8,a)') path_restart//'vel_tt', nums, '.out'

disp0 = 0
call read_file_mpi_3d(u(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
disp0 = disp0+disp3d
call read_file_mpi_3d(v(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
disp0 = disp0+disp3d
call read_file_mpi_3d(w(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
disp0 = disp0+disp3d
call read_file_mpi_3d(RHSx(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
disp0 = disp0+disp3d
call read_file_mpi_3d(RHSy(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
disp0 = disp0+disp3d
call read_file_mpi_3d(RHSz(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
disp0 = disp0+disp3d
call read_file_mpi_3d(Cs_opt2(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
disp0 = disp0+disp3d
call read_file_mpi_3d(F_LM(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
disp0 = disp0+disp3d
call read_file_mpi_3d(F_MM(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
disp0 = disp0+disp3d
call read_file_mpi_3d(F_QN(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
disp0 = disp0+disp3d
call read_file_mpi_3d(F_NN(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
! ---
end subroutine read_flow_field
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine read_precursor_field
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
use types, only:rprec
use param
use sim_param
use sgsmodule
use io, only:fcumulative_time
implicit none
! ---
real(rprec), dimension(nxt/2, nynpy, nz-1) :: u_tmp, v_tmp, w_tmp, RHSx_tmp,    &
                                              RHSy_tmp, RHSz_tmp, Cs_opt2_tmp,   &
                                              F_LM_tmp, F_MM_tmp, F_QN_tmp, F_NN_tmp

character(80) :: fname
logical :: exst
integer(MPI_OFFSET_KIND) :: disp0, disp3d 
integer :: jx, jy, jz
! ---
inquire (file=fcumulative_time, exist=exst)
if (exst) then
    open(1, file=fcumulative_time)
    read(1, *) nums
    close(1)
else
    write(*, *) "Restart file NOT found"
    nums = 0
end if

disp3d = np  * sizeof(u_tmp(1:nxt/2, 1:nynpy, 1:nz-1))

write (fname, '(a,i8.8,a)') path_restart//'vel_tt', nums, '.out'

disp0 = 0
call read_file_mpi_3d_inflow(u_tmp(1:nxt/2, 1:nynpy, 1:nz-1), fname, disp0)
disp0 = disp0+disp3d
call read_file_mpi_3d_inflow(v_tmp(1:nxt/2, 1:nynpy, 1:nz-1), fname, disp0)
disp0 = disp0+disp3d
call read_file_mpi_3d_inflow(w_tmp(1:nxt/2, 1:nynpy, 1:nz-1), fname, disp0)
disp0 = disp0+disp3d
call read_file_mpi_3d_inflow(RHSx_tmp(1:nxt/2, 1:nynpy, 1:nz-1), fname, disp0)
disp0 = disp0+disp3d
call read_file_mpi_3d_inflow(RHSy_tmp(1:nxt/2, 1:nynpy, 1:nz-1), fname, disp0)
disp0 = disp0+disp3d
call read_file_mpi_3d_inflow(RHSz_tmp(1:nxt/2, 1:nynpy, 1:nz-1), fname, disp0)
disp0 = disp0+disp3d
call read_file_mpi_3d_inflow(Cs_opt2_tmp(1:nxt/2, 1:nynpy, 1:nz-1), fname, disp0)
disp0 = disp0+disp3d
call read_file_mpi_3d_inflow(F_LM_tmp(1:nxt/2, 1:nynpy, 1:nz-1), fname, disp0)
disp0 = disp0+disp3d
call read_file_mpi_3d_inflow(F_MM_tmp(1:nxt/2, 1:nynpy, 1:nz-1), fname, disp0)
disp0 = disp0+disp3d
call read_file_mpi_3d_inflow(F_QN_tmp(1:nxt/2, 1:nynpy, 1:nz-1), fname, disp0)
disp0 = disp0+disp3d
call read_file_mpi_3d_inflow(F_NN_tmp(1:nxt/2, 1:nynpy, 1:nz-1), fname, disp0)
! ---
u(1:nxt/2, 1:nynpy, 1:nz-1) = u_tmp(1:nxt/2, 1:nynpy, 1:nz-1)
v(1:nxt/2, 1:nynpy, 1:nz-1) = v_tmp(1:nxt/2, 1:nynpy, 1:nz-1)
w(1:nxt/2, 1:nynpy, 1:nz-1) = w_tmp(1:nxt/2, 1:nynpy, 1:nz-1)
RHSx(1:nxt/2, 1:nynpy, 1:nz-1) = RHSx_tmp(1:nxt/2, 1:nynpy, 1:nz-1)
RHSy(1:nxt/2, 1:nynpy, 1:nz-1) = RHSy_tmp(1:nxt/2, 1:nynpy, 1:nz-1)
RHSz(1:nxt/2, 1:nynpy, 1:nz-1) = RHSz_tmp(1:nxt/2, 1:nynpy, 1:nz-1)
Cs_opt2(1:nxt/2, 1:nynpy, 1:nz-1) = Cs_opt2_tmp(1:nxt/2, 1:nynpy, 1:nz-1)
F_LM(1:nxt/2, 1:nynpy, 1:nz-1) = F_LM_tmp(1:nxt/2, 1:nynpy, 1:nz-1)
F_MM(1:nxt/2, 1:nynpy, 1:nz-1) = F_MM_tmp(1:nxt/2, 1:nynpy, 1:nz-1)
F_QN(1:nxt/2, 1:nynpy, 1:nz-1) = F_QN_tmp(1:nxt/2, 1:nynpy, 1:nz-1)
F_NN(1:nxt/2, 1:nynpy, 1:nz-1) = F_NN_tmp(1:nxt/2, 1:nynpy, 1:nz-1)

u(nxt/2+1:nxt, 1:nynpy, 1:nz-1) = u_tmp(1:nxt/2, 1:nynpy, 1:nz-1)
v(nxt/2+1:nxt, 1:nynpy, 1:nz-1) = v_tmp(1:nxt/2, 1:nynpy, 1:nz-1)
w(nxt/2+1:nxt, 1:nynpy, 1:nz-1) = w_tmp(1:nxt/2, 1:nynpy, 1:nz-1)
RHSx(nxt/2+1:nxt, 1:nynpy, 1:nz-1) = RHSx_tmp(1:nxt/2, 1:nynpy, 1:nz-1)
RHSy(nxt/2+1:nxt, 1:nynpy, 1:nz-1) = RHSy_tmp(1:nxt/2, 1:nynpy, 1:nz-1)
RHSz(nxt/2+1:nxt, 1:nynpy, 1:nz-1) = RHSz_tmp(1:nxt/2, 1:nynpy, 1:nz-1)
Cs_opt2(nxt/2+1:nxt, 1:nynpy, 1:nz-1) = Cs_opt2_tmp(1:nxt/2, 1:nynpy, 1:nz-1)
F_LM(nxt/2+1:nxt, 1:nynpy, 1:nz-1) = F_LM_tmp(1:nxt/2, 1:nynpy, 1:nz-1)
F_MM(nxt/2+1:nxt, 1:nynpy, 1:nz-1) = F_MM_tmp(1:nxt/2, 1:nynpy, 1:nz-1)
F_QN(nxt/2+1:nxt, 1:nynpy, 1:nz-1) = F_QN_tmp(1:nxt/2, 1:nynpy, 1:nz-1)
F_NN(nxt/2+1:nxt, 1:nynpy, 1:nz-1) = F_NN_tmp(1:nxt/2, 1:nynpy, 1:nz-1)
! ---
end subroutine read_precursor_field
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine ekman_layer(ubar, vbar)
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------    
use param
use stokes_drift, only: U_stokes, wavenm_w
implicit none
! ---
integer :: jz
real(kind=rprec):: z, nu
real(kind=rprec), dimension(nz), intent(inout) :: ubar, vbar
complex(kind=rprec) :: gamma, Ubarcomp
! ---
nu = nu_ek/(u_scale*z_i)

do jz = 1, nz
    z = zuvp(jz)
    
    gamma = cmplx(0._rprec, 1._rprec)*coriol*U_stokes/            &
            (4._rprec*wavenm_w**2*nu - cmplx(0._rprec, 1._rprec)*coriol)

    Ubarcomp = cmplx(1._rprec, -1._rprec) / sqrt(2._rprec*coriol*nu) &
             * (u_star - 2._rprec*wavenm_w*nu*gamma &
             * cmplx(cos(agl_Ust), sin(agl_Ust)))*exp(cmplx(1._rprec, 1._rprec) &
             / sqrt(2._rprec)*sqrt(coriol/nu)*(-z)) &
             + cmplx(cos(agl_Ust), sin(agl_Ust))*gamma*exp(2._rprec*wavenm_w*(-z))
    ubar(jz) = real(Ubarcomp) + ug
    vbar(jz) = aimag(Ubarcomp) + vg
end do
! ---
end subroutine ekman_layer
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine init_profile(ubar, vbar, z_inv)
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------    
use param
use stokes_drift, only: U_stokes, wavenm_w
implicit none
! ---
integer :: jz
real(kind=rprec):: z, nu, Kvis
real(kind=rprec), intent(in) :: z_inv
real(kind=rprec), dimension(nz), intent(inout) :: ubar, vbar
complex(kind=rprec) :: gamma, Ubarcomp
! ---
nu = nu_ek/(u_scale*z_i)
Kvis = 1.e-4_rprec

do jz = 1, nz
    z = zuvp(jz)
    
    if ( z .le. z_inv/z_i ) then
        gamma = cmplx(0._rprec, 1._rprec)*coriol*U_stokes/            &
                (4._rprec*wavenm_w**2*nu - cmplx(0._rprec, 1._rprec)*coriol)

        Ubarcomp = cmplx(1._rprec, -1._rprec) / sqrt(2._rprec*coriol*nu) &
                    * (u_star - 2._rprec*wavenm_w*nu*gamma &
                    * cmplx(cos(agl_Ust), sin(agl_Ust)))*exp(cmplx(1._rprec, 1._rprec) &
                    / sqrt(2._rprec)*sqrt(coriol/nu)*(-z)) &
                    + cmplx(cos(agl_Ust), sin(agl_Ust))*gamma*exp(2._rprec*wavenm_w*(-z))
        ubar(jz) = real(Ubarcomp) + ug
        vbar(jz) = aimag(Ubarcomp) + vg
    else
        gamma = sqrt(freq_coriolis/2._rprec/Kvis)
        ubar(jz) = 0.01_rprec + ug*(1-exp(-gamma*(lz_tot-z)*z_i)*cos(gamma*(lz_tot-z)*z_i))
        vbar(jz) = ug*exp(-gamma*(lz_tot-z)*z_i)*sin(gamma*(lz_tot-z)*z_i)
    end if
end do
! ---
end subroutine init_profile
