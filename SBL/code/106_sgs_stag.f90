!--------------------------------------------------------------------!
! put everything onto w-nodes, follow original version
!--provides txx, txy, tyy, tzz for jz=1:nz-1; txz, tyz for 1:nz
!--------------------------------------------------------------------!  
subroutine sgs_stag()
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!  
use types, only: rprec
use param
use sim_param, only: u, v, w, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx,  &
                     dwdy, dwdz, txx, txy, txz, tyy, tyz, tzz
use sgsmodule, only: u_lag, v_lag, w_lag, Cs_opt2, Nu_t, dissip, magS
!------------------------- Vij Comment begins---------------
! 04/14/2004 - Added Nu_t to the list of variables from sgsmodule
! 05/19/2004 - Replaced all references to visc by Nu_t; deleted local
!              declaration of visc
! 05/19/2004 - Replace all .5 by 0.5
!-------------------------Vij Comment ends ------------------
use test_filtermodule, only: filter_size
implicit none
! ---  
real(kind=rprec), dimension(nz) :: l, ziko, zz
real(kind=rprec), dimension(ldx, nynpy, nz) :: S11, S12, S22, S33, S13, S23
!real(kind=rprec),dimension(ldx,nynpy,nz):: dissip
real(kind=rprec) :: ux, uy, uz, vx, vy, vz, wx, wy, wz
real(kind=rprec), dimension(ldx, nynpy) :: txzp, tyzp, S
real(kind=rprec) :: delta, nu, const
integer :: jx, jy, jz, jz_min, jz_max
! ---
character(80) :: fname
! ---
delta = filter_size*(dx*dy*dz)**(1._rprec/3._rprec) ! nondimensional
! Cs is Smagorinsky's constant. l is a filter size (non-dim.)
nu = 0._rprec
if (molec) then
    nu = (nu_molec/(u_scale*z_i)) ! dim/less
end if

if ( coordz == 0 ) then
    ! save the wall level
    ! txzp(:, :) = txz(:, :, 1)
    ! tyzp(:, :) = tyz(:, :, 1)

    ! calculate Sij on w-nodes
    ! at level z=dz/2.  (Where dudz and dvdz are on UVP nodes)
    do jy = 1, nynpy
    do jx = 1, nxt
        ux = dudx(jx, jy, 1) ! uvp-node
        uy = dudy(jx, jy, 1) ! uvp-node
        uz = dudz(jx, jy, 1) ! uvp-node
        vx = dvdx(jx, jy, 1) ! uvp-node
        vy = dvdy(jx, jy, 1) ! uvp-node
        vz = dvdz(jx, jy, 1) ! uvp-node
        ! special case
        wx = 0.5_rprec*(dwdx(jx, jy, 1) + dwdx(jx, jy, 2)) ! uvp-node
        wy = 0.5_rprec*(dwdy(jx, jy, 1) + dwdy(jx, jy, 2)) ! uvp-node
        wz = dwdz(jx, jy, 1) ! uvp-node
        S11(jx, jy, 1) = ux ! uvp-node
        S12(jx, jy, 1) = 0.5_rprec*(uy + vx) ! uvp-node
        ! taken care of with wall stress routine
        S13(jx, jy, 1) = 0.5_rprec*(uz + wx) ! uvp
        S22(jx, jy, 1) = vy ! uvp-node
        ! taken care of with wall stress routine
        S23(jx, jy, 1) = 0.5_rprec*(vz + wy) ! uvp
        S33(jx, jy, 1) = wz ! uvp-node
    end do
    end do

    jz_min = 2
else
    jz_min = 1
end if


if ( coordz == npz - 1 ) then
    ! save the wall level
    ! txzp(:, :) = txz(:, :, nz)
    ! tyzp(:, :) = tyz(:, :, nz)

    ! calculate Sij on w-nodes
    ! at level z=dz/2.  (Where dudz and dvdz are on UVP nodes)
    do jy = 1, nynpy
    do jx = 1, nxt
        ux = dudx(jx, jy, nz-1) ! uvp-node
        uy = dudy(jx, jy, nz-1) ! uvp-node
        uz = dudz(jx, jy, nz) ! uvp-node
        vx = dvdx(jx, jy, nz-1) ! uvp-node
        vy = dvdy(jx, jy, nz-1) ! uvp-node
        vz = dvdz(jx, jy, nz) ! uvp-node
        ! special case
        wx = 0.5_rprec*(dwdx(jx, jy, nz-1) + dwdx(jx, jy, nz)) ! uvp-node
        wy = 0.5_rprec*(dwdy(jx, jy, nz-1) + dwdy(jx, jy, nz)) ! uvp-node
        wz = dwdz(jx, jy, nz-1) ! uvp-node
        S11(jx, jy, nz) = ux ! uvp-node
        S12(jx, jy, nz) = 0.5_rprec*(uy + vx) ! uvp-node
        ! taken care of with wall stress routine
        S13(jx, jy, nz) = 0.5_rprec*(uz + wx) ! uvp
        S22(jx, jy, nz) = vy ! uvp-node
        ! taken care of with wall stress routine
        S23(jx, jy, nz) = 0.5_rprec*(vz + wy) ! uvp
        S33(jx, jy, nz) = wz ! uvp-node
    end do
    end do

    jz_max = nz - 1
else
    jz_max = nz
end if

!
!--this is only required b/c of the unnatural placement of all strains
!  onto w-nodes, be careful not to overwrite nz on top process with garbage
!--dwdz(jz=0) is already known, except at bottom process (OK)
!call mpi_sendrecv (dwdz(1, 1, nz-1), ldx*nynpy, MPI_RPREC, up,   1,    &
!                   dwdz(1, 1, 0),    ldx*nynpy, MPI_RPREC, down, 1,    &
!                   comm, status, ierr)
call mpi_sendrecv(dwdz(1, 1, 1),  ldx*nynpy, MPI_RPREC, down, 2,   &
                  dwdz(1, 1, nz), ldx*nynpy, MPI_RPREC, up,   2,   &
                  comm, status, ierr)
!
! calculate derivatives/strain on w-nodes
!--in MPI version, calculating up to nz saves some interprocess exchange
!  later but note that dwdz is not provided w/o some communication
!  (unless its the top process)
do jz = jz_min, jz_max
do jy = 1, nynpy
do jx = 1, nxt
    ux = 0.5_rprec*(dudx(jx, jy, jz) + dudx(jx, jy, jz - 1)) ! w-node
    uy = 0.5_rprec*(dudy(jx, jy, jz) + dudy(jx, jy, jz - 1)) ! w-node
    uz = dudz(jx, jy, jz) ! w-node
    vx = 0.5_rprec*(dvdx(jx, jy, jz) + dvdx(jx, jy, jz - 1)) ! w-node
    vy = 0.5_rprec*(dvdy(jx, jy, jz) + dvdy(jx, jy, jz - 1)) ! w-node
    vz = dvdz(jx, jy, jz) ! w-node
    wx = dwdx(jx, jy, jz) ! w-node
    wy = dwdy(jx, jy, jz) ! w-node
    wz = 0.5_rprec*(dwdz(jx, jy, jz) + dwdz(jx, jy, jz - 1)) ! w-node
    S11(jx, jy, jz) = ux ! w-node
    S12(jx, jy, jz) = 0.5_rprec*(uy + vx) ! w-node
    S13(jx, jy, jz) = 0.5_rprec*(uz + wx) ! w-node
    S22(jx, jy, jz) = vy ! w-node
    S23(jx, jy, jz) = 0.5_rprec*(vz + wy) ! w-node
    S33(jx, jy, jz) = wz ! w-node
end do
end do
end do


! This part computes the average velocity during cs_count times steps
! This is used with the lagrangian model only
if (model == 4 .or. model == 5) then
    u_lag = u_lag + u
    v_lag = v_lag + v
    w_lag = w_lag + w
end if

if (sgs) then
    if (model == 1) then !For traditional Smagorinsky   ref02

        ! Define parameters (Co and nn) for wallfunction
        Cs_opt2 = Co**2 ! constant coefficient

        do jz = 1, nz
            if ( coordz == 0 .and. jz == 1 ) then
                zz(jz) = (jz - 0.5_rprec)*dz
            else if ( coordz == npz-1 .and. jz == nz ) then
                zz(jz) = ((jz-1-0.5_rprec) + coordz*(nz-1))*dz
            else
                zz(jz) = ((jz-1) + coordz*(nz-1))*dz
            end if

            if ( ubc_mom == 'wall' .and. zz(jz) .ge. lz_tot/2. ) then
                zz(jz) = lz_tot - zz(jz)
            end if

            l(jz) = ( Co**(nnn)*(vonk*zz(jz))**(-nnn) + &
                        (delta)**(-nnn) )**(-1._rprec/nnn)
        end do

    else ! for dynamic procedures: initialization  !cont ref02
        l = delta ! constant equal to delta
        if ( jt == 1 .and. inilag ) then
            Cs_opt2 = 0.03_rprec
        ! --- make sure that the "node conventions" in these dynamic models match with those done here
        else if ((jt .ge. DYN_init .or. restart_vel) .and. (mod(jt, cs_count) == 0)) then
            if (model == 2) then
                !--provides ziko 1:nz
                call std_dynamic(ziko, S11, S12, S13, S22, S23, S33)
                forall (jz=1:nz) Cs_opt2(:, :, jz) = ziko(jz)
            else if (model == 3) then
                ! Plane average dynamic  continue ref 05
                call scaledep_dynamic(ziko, u, v, w, S11, S12, S13, S22, S23, S33)
                forall (jz=1:nz) Cs_opt2(:, :, jz) = ziko(jz)
            else if (model == 4) then
                call lagrange_Ssim(S11, S12, S13, S22, S23, S33)
            else if (model == 5) then
                call lagrange_Sdep(S11, S12, S13, S22, S23, S33)
            end if
        end if
    end if
end if


! define |S| and viscosity on w-nodes (on uvp node for jz=1)
!--MPI: going up to nz saves communication
do jz = 1, nz
do jy = 1, nynpy
do jx = 1, nxt
    S(jx, jy) = sqrt(2._rprec*(S11(jx, jy, jz)**2 + S22(jx, jy, jz)**2 + S33(jx, jy, jz)**2 + &
                     2._rprec*(S12(jx, jy, jz)**2 + S13(jx, jy, jz)**2 + S23(jx, jy, jz)**2)))
    Nu_t(jx, jy, jz) = S(jx, jy)*Cs_opt2(jx, jy, jz)*(l(jz)**2)     ! evaluate nu_t= c_s^2 l^2 |S|
    dissip(jx, jy, jz) = -(S(jx, jy)**3)*Cs_opt2(jx, jy, jz)*(l(jz)**2)    
    magS(jx, jy, jz) = S(jx, jy)        ! Store |S|
end do
end do
end do

! ---
if ( coordz == 0 ) then
    ! bottom
    do jy = 1, nynpy
    do jx = 1, nxt
        const = 0._rprec
        if (sgs) then
            const = Nu_t(jx, jy, 1)
        end if
        ! everything on right node here
        txx(jx, jy, 1) = -2._rprec*(const + nu)*S11(jx, jy, 1)
        txy(jx, jy, 1) = -2._rprec*(const + nu)*S12(jx, jy, 1)
        tyy(jx, jy, 1) = -2._rprec*(const + nu)*S22(jx, jy, 1)
        tzz(jx, jy, 1) = -2._rprec*(const + nu)*S33(jx, jy, 1)
    end do
    end do

    jz_min = 2 !--wall level set by wallstress routine
else
    jz_min = 1
end if

if ( coordz == npz - 1 ) then
    ! bottom
    do jy = 1, nynpy
    do jx = 1, nxt
        const = 0._rprec
        if (sgs) then
            const = Nu_t(jx, jy, nz)
        end if
        ! everything on right node here
        txx(jx, jy, nz-1) = -2._rprec*(const + nu)*S11(jx, jy, nz)
        txy(jx, jy, nz-1) = -2._rprec*(const + nu)*S12(jx, jy, nz)
        tyy(jx, jy, nz-1) = -2._rprec*(const + nu)*S22(jx, jy, nz)
        tzz(jx, jy, nz-1) = -2._rprec*(const + nu)*S33(jx, jy, nz)
    end do
    end do

    jz_max = nz - 2 !--wall level set by wallstress routine
else
    jz_max = nz - 1
end if

!--note that nz level is available for the interpolations at (jz+1) here
!--txx, tyy, tzz, txy not needed at nz
do jz = jz_min, jz_max
do jy = 1, nynpy
do jx = 1, nxt
    const = 0._rprec
    if (sgs) then
        const = 0.5_rprec*(Nu_t(jx, jy, jz) + Nu_t(jx, jy, jz + 1))
    end if
    txx(jx, jy, jz) = -2._rprec*(const + nu)* &
           0.5_rprec*(S11(jx, jy, jz) + S11(jx, jy, jz + 1))
    txy(jx, jy, jz) = -2._rprec*(const + nu)* &
           0.5_rprec*(S12(jx, jy, jz) + S12(jx, jy, jz + 1))
    tyy(jx, jy, jz) = -2._rprec*(const + nu)* &
           0.5_rprec*(S22(jx, jy, jz) + S22(jx, jy, jz + 1))
    tzz(jx, jy, jz) = -2._rprec*(const + nu)* &
           0.5_rprec*(S33(jx, jy, jz) + S33(jx, jy, jz + 1))
end do
end do
end do


do jz = jz_min, nz - 1
do jy = 1, nynpy
do jx = 1, nxt
    const = 0._rprec
    if (sgs) const = Nu_t(jx, jy, jz)
    txz(jx, jy, jz) = -2._rprec*(const + nu) * S13(jx, jy, jz)
    tyz(jx, jy, jz) = -2._rprec*(const + nu) * S23(jx, jy, jz)
end do
end do
end do
!
!--recv information for top nodes: txy, txz only
!--other components not needed, since no equation is solved there
!  (values are just copied)
call mpi_sendrecv(txz(1, 1, 1),  ldx*nynpy, MPI_RPREC, down, 3,    &
                  txz(1, 1, nz), ldx*nynpy, MPI_RPREC, up,   3,    &
                  comm, status, ierr)
call mpi_sendrecv(tyz(1, 1, 1),  ldx*nynpy, MPI_RPREC, down, 4,    &
                  tyz(1, 1, nz), ldx*nynpy, MPI_RPREC, up,   4,    &
                  comm, status, ierr)

!--also set 0-layer to bogus values
txx(:, :, 0) = BOGUS
txy(:, :, 0) = BOGUS
txz(:, :, 0) = BOGUS
tyy(:, :, 0) = BOGUS
tyz(:, :, 0) = BOGUS
tzz(:, :, 0) = BOGUS !--tzz is updated in main (move to here)

txx(:, :, nz) = BOGUS
txy(:, :, nz) = BOGUS
tyy(:, :, nz) = BOGUS
tzz(:, :, nz) = BOGUS

! ---
end subroutine sgs_stag
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!  
