subroutine output_field(num)
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
use types, only: rprec
use param
use sim_param, only: u, v, w, p, theta, PCon
use stokes_drift, only: ust, vst
use sgsmodule, only: dissip
implicit none
! ---
integer, intent(in) :: num

real(rprec), dimension(nxt, nynpy, 0:nz) :: p_org
character(80) :: fname
integer(MPI_OFFSET_KIND) :: disp0, disp3d
integer :: jz, ipcon
! ---
do jz = 1, nz-1
    p_org(:, :,jz) = p(:, :,jz) - 0.5*(u(:, :, jz)*u(:, :, jz) +   &
                                       v(:, :, jz)*v(:, :, jz) +   &
                                       w(:, :, jz)*w(:, :, jz))    &
                                -     (u(:, :, jz)*ust(jz)+v(:, :, jz)*vst(jz))   &
                                - 0.5*(ust(jz)**2 + vst(jz)**2)
end do
! ---
disp3d = np * sizeof(real(u(1:nxt, 1:nynpy, 1:nz-1)))

write(fname, '(a,i8.8,a)') '../output/out/field_', num, '.out'
disp0 = 0
call write_file_mpi_3d_sp(real(u(1:nxt, 1:nynpy, 1:nz-1)), fname, disp0)
disp0 = disp0+disp3d
call write_file_mpi_3d_sp(real(v(1:nxt, 1:nynpy, 1:nz-1)), fname, disp0)
disp0 = disp0+disp3d
call write_file_mpi_3d_sp(real(w(1:nxt, 1:nynpy, 1:nz-1)), fname, disp0)
disp0 = disp0+disp3d

if (flag_opt_pre) then
    call write_file_mpi_3d_sp(real(p_org(1:nxt, 1:nynpy, 1:nz-1)), fname, disp0)
    disp0 = disp0+disp3d
end if

if (flag_opt_eps) then
    call write_file_mpi_3d_sp(real(dissip(1:nxt, 1:nynpy, 1:nz-1)), fname, disp0)
    disp0 = disp0+disp3d
end if

if ( pcon_flag ) then
    do ipcon = 1, npcon 
        call write_file_mpi_3d_sp(real(PCon(1:nxt, 1:nynpy, 1:nz-1, ipcon)), fname, disp0)
        disp0 = disp0+disp3d
    end do
end if

if ( theta_flag ) then
    call write_file_mpi_3d_sp(real(theta(1:nxt, 1:nynpy, 1:nz-1)), fname, disp0)
end if
! ---
end subroutine output_field
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------!
subroutine output_field_fraction(num)
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
use types, only: rprec
use param
use sim_param
use stokes_drift, only: ust, vst
use sgsmodule, only: Nu_t, dissip
use scalars_module, only:sgs_t3, sgs_PCon1, sgs_PCon3, res_PCon3
implicit none
! ---
integer, intent(in) :: num

real(rprec), dimension(nxt, nynpy, 0:nz) :: p_org
character(80) :: fname
integer(MPI_OFFSET_KIND) :: disp0, disp3d
integer :: jz, ipcon
! ---
disp3d = nyout * nzout * sizeof(u(1:nxout, 1, 1))

! --- flow field
write(fname, '(a,i8.8,a)') '../output/flow/flow_', num, '.out'
disp0 = 0
call write_file_mpi_3dFraction(u(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1), fname, disp0, bufsize)
disp0 = disp0+disp3d
call write_file_mpi_3dFraction(v(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1), fname, disp0, bufsize)
disp0 = disp0+disp3d
call write_file_mpi_3dFraction(w(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1), fname, disp0, bufsize)
disp0 = disp0+disp3d

if (flag_opt_pre) then
    do jz = 1, nz-1
        p_org(:, :,jz) = p(:, :,jz) - 0.5*(u(:, :, jz)*u(:, :, jz) +   &
                                           v(:, :, jz)*v(:, :, jz) +   &
                                           w(:, :, jz)*w(:, :, jz))    &
                                    -     (u(:, :, jz)*ust(jz)+v(:, :, jz)*vst(jz))   &
                                    - 0.5*(ust(jz)**2 + vst(jz)**2)
    end do
    ! ---
    call write_file_mpi_3dFraction(p_org(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1), fname, disp0, bufsize)    
    disp0 = disp0+disp3d
end if

if (flag_opt_eps) then
    call write_file_mpi_3dFraction(dissip(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1), fname, disp0, bufsize)    
    disp0 = disp0+disp3d
end if

if ( theta_flag ) then
    call write_file_mpi_3dFraction(theta(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1), fname, disp0, bufsize)    
    ! disp0 = disp0+disp3d
end if

! --- shear gradient
if (flag_opt_shear) then
    write(fname, '(a,i8.8,a)') '../output/shear/shear_', num, '.out'
    disp0 = 0
    call write_file_mpi_3dFraction(dudx(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1), fname, disp0, bufsize)
    disp0 = disp0+disp3d
    call write_file_mpi_3dFraction(dudy(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1), fname, disp0, bufsize)
    disp0 = disp0+disp3d
    call write_file_mpi_3dFraction(dvdx(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1), fname, disp0, bufsize)
    disp0 = disp0+disp3d
    call write_file_mpi_3dFraction(dvdy(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1), fname, disp0, bufsize)
    disp0 = disp0+disp3d
    call write_file_mpi_3dFraction(dwdx(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1), fname, disp0, bufsize)
    disp0 = disp0+disp3d
    call write_file_mpi_3dFraction(dwdy(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1), fname, disp0, bufsize)
    ! disp0 = disp0+disp3d
end if

! --- sgs field
if (flag_opt_sgs) then
    write(fname, '(a,i8.8,a)') '../output/sgs/sgs_', num, '.out'
    disp0 = 0
    call write_file_mpi_3dFraction(txx(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1), fname, disp0, bufsize)
    disp0 = disp0+disp3d
    call write_file_mpi_3dFraction(txy(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1), fname, disp0, bufsize)
    disp0 = disp0+disp3d
    call write_file_mpi_3dFraction(txz(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1), fname, disp0, bufsize)
    disp0 = disp0+disp3d
    call write_file_mpi_3dFraction(tyy(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1), fname, disp0, bufsize)
    disp0 = disp0+disp3d
    call write_file_mpi_3dFraction(tyz(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1), fname, disp0, bufsize)
    disp0 = disp0+disp3d
    call write_file_mpi_3dFraction(tzz(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1), fname, disp0, bufsize)
    disp0 = disp0+disp3d
    
    if ( theta_flag ) then
        call write_file_mpi_3dFraction(sgs_t3(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1), fname, disp0, bufsize)
        disp0 = disp0+disp3d
    end if   

    call write_file_mpi_3dFraction(Nu_t(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1), fname, disp0, bufsize)
    ! disp0 = disp0+disp3d
end if
! --- PCon field
if ( pcon_flag ) then
    write(fname, '(a,i8.8,a)') '../output/pcon/pcon_', num, '.out'
    disp0 = 0
    do ipcon = 1, npcon 
        call write_file_mpi_3dFraction(PCon(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1, ipcon), fname, disp0, bufsize)
        disp0 = disp0+disp3d
        call write_file_mpi_3dFraction(res_PCon3(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1, ipcon), fname, disp0, bufsize)
        disp0 = disp0+disp3d

        if (flag_opt_sgs) then
            call write_file_mpi_3dFraction(sgs_PCon1(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1, ipcon), fname, disp0, bufsize)
            disp0 = disp0+disp3d
            call write_file_mpi_3dFraction(sgs_PCon3(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1, ipcon), fname, disp0, bufsize)
            disp0 = disp0+disp3d
        end if 
    end do
end if
! ---
end subroutine output_field_fraction
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------!
subroutine output_field_fraction_sp(num)
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
use types, only: rprec
use param
use sim_param
use stokes_drift, only: ust, vst
use sgsmodule, only: Nu_t, dissip
use scalars_module, only:sgs_t3, sgs_PCon1, sgs_PCon3, res_PCon3
implicit none
! ---
integer, intent(in) :: num

real(rprec), dimension(nxt, nynpy, 0:nz) :: p_org
character(80) :: fname
integer(MPI_OFFSET_KIND) :: disp0, disp3d
integer :: jz, ipcon
! ---
disp3d = nyout * nzout * sizeof(real(u(1:nxout, 1, 1)))

! --- flow field
write(fname, '(a,i8.8,a)') '../output/flow/flow_', num, '.out'
disp0 = 0
call write_file_mpi_3dFraction_sp(real(u(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1)), fname, disp0, bufsize)
disp0 = disp0+disp3d
call write_file_mpi_3dFraction_sp(real(v(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1)), fname, disp0, bufsize)
disp0 = disp0+disp3d
call write_file_mpi_3dFraction_sp(real(w(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1)), fname, disp0, bufsize)
disp0 = disp0+disp3d

if (flag_opt_pre) then
    do jz = 1, nz-1
        p_org(:, :,jz) = p(:, :,jz) - 0.5*(u(:, :, jz)*u(:, :, jz) +   &
                                           v(:, :, jz)*v(:, :, jz) +   &
                                           w(:, :, jz)*w(:, :, jz))    &
                                    -     (u(:, :, jz)*ust(jz)+v(:, :, jz)*vst(jz))   &
                                    - 0.5*(ust(jz)**2 + vst(jz)**2)
    end do
    ! ---
    call write_file_mpi_3dFraction_sp(real(p_org(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1)), fname, disp0, bufsize)    
    disp0 = disp0+disp3d
end if

if (flag_opt_eps) then
    call write_file_mpi_3dFraction_sp(real(dissip(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1)), fname, disp0, bufsize)    
    disp0 = disp0+disp3d
end if

if ( theta_flag ) then
    call write_file_mpi_3dFraction_sp(real(theta(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1)), fname, disp0, bufsize)    
    ! disp0 = disp0+disp3d
end if

! --- shear gradient
if (flag_opt_shear) then
    write(fname, '(a,i8.8,a)') '../output/shear/shear_', num, '.out'
    disp0 = 0
    call write_file_mpi_3dFraction_sp(real(dudx(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1)), fname, disp0, bufsize)
    disp0 = disp0+disp3d
    call write_file_mpi_3dFraction_sp(real(dudy(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1)), fname, disp0, bufsize)
    disp0 = disp0+disp3d
    call write_file_mpi_3dFraction_sp(real(dvdx(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1)), fname, disp0, bufsize)
    disp0 = disp0+disp3d
    call write_file_mpi_3dFraction_sp(real(dvdy(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1)), fname, disp0, bufsize)
    disp0 = disp0+disp3d
    call write_file_mpi_3dFraction_sp(real(dwdx(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1)), fname, disp0, bufsize)
    disp0 = disp0+disp3d
    call write_file_mpi_3dFraction_sp(real(dwdy(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1)), fname, disp0, bufsize)
    ! disp0 = disp0+disp3d
end if

! --- sgs field
if (flag_opt_sgs) then
    write(fname, '(a,i8.8,a)') '../output/sgs/sgs_', num, '.out'
    disp0 = 0
    call write_file_mpi_3dFraction_sp(real(txx(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1)), fname, disp0, bufsize)
    disp0 = disp0+disp3d
    call write_file_mpi_3dFraction_sp(real(txy(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1)), fname, disp0, bufsize)
    disp0 = disp0+disp3d
    call write_file_mpi_3dFraction_sp(real(txz(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1)), fname, disp0, bufsize)
    disp0 = disp0+disp3d
    call write_file_mpi_3dFraction_sp(real(tyy(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1)), fname, disp0, bufsize)
    disp0 = disp0+disp3d
    call write_file_mpi_3dFraction_sp(real(tyz(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1)), fname, disp0, bufsize)
    disp0 = disp0+disp3d
    call write_file_mpi_3dFraction_sp(real(tzz(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1)), fname, disp0, bufsize)
    disp0 = disp0+disp3d
    
    if ( theta_flag ) then
        call write_file_mpi_3dFraction_sp(real(sgs_t3(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1)), fname, disp0, bufsize)
        disp0 = disp0+disp3d
    end if   

    call write_file_mpi_3dFraction_sp(real(Nu_t(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1)), fname, disp0, bufsize)
    ! disp0 = disp0+disp3d
end if
! --- PCon field
if ( pcon_flag ) then
    write(fname, '(a,i8.8,a)') '../output/pcon/pcon_', num, '.out'
    disp0 = 0
    do ipcon = 1, npcon 
        call write_file_mpi_3dFraction_sp(real(PCon(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1, ipcon)), fname, disp0, bufsize)
        disp0 = disp0+disp3d
        call write_file_mpi_3dFraction_sp(real(res_PCon3(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1, ipcon)), fname, disp0, bufsize)
        disp0 = disp0+disp3d

        if (flag_opt_sgs) then
            call write_file_mpi_3dFraction_sp(real(sgs_PCon1(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1, ipcon)), fname, disp0, bufsize)
            disp0 = disp0+disp3d
            call write_file_mpi_3dFraction_sp(real(sgs_PCon3(jx_opt_start:jx_opt_start-1+nxout, 1:nynpy, 1:nz-1, ipcon)), fname, disp0, bufsize)
            disp0 = disp0+disp3d
        end if 
    end do
end if
! ---
end subroutine output_field_fraction_sp
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------!
subroutine output_fringe(num)
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
use types, only: rprec
use param
use sim_param, only: u, v, w, theta, PCon   ! , p, RHSx, RHSy, RHSz
! use sgsmodule, only : Cs_opt2, F_LM, F_MM, F_QN, F_NN, F_KX, F_XX,F_KX2,F_XX2
! use scalars_module, only : RHS_T
implicit none
! ---
integer, intent(in) :: num

character(80) :: fname
integer(MPI_OFFSET_KIND) :: disp0, disp3d
integer :: jz, ipcon
! ---
disp3d = npy * npz * sizeof(u(jx_fringe+1:nxt, 1:nynpy, 1:nz-1))
! --- output velocity field at the outlet
write(fname, '(a,i8.8,a)') '../output/inflow/vel_tt', num, '.out'
disp0 = 0
call write_file_mpi_fringe(u(jx_fringe+1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
disp0 = disp0 + disp3d
call write_file_mpi_fringe(v(jx_fringe+1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
disp0 = disp0 + disp3d
call write_file_mpi_fringe(w(jx_fringe+1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
! disp0 = disp0 + disp2d
! call write_file_mpi_yz(RHSx(nxt, 1:nynpy, 1:nz-1), fname, disp0)
! disp0 = disp0 + disp2d
! call write_file_mpi_yz(RHSy(nxt, 1:nynpy, 1:nz-1), fname, disp0)
! disp0 = disp0 + disp2d
! call write_file_mpi_yz(RHSz(nxt, 1:nynpy, 1:nz-1), fname, disp0)
! disp0 = disp0 + disp2d
! call write_file_mpi_yz(Cs_opt2(nxt, 1:nynpy, 1:nz-1), fname, disp0)
! disp0 = disp0 + disp2d
! call write_file_mpi_yz(F_LM(nxt, 1:nynpy, 1:nz-1), fname, disp0)
! disp0 = disp0 + disp2d
! call write_file_mpi_yz(F_MM(nxt, 1:nynpy, 1:nz-1), fname, disp0)
! disp0 = disp0 + disp2d
! call write_file_mpi_yz(F_QN(nxt, 1:nynpy, 1:nz-1), fname, disp0)
! disp0 = disp0 + disp2d
! call write_file_mpi_yz(F_NN(nxt, 1:nynpy, 1:nz-1), fname, disp0)
! --- output the temperature field
if (theta_flag) then
    write(fname, '(a,i8.8,a)') '../output/inflow/temp_tt', num, '.out'
    disp0 = 0
    call write_file_mpi_fringe(theta(jx_fringe+1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
    ! disp0 = disp0 + disp2d
    ! call write_file_mpi_yz(RHS_T(nxt, 1:nynpy, 1:nz-1), fname, disp0)
end if

if (PCon_flag) then		
    write(fname, '(a,i8.8,a)') '../output/tracer/pcon_tt', num, '.out'		
    disp0 = 0		
    do ipcon = 1, npcon		
        call write_file_mpi_fringe(PCon(jx_fringe+1:nxt, 1:nynpy, 1:nz-1, ipcon), fname, disp0)		
        disp0 = disp0 + disp3d		
    end do		
end if
! ---
end subroutine output_fringe
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------!
!subroutine output_loop
!!--------------------------------------------------------------------!
!!BC revised by Bicheng Chen to add output sample
!!-BC add variable flagVSpl, ttVSpl, flagCSpl, ttCSpl
!!-BC add new subroutine checkVSpl and checkCSpl                                               
!!--------------------------------------------------------------------! 
!use param,only:path, output, c_count, theta_flag, theta_init_time, jt, jt_total,  &
!               jan_diurnal_run, flagVSpl, ttVSpl, flagCSpl,             &
!               ttCSpl, flag_srfV, tt_srfV, PCon_FLAG, PCon_init,        &
!               use_avgslice, base_time
!use scalars_module2,only:theta_slice, pcon_slice, budget_TKE_scalar
!use io, only:calc_mean, avgslice, checkpoint, checkVSpl, checkCSpl,     &
!             check_srfV, post_spec, io_spec, spec_write_start,          &
!             spec_write_end, spec_write_freqz
!implicit none
!! ---
!jt_total = jt_total + 1
!!call calc_mean()
!
!!cyan if (output) then
!!    if (mod(jt_total, base_time)==0) then
!!    !-- move all stuff into checkpoint (Bicheng Chen 06/26/2015)
!!        call checkpoint()
!!    end if
!!end if
!! --- BC added by Bicheng Chen for output the sample data
!if (flagVSpl) then
!    if (mod(jt_total, ttVSpl)==0) then
!        call checkVSpl()
!    end if
!end if
!
!if (flagCSpl) then
!    if (mod(jt_total, ttCSpl)==0) then
!        call checkCSpl()
!    end if
!end if
!!BC END
!! ---
!if (flag_srfV) then
!    if (mod(jt_total, tt_srfV)==0) then
!        call check_srfV()
!    end if
!end if
!
!! ---
!!cyanif ((io_spec) .and. (jt_total .gt. spec_write_start .and. jt_total .le. spec_write_end)) then
!!!  if ((io_spec) .and. mod(jt_total,spec_write_freqz)==0) then
!!    if (mod(jt_total,spec_write_freqz)==0) call post_spec(jt_total)
!!end if
!
!!cyan if (time_spec.gt.0) call timeseries_spec
!! ---
!end subroutine output_loop
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------