!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------!    
subroutine screen_diagnose
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------!    
use types, only:rprec
use param
implicit none
! ---
integer, parameter :: wbase = 1 ! controls the frequency of screen diagnostics
! ---
if (modulo(jt, wbase) == 0) then
    call check_rmsdiv()             ! calculate the divergence of the flow field
    !call output_loop()
    !  call output_slice_loop()
end if
! ---    
end subroutine screen_diagnose
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------!  
subroutine check_rmsdiv()
!--------------------------------------------------------------------!
! actually, this is NOT the rms divergence of velocity, its like an
! l_1 norm or something.
!--------------------------------------------------------------------!  
use types, only: rprec
use param
use sim_param, only: du=>dudx, dv=>dvdy, dw=>dwdz
implicit none
! ---  
integer :: jx, jy, jz, jz_max
real(kind=rprec) :: rms, rms_global
! ---
if ( coordz == npz-1 ) then
    jz_max = nz-1
else
    !jz_max = nz
    jz_max = nz-1
end if

rms = 0._rprec
do jz = 1, jz_max
do jy = 1, nynpy
do jx = 1, nxt
    rms = rms + abs(du(jx,jy,jz)+dv(jx,jy,jz)+dw(jx,jy,jz))
end do
end do
end do
rms = rms / (nxt*nyt*(jz_max))  ! should be nyt here

call mpi_reduce(rms, rms_global, 1, MPI_RPREC, MPI_SUM, 0, comm, ierr)
if (rank == 0) then
    ! rms = rms_global / np
    open(13,file=path//'output/check_div.out',status="unknown",position="append")
    write(13, *) jt, rms
    close(13)
end if
! ---
end subroutine check_rmsdiv
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------! 
subroutine check_cfl(CFL, visc_stab)
!--------------------------------------------------------------------!
! This subroutine computes CFl and viscous stability and is called 
! every wbase timesetps from main.f90
!--------------------------------------------------------------------!
use types, only: rprec
use param
use sim_param, only: u, v, w
use sgsmodule, only: Nu_t
implicit none
! ---
real(kind=rprec), parameter :: CFL_max = 0.6_rprec     ! threshold value of CFL
real(kind=rprec), parameter :: CFL_min = 0.001_rprec
! ---
real(kind=rprec), intent(out) :: CFL, visc_stab
real(kind=rprec) :: u_max, v_max, w_max
real(kind=rprec) :: temp, nu_max
real(kind=rprec) :: cflx, cfly, cflz
real(kind=rprec) :: cflt, cflp, nu_max0

integer :: i, j, k
! ---
cflp = 0._rprec
u_max  = maxval(u(1:nxt, 1:nynpy, 1:nz - 1))
v_max  = maxval(v(1:nxt, 1:nynpy, 1:nz - 1))
w_max  = maxval(w(1:nxt, 1:nynpy, 1:nz - 1))
nu_max = maxval(Nu_t(1:nxt, 1:nynpy, 1:nz - 1))

cflx = u_max * dt / dx
cfly = v_max * dt / dy
cflz = w_max * dt / dz
nu_max0 = nu_max * dt / min(dx, dy, dz)**2

cflt = cflx + cfly + cflz

! do k = 1, nz-1
! do j = 1, nynpy
! do i = 1, nxt
!     cflx = abs(u(i,j,k) * dt / dx)
!     cfly = abs(v(i,j,k) * dt / dy)
!     cflz = abs(w(i,j,k) * dt / dz)
!     cflt = cflx + cfly + cflz
!     cflp = max(cflt, cflp)

!     if ( cflp .ge. CFL_max ) then
!         write(*,'(2i, 3e11.4, 3i)') coordy, coordz, cflx, cfly, cflz, i, j, k
!     end if
! !     ! if ( cflp .ge. CFL_max ) then
! !     !     write(*,'(i4.4, 6e11.4)') rank, cflx, cfly, cflz, u_max, v_max, w_max
! !     ! end if
! end do
! end do
! end do

call mpi_allreduce(cflt, CFL, 1, mpi_rprec, mpi_max, comm, ierr)
call mpi_allreduce(nu_max0, visc_stab, 1, mpi_rprec, mpi_max, comm, ierr)
! ---
if ( CFL .gt. CFL_max ) then
    if ( rank == 0 ) then
        write(*,*) "CFL = ", CFL, "IS TOO LARGE, PLEASE REDUCE THE TIME STEP"
    end if
    call mpi_finalize(ierr)
    stop
!else if ( CFL .lt. CFL_min ) then
!    if ( rank == 0 ) then
!        write(*,*) "CFL = ", CFL, "IS TOO SMALL, PLEASE INCREASE THE TIME STEP"
!    end if
!    call mpi_finalize(ierr)
!    stop
end if

!if ( rank == 0 ) then
!    open(13,file=path//'output/check_cfl.out',status="unknown",position="append")
!    write (13, *) cfl, visc_stab
!    close(13)
!end if
! ---
end subroutine check_cfl
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
!use scalars_module2,only:scalar_slice, pollen_slice, budget_TKE_scalar
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
subroutine check_field(num)
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
use types, only: rprec
use param
use sim_param, only: u, v, w, PCon
implicit none
! ---
integer, intent(in) :: num
real, dimension(nxout, nyt)   :: u_tot_xy, v_tot_xy, w_tot_xy
real, dimension(nxout, nzt-1) :: u_tot_xz, v_tot_xz, w_tot_xz
real, dimension(nxout, nyt, npcon) :: PCon_tot_xy
real, dimension(nxout, nzt-1, npcon) :: PCon_tot_xz

integer, parameter :: jy_index = 200
integer, parameter :: jz_index = 20
character(80) :: fname
integer :: jx, jy, jz, ipcon
real(rprec) :: xx, yy, zz
! ---
if (coordz == (jz_index-1)/(nz-1)) then
    jz = jz_index - coordz*(nz-1)
    do jx = jx_opt_start, jx_opt_start+nxout-1
        call collect_data_y_sp(u(jx,:,jz), u_tot_xy(jx-jx_opt_start+1,:))
        call collect_data_y_sp(v(jx,:,jz), v_tot_xy(jx-jx_opt_start+1,:))
        call collect_data_y_sp(w(jx,:,jz), w_tot_xy(jx-jx_opt_start+1,:))
        
        if (PCon_flag) then
            do ipcon = 1, npcon
                call collect_data_y_sp(PCon(jx,1:nynpy,jz,ipcon), PCon_tot_xy(jx-jx_opt_start+1,:,ipcon))
            end do
        end if
    end do

    if (coordy == 0) then
        write(fname, '(a,i8.8,a)') '../output/anim/flow/xy/vel_xy_check_', num, '.csv'
        open(9,file=fname)
        write(9, '(5a)') 'u', ',', 'v', ',', 'w'
        do jy = jy_opt_start, jy_opt_start+nyout-1
        do jx = 1, nxout
            write(9, '(e15.7, a, e15.7, a, e15.7)')      &
                u_tot_xy(jx, jy), ',', v_tot_xy(jx, jy), ',', w_tot_xy(jx, jy)
        end do
        end do
        close(9)

        ! write(fname, '(a,i8.8,a)') '../output/anim/flow/xy/vel_xy_check_', num, '.dat'
        ! open(9,file=fname)
        ! write(9,*) 'variables = "x","y","u","v","w"'
        ! write(9,*) 'zone i=', nxout, 'j=', nyout
        ! do jy = jy_opt_start, jy_opt_start+nyout-1
        ! do jx = 1, nxout
        !     xx = (jx-1) * lx_tot * z_i / nxt
        !     yy = (jy-1) * ly_tot * z_i / nyt
        !     write(9, '(5e15.7)') xx, yy, u_tot_xy(jx, jy), v_tot_xy(jx, jy), w_tot_xy(jx, jy)
        ! end do
        ! end do
        ! close(9)

        if (PCon_flag) then
            write(fname, '(a,i8.8,a)') '../output/anim/pcon/xy/tracer_xy_check_', num, '.csv'
            open(9,file=fname)
            write(9, '(5a)') 'tracer1', ',', 'tracer2', ',', 'tracer3'
            do jy = jy_opt_start, jy_opt_start+nyout-1
            do jx = 1, nxout
                write(9, '(e15.7, a, e15.7, a, e15.7)')      &
                PCon_tot_xy(jx, jy, 1), ',', PCon_tot_xy(jx, jy, 2), ',', PCon_tot_xy(jx, jy, 3)
            end do
            end do
            close(9)
        end if
    end if
end if


if (coordy == (jy_index-1)/nynpy) then
    jy = jy_index - coordy*nynpy
    do jx = jx_opt_start, jx_opt_start+nxout-1
        call collect_data_z_sp(u(jx,jy,1:nz-1), u_tot_xz(jx-jx_opt_start+1,:))
        call collect_data_z_sp(v(jx,jy,1:nz-1), v_tot_xz(jx-jx_opt_start+1,:))
        call collect_data_z_sp(w(jx,jy,1:nz-1), w_tot_xz(jx-jx_opt_start+1,:))
        
        if (PCon_flag) then
            do ipcon = 1, npcon
                call collect_data_z_sp(PCon(jx,jy,1:nz-1,ipcon), PCon_tot_xz(jx-jx_opt_start+1,:,ipcon))
            end do
        end if
    end do

    if (coordz == 0) then 
        write(fname, '(a,i8.8,a)') '../output/anim/flow/xz/vel_xz_check_', num, '.csv'
        open(9,file=fname)
        write(9, '(5a)') 'u', ',', 'v', ',', 'w'
        do jz = jz_opt_start, jz_opt_start+nzout-1
        do jx = 1, nxout
            write(9, '(e15.7, a, e15.7, a, e15.7)')      &
                u_tot_xz(jx, jz), ',', v_tot_xz(jx, jz), ',', w_tot_xz(jx, jz)
        end do
        end do 
        close(9)

        ! write(fname, '(a,i8.8,a)') '../output/anim/flow/xz/vel_xz_check_', num, '.dat'
        ! open(9,file=fname)
        ! write(9,*) 'variables = "x","z","u","v","w"'
        ! write(9,*) 'zone i=', nxout, 'j=', nzout
        ! do jz = jz_opt_start, jz_opt_start+nzout-1
        ! do jx = 1, nxout
        !     xx = (jx-1) * lx_tot * z_i / nxt
        !     zz = (jz-1) * z_i / (nzt-1)
        !     write(9, '(5e15.7)') xx, zz, u_tot_xz(jx, jz), v_tot_xz(jx, jz), w_tot_xz(jx, jz)
        ! end do
        ! end do
        ! close(9)
        
        if (PCon_flag) then
            write(fname, '(a,i8.8,a)') '../output/anim/pcon/xz/tracer_xz_check_', num, '.csv'
            open(9,file=fname)
            write(9, '(5a)') 'tracer1', ',', 'tracer2', ',', 'tracer3'
            do jz = jz_opt_start, jz_opt_start+nzout-1
            do jx = 1, nxout
                write(9, '(e15.7, a, e15.7, a, e15.7)')      &
                PCon_tot_xz(jx, jz, 1), ',', PCon_tot_xz(jx, jz, 2), ',', PCon_tot_xz(jx, jz, 3)
            end do
            end do
            close(9)
        end if
    end if
end if 
! ---
end subroutine check_field
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------!
subroutine check_field_bkp(num)
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
use types, only: rprec
use param
use sim_param, only: u, v, w, p, theta, PCon
use stokes_drift, only: ust, vst
implicit none
! ---
integer, intent(in) :: num
real, dimension(nxt, nyt, nzt) :: u_tot
real, dimension(nxt, nyt, nzt) :: v_tot
real, dimension(nxt, nyt, nzt) :: w_tot
! real, dimension(nxt, nyt, nzt) :: p_tot
! real, dimension(nxt, nyt, nzt) :: p_org
! real, dimension(nxt, nyt, nzt) :: theta_tot
real, dimension(nxt, nyt, nzt, npcon) :: PCon_tot

! real, dimension(nzt-1) :: ust_tot, vst_tot

character(80) :: fname
integer :: jx, jy, jz, ipcon
real(rprec) :: xx, yy, zz

! ---
call collect_data_3d_sp(u(1:nxt,:,:), u_tot)
call collect_data_3d_sp(v(1:nxt,:,:), v_tot)
call collect_data_3d_sp(w(1:nxt,:,:), w_tot)
! call collect_data_3d_sp(p(1:nxt,:,:), p_tot)
! call collect_data_3d_sp(theta(1:nxt,:,:), theta_tot)  

if (PCon_flag) then
    do ipcon = 1, npcon
        call collect_data_3d_sp(PCon(1:nxt, 1:nynpy, 0:nz, ipcon), PCon_tot(:,:,:,ipcon))
    end do
end if

! call collect_data_z_sp(ust(1:nz-1), ust_tot)
! call collect_data_z_sp(vst(1:nz-1), vst_tot)

if ( rank == 0 ) then
    ! do jz = 1, nzt-1
    !     p_org(:,:,jz) = p_tot(:,:,jz) - 0.5*(u_tot(:, :, jz)*u_tot(:, :, jz) +  &
    !                                          v_tot(:, :, jz)*v_tot(:, :, jz) +  &
    !                                          w_tot(:, :, jz)*w_tot(:, :, jz))   &
    !                                   -     (u_tot(:, :, jz)*ust_tot(jz)+       &
    !                                          v_tot(:, :, jz)*vst_tot(jz))       &
    !                                   - 0.5*(ust_tot(jz)*ust_tot(jz) + vst_tot(jz)*vst_tot(jz))     
    ! end do

    ! write(fname, '(a,i8.8,a)') '../output/z/vel_1D_check_', num, '.csv'
    ! open(9,file=fname)
    ! write(9, '(11a)') 'u', ',', 'v', ',', 'w',',', 'p', ',', 'theta', ',', 'PCon' 
    ! do jz = 1, nzt-1
    !    !zz = (jz-1) * z_i / (nzt-1)
    !     write(9, '(e15.7, a, e15.7, a, e15.7, a, e15.7, a, e15.7, a, e15.7)')                      & 
    !                        sum(u_tot(:, :, jz))/float(nxt*nyt)*u_scale,     & 
    !                   ',', sum(v_tot(:, :, jz))/float(nxt*nyt)*u_scale,     &
    !                   ',', sum(w_tot(:, :, jz))/float(nxt*nyt)*u_scale,     &
    !                   ',', sum(p_org(:, :, jz))/float(nxt*nyt),             &
    !                    ',', sum(theta_tot(:, :, jz))/float(nxt*nyt)*T_scale, &
    !                   ',', sum(PCon_tot(:, :, jz, 1))/float(nxt*nyt)
    ! end do
    ! close(9)

    write(fname, '(a,i8.8,a)') '../output/anim/flow/xz/vel_xz_check_', num, '.csv'
    open(9,file=fname)
    write(9, '(5a)') 'u', ',', 'v', ',', 'w'
    jy = nyt/2
    do jz = jz_opt_start, jz_opt_start+nzout-1
    do jx = jx_opt_start, jx_opt_start+nxout-1
        write(9, '(e15.7, a, e15.7, a, e15.7)')      &
            u_tot(jx, jy, jz), ',', v_tot(jx, jy, jz), ',', w_tot(jx, jy, jz) 
    end do
    end do
    close(9)

    write(fname, '(a,i8.8,a)') '../output/anim/flow/xy/vel_xy_check_', num, '.csv'
    open(9,file=fname)
    write(9, '(5a)') 'u', ',', 'v', ',', 'w'
    jz = 20
    do jy = jy_opt_start, jy_opt_start+nyout-1
    do jx = jx_opt_start, jx_opt_start+nxout-1
        write(9, '(e15.7, a, e15.7, a, e15.7)')      &
            u_tot(jx, jy, jz), ',', v_tot(jx, jy, jz), ',', w_tot(jx, jy, jz) 
    end do
    end do
    close(9)

    if (PCon_flag) then
        write(fname, '(a,i8.8,a)') '../output/anim/pcon/xz/tracer_xz_check_', num, '.csv'
        open(9,file=fname)
        write(9, '(5a)') 'tracer1', ',', 'tracer2', ',', 'tracer3'
        jy = nyt/2
        do jz = jz_opt_start, jz_opt_start+nzout-1
        do jx = jx_opt_start, jx_opt_start+nxout-1
            write(9, '(e15.7, a, e15.7, a, e15.7)')      &
                PCon_tot(jx, jy, jz, 1), ',', PCon_tot(jx, jy, jz, 2), ',', PCon_tot(jx, jy, jz, 3) 
        end do
        end do
        close(9)
    end if
    
    ! write(fname, '(a,i8.8,a)') '../output/xy/vel_xy_s_check_', num, '.csv'
    ! open(9,file=fname)
    ! write(9, '(7a)') 'u', ',', 'v', ',', 'w', ',', 'theta'
    ! jz = 2
    ! do jy = 1, nyt
    ! do jx = 1, nxt
    !     write(9, '(e15.7, a, e15.7, a, e15.7, a, e15.7)')       &
    !         u_tot(jx, jy, jz)*u_scale, ',', v_tot(jx, jy, jz)*u_scale, ',',     &
    !         w_tot(jx, jy, jz)*u_scale, ',', theta_tot(jx, jy, jz)*T_scale 
    ! end do
    ! end do
    ! close(9)

    ! write(fname, '(a,i8.8,a)') '../output/xy/vel_xy_check_', num, '.csv'
    ! open(9,file=fname)
    ! write(9, '(9a)') 'u', ',', 'v', ',', 'w', ',', 'theta'  !, ',', 'PCon'
    ! jz = 2
    ! do jy = 1, nyt
    ! do jx = 1, nxt
    !     write(9, '(e15.7, a, e15.7, a, e15.7, a, e15.7)')       &
    !         u_tot(jx, jy, jz)*u_scale, ',', v_tot(jx, jy, jz)*u_scale, ',',     &
    !         w_tot(jx, jy, jz)*u_scale, ',', theta_tot(jx, jy, jz)*T_scale
    !         ! , ',', &PCon_tot(jx, jy, jz, 1) 
    ! end do
    ! end do
    ! close(9)

    ! write(fname, '(a,i8.8,a)') '../output/xy/vel_xy_m_check_', num, '.csv'
    ! open(9,file=fname)
    ! write(9, '(7a)') 'u', ',', 'v', ',', 'w', ',', 'theta'
    ! jz = 84
    ! do jy = 1, nyt
    ! do jx = 1, nxt
    !     write(9, '(e15.7, a, e15.7, a, e15.7, a, e15.7)')       &
    !         u_tot(jx, jy, jz)*u_scale, ',', v_tot(jx, jy, jz)*u_scale, ',',     &
    !         w_tot(jx, jy, jz)*u_scale, ',', theta_tot(jx, jy, jz)*T_scale
    ! end do
    ! end do
    ! close(9)

    ! write(fname, '(a,i8.8,a)') '../output/xy/vel_xy_b_check_', num, '.csv'
    ! open(9,file=fname)
    ! write(9, '(7a)') 'u', ',', 'v', ',', 'w', ',', 'theta'
    ! jz = 114
    ! do jy = 1, nyt
    ! do jx = 1, nxt
    !     write(9, '(e15.7, a, e15.7, a, e15.7, a, e15.7)')       &
    !         u_tot(jx, jy, jz)*u_scale, ',', v_tot(jx, jy, jz)*u_scale, ',',     &
    !         w_tot(jx, jy, jz)*u_scale, ',', theta_tot(jx, jy, jz)*T_scale 
    ! end do
    ! end do
    ! close(9)

    ! write(fname, '(a,i8.8,a)') '../output/out/field_', num, '.out'
    ! open(9, file=fname, form='unformatted')
    ! write(9) u_tot(:,:,1:nzt), v_tot(:,:,1:nzt), w_tot(:,:,1:nzt),  &
    !          p_org(:,:,1:nzt), theta_tot(:,:,1:nzt)
    ! close(9)

    ! write(fname, '(a,i8.8,a)') '../output/vel_3D_check_', num, '.csv'
    ! open(9,file=fname)
    ! write(9, '(9a)') 'u', ',', 'v', ',', 'w', ',', 'theta', ',', 'PCon'
    ! do jz = 1, nzt-1
    ! do jy = 1, nyt
    ! do jx = 1, nxt
    !     write(9, '(e15.7, a, e15.7, a, e15.7, a, e15.7, a, e15.7)')       &
    !         u_tot(jx, jy, jz), ',', v_tot(jx, jy, jz), ',',     &
    !         w_tot(jx, jy, jz), ',', theta_tot(jx, jy, jz), ',', &
    !         PCon_tot(jx, jy, jz, 1) 
    ! end do
    ! end do
    ! end do
    ! close(9)

    !  write(fname, '(a,i8.8,a)') '../output/vel_1D_check', num, '.dat'
    !  open(9,file=fname)
    !  write(9,*) 'variables = "z","u","v","w","p","theta"'
    !  write(9,*) 'zone i=', nzt-1, 'f=point'
    !  do jz = 1, nzt-1
    !      zz = (jz-1) * z_i / (nzt-1)
    !      write(9, '(4e15.7)') -zz, sum(u_tot(:, :, jz))/float(nxt*nyt)*u_scale,         &      
    !                                sum(v_tot(:, :, jz))/float(nxt*nyt)*u_scale,         &
    !                                sum(w_tot(:, :, jz))/float(nxt*nyt)*u_scale,         &
	! 			                   sum(p_tot(:, :, jz))/float(nxt*nyt),                 &
    !                                sum(theta_tot(:, :, jz))/float(nxt*nyt)*T_scale  
    !  end do
    !  close(9)
    
    !  write(fname, '(a,i8.8,a)') '../output/xz/vel_xz_check', num, '.dat'
    !  open(9,file=fname)
    !  write(9,*) 'variables = "x","z","u","v","w","theta"'
    !  write(9,*) 'zone i=', nxt, 'j=', nzt-1, 'f=point'
    !  jy = nyt/2
    !  do jz = 1, nzt-1
    !  do jx = 1, nxt
    !      xx = (jx-1) * lx_tot * z_i / nxt
    !      zz = (jz-1) * z_i / (nzt-1)
    !      write(9, '(6e15.7)') xx, -zz, u_tot(jx, jy, jz)*u_scale,     &
    !                           v_tot(jx, jy, jz)*u_scale, w_tot(jx, jy, jz)*u_scale,                     &
    !                           theta_tot(jx, jy, jz)*T_scale
    !  end do
    !  end do
    !  close(9)
    
    !  write(fname, '(a,i8.8,a)') '../output/vel_xy_check', num, '.dat'
    !  open(9,file=fname)
    !  write(9,*) 'variables = "x","y","u","v","w","theta"'
    !  write(9,*) 'zone i=', nxt, 'j=', nyt, 'f=point'
    !  jz = 2
    !  do jy = 1, nyt
    !  do jx = 1, nxt
    !      xx = (jx-1) * lx_tot * z_i / nxt
    !      yy = (jy-1) * ly_tot * z_i / nyt
    !      write(9, *) xx, yy, u_tot(jx, jy, jz)*u_scale,   &
    !                  v_tot(jx, jy, jz)*u_scale, w_tot(jx, jy, jz)*u_scale,                  &
    !                  theta_tot(jx, jy, jz)*T_scale 
    !  end do
    !  end do
    !  close(9)
    
    ! write(fname, '(a,i8.8,a)') '../output/vel_3D_check', num, '.dat'
    ! open(9,file=fname)
    ! write(9,*) 'variables = "x","y","z","u","v","w"'
    ! write(9,*) 'zone i=', nxt, 'j=', nyt, 'k=', nzt-1, 'f=point'
    ! do jz = 1, nzt-1
    ! do jy = 1, nyt
    ! do jx = 1, nxt
    !     xx = (jx-1) * lx_tot * z_i / nxt
    !     yy = (jy-1) * ly_tot * z_i / nyt
    !     zz = (jz-1) * z_i / (nzt-1)
    !     write(9, '(6e15.7)') xx, yy, zz, u_tot(jx, jy, jz), v_tot(jx, jy, jz), w_tot(jx, jy, jz) 
    ! end do
    ! end do
    ! end do
    ! close(9)
end if
! ---
end subroutine check_field_bkp
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------!
! subroutine check_betascal(num)
! !--------------------------------------------------------------------!
! !                                               
! !--------------------------------------------------------------------!
! use types, only: rprec
! use param
! use scalars_module, only: beta_scal
! implicit none
! ! ---
! integer, intent(in) :: num
! real, dimension(nxt, nyt, nzt) :: beta_scal_tot

! character(80) :: fname
! integer :: jx, jy, jz
! real(rprec) :: xx, yy, zz

! ! ---
! call collect_data_3d(beta_scal, beta_scal_tot)

! if ( rank == 0 ) then
!      write(fname, '(a,i8.8,a)') '../output/beta_1D_check', num, '.dat'
!      open(9,file=fname)
!      write(9,*) 'variables = "z","beta"'
!      write(9,*) 'zone i=', nzt-1, 'f=point'
!      do jz = 1, nzt-1
!          zz = (jz-1) * z_i / (nzt-1)
!          write(9, '(4e15.4)') -zz, sum(beta_scal_tot(:, :, jz))/float(nxt*nyt)  
!      end do
!      close(9)
    
!      write(fname, '(a,i8.8,a)') '../output/beta_xz_check', num, '.dat'
!      open(9,file=fname)
!      write(9,*) 'variables = "x","z","beta"'
!      write(9,*) 'zone i=', nxt, 'j=', nzt-1, 'f=point'
!      jy = nyt/2
!      do jz = 1, nzt-1
!      do jx = 1, nxt
!          xx = (jx-1) * lx_tot * z_i / nxt
!          zz = (jz-1) * z_i / (nzt-1)
!          write(9, '(5e15.4)') xx, -zz, beta_scal_tot(jx, jy, jz) 
!      end do
!      end do
!      close(9)
    
!      write(fname, '(a,i8.8,a)') '../output/beta_xy_check', num, '.dat'
!      open(9,file=fname)
!      write(9,*) 'variables = "x","y","beta"'
!      write(9,*) 'zone i=', nxt, 'j=', nyt, 'f=point'
!      jz = 4
!      do jy = 1, nyt
!      do jx = 1, nxt
!          xx = (jx-1) * lx_tot * z_i / nxt
!          yy = (jy-1) * ly_tot * z_i / nyt
!          write(9, *) xx, yy, beta_scal_tot(jx, jy, jz)
!      end do
!      end do
!      close(9)    
! end if
! ! ---
! end subroutine check_betascal
! !--------------------------------------------------------------------!
! !                                               
! !--------------------------------------------------------------------!
subroutine check_pre(num)
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
use types, only: rprec
use param
use sim_param, only: p, dpdx, dpdy
implicit none
! ---
integer, intent(in) :: num
real, dimension(nxt, nyt, 1:nzt) :: p_tot, dpdx_tot, dpdy_tot

character(80) :: fname
integer :: jx, jy, jz
real(rprec) :: xx, yy, zz

! ---
call collect_data_3d_sp(p(1:nxt,:,:), p_tot)
call collect_data_3d_sp(dpdx(1:nxt,:,:), dpdx_tot)
call collect_data_3d_sp(dpdy(1:nxt,:,:), dpdy_tot)
   
if ( rank == 0 ) then
!    write(fname, '(a,i8.8,a)') '../output/debug/pre_1D_check', num, '.dat'
!    open(9,file=fname)
!    write(9,*) 'variables = "z","p","dpdx","dpdy"'
!    write(9,*) 'zone i=', nzt-1, 'f=point'
!    do jz = 1, nzt-1
!        zz = (jz-1) * z_i / (nzt-1)
!        write(9, '(4e15.4)') zz, sum(p_tot(:, :, jz))/float(nxt*nyt),    &
!                                 sum(dpdx_tot(:, :, jz))/float(nxt*nyt),    &
!                                 sum(dpdy_tot(:, :, jz))/float(nxt*nyt)
!    end do
!    close(9)
   
   write(fname, '(a,i8.8,a)') '../output/debug/pre_xz_check', num, '.dat'
   open(9,file=fname)
   write(9,*) 'variables = "x","z","p","dpdx","dpdy"'
   write(9,*) 'zone i=', nxt, 'j=', nzt-1, 'f=point'
   jy = nyt/2
   do jz = 1, nzt-1
   do jx = 1, nxt
       xx = (jx-1) * lx_tot * z_i / nxt
       zz = (jz-1) * z_i / (nzt-1)
       write(9, '(3e15.4)') xx, zz, p_tot(jx, jy, jz), dpdx_tot(jx, jy, jz), dpdy_tot(jx, jy, jz)
   end do
   end do
   close(9)

   write(fname, '(a,i8.8,a)') '../output/debug/pre_xy_check', num, '.dat'
   open(9,file=fname)
   write(9,*) 'variables = "x","y","p","dpdx","dpdy"'
   write(9,*) 'zone i=', nxt, 'j=', nyt, 'f=point'
   jz = 4
   do jy = 1, nyt
   do jx = 1, nxt
       xx = (jx-1) * lx_tot * z_i / nxt
       yy = (jy-1) * ly_tot * z_i / nyt
       write(9, *) xx, yy, p_tot(jx, jy, jz), dpdx_tot(jx, jy, jz), dpdy_tot(jx, jy, jz) 
   end do
   end do
   close(9)
end if
! ---
end subroutine check_pre
!--------------------------------------------------------------------!
!                                            
!--------------------------------------------------------------------!
!subroutine check_pg(num)
!!-----------------------------------------------------------------------
!!   
!!-----------------------------------------------------------------------
!use types, only: rprec
!use param
!use sim_param, only: dpdx, dpdy, dpdz
!implicit none
!! ---
!integer, intent(in) :: num
!real(rprec), dimension(nxt, nyt, 1:nzt) :: dpdx_tot
!real(rprec), dimension(nxt, nyt, 1:nzt) :: dpdy_tot
!real(rprec), dimension(nxt, nyt, 1:nzt) :: dpdz_tot
!
!integer, dimension(npz) :: sendcounts, recvcounts, displs
!
!character(80) :: fname
!integer :: jx, jy, jz
!real(rprec) :: x, y, z
!
!! ---
!sendcounts = size (dpdx(:,:,1:nz))
!recvcounts = size (dpdx(:,:,1:nz-1))
!displs = coord_of_rank * recvcounts
!
!call mpi_sendrecv(dpdx(1, 1, nz - 1), nxt*nyt, MPI_RPREC, up,   1, &
!                  dpdx(1, 1, 0),      nxt*nyt, MPI_RPREC, down, 1, &
!                  comm, status, ierr)
!call mpi_sendrecv(dpdx(1, 1, 1),  nxt*nyt, MPI_RPREC, down, 2, &
!                  dpdx(1, 1, nz), nxt*nyt, MPI_RPREC, up,   2, &
!                  comm, status, ierr)
!call mpi_sendrecv(dpdx(1, 1, nz - 1), nxt*nyt, MPI_RPREC, up,   3, &
!                  dpdx(1, 1, 0),      nxt*nyt, MPI_RPREC, down, 3, &
!                  comm, status, ierr)
!call mpi_sendrecv(dpdx(1, 1, 1),  nxt*nyt, MPI_RPREC, down, 4, &
!                  dpdx(1, 1, nz), nxt*nyt, MPI_RPREC, up,   4, &
!                  comm, status, ierr)
!call mpi_sendrecv(dpdx(1, 1, nz - 1), nxt*nyt, MPI_RPREC, up,   5, &
!                  dpdx(1, 1, 0),      nxt*nyt, MPI_RPREC, down, 5, &
!                  comm, status, ierr)
!call mpi_sendrecv(dpdx(1, 1, 1),  nxt*nyt, MPI_RPREC, down, 6, &
!                  dpdx(1, 1, nz), nxt*nyt, MPI_RPREC, up,   6, &
!                  comm, status, ierr)
!
!call mpi_gatherv(dpdx(1,1,1), sendcounts, MPI_RPREC,           &
!                 dpdx_tot(1,1,1), sendcounts, displs,          &
!                 MPI_RPREC, rank_of_coord(0), comm, ierr)
!call mpi_gatherv(dpdy(1,1,1), sendcounts, MPI_RPREC,           &
!                 dpdy_tot(1,1,1), sendcounts, displs,          &
!                 MPI_RPREC, rank_of_coord(0), comm, ierr)
!call mpi_gatherv(dpdz(1,1,1), sendcounts, MPI_RPREC,           &
!                 dpdz_tot(1,1,1), sendcounts, displs,          &
!                 MPI_RPREC, rank_of_coord(0), comm, ierr)                 
!    
!if ( rank == 0 ) then
!    write(fname, '(a,i8.8,a)') '../output/pg_1D_check', num, '.dat'
!    open(9,file=fname)
!    write(9,*) 'variables = "z","dpdx","dpdy","dpdz"'
!    write(9,*) 'zone i=', nzt-1, 'f=point'
!    do jz = 1, nzt-1
!        z = (jz-1) * z_i / (nzt-1)
!        write(9, '(4e15.4)') z, sum(dpdx_tot(:, :, jz))/float(nxt*nyt),   &
!                                sum(dpdy_tot(:, :, jz))/float(nxt*nyt),   &
!                                sum(dpdz_tot(:, :, jz))/float(nxt*nyt)
!    end do
!    close(9)
!    
!    write(fname, '(a,i8.8,a)') '../output/pg_xz_check', num, '.dat'
!    open(9,file=fname)
!    write(9,*) 'variables = "x","z","dpdx","dpdy","dpdz"'
!    write(9,*) 'zone i=', nxt, 'j=', nzt-1, 'f=point'
!    jy = nyt/2
!    do jz = 1, nzt-1
!    do jx = 1, nxt
!        x = (jx-1) * lx_tot * z_i / nxt
!        z = (jz-1) * z_i / (nzt-1)
!        write(9, '(3e15.4)') x, z, dpdx_tot(jx, jy, jz), dpdy_tot(jx, jy, jz),  &
!                                   dpdz_tot(jx, jy, jz)
!    end do
!    end do
!    close(9)
!
!    write(fname, '(a,i8.8,a)') '../output/pg_xy_check', num, '.dat'
!    open(9,file=fname)
!    write(9,*) 'variables = "x","y","dpdx","dpdy","dpdz"'
!    write(9,*) 'zone i=', nxt, 'j=', nyt, 'f=point'
!    jz = nz - 5
!    do jy = 1, nyt
!    do jx = 1, nxt
!        x = (jx-1) * lx_tot * z_i / nxt
!        y = (jy-1) * ly_tot * z_i / nyt
!        write(9, *) x, y, dpdx_tot(jx, jy, jz), dpdy_tot(jx, jy, jz),  &
!                          dpdz_tot(jx, jy, jz) 
!    end do
!    end do
!    close(9)
!end if
!! ---
!end subroutine check_pg
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------! 
subroutine check_RHS(num)
!--------------------------------------------------------------------!
! 
!--------------------------------------------------------------------!
use types, only:rprec
use param
use sim_param
implicit none
! ---
integer, intent(in) :: num

real, dimension(nxt, nyt, 1:nzt) :: RHSx_tot, RHSy_tot, RHSz_tot
! ---
character(80) :: fname
integer :: jx, jy, jz
real(rprec) :: xx, yy, zz
! ---
call collect_data_3d_sp(RHSx(1:nxt,:,:), RHSx_tot)
call collect_data_3d_sp(RHSy(1:nxt,:,:), RHSy_tot)
call collect_data_3d_sp(RHSz(1:nxt,:,:), RHSz_tot)

if ( rank == 0 ) then
    write(fname, '(a,i8.8,a)') '../output/rhs_1D_check_', num, '.dat'
    open(9,file=fname)
    write(9,*) 'variables = "z","rhsx","rhsy","rhsz"'
    do jz = 1, nzt-1
        zz = (jz-1) * z_i / (nzt-1)
        write(9, '(4e15.4)') -zz, sum(RHSx_tot(:, :, jz))/float(nxt*nyt),   & 
                                sum(RHSy_tot(:, :, jz))/float(nxt*nyt),   &
                                sum(RHSz_tot(:, :, jz))/float(nxt*nyt)
    end do
    close(9)

    write(fname, '(a,i8.8,a)') '../output/rhs_xz_check', num, '.dat'
    open(9,file=fname)
    write(9,*) 'variables = "x","z","rhsx","rhsy","rhsz"'
    write(9,*) 'zone i=', nxt, 'j=', nzt-1, 'f=point'
    jy = nyt/2
    do jz = 1, nzt-1
    do jx = 1, nxt
        xx = (jx-1) * lx_tot * z_i / nxt
        zz = (jz-1) * z_i / (nzt-1)
        write(9, '(5e15.4)') xx, -zz, RHSx_tot(jx, jy, jz), RHSy_tot(jx, jy, jz),     &
                                RHSz_tot(jx, jy, jz)
    end do
    end do
    close(9)

    write(fname, '(a,i8.8,a)') '../output/rhs_xy_check', num, '.dat'
    open(9,file=fname)
    write(9,*) 'variables = "x","y","rhsx","rhsy","rhsz"'
    write(9,*) 'zone i=', nxt, 'j=', nyt, 'f=point'
    jz = 4
    do jy = 1, nyt
    do jx = 1, nxt
        xx = (jx-1) * lx_tot * z_i / nxt
        yy = (jy-1) * ly_tot * z_i / nyt
        write(9, *) xx, yy, RHSx_tot(jx, jy, jz), RHSy_tot(jx, jy, jz),     &
                    RHSz_tot(jx, jy, jz)
    end do
    end do
    close(9)
end if
! ---
end subroutine check_RHS
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------! 
! subroutine check_RHST(num)
! !--------------------------------------------------------------------!
! ! 
! !--------------------------------------------------------------------!
! use types, only:rprec
! use param
! use scalars_module, only: RHS_T
! implicit none
! ! ---
! integer, intent(in) :: num

! real, dimension(nxt, nyt, 1:nzt) :: RHST_tot
! ! ---
! character(80) :: fname
! integer :: jx, jy, jz
! real(rprec) :: xx, yy, zz
! ! ---
! call collect_data_3d(RHS_T, RHST_tot)

! if ( rank == 0 ) then
!     write(fname, '(a,i8.8,a)') '../output/rhsT_1D_check_', num, '.dat'
!     open(9,file=fname)
!     write(9,*) 'variables = "z","rhsT"'
!     do jz = 1, nzt-1
!         zz = (jz-1) * z_i / (nzt-1)
!         write(9, '(4e15.4)') -zz, sum(RHST_tot(:, :, jz))/float(nxt*nyt)
!     end do
!     close(9)

!     write(fname, '(a,i8.8,a)') '../output/rhsT_xz_check', num, '.dat'
!     open(9,file=fname)
!     write(9,*) 'variables = "x","z","rhsT"'
!     write(9,*) 'zone i=', nxt, 'j=', nzt-1, 'f=point'
!     jy = nyt/2
!     do jz = 1, nzt-1
!     do jx = 1, nxt
!         xx = (jx-1) * lx_tot * z_i / nxt
!         zz = (jz-1) * z_i / (nzt-1)
!         write(9, '(5e15.4)') xx, -zz, RHST_tot(jx, jy, jz)
!     end do
!     end do
!     close(9)

!     write(fname, '(a,i8.8,a)') '../output/rhsT_xy_check', num, '.dat'
!     open(9,file=fname)
!     write(9,*) 'variables = "x","y","rhsT"'
!     write(9,*) 'zone i=', nxt, 'j=', nyt, 'f=point'
!     jz = 4
!     do jy = 1, nyt
!     do jx = 1, nxt
!         xx = (jx-1) * lx_tot * z_i / nxt
!         yy = (jy-1) * ly_tot * z_i / nyt
!         write(9, *) xx, yy, RHST_tot(jx, jy, jz)
!     end do
!     end do
!     close(9)
! end if
! ! ---
! end subroutine check_RHST
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------!
! subroutine check_divstress(num)
! !--------------------------------------------------------------------!
! ! 
! !--------------------------------------------------------------------!
! use types, only:rprec
! use param
! use sim_param
! implicit none
! ! ---
! integer, intent(in) :: num

! real, dimension(nxt, nyt, 1:nzt) :: divtx_tot, divty_tot, divtz_tot
! integer, dimension(npz) :: sendcounts, recvcounts, displs
! ! ---
! character(*), parameter :: fmt = '(e15.4, a, e15.4, a, e15.4)'
! character(80) :: fname
! integer :: jx, jy, jz
! real(rprec) :: xx, yy, zz
! ! ---
! call collect_data_3d_sp(divtx, divtx_tot)
! call collect_data_3d_sp(divty, divty_tot)
! call collect_data_3d_sp(divtz, divtz_tot)

! if ( rank == 0 ) then
!    write(fname, '(a,i8.8,a)') '../output/divstress_1D_check_', num, '.dat'
!    open(9,file=fname)
!    write(9,*) 'variables = "z","divtx","divty","divtz"'
!    do jz = 1, nzt-1
!        zz = (jz-1) * z_i / (nzt-1)
!        write(9, '(4e15.4)') zz, sum(divtx_tot(:, :, jz))/float(nxt*nyt),   & 
!                                sum(divty_tot(:, :, jz))/float(nxt*nyt),   &
!                                sum(divtz_tot(:, :, jz))/float(nxt*nyt)
!    end do
!    close(9)

!    write(fname, '(a,i8.8,a)') '../output/divstress_xz_check', num, '.dat'
!    open(9,file=fname)
!    write(9,*) 'variables = "x","z","divtx","divty","divtz"'
!    write(9,*) 'zone i=', nxt, 'j=', nzt-1, 'f=point'
!    jy = nyt/2
!    do jz = 1, nzt-1
!    do jx = 1, nxt
!        xx = (jx-1) * lx_tot * z_i / nxt
!        zz = (jz-1) * z_i / (nzt-1)
!        write(9, '(5e15.4)') xx, zz, divtx_tot(jx, jy, jz), divty_tot(jx, jy, jz),    &
!                                   divtz_tot(jx, jy, jz)
!    end do
!    end do
!    close(9)
   
!    write(fname, '(a,i8.8,a)') '../output/divstress_xy_check', num, '.dat'
!    open(9,file=fname)
!    write(9,*) 'variables = "x","y","txx","txy","txz","tyy","tyz","tzz"'
!    write(9,*) 'zone i=', nxt, 'j=', nyt, 'f=point'
!    jz = 4
!    do jy = 1, nyt
!    do jx = 1, nxt
!        xx = (jx-1) * lx_tot * z_i / nxt
!        yy = (jy-1) * ly_tot * z_i / nyt
!        write(9, '(5e15.4)') xx, yy, divtx_tot(jx, jy, jz), divty_tot(jx, jy, jz),    &
!                                   divtz_tot(jx, jy, jz)
!    end do
!    end do
!    close(9)
! end if
! ! ---
! end subroutine check_divstress
! !--------------------------------------------------------------------!
! !                                               
! !--------------------------------------------------------------------! 
! subroutine check_stress(num)
! !--------------------------------------------------------------------!
! ! 
! !--------------------------------------------------------------------!
! use types, only:rprec
! use param
! use sim_param
! implicit none
! ! ---
! integer, intent(in) :: num

! real, dimension(nxt, nyt, 1:nzt) :: txx_tot, txy_tot, txz_tot,     &
!                                     tyy_tot, tyz_tot, tzz_tot
! integer, dimension(npz) :: sendcounts, recvcounts, displs
! ! ---
! character(*), parameter :: fmt = '(e15.4, a, e15.4, a, e15.4)'
! character(80) :: fname
! integer :: jx, jy, jz
! real(rprec) :: xx, yy, zz
! ! ---
! call collect_data_3d_sp(txx, txx_tot)
! call collect_data_3d_sp(txy, txy_tot)
! call collect_data_3d_sp(txz, txz_tot)
! call collect_data_3d_sp(tyy, tyy_tot)
! call collect_data_3d_sp(tyz, tyz_tot)
! call collect_data_3d_sp(tzz, tzz_tot)

! if ( rank == 0 ) then
!    write(fname, '(a,i8.8,a)') '../output/stress_1D_check_', num, '.dat'
!    open(9,file=fname)
!    write(9,*) 'variables = "z","txx","txy","txz","tyy","tyz","tzz"'
!    do jz = 1, nzt-1
!        zz = (jz-1) * z_i / (nzt-1)
!        write(9, '(7e15.4)') zz, sum(txx_tot(:, :, jz))/float(nxt*nyt),   & 
!                                sum(txy_tot(:, :, jz))/float(nxt*nyt),   &
!                                sum(txz_tot(:, :, jz))/float(nxt*nyt),   &
!                                sum(tyy_tot(:, :, jz))/float(nxt*nyt),   & 
!                                sum(tyz_tot(:, :, jz))/float(nxt*nyt),   &
!                                sum(tzz_tot(:, :, jz))/float(nxt*nyt)
!    end do
!    close(9)

!    write(fname, '(a,i8.8,a)') '../output/stress_xz_check', num, '.dat'
!    open(9,file=fname)
!    write(9,*) 'variables = "x","z","txx","txy","txz","tyy","tyz","tzz"'
!    write(9,*) 'zone i=', nxt, 'j=', nzt-1, 'f=point'
!    jy = nyt/2
!    do jz = 1, nzt-1
!    do jx = 1, nxt
!        xx = (jx-1) * lx_tot * z_i / nxt
!        zz = (jz-1) * z_i / (nzt-1)
!        write(9, '(8e15.4)') xx, zz, txx_tot(jx, jy, jz), txy_tot(jx, jy, jz),    &
!                                   txz_tot(jx, jy, jz), tyy_tot(jx, jy, jz),    &
!                                   tyz_tot(jx, jy, jz), tzz_tot(jx, jy, jz)
!    end do
!    end do
!    close(9)
   
!    write(fname, '(a,i8.8,a)') '../output/stress_xy_check', num, '.dat'
!    open(9,file=fname)
!    write(9,*) 'variables = "x","y","txx","txy","txz","tyy","tyz","tzz"'
!    write(9,*) 'zone i=', nxt, 'j=', nyt, 'f=point'
!    jz = 4
!    do jy = 1, nyt
!    do jx = 1, nxt
!        xx = (jx-1) * lx_tot * z_i / nxt
!        yy = (jy-1) * ly_tot * z_i / nyt
!        write(9, '(8e15.4)') xx, yy, txx_tot(jx, jy, jz), txy_tot(jx, jy, jz),    &
!                                   txz_tot(jx, jy, jz), tyy_tot(jx, jy, jz),    &
!                                   tyz_tot(jx, jy, jz), tzz_tot(jx, jy, jz) 
!    end do
!    end do
!    close(9)
! end if
! ! ---
! end subroutine check_stress
! !--------------------------------------------------------------------!
! !                                               
! !--------------------------------------------------------------------! 
! subroutine check_derivative(num)
! !--------------------------------------------------------------------!
! ! 
! !--------------------------------------------------------------------!
! use types, only:rprec
! use param
! use sim_param
! implicit none
! ! ---
! integer, intent(in) :: num

! real, dimension(nxt, nyt, 1:nzt) :: dudz_tot, dvdz_tot, dwdz_tot
! integer, dimension(npz) :: sendcounts, recvcounts, displs
! ! ---
! character(*), parameter :: fmt = '(e15.4, a, e15.4, a, e15.4)'
! character(80) :: fname
! integer :: jx, jy, jz
! real(rprec) :: xx, yy, zz
! ! ---
! call collect_data_3d_sp(dudz, dudz_tot)
! call collect_data_3d_sp(dvdz, dvdz_tot)
! call collect_data_3d_sp(dwdz, dwdz_tot)

! if ( rank == 0 ) then
!    write(fname, '(a,i8.8,a)') '../output/dwdz_1D_check_', num, '.dat'
!    open(9,file=fname)
!    write(9,*) 'variables = "z","dudz","dvdz","dwdz"'
!    do jz = 1, nzt-1
!        zz = (jz-1) * z_i / (nzt-1)
!        write(9, '(2e15.4)') zz, sum(dudz_tot(:, :, jz))/float(nxt*nyt),  &
!                                sum(dvdz_tot(:, :, jz))/float(nxt*nyt),  &
!                                sum(dwdz_tot(:, :, jz))/float(nxt*nyt)
!    end do
!    close(9)

!    write(fname, '(a,i8.8,a)') '../output/dwdz_xz_check', num, '.dat'
!    open(9,file=fname)
!    write(9,*) 'variables = "x","z","dudz","dvdz","dwdz"'
!    write(9,*) 'zone i=', nxt, 'j=', nzt-1, 'f=point'
!    jy = nyt/2
!    do jz = 1, nzt-1
!    do jx = 1, nxt
!        xx = (jx-1) * lx_tot * z_i / nxt
!        zz = (jz-1) * z_i / (nzt-1)
!        write(9, '(3e15.4)') xx, zz, dudz_tot(jx, jy, jz), dvdz_tot(jx, jy, jz), dwdz_tot(jx, jy, jz)
!    end do
!    end do
!    close(9)
   
!    write(fname, '(a,i8.8,a)') '../output/dwdz_xy_check', num, '.dat'
!    open(9,file=fname)
!    write(9,*) 'variables = "x","y","dudz","dvdz","dwdz"'
!    write(9,*) 'zone i=', nxt, 'j=', nyt, 'f=point'
!    jz = 4
!    do jy = 1, nyt
!    do jx = 1, nxt
!        xx = (jx-1) * lx_tot * z_i / nxt
!        yy = (jy-1) * ly_tot * z_i / nyt
!        write(9, '(8e15.4)') xx, yy, dudz_tot(jx, jy, jz), dvdz_tot(jx, jy, jz), dwdz_tot(jx, jy, jz) 
!    end do
!    end do
!    close(9)
! end if
! ! ---
! end subroutine check_derivative
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------!
! subroutine check_divt(num)
! !--------------------------------------------------------------------!
! ! 
! !--------------------------------------------------------------------!
! use types, only:rprec
! use param
! use sim_param
! implicit none
! ! ---
! integer, intent(in) :: num
! ! ---
! character(*), parameter :: fmt = '(e15.4, a, e15.4, a, e15.4)'
! character(80) :: fname
! integer :: jx, jy, jz
! ! ---
! if ( rank == 0 ) then
!     write(fname, '(a,i8.8,a)') '../output/divt_check_', num, '.csv'
!     open(9,file=fname)
!     write(9, '(5a)') 'divtx', ',', 'divty', ',', 'divtxz'
!     do jz = 1, nz-1
!     do jy = 1, nyt
!     do jx = 1, nxt
!         write(9, fmt) divtx(jx, jy, jz), ',', divty(jx, jy, jz), ',', divtz(jx, jy, jz)
!     end do
!     end do
!     end do
!     close(9)
! end if
! ! ---
! end subroutine check_divt
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------! 
! subroutine check_coriolis(num)
! !--------------------------------------------------------------------!
! ! 
! !--------------------------------------------------------------------!
! use types, only:rprec
! use param
! use sim_param
! implicit none
! ! ---
! integer, intent(in) :: num
! ! ---
! character(*), parameter :: fmt = '(e15.4, a, e15.4, a, e15.4)'
! character(80) :: fname
! integer :: jx, jy, jz
! real(rprec) :: xx, zz
! ! ---
! if ( rank == 0 ) then
!     ! write(fname, '(a,i8.8,a)') '../output/coriolis_check_', num, '.csv'
!     ! open(9,file=fname)
!     ! write(9, '(5a)') 'coriol_v', ',', 'coriol_u'
!     ! do jz = 1, nz-1
!     ! do jy = 1, nyt
!     ! do jx = 1, nxt
!     !     write(9, fmt) coriol*v(jx, jy, jz), ',', coriol*u(jx, jy, jz)
!     ! end do
!     ! end do
!     ! end do
!     ! close(9)

!     write(fname, '(a,i8.8,a)') '../output/coriolis_check', num, '.dat'
!     open(9,file=fname)
!     write(9,*) 'variables = "x","z","fv","fu"'
!     write(9,*) 'zone i=', nxt, 'j=', nzt-1, 'f=point'
!     jy = nyt/2
!     do jz = 1, nzt-1
!     do jx = 1, nxt
!         xx = (jx-1) * lx_tot * z_i / nxt
!         zz = (jz-1) * z_i / (nzt-1)
!         write(9, '(4e15.4)') xx, zz, coriol*v(jx, jy, jz), coriol*u(jx, jy, jz) 
!     end do
!     end do
!     close(9)
! end if
! ! ---
! end subroutine check_coriolis
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------! 
! subroutine check_Nut(num)
! !--------------------------------------------------------------------!
! ! 
! !--------------------------------------------------------------------!
! use param
! use sgsmodule, only: Cs_opt2, Nu_t
! implicit none
! ! ---
! integer, intent(in) :: num
! real(rprec), dimension(nxt, nyt, 1:nzt) :: Nu_t_tot
! integer, dimension(npz) :: sendcounts, recvcounts, displs
! ! ---
! character(*), parameter :: fmt = '(e15.4, a, e15.4, a, e15.4, a, e15.4)'
! character(80) :: fname
! integer :: jx, jy, jz
! real(rprec) :: x, y, z
! ! ---
! sendcounts = size (Nu_t(:,:,1:nz))
! recvcounts = size (Nu_t(:,:,1:nz-1))
! displs = coord_of_rank * recvcounts

! call mpi_gatherv(Nu_t(1,1,1), sendcounts, MPI_RPREC,           &
!                  Nu_t_tot(1,1,1), sendcounts, displs,          &
!                  MPI_RPREC, rank_of_coord(0), comm, ierr)

! if ( rank == 0 ) then
!     write(fname, '(a,i8.8,a)') '../output/Nut_1D_check_', num, '.dat'
!     open(9,file=fname)
!     write(9,*) 'variables = "z","Nu_t"'
!     write(9,*) 'zone i=', nzt-1, 'f=point'
!     do jz = 1, nzt-1
!         z = (jz-1) * z_i / (nzt-1)
!         write(9, '(2e15.4)') z, sum(Nu_t_tot(:, :, jz))/float(nxt*nyt)
!     end do
!     close(9)
    
!     write(fname, '(a,i8.8,a)') '../output/Nut_xz_check', num, '.dat'
!     open(9,file=fname)
!     write(9,*) 'variables = "x","z","Nu_t"'
!     write(9,*) 'zone i=', nxt, 'j=', nzt-1, 'f=point'
!     jy = nyt/2
!     do jz = 1, nzt-1
!     do jx = 1, nxt
!         x = (jx-1) * lx_tot * z_i / nxt
!         z = (jz-1) * z_i / (nzt-1)
!         write(9, '(3e15.4)') x, z, Nu_t_tot(jx, jy, jz)
!     end do
!     end do
!     close(9)
! end if
! ! ---
! end subroutine check_Nut
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------! 
subroutine check_lad
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
use param
implicit none
! ---
real(kind=rprec), dimension(nxt, nyt, nzt) :: lad_uv_tot, lad_w_tot

character(80) :: fname
integer :: jx, jy, jz
real(rprec) :: xx, yy, zz
! ---
call collect_data3_3d(a_leafz_uv, lad_uv_tot)
call collect_data3_3d(a_leafz_w, lad_w_tot)

if ( rank == 0 ) then
    write(fname, '(a,i8.8,a)') '../output/lad_check.dat'
    open(9,file=fname)
    write(9,*) 'variables = "x","y","z","lad_uv","lad_w"'
    write(9,*) 'zone i=', nxt, 'j=', nyt, 'k=', nzt-1, 'f=point'
    do jz = 1, nzt-1
    do jy = 1, nyt
    do jx = 1, nxt
        xx = (jx-1) * lx_tot * z_i / nxt
        yy = (jy-1) * ly_tot * z_i / nyt
        zz = (jz-1) * z_i / (nzt-1)
        write(9, '(5e15.7)') xx, yy, zz, lad_uv_tot(jx, jy, jz), lad_w_tot(jx, jy, jz)
    end do
    end do
    end do
    close(9)
end if
! ---
end subroutine check_lad
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine check_sponge()
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
use types, only: rprec
use param
use topbc, only: sponge
implicit none
! ---
real(rprec), dimension(nzt-1) :: sponge_tot

character(80) :: fname
integer :: jz
real(rprec) :: zz
! ---

call collect_data_z(sponge(1:nz-1), sponge_tot)

if ( rank == 0 ) then
    write(fname, '(a,i8.8,a)') '../output/sponge.csv'
    open(9,file=fname)
    write(9, '(3a)') 'z', ',', 'sponge' 
    do jz = 1, nzt-1
        zz = (jz-1) * z_i / (nzt-1)
        write(9, '(e15.7, a, e15.7)')  zz, ',', sponge_tot(jz)
    end do
    close(9)
end if
! ---
end subroutine check_sponge
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------!