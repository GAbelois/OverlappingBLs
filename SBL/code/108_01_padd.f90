!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine padd(u_big, u)
!--------------------------------------------------------------------!
! puts arrays into larger, zero-padded arrays
! automatically zeroes the oddballs                                              
!--------------------------------------------------------------------!
use types, only: rprec
use param, only: lhx, lhx2, nxt, nyt, nyt2, nxhnpy, coordy
implicit none
! --- note we're calling with 2D arrays
complex(kind=rprec), dimension(nxhnpy, nyt), intent(in)  :: u
complex(kind=rprec), dimension(lhx2, nyt2), intent(out) :: u_big
! ---
complex(kind=rprec), dimension(lhx, nyt) :: tmp
integer :: jx, jy
! make sure the big array is zeroed!
call collect_data_xh(u, tmp)

if ( coordy == 0 ) then
    u_big = 0._rprec
    do jy = 1, nyt/2
    do jx = 1, nxt/2
        ! skip the Nyquist frequency since it should be zero anyway
        u_big(jx, jy) = tmp(jx, jy)
    end do
    end do

    do jy = 1, nyt/2 - 1
    do jx = 1, nxt/2
        u_big(jx, jy+nyt2-nyt/2+1) = tmp(jx, jy+nyt/2+1)
    end do
    end do
end if
! ---
end subroutine padd
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine unpadd(cc, cc_big)
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
use types, only: rprec
use param, only: nxt, nyt, nyt2, lhx2, lhx
implicit none
! ---  
integer :: jx, jy
! note using 2d definitions here!
complex(kind=rprec), dimension(lhx2, nyt2), intent(in)  :: cc_big
complex(kind=rprec), dimension(lhx, nyt),   intent(out) :: cc

do jy = 1, nyt/2
do jx = 1, nxt/2
    cc(jx, jy) = cc_big(jx, jy)
end do
end do
! oddballs
cc(nxt/2 + 1, :) = 0._rprec
cc(:, nyt/2 + 1) = 0._rprec
  
do jy = 1, nyt/2 - 1
do jx = 1, nxt/2
    cc(jx, jy + nyt/2 + 1) = cc_big(jx, jy + nyt2 - nyt/2 + 1)
end do
end do
! ---
end subroutine unpadd
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine zero_padding(var, var_big)
!--------------------------------------------------------------------!
! puts arrays into larger, zero-padded arrays
! automatically zeroes the oddballs                                              
!--------------------------------------------------------------------!
use types, only: rprec
use param, only: ldx, nxt, nxt2, nynpy, ny2npy, lhx, lhx2, lhy, lhy2, nyt, nyt2, nx2npy
use fft 
implicit none
! --- note we're calling with 2D arrays
real(kind=rprec), dimension(ldx, nynpy), intent(in)  :: var
real(kind=rprec), dimension(nxt2, ny2npy), intent(out) :: var_big
! ---
complex(kind=rprec), dimension(lhx, nynpy) :: var_hat_x
complex(kind=rprec), dimension(nx2npy, lhy) :: var_hat_y
complex(kind=rprec), dimension(lhx2, nynpy) :: var_hat_x_big
complex(kind=rprec), dimension(nx2npy, lhy2) :: var_hat_y_big
real(kind=rprec), dimension(nxt2, nynpy) :: var_x_big
real(kind=rprec), dimension(nx2npy, nyt) :: var_y
real(kind=rprec), dimension(nx2npy, nyt2) :: var_y_big

real(rprec), dimension(nxt, nynpy) :: tmp
integer :: jx, jy
! --- make sure the big array is zeroed!
var_hat_x_big = 0._rprec
var_hat_y_big = 0._rprec

tmp(:, :) = var(1:nxt, :) / nxt

do jy = 1, nynpy
    call dfftw_execute_dft_r2c(plan_xf, tmp(:,jy), var_hat_x(:,jy))       ! 1D real to complex FFT in x direction
end do

do jx = 1, nxt/2
    var_hat_x_big(jx, :) = var_hat_x(jx, :)
end do

do jy = 1, nynpy
    call dfftw_execute_dft_c2r(plan_x2b, var_hat_x_big(:,jy), var_x_big(:,jy))
end do

call transpose_y_to_x2(var_x_big, var_y)

var_y = var_y / nyt
do jx = 1, nx2npy
    call dfftw_execute_dft_r2c(plan_yf, var_y(jx,:), var_hat_y(jx,:))
end do

do jy = 1, nyt/2
    var_hat_y_big(:, jy) = var_hat_y(:, jy)
end do

do jx = 1, nx2npy
    call dfftw_execute_dft_c2r(plan_y2b, var_hat_y_big(jx,:), var_y_big(jx,:)) 
end do

call transpose_x2_to_y(var_y_big, var_big)

! ---
end subroutine zero_padding
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine zero_unpadding(var_big, var)
!--------------------------------------------------------------------!
! puts arrays into larger, zero-padded arrays
! automatically zeroes the oddballs                                              
!--------------------------------------------------------------------!
use types, only: rprec
use param, only: ldx, nxt, nxt2, nxnpy, nynpy, ny2npy, lhx, lhx2, lhy, lhy2, nyt, nyt2
use fft 
implicit none
! --- note we're calling with 2D arrays
real(kind=rprec), dimension(nxt2, ny2npy), intent (in)  :: var_big
real(kind=rprec), dimension(ldx, nynpy), intent(out) :: var
! ---
complex(kind=rprec), dimension(lhx2, ny2npy) :: var_hat_x_big
complex(kind=rprec), dimension(nxnpy, lhy2) :: var_hat_y_big
complex(kind=rprec), dimension(lhx, ny2npy) :: var_hat_x
complex(kind=rprec), dimension(nxnpy, lhy) :: var_hat_y
real(kind=rprec), dimension(nxt, ny2npy) :: var_x
real(kind=rprec), dimension(nxnpy, nyt) :: var_y
real(kind=rprec), dimension(nxnpy, nyt2) :: var_y_big

real(rprec), dimension(nxt2, ny2npy) :: tmp
real(rprec), dimension(nxt, nynpy) :: var_tmp
integer :: jx, jy
! ---
tmp = var_big / nxt2

do jy = 1, ny2npy
    call dfftw_execute_dft_r2c(plan_x2f, tmp(:,jy), var_hat_x_big(:,jy))
end do

do jx = 1, nxt/2
    var_hat_x(jx, :) = var_hat_x_big(jx, :)
end do
var_hat_x(nxt/2+1, :) = 0._rprec

do jy = 1, ny2npy
    call dfftw_execute_dft_c2r(plan_xb, var_hat_x(:,jy), var_x(:,jy))
end do

call transpose_y2_to_x(var_x, var_y_big)

var_y_big = var_y_big / nyt2 
do jx = 1, nxnpy
    call dfftw_execute_dft_r2c(plan_y2f, var_y_big(jx,:), var_hat_y_big(jx,:))
end do

do jy = 1, nyt/2
    var_hat_y(:, jy) = var_hat_y_big(:, jy)
end do
var_hat_y(:, nyt/2+1) = 0._rprec

do jx = 1, nxnpy
    call dfftw_execute_dft_c2r(plan_yb, var_hat_y(jx,:), var_y(jx,:))
end do

call transpose_phy_x_to_y(var_y, var_tmp)

var(1:nxt, :) = var_tmp(1:nxt, :)
var(nxt+1:ldx, :) = var(1:2, :)
! ---
end subroutine zero_unpadding
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!