!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine dealias_operation(var, var_d)
!--------------------------------------------------------------------!
!   2/3 rule                                             
!--------------------------------------------------------------------!
use types, only:rprec
use param, only:ldx, nxt, nyt, nxhnpy, nynpy, nz
use fft
implicit none
! ---
real(rprec), dimension(ldx, nynpy, 0:nz), intent(in)  :: var
real(rprec), dimension(ldx, nynpy, 0:nz), intent(out) :: var_d
! ---
real(rprec), dimension(ldx, nynpy, 0:nz) :: var_tmp
complex(rprec), dimension(nxhnpy, nyt) :: var_hat
real(rprec) :: const
integer :: jz
! ---
const = 1._rprec / (nxt*nyt)
var_tmp = const * var

do jz = 0, nz
    call fftn_mpi(var_tmp(:,:,jz), var_hat)
    var_hat(:,:) = var_hat(:,:) * dealias(:,:)
    call ifftn_mpi(var_hat, var_d(:,:,jz))
end do
! ---
end subroutine dealias_operation
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine dealias1(u, u_big)
! puts an array into the 'big' size
! doesn't dealias anything
! note: this sort of trashes u
  !!! 01/18/2018 - Update to FFTW3 by Bicheng Chen (chabby@ucla.edu)
  use types, only:rprec
  use param, only:ldx, ldx2, nxt, nyt, nz, nyt2
  use fft
  implicit none

  integer::jz
  real(rprec), dimension(ldx, nyt, 0:nz), intent(in):: u
  real(rprec), dimension(ldx2, nyt2, 0:nz), intent(inout)::u_big
  real(rprec) :: const
  ! ahh get rid of this
  real(rprec), dimension(ldx, nyt, 0:nz)::temp
  ! be careful using u after calling this subroutine!
  const = 1._rprec/(nxt*nyt)
  temp = const*u

  do jz = 0, nz !LV1
    ! still need to normalize
    call dfftw_execute_dft_r2c(forw, temp(:, :, jz), temp(:, :, jz))
    ! check padd syntax
    call padd(u_big(:, :, jz), temp(:, :, jz))
    call dfftw_execute_dft_c2r(back2, u_big(:, :, jz), u_big(:, :, jz))
  end do !LV1
end subroutine dealias1
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine dealias2(u, u_big)
! puts back into small array
  use types, only:rprec
  use param, only:ldx, ldx2, nxt, nyt, nz, nxt2, nyt2
  use fft
  implicit none

  integer::jz
  real(rprec), dimension(ldx, nyt, 0:nz), intent(in):: u
  real(rprec), dimension(ldx2, nyt2, 0:nz), intent(inout)::u_big
  real(rprec) :: const

  ! normalize
  const = 1._rprec/(nxt2*nyt2)
  u_big = const*u_big
  ! Loop through horizontal slices
  do jz = 0, nz
    ! perform forward FFT
    call dfftw_execute_dft_r2c(forw2, u_big(:, :, jz), u_big(:, :, jz))
    call unpadd(u(:, :, jz), u_big(:, :, jz))
    ! Back to physical space
    call dfftw_execute_dft_c2r(back, u(:, :, jz), u(:, :, jz))
  end do

! sc: do we need this?
!.....Making filter circular !!!!!!!!!!!
!        call filter_data(uu)
!.......................................
end subroutine dealias2
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
