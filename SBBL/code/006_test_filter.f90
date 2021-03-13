!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!    
module test_filtermodule
!--------------------------------------------------------------------!
! Modified to fit namelist feature, use allocatable attribute 
! (Bicheng Chen 06/12/2016)
!--------------------------------------------------------------------!
use types, only:rprec
!TS Truely grid refinement test needs to keep the filter_size
!TS the same as that in coarse grid (Double grid size: filter_size=2. etc)
real(kind=rprec), parameter :: filter_size = 1._rprec
real(kind=rprec), dimension(:,:), allocatable :: G_test, G_test_test
! ---
end module test_filtermodule
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------     
subroutine test_filter(f, G_test)
!--------------------------------------------------------------------!
!                                             
!--------------------------------------------------------------------!  
use types, only: rprec
use param, only: ldx, nyt, nynpy, nxhnpy
use fft
implicit none
! --- 
real(kind=rprec), dimension(ldx, nynpy), intent(inout) :: f
real(kind=rprec), dimension(nxhnpy, nyt), intent(in)   :: G_test
! ---
complex(kind=rprec), dimension(nxhnpy, nyt) :: tmp
! ---
call fftn_mpi(f, tmp)
tmp = G_test * tmp      ! the normalization and Nyquist taken care of with G_test
call ifftn_mpi(tmp, f)
! ---
end subroutine test_filter
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
subroutine init_test_filter(alpha, G_test)
!-----------------------------------------------------------------------
! spectral cutoff filter at width alpha*delta
! note the normalization for FFT's is already in G! (see 1/(nxt*nyt))
!-----------------------------------------------------------------------
use types, only: rprec
use param, only: lhx, nxt, nyt, dx, dy, pi, ifilter, model, coordy, nxhnpy, npy
use fft
implicit none
! ---  
real(kind=rprec), intent(in) :: alpha            ! ratio of test filter to grid filter widths
real(kind=rprec), dimension(nxhnpy, nyt), intent(out) :: G_test

real(kind=rprec) :: delta, kc2
! ---
G_test = 1._rprec/(nxt*nyt)       ! normalization for the forward FFT
delta = alpha*sqrt(dx*dy)       ! "2d-delta", not full 3d one

if ( ifilter == 1 ) then        ! spectral cutoff filter
    kc2 = (pi/(delta))**2
    where (real(k2_mpi) >= kc2) G_test = 0._rprec
    
elseif ( ifilter == 2 ) then    ! Gaussian filter
    G_test = exp(-delta**2*k2_mpi/(4._rprec*6._rprec))*G_test
    
elseif ( ifilter == 3 ) then    ! Top-hat (Box) filter
    G_test = (sin(kx_2d_mpi*delta/2._rprec)*sin(ky_2d_mpi*delta/2._rprec)+1E-8)/  &
            (kx_2d_mpi*delta/2._rprec*ky_2d_mpi*delta/2._rprec+1E-8)*G_test
end if
! since our k2 has zero at Nyquist, we have to do this by hand
!if ( coordy == npy - 1 ) then
!!if ( coordy == 0 ) then
!    G_test(nxhnpy, :) = 0._rprec
!end if
!G_test(:, nyt/2+1) = 0._rprec
! ---
end subroutine init_test_filter
!!--------------------------------------------------------------------!
!! kills the oddball components in f, and calculates in plane derivatives
!!--supplies results for jz=$lbz:nz
!!--------------------------------------------------------------------!
subroutine filter_data(f_c)
!--------------------------------------------------------------------!
! Modified it to only filter data
!--------------------------------------------------------------------!
use types, only: rprec
use param, only: ldx, nxt, nyt, nynpy, nxhnpy, nz, coordy, npy, comm, ierr
use fft
implicit none
! ---
real(rprec), dimension(ldx, nynpy, 0:nz), intent(inout) :: f_c
! ---
complex(kind=rprec), dimension(nxhnpy, nyt) :: tmp
real(rprec) :: const
integer :: jz

! --- loop through horizontal slices
const = 1._rprec/(nxt*nyt)
do jz = 0, nz
    f_c(:, :, jz) = const * f_c(:, :, jz) !normalize
    call fftn_mpi(f_c(:, :, jz), tmp)
    call ifftn_mpi(tmp, f_c(:, :, jz))
end do
! ---
end subroutine filter_data
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------    
