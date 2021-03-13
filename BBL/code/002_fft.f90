!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!    
module fft
!--------------------------------------------------------------------!
! Fast Fourier Transformation (FFT) for fftw3.3.X. 
! The module is based on both 1D and 2D FFT.
!--------------------------------------------------------------------!    
use types, only: rprec
use param, only: nxt, nyt, ldx, ldy, lhx, lhy, nxt2, nyt2, lhx2, ldx2,    &
                 lhy2, ldy2, lx_tot, ly_tot
implicit none
include 'fftw3.f'
save
public
! ---
real(rprec), dimension(:), allocatable :: kx, ky
real(rprec), dimension(:,:), allocatable :: kx_2d, ky_2d, k2
real(rprec), dimension(:,:), allocatable :: kx_2d_mpi, ky_2d_mpi, k2_mpi
integer, dimension(:,:), allocatable :: dealias
logical, dimension(:,:), allocatable :: notNyquist, notNyquist_mpi
! --- 2/3 of the Nyquist mode
integer :: kx_max_dealias, ky_max_dealias

integer(kind=8) :: xforw, yforw, forw, forw2, forw_p,   &
                   xback, yback, back, back2
integer(kind=8) :: plan_xb, plan_yb, plan_xf, plan_yf,  &       ! real to complex transform
                   plan_xcb, plan_ycb, plan_xcf, plan_ycf       ! complex to complex transform
integer(kind=8) :: plan_x2f, plan_x2b, plan_y2f, plan_y2b
complex(rprec), parameter :: eye = (0._rprec, 1._rprec)

contains

!    subroutine init_fft()
!!--------------------------------------------------------------------!
!!!! Initialize the 1D FFT environment
!!!! 01/11/2019 - First created by Bicheng Chen (chabby@ucla.edu)
!!!! Variables:
!!!! arr_nx, arr_ny, arr_nx2, arr_ny2 - arrays used to create plan                                          
!!--------------------------------------------------------------------!    
!    implicit none
!    ! ---
!    real(rprec), dimension(nxt) :: arr_nx
!    real(rprec), dimension(nyt) :: arr_ny
!    real(rprec), dimension(ldx) :: arr_ldx
!    real(rprec), dimension(ldy) :: arr_ldy
!    real(rprec), dimension(ldx, nyt) :: mat
!    real(rprec), dimension(nxt, nyt) :: mat_p
!    real(rprec), dimension(ldx2, nyt2) :: mat2
!    ! Create plans
!    ! forward
!    call dfftw_plan_dft_r2c_1d(xforw, nxt, arr_nx, arr_ldx,&
!                               FFTW_PATIENT, FFTW_UNALIGNED)
!    call dfftw_plan_dft_r2c_1d(yforw, nyt, arr_ny, arr_ldy,&
!                               FFTW_PATIENT, FFTW_UNALIGNED)
!    call dfftw_plan_dft_r2c_2d(forw, nxt, nyt, mat, mat,&
!                               FFTW_PATIENT, FFTW_UNALIGNED)
!    call dfftw_plan_dft_r2c_2d(forw2, nxt2, nyt2, mat2, mat2,&
!                               FFTW_PATIENT, FFTW_UNALIGNED)
!    call dfftw_plan_dft_r2c_2d(forw_p, nxt, nyt, mat_p, mat,&
!                               FFTW_PATIENT, FFTW_UNALIGNED)
!    ! backward
!    call dfftw_plan_dft_c2r_1d(xback, nxt, arr_ldx, arr_nx,&
!                               FFTW_PATIENT, FFTW_UNALIGNED)
!    call dfftw_plan_dft_c2r_1d(yback, nyt, arr_ldy, arr_ny,&
!                               FFTW_PATIENT, FFTW_UNALIGNED)
!    call dfftw_plan_dft_c2r_2d(back, nxt, nyt, mat, mat,&
!                               FFTW_PATIENT, FFTW_UNALIGNED)
!    call dfftw_plan_dft_c2r_2d(back2, nxt2, nyt2, mat2, mat2,&
!                               FFTW_PATIENT, FFTW_UNALIGNED)
!    ! Initilized wave number
!    call init_wavenumber()
!    ! ---
!    end subroutine init_fft
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!    
    subroutine init_fft_plan()
!--------------------------------------------------------------------!
!   create fft plans                                          
!--------------------------------------------------------------------!    
    implicit none
    ! ---
    real(rprec), dimension(nxt) :: arr_nxt
    real(rprec), dimension(nyt) :: arr_nyt
    real(rprec), dimension(nxt2) :: arr_nxt2
    real(rprec), dimension(nyt2) :: arr_nyt2
    complex(rprec), dimension(lhx) :: arr_lhx
    complex(rprec), dimension(lhy) :: arr_lhy
    complex(rprec), dimension(lhx2) :: arr_lhx2
    complex(rprec), dimension(lhy2) :: arr_lhy2 
    !complex(rprec), dimension(nxt) :: arr_x_in
    !complex(rprec), dimension(nxt) :: arr_x_out
    complex(rprec), dimension(nyt) :: arr_y_in
    complex(rprec), dimension(nyt) :: arr_y_out
    ! Create plans
    call dfftw_plan_dft_r2c_1d(plan_xf, nxt, arr_nxt, arr_lhx, FFTW_PATIENT, FFTW_UNALIGNED)
    call dfftw_plan_dft_c2r_1d(plan_xb, nxt, arr_lhx, arr_nxt, FFTW_PATIENT, FFTW_UNALIGNED)
    
    call dfftw_plan_dft_r2c_1d(plan_x2f, nxt2, arr_nxt2, arr_lhx2, FFTW_PATIENT, FFTW_UNALIGNED)
    call dfftw_plan_dft_c2r_1d(plan_x2b, nxt2, arr_lhx2, arr_nxt2, FFTW_PATIENT, FFTW_UNALIGNED)
    
    call dfftw_plan_dft_r2c_1d(plan_yf, nyt, arr_nyt, arr_lhy, FFTW_PATIENT, FFTW_UNALIGNED)
    call dfftw_plan_dft_c2r_1d(plan_yb, nyt, arr_lhy, arr_nyt, FFTW_PATIENT, FFTW_UNALIGNED)
    
    call dfftw_plan_dft_r2c_1d(plan_y2f, nyt2, arr_nyt2, arr_lhy2, FFTW_PATIENT, FFTW_UNALIGNED)
    call dfftw_plan_dft_c2r_1d(plan_y2b, nyt2, arr_lhy2, arr_nyt2, FFTW_PATIENT, FFTW_UNALIGNED)
    
    !call dfftw_plan_dft_1d(plan_xcf, nxt, arr_x_in, arr_x_out, FFTW_FORWARD,  FFTW_ESTIMATE)
    !call dfftw_plan_dft_1d(plan_xcb, nxt, arr_x_out, arr_x_in, FFTW_BACKWARD, FFTW_ESTIMATE)
    
    call dfftw_plan_dft_1d(plan_ycf, nyt, arr_y_in, arr_y_out, FFTW_FORWARD,  FFTW_ESTIMATE)
    call dfftw_plan_dft_1d(plan_ycb, nyt, arr_y_out, arr_y_in, FFTW_BACKWARD, FFTW_ESTIMATE)
    ! Initilized wave number
    call init_wavenumber()
    ! ---
    end subroutine init_fft_plan
!--------------------------------------------------------------------!
!                              
!--------------------------------------------------------------------!    
    subroutine init_wavenumber()
!--------------------------------------------------------------------!
!   Initialize the decomposed wavenumber meshes                                           
!--------------------------------------------------------------------!
    use param, only: pi, coordy, nxhnpy    
    implicit none
    ! ---
    integer :: ind, st, en
    integer :: jx, jy
    ! 1D FFT
    kx = [(real(ind-1, kind=rprec), ind=1, lhx)] * 2._rprec*pi/lx_tot
    ky = [(real(ind-1, kind=rprec), ind=1, lhy)] * 2._rprec*pi/ly_tot

    ! Remove Nyquist frequency if it exists
    if (mod(nxt, 2)==0) then
        kx(lhx) = 0._rprec
    end if
    if (mod(nyt, 2)==0) then
        ky(lhy) = 0._rprec
    end if

    ! 2D FFT
    do ind = 1, nyt
        kx_2d(:, ind) = kx
    end do
    
    ky_2d(1, :) = [(real(mod(ind-1+nyt/2, nyt)-nyt/2, kind=rprec), ind=1, nyt)]&
                    *2._rprec*pi/ly_tot
    if (lhx>=2) then 
        do ind=2, lhx
            ky_2d(ind, :) = ky_2d(1, :)
        end do
    end if

    ! Remove Nyquist frequency if it exists
    notNyquist = .true.
    if (mod(nxt, 2)==0) then
        ky_2d(lhx, :) = 0._rprec
        notNyquist(lhx, :) = .false.
    end if
    if (mod(nyt, 2)==0) then
        kx_2d(:, nyt/2+1) = 0._rprec
        ky_2d(:, nyt/2+1) = 0._rprec
        notNyquist(:, nyt/2+1) = .false.
    end if

    k2 = kx_2d*kx_2d + ky_2d*ky_2d
    
    st = coordy * (nxhnpy-1) + 1
    en = (coordy+1) * (nxhnpy-1) + 1
    notNyquist_mpi = notNyquist(st:en,:)
    kx_2d_mpi = kx_2d(st:en,:)
    ky_2d_mpi = ky_2d(st:en,:)
    k2_mpi = kx_2d_mpi*kx_2d_mpi + ky_2d_mpi*ky_2d_mpi
    ! --- 2/3 of the Nyquist mode
    kx_max_dealias = nxt / 3
    ky_max_dealias = nyt / 3
    
    dealias = 0
    do jy = 1, nyt
    do jx = 1, nxhnpy    
        if ( abs(kx_2d_mpi(jx, jy)) < kx_max_dealias .and.  &
             abs(ky_2d_mpi(jx, jy)) < ky_max_dealias ) then
            dealias = 1
        end if
    end do
    end do
    ! ---
    end subroutine init_wavenumber
! ---
end module fft

