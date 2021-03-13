!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------     
program main
!-----------------------------------------------------------------------
!  Purpose:
!    LES of canopy flow in atmospheric/oceanic boundary layer
!
!  Record of revisions:
!		date		programmer			description of change
!		====		==========			=====================
!     02/23/18		 Chao Yan				  
!     04/21/18		 Chao Yan             shallow ocean flow
!     05/12/18		 Chao Yan                2D parallel
!
!  Email: yanchao@ucla.edu
!-----------------------------------------------------------------------
use types, only:rprec
use param
use mpi
implicit none
! ---
integer :: ip, ipy, ipz, nproc, coords(2)
! --- Initialize computational environment (MPI and number of processors)
open(fid_param, file=fn_param, form='formatted', status='old')
read(fid_param, nml=mpi_param)
close(fid_param)
    
! --- mpi environment initialization
call mpi_init(ierr)
call mpi_comm_size(mpi_comm_world, nproc, ierr)
call mpi_comm_rank(mpi_comm_world, global_rank, ierr)

! --- check if run-time number of processors agrees with npz parameter
np = npy * npz
if (np /= nproc) then
    write(*, *) 'running processors .ne. the processors required'
    write(*, *) 'number of procs = ', nproc, 'value of processors = ', np
    call mpi_finalize(ierr)
    stop
end if

! --- set up a 2d cartesian topology
dims(1) = npy; dims(2) = npz
periodic(1) = .true.; periodic(2) = .false.

call mpi_cart_create(mpi_comm_world, 2, dims, periodic, .true., comm, ierr)
call mpi_comm_rank(comm, rank, ierr)
call mpi_cart_coords(comm, rank, 2, coords, ierr)
coordy = coords(1); coordz = coords(2)

call mpi_cart_shift(comm, 0, 1, south, north, ierr)
call mpi_cart_shift(comm, 1, 1, down, up, ierr)
call mpi_cart_sub(comm, (/.true., .false./), comm_ver, ierr)
call mpi_cart_sub(comm, (/.false., .true./), comm_col, ierr)

call mpi_comm_rank(comm_ver, ver_rank, ierr)
! --- rank->coord and coord->rank conversions
allocate(rank_of_coord(0:npy-1, 0:npz-1))
allocate(rank_ver(0:npy-1), rank_col(0:npz-1))
allocate(idy(0:np-1), idz(0:np-1))
allocate(idy_ver(0:npy-1), idz_col(0:npz-1))

do ipz = 0, npz - 1
do ipy = 0, npy - 1
    call mpi_cart_rank(comm, (/ipy, ipz/), rank_of_coord(ipy, ipz), ierr)
end do
end do

do ip = 0, np - 1
    call mpi_cart_coords(comm, ip, 2, coords, ierr)
    idy(ip) = coords(1)
    idz(ip) = coords(2)
end do

do ipz = 0, npz - 1
    call mpi_cart_rank(comm_col, (/ipz/), rank_col(ipz), ierr)
    call mpi_cart_coords(comm_col, ipz, 1, coord_col, ierr)
    idz_col(ipz) = coord_col
end do

do ipy = 0, npy - 1
    call mpi_cart_rank(comm_ver, (/ipy/), rank_ver(ipy), ierr)
    call mpi_cart_coords(comm_ver, ipy, 1, coord_ver, ierr)
    idy_ver(ipy) = coord_ver
end do

! --- set the MPI_RPREC variable
if (rprec == kind(1.e0)) then
    mpi_rprec = mpi_real
    mpi_cprec = mpi_complex
else if (rprec == kind(1.d0)) then
    mpi_rprec = mpi_double_precision
    mpi_cprec = mpi_double_complex
else
    write(*, *) 'error defining MPI_RPREC/MPI_CPREC'
    stop
end if

tt = 0.0_rprec
! --- I.C.s
call init()
! --- time Loop
do jt = 1, nsteps
    call flosol()
end do
! ---
if (rank == 0) write(*,*) "Max Iteration Reached, Job Finished."
call mpi_finalize(ierr)
! --- 		                                                      
end program main
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
