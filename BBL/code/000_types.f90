!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------!    
module types
!--------------------------------------------------------------------!
!   purpose: to declare data to share                                             
!--------------------------------------------------------------------!  
implicit none
save  
public
! --- data dictionary: declare constants
integer, parameter :: rprec = kind(1.d0)
!integer, parameter :: rprec = kind (1.e0)

! --- the maximum dimension and the default value for derived type read from namelist (Bicheng 05/25/2015)
integer, parameter :: nMax_dim = 256
real(rprec), parameter :: float_missing = -1234567890._rprec
integer, parameter :: int_missing = -9999

! --- Type of information for domain pararmeters (Bicheng 07/02/2015)
type type_infoDm
    integer :: nxt, nyt           ! the grid number of horizontal direction
    integer :: nz               ! the grid number of vertical direction in each processor
    integer :: nzt           ! the total grid number of vertical direction
    real(rprec) :: lx_tot, ly_tot, lz, lz_tot
    real(rprec) :: dx, dy, dz   ! the grid spacing in three directions
    real(rprec) :: z_i
    real(rprec), dimension(:), allocatable :: z_uvp     ! the z-coordinate of u, v, p, concentration
    real(rprec), dimension(:), allocatable :: z_w       ! the z-coordinate of w
    real(rprec), dimension(:), allocatable :: z_tau     ! the z-coordinate of Reynolds stress
    character(20), private :: unit_l = 'meter'
contains
    procedure :: getunit_l
end type type_infoDm

! --- Type of information for time parameters (Bicheng 06/03/2016)
type type_infoTime
    integer :: nsteps
    real(rprec) :: dt
    character(20), private :: unit_t = 'second'
    logical :: flag_restart = .false.
    integer :: nt_restart = int_missing
contains
    procedure :: getunit_t
end type type_infoTime

! --- Type of information for scalars(con), BC added by Bicheng 06/24/2015 
type type_infoCon
    integer, pointer :: n_con       ! the number of scalars(con) (represent the npcon in old version)
    real(rprec), dimension(:), pointer :: vel_settling      ! settling velocity or rise velocity
    real(rprec), dimension(:), pointer :: ratio_dens        ! the density ratio bewteen scalar and fluid
    real(rprec), dimension(:), pointer :: vol_spec          ! specific volumn of scalars
    real(rprec), dimension(:), pointer :: pcon_sfc_flux
    !real(rprec), dimension(nMax_dim) :: vel_settling
    !real(rprec), dimension(nMax_dim) :: ratio_dens = real_default
    !real(rprec), dimension(nMax_dim) :: vol_spec = real_default
    character(20), private :: unit_v = 'meter per second'
contains
    procedure :: getunit_v
end type type_infoCon

! Type of replicated domain
! n - total number of the replicated domain
type rp_dm
   integer :: row, col
   integer :: n
end type rp_dm

! Type of source release
! BC added by Bicheng for maximum number of single source
! BC release. Can be improved in ifort v14
! BC (no on this machine yet) 05/22/2015
! nmax_src - maximum number of single source release
! con_2D - for the planar source
integer, parameter :: nmax_src = 10
type src_rls
   integer :: n
   integer, dimension(nmax_src) :: ix, iy, iz
   real(rprec), dimension(nmax_src) :: xs, ys, zs
   integer, dimension(nmax_src) :: row, col
   real(rprec), dimension(nmax_src) :: rls
   real(rprec), dimension(:, :), allocatable :: con_2D
   real(rprec), dimension(:, :), allocatable :: dcon_2D
   !- For MPI only
   !- icpu - the number of the runing cpu
   !- iz_cpu - the local iz in the runing cpu
   integer, dimension(nmax_src) :: icpu_y
   integer, dimension(nmax_src) :: icpu_z
   integer, dimension(nmax_src) :: iy_cpu
   integer, dimension(nmax_src) :: iz_cpu
end type src_rls

! Type of information for ENDLESS method
! BC added by Bicheng for ENDLESS 06/24/2015
! flag_dynDm - flage for dynamicly ajusting of scalar(con) field domain
! nMax_dynDm - the maximum domain number for each scalar(con)
! con_thd - concentraiong threshold value to add or remove dynamic domain
! n_used - the number of used dynamic domain for each scalar(con)
type type_infoEndless
   logical :: flag_dynDm = .false.
   integer :: nMax_dynDm = int_missing
   real(rprec) :: con_thd = 1.0e-50_rprec
   integer, dimension(nMax_dim) :: n_used = int_missing
end type type_infoEndless

! Type to save parameters of for cellular flow
! BC added by Bicheng for ENDLESS 06/06/2016
type type_paramCell
   integer :: n_mod = int_missing
   real(rprec), dimension(:), allocatable :: u0
   real(rprec), dimension(:), allocatable :: u_adv, v_adv
   real(rprec) :: tt0
   real(rprec), dimension(:), allocatable :: depth
   real(rprec), dimension(:), allocatable :: diameter
   real(rprec), dimension(:), allocatable :: piPhi_x, piPhi_y
end type type_paramCell

! Type of mapping between domain of scalar(con) and its location
! BC added by Bicheng for ENDLESS and adaptive domain 06/24/2015
! flag_used - flag of used domain or not
! flag_src - flag of domain with source or not
! row, col - location(row and column) of the domain
! icon - the index of scalar(con)
! x0, xL, y0, yL - the nearby domains at each of four boundary
! no_nb - the index repsenting no nearby domain existing
type type_map
   logical, dimension(:), allocatable :: flag_used
   logical, dimension(:), allocatable :: flag_src
   integer, dimension(:), allocatable :: row, col
   integer, dimension(:), allocatable :: icon
   integer, dimension(:), allocatable :: x0, xL, y0, yL
   integer :: no_nb = -9999
end type type_map

! Type of linked list for usage of scalar(con) fields (pcon(:,:,:,:))
! BC added by Bicheng for ENDLESS and adaptive domain 06/24/2015
! list - store the index of used scalar fields and free scalar fields
!   the structure is list=[ind-of-used, ind-of-free]
! n_used - the number of used scalar(con) fields
! n_free - the number of free scalar(con) fields
!  n_used + n_free = nMax_dynDm in type_infoEndless
type type_usageCon
   integer :: n_used, n_free
   integer, dimension(:), allocatable :: list
end type type_usageCon

! Type of maximum concentration for different location of domain
! BC added by Bicheng for ENDLESS and adaptive domain 06/24/2015
! x0, xL, y0, yL - the maximum concentration at four horizontal sides of
!   domain
! dm - the maximum concentration of whold domain
type type_conMax
   real(rprec), dimension(:,:), allocatable :: x0, xL, y0, yL
   real(rprec), dimension(:,:), allocatable :: dm
end type type_conMax

contains

function getunit_l(info)
    class(type_infoDm), intent(in) :: info
    character(20) :: getunit_l
    getunit_l = info%unit_l
end function getunit_l

function getunit_t(info)
    class(type_infoTime), intent(in) :: info
    character(20) :: getunit_t
    getunit_t = info%unit_t
end function getunit_t

function getunit_v(info)
    class(type_infoCon), intent(in) :: info
    character(20) :: getunit_v
    getunit_v = info%unit_v
end function getunit_v
! ---
end module types
!--------------------------------------------------------------------!
!                                                
!--------------------------------------------------------------------!


    