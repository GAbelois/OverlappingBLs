!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
module bottombc
!--------------------------------------------------------------------!
! Bicheng Chen edit some variables which need nxt and nyt to 
! have attribution allocatable                                              
!--------------------------------------------------------------------!
use types, only:rprec
use param, only:nxt, nyt, ldx, rag06, fieldsize, fieldsize_stripe,    &
                field_xs, field_xe, field_ys, field_ye, PCon_sfc
implicit none
save
public
! ---
integer, parameter :: num_patch = 1     ! numbr of patches, zo?=surface roughness for the patches types
integer, parameter :: ptypes = 2        ! number of patches types to be used, usually we use 2
real(kind=rprec), parameter :: zo1=0.01_rprec, zo2=0.125_rprec
real(kind=rprec), parameter :: do1=0.00_rprec, do2=0.75_rprec
integer, parameter :: square = 0
real(kind=rprec), parameter :: zo_out=.15_rprec, zo_square=.15_rprec
integer, parameter :: Nx_sq = 5, Ny_sq = 3
!real(kind=rprec),parameter::t_scale=300._rprec
real(kind=rprec), parameter :: theta_s1=288._rprec, theta_s2=320.15_rprec
real(kind=rprec), parameter :: B = 5.2_rprec           ! the constant in the log law of the smooth wall

real(kind=rprec),dimension(:,:), allocatable :: ustar_avg

! This is the common local ustar computed in obukhov and used in wallstress.f90
real(kind=rprec),dimension(:,:), allocatable :: zo, T_s, q_s
real(kind=rprec),dimension(:,:), allocatable :: z_os
!TS add for non-neutral case
real(kind=rprec),dimension(:,:), allocatable :: phi_m,psi_m,phi_h,psi_h
!VK The obukhov similarity functions are computed using obukhov(scalars_module.f90)
!VK for non-neutral scenario
integer,dimension(:,:), allocatable :: patch
integer, dimension(num_patch) :: patchnum

! Variables added for pollen Chamecki - 08/04/2006
real(kind=rprec),dimension(:,:), allocatable :: zo_PCon     ! Roughness for pollen concentration
real(kind=rprec),dimension(:,:), allocatable :: PCon_s      ! Surface concentration of pollen
real(kind=rprec),dimension(:,:), allocatable :: d0          ! Displacement height
! ---
real(kind=rprec):: T_s_min, T_s_max, z_o_min, z_o_max

contains

    subroutine patches()
    !VK This assigns momentum roughness, temperature and wetness
    !VK for the different patches
    !VK and fills the lookup tables patch and patchnum
    !VK This is called from routines patch_or_remote.f90 depending
    !VK whether to use remotely-sensed data or patches
    use param
    implicit none
    ! ---
    !integer :: i, j, patchtype, begini, endi, type1
    !integer,dimension(ptypes),intent(inout)::patchnum
    !integer,dimension(nxt,nyt),intent(inout)::patch

    ! sc: without this, some compilers may give total junk
    ! this was cause the u_avg to be WAY off, and so txz was WAY off, etc.
    ! patchnum = 0
    
    ! do j = 1, nynpy
    !     endi = 0
    !     type1 = 1
    !     do patchtype = 1, num_patch
    !         begini = endi + 1
    !         endi = (nxt*patchtype)/num_patch
    !         do i = begini, endi
    !             if ( type1 .eq. 1 ) then
    !                 zo(i,j) = zo1/z_i
    !                 if (theta_flag) then
    !                     T_s(i,j) = theta_s1/t_scale
    !                     q_s(i,j) = q_s1
    !                 end if
    !                 ! Added for pollen
    !                 if (PCon_flag) then
    !                     PCon_s(i,j) = PCon_sfc
    !                 end if
        
    !                 patch(i,j) = type1
    !             end if
        
    !             if ( type1 .eq. 2 ) then
    !                 zo(i,j) = zo2/z_i
        
    !                 if (theta_flag) then
    !                     T_s(i,j) = theta_s2/t_scale
    !                     q_s(i,j) = q_s2
    !                 end if
    !                 patch(i,j) = type1
    !             end if
        
    !             patchnum(type1) = patchnum(type1)+1
    !         end do
    
    !         if ( type1 .eq. 1 ) then
    !             type1 = 2
    !         else if ( type1 .eq. 2 ) then
    !             type1 = 1
    !         end if
    !     end do
    ! end do

    patch(:,:) = 1

    if (theta_flag) then
        T_s(:,:) = theta_s1/t_scale
        q_s(:,:) = q_s1
    end if
    ! Added for pollen
    if (PCon_flag) then
        PCon_s(:,:) = PCon_sfc
    end if

    ! Homogeneous
    zo(:,:) = zo1/z_i
    d0(:,:) = do1/z_i
    PCon_s(:,:) = PCon_sfc
    
    ! Set reference height for pollen (usually set in scalarr_module2.f90)
    zo_PCon(:,:) = zo(:,:) + d0(:,:)
        
    ! z_os already defined in scalars_module
    z_os(:, :) = (1._rprec/10._rprec)*zo(:, :)
    
    ! Calculate minimum and maximum values of T_s and zo for use in ic_scal()
    if (lbc .eq. 0) then
        T_s_min = minval(T_s)
        !print *, 'T_s_min', T_s_min
    end if

    ! if ( rank == 0 ) then
    !     open(1,file=path//'output/check_patches.out',status="unknown")
    !     do i = 1, nxt
    !         write(1,'(201f18.8)') i*dx, (zo(i,j),j=1,nyt)
    !     end do
    !     close(1)
    ! end if
    ! ---
    end subroutine patches
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
    subroutine remote_to_patch(T_s_in,T_or_z)
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
    implicit none
    ! ---
    integer :: T_or_z
    real(kind=rprec), dimension(nxt,nyt) :: T_s_in, dummy
    real(kind=rprec) :: sigma_multiplier, crap1, crap2, crap3
    real(kind=rprec) :: patch_cold, patch_hot
    ! ---
    sigma_multiplier = 1._rprec
    dummy = 0._rprec
    if (T_or_z == 1) then
        crap1=sum(T_s_in)/float(nxt*nyt)
    elseif (T_or_z == 2) then
        crap1=exp(sum(dlog(T_s_in))/float(nxt*nyt))
    else
        print *,'Wrong choice of T_or_z in remote_to_patch().. STOPPING'
        stop
    end if

    crap2 = sqrt(sum((T_s_in-crap1)**2)/float(nxt*nyt))
    !crap2=sqrt(sum((T_s_in-crap1)**2._rprec)/(real(nxt)*real(nyt)))

    patch_hot=crap1+sigma_multiplier*crap2;
    patch_cold=crap1-sigma_multiplier*crap2;

    print *,'mean, std, patch_hot, patch_cold = ',crap1,crap2,patch_hot,patch_cold
    !First do the patch business for temperature
    if ((patch_hot .lt. 0._rprec) .OR. (patch_cold .lt. 0._rprec)) then
        print *,'Hot & Cold patch calculation yields negative T_s'
        !print *,'mean, std, patch_hot, patch_cold = ',crap1,crap2,patch_hot,patch_cold
        print *,'Trying sigma_multiplier = 0.75'
        if (patch_cold < 0._rprec) then
            sigma_multiplier=0.75_rprec
            patch_cold=crap1-sigma_multiplier*crap2;
            patch_hot=crap1+sigma_multiplier*crap2;
            crap3=0.5_rprec*(patch_cold+patch_hot)
            print *,'NEW:mean, patch_hot, patch_cold = ',crap3,patch_hot,patch_cold
            if (patch_cold < 0._rprec) then
                print *,'sigma = 0.75 FAILED, Trying sigma_= 0.5'
                sigma_multiplier=0.5_rprec
                patch_cold=crap1-sigma_multiplier*crap2;
                patch_hot=crap1+sigma_multiplier*crap2;
                crap3=0.5_rprec*(patch_cold+patch_hot)
                print *,'NEW:mean, patch_hot, patch_cold = ',crap3,patch_hot,patch_cold
                if (patch_cold < 0._rprec) then
                    print *,'sigma = 0.5 FAILED, STOPPING NOW...'
                    print *,'This message is from the subroutine patch_to_remote in scalars_module2.f90.'
                end if
            end if
        end if
    end if

    ! Assign to patches
    if (T_or_z .eq. 1) then
        dummy(1:nxt/2,:)=patch_hot
        dummy(nxt/2+1:nxt,:)=patch_cold
        print *,'2 patch T_s field generated: patch_hot,patch_cold = ',patch_hot,patch_cold
    elseif (T_or_z .eq. 2) then
    ! Assign cold patch to the first patch and hot patch to the second one
    ! Note that this is exactly opposite to what we do for temperature as
    ! in REALITY, hot surfaces have low roughness and vice-versa
        dummy(1:nxt/2,:)=patch_cold
        dummy(nxt/2+1:nxt,:)=patch_hot
        print *,'2 patch roughnesss field generated: patch_smooth,patch_rough = ',patch_cold,patch_hot
    end if

    T_s_in = dummy
    ! ---
    end subroutine remote_to_patch
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
!subroutine avgpatch(u_avg)
!! computes the averaged value of a variable (at the wall) over a patch
!! and assigns it to an nxt X nyt array
!
!! sc:
!! note: this is inefficient: should calculate directly to u_avg
!! --no maybe its good
!!use types,only:rprec
!  use sim_param,only:u
!!use bottombc,only:ptypes
!  implicit none
!  real(kind=rprec),dimension(nxt,nyt),intent(inout)::u_avg
!!integer,dimension(ptypes),intent(in)::patchnum
!!integer,dimension(nxt,nyt),intent(in)::patch
!  integer::i,j,k
!  real(kind=rprec),dimension(ptypes)::temp
!!real(kind=rprec),dimension(ld,nyt,1),intent(in)::u
!
!  temp=0._rprec
!  do j=1,nyt
!    do i=1,nxt
!      do k=1,ptypes
!        if (patch(i,j).eq.k) then
!          temp(patch(i,j))=temp(patch(i,j))+u(i,j,1)/real(patchnum(patch(i,j)))
!        end if
!      end do
!    end do
!  end do
!
!  do j=1,nyt
!    do i=1,nxt
!      do k=1,ptypes
!        if (patch(i,j).eq.k) then
!          u_avg(i,j)=real(temp(k))
!        end if
!      end do
!    end do
!  end do
!end subroutine avgpatch

! ---
end module bottombc
