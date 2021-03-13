!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
subroutine divstress_uv(divt, tx, ty, tz)
!-----------------------------------------------------------------------
! provides divt, jz=1:nz-1
!----------------------------------------------------------------------- 
use types, only: rprec
use param, only: ldx, nxt, nynpy, nz, BOGUS, coordz, npz
use derivatives, only: ddx, ddy, ddz_w
implicit none
! ---
real(kind=rprec), dimension(ldx, nynpy, 0:nz), intent(out) :: divt
real(rprec), dimension(ldx, nynpy, 0:nz), intent(in) :: tx, ty, tz
! sc: we shouldx be able to save some memory here!
! do this by writing a divergence subroutine--then do not store derivs
real(kind=rprec), dimension(ldx, nynpy, 0:nz) :: dtxdx, dtydy, dtzdz

! compute stress gradients
!--MPI: tx 1:nz-1 => dtxdx 1:nz-1
call ddx(dtxdx, tx) !--really shouldx replace with ddxy (save an fft)

!--MPI: ty 1:nz-1 => dtdy 1:nz-1
call ddy(dtydy, ty)

!--MPI: tz 1:nz => ddz_w limits dtzdz to 1:nz-1, except top process 1:nz
call ddz_w(dtzdz, tz)

!--MPI following comment only true at bottom process
! the following gives bad results...but it seems like i the
! issue shouldx be taken care of somewhere
! need to correct wall level, since tz will be on uv-node there
!      dtzdz(:,:,1) = (tz(:,:,2)-tz(:,:,1))/(0.5*dz)

!--only 1:nz-1 are valid
divt(:, :, 1:nz - 1) = dtxdx(:, :, 1:nz - 1) + dtydy(:, :, 1:nz - 1) + &
                       dtzdz(:, :, 1:nz - 1)

!--Set ldx-1, ldx to 0 (or couldx do BOGUS)
divt(ldx - 1:ldx, :, 1:nz - 1) = 0._rprec
divt(:, :, 0) = BOGUS
divt(:, :, nz) = BOGUS
! ---
end subroutine divstress_uv
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
subroutine divstress_w(divt, tx, ty, tz)
!-----------------------------------------------------------------------
!--provides divt 1:nz
!--nz may not be used anyway (BC is used instead)
!--MPI: provides 1:nz-1, except at top 1:nz
! 01/19/2018 - Replace do loop with vectorized operation by Bicheng Chen
! >>(chabby@ucla.edu)
!----------------------------------------------------------------------- 
use types, only: rprec
use param, only: ldx, nxt, nynpy, nz, npz, coordz, BOGUS
use derivatives, only: ddx, ddy, ddz_uv
implicit none
! ---
real(kind=rprec), dimension(ldx, nynpy, 0:nz), intent(out) :: divt
real(kind=rprec), dimension(ldx, nynpy, 0:nz), intent(in) :: tx, ty, tz
real(kind=rprec), dimension(ldx, nynpy, 0:nz) :: dtxdx, dtydy, dtzdz

! compute stress gradients
!--tx 1:nz => dtxdx 1:nz
call ddx(dtxdx, tx)
dtxdx(:, :, 0) = BOGUS

!--ty 1:nz => dtydy 1:nz
call ddy(dtydy, ty)
dtydy(:, :, 0) = BOGUS

!--tz 0:nz-1 (special case) => dtzdz 1:nz-1 (default), 2:nz-1 (bottom),
!                                    1:nz (top)
call ddz_uv(dtzdz, tz)
dtzdz(:, :, 0) = BOGUS

divt(:, :, 0) = BOGUS

if ( coordz == 0 ) then
    ! at wall we have to assume that dz(tzz)=0.0.  Any better ideas?
    ! in old version, it looks like some people tried some stuff with dwdz here
    ! but they were zeroed out, so they were left out of this version
    divt(1:nxt, 1:nynpy, 1) = dtxdx(1:nxt, 1:nynpy, 1) + dtydy(1:nxt, 1:nynpy, 1)
else
    divt(1:nxt, 1:nynpy, 1) = dtxdx(1:nxt, 1:nynpy, 1) + dtydy(1:nxt, 1:nynpy, 1)   &
                            + dtzdz(1:nxt, 1:nynpy, 1)
end if

divt(1:nxt, 1:nynpy, 2:nz-1) = dtxdx(1:nxt, 1:nynpy, 2:nz-1)    &
     + dtydy(1:nxt, 1:nynpy, 2:nz-1) + dtzdz(1:nxt, 1:nynpy, 2:nz-1)

!--set ldx-1, ldx to 0 (could maybe do BOGUS)
divt(ldx - 1:ldx, :, 1:nz - 1) = 0._rprec

if ( coordz == npz - 1 ) then
    divt(1:nxt, 1:nynpy, nz) = dtxdx(1:nxt, 1:nynpy, nz) + dtydy(1:nxt, 1:nynpy, nz)
    divt(ldx - 1:ldx, :, nz) = 0._rprec
else
    divt(:, :, nz) = BOGUS
end if
! ---
end subroutine divstress_w
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
