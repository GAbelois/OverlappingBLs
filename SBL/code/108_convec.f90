!-----------------------------------------------------------------------
! c = - (u X vort)
!...Compute rotational convective term in physical space  (X-mom eq.)
!...For staggered grid
!...v, dv, and du4 are on W nodes, rest on UVP nodes
!...(Fixed k=Nz plane on 1/24)
!...Corrected near wall on 4/14/96 {w(DZ/2)=0.5*w(DZ)}
! sc: better to have 0.25*w(dz), since really parabolic.
!--provides cx, cy, cz at 1:nz-1
!-----------------------------------------------------------------------
subroutine convec (cx,cy,cz)
!-----------------------------------------------------------------------
! Calculate the convection terms in NS equation
! 01/18/2018 - Update to FFTW3 by Bicheng Chen (chabby@ucla.edu)
!-----------------------------------------------------------------------
use types, only: rprec
use param
use sim_param, only: u, v, w, dudy, dudz, dvdx, dvdz, dwdx, dwdy
use fft
use intermediate, only: cross_x_big, cross_y_big, cross_z_big, vort1, vort2, vort3,	& 
                        u_big, v_big, w_big, vort1_big, vort2_big, vort3_big
implicit none
! ---
real(rprec), dimension(ldx, nynpy, 0:nz), intent(out) :: cx, cy, cz
! ---
!complex(rprec), dimension(nxhnpy, nyt) :: cx_hat, cy_hat, cz_hat
!real(rprec) :: const
integer :: jz

!...Recall dudz, and dvdz are on UVP node for k=1 only
!...So dv does not vary from arg2a to arg2b in 1st plane (k=1)

! sc: it would be nice to NOT have to loop through slices here
! Loop through horizontal slices

! const = 1._rprec / (nxt*nyt)

do jz = 0, nz
    call zero_padding(u(:,:,jz), u_big(:,:,jz))
    call zero_padding(v(:,:,jz), v_big(:,:,jz))
    call zero_padding(w(:,:,jz), w_big(:,:,jz))
end do

do jz = 1, nz
    !--if dudz, dvdz are on u-nodes for jz=1 and jz=nz, then we need a special
    !  definition of the vorticity in that case which also interpolates
    !  dwdx, dwdy to the u-node at jz=1 and jz=nz
    if ( coordz == 0 .and. jz == 1 ) then
        !--dwdy(jz=1) should be 0, so we could use this
        vort1(:, :, 1) = 0.5_rprec * (dwdy(:, :, 1)  +   &
                                      dwdy(:, :, 2)) -   &
                                      dvdz(:, :, 1)
        !--dwdx(jz=1) should be 0, so we could use this
        vort2(:, :, 1) =              dudz(:, :, 1) -    &
                         0.5_rprec * (dwdx(:, :, 1) +    &
                                      dwdx(:, :, 2))
    else if ( coordz == npz - 1 .and. jz == nz ) then
        vort1(:, :, nz) = 0.5_rprec * (dwdy(:, :, nz)    +     &
                                       dwdy(:, :, nz-1)) -     &
                                       dvdz(:, :, nz)
        vort2(:, :, nz) =              dudz(:, :, nz) -        &
                          0.5_rprec * (dwdx(:, :, nz) +        &
                                       dwdx(:, :, nz-1))
    else
        vort1(:, :, jz) = dwdy(:, :, jz) - dvdz(:, :, jz)
        vort2(:, :, jz) = dudz(:, :, jz) - dwdx(:, :, jz)
    end if
    vort3(:, :, jz) = dvdx(:, :, jz) - dudy(:, :, jz)
end do

do jz = 1, nz
    call zero_padding(vort1(:,:,jz), vort1_big(:,:,jz))
    call zero_padding(vort2(:,:,jz), vort2_big(:,:,jz))
    call zero_padding(vort3(:,:,jz), vort3_big(:,:,jz))
end do

! --- cross product of u X vort
!const = 1._rprec / (nxt2*nyt2)
call cross_product()
! cross_x_big = const * cross_x_big
! cross_y_big = const * cross_y_big
! cross_z_big = const * cross_z_big

do jz = 1, nz-1
    call zero_unpadding(cross_x_big(:,:,jz), cx(:,:,jz))
    call zero_unpadding(cross_y_big(:,:,jz), cy(:,:,jz))
    call zero_unpadding(cross_z_big(:,:,jz), cz(:,:,jz))
end do
! ---
cx(:, :, 0) = BOGUS
cy(:, :, 0) = BOGUS
cz(: ,:, 0) = BOGUS

!--top level is not valid
cx(:, :, nz) = BOGUS
cy(:, :, nz) = BOGUS
! cz(:, :, nz) = BOGUS
! ---
end subroutine convec
!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
subroutine cross_product()
!-----------------------------------------------------------------------
! calculate cross product of u and vort in physical space
!-----------------------------------------------------------------------
use param
use sim_param, only: u, v, w
use intermediate, only: cross_x_big, cross_y_big, cross_z_big, vort1_big,   &
                        u_big, v_big, w_big, vort2_big, vort3_big
use stokes_drift
implicit none
! ---
integer :: jz, jz_min, jz_max
! ---
if ( coordz == 0 ) then
    ! the cc's contain the normalization factor for the upcoming fft's
    cross_x_big(:,:,1) = (v_big(:,:,1)+vst(1))*(-vort3_big(:,:,1)) &
                       + 0.5_rprec*w_big(:,:,2)*vort2_big(:,:,2)
    !--vort2(jz=1) is located on uvp-node        ^  try with 1 (experimental)
    !--the 0.5 * w(:,:,2) is the interpolation of w to the first uvp node
    !  above the wall (could arguably be 0.25 * w(:,:,2))
    jz_min = 2
else
    jz_min = 1
end if

! if ( coordz == npz - 1 ) then
!     if (ocean_flag .and. stokes_flag) then
!         cross_x_big(:,:,nz-1) = (v_big(:,:,nz-1)+vst(nz-1))*(-vort3_big(:,:,nz-1))  &
!                       + 0.5_rprec*w_big(:,:,nz-1)*(vort2_big(:,:,nz-1))
!     else
!         cross_x_big(:,:,nz-1) = v_big(:,:,nz-1)*(-vort3_big(:,:,nz-1))          &
!                       + 0.5_rprec*w_big(:,:,nz-1)*(vort2_big(:,:,nz-1))
!     end if

!     jz_max = nz - 2
! else
!     jz_max = nz - 1
! end if

do jz = jz_min, nz - 1
    cross_x_big(:,:,jz) = (v_big(:,:,jz)+vst(jz))*(-vort3_big(:,:,jz))          &
                        + 0.5_rprec*(w_big(:,:,jz+1)*(vort2_big(:,:,jz+1))      &
                            + w_big(:,:,jz)  *(vort2_big(:,:,jz)))
end do
! ---
if ( coordz == 0 ) then
    ! the cc's contain the normalization factor for the upcoming fft's
    cross_y_big(:,:,1) = (u_big(:,:,1)+ust(1))*(vort3_big(:,:,1)) &
                       + 0.5_rprec*w_big(:,:,2)*(-vort1_big(:,:,2))
    !--vort1(jz=1) is uvp-node                    ^ try with 1 (experimental)
    !--the 0.5 * w(:,:,2) is the interpolation of w to the first uvp node
    !  above the wall (could arguably be 0.25 * w(:,:,2))

    jz_min = 2
else
    jz_min = 1
end if

! if ( coordz == npz - 1 ) then
!     if (ocean_flag .and. stokes_flag) then
!         cross_y_big(:,:,nz-1) = (u_big(:,:,nz-1)+ust(nz-1))*(vort3_big(:,:,nz-1))  &
!                          + 0.5_rprec*w_big(:,:,nz-1)*(-vort1_big(:,:,nz-1))
!     else
!         cross_y_big(:,:,nz-1) = u_big(:,:,nz-1)*(vort3_big(:,:,nz-1))             &
!                       + 0.5_rprec*   w_big(:,:,nz-1)*(-vort1_big(:,:,nz-1))
!     endif

!     jz_max = nz - 2
! else
!     jz_max = nz - 1
! end if

do jz = jz_min, nz - 1
    cross_y_big(:,:,jz) = (u_big(:,:,jz)+ust(jz))*(vort3_big(:,:,jz))&
                        + 0.5_rprec*(w_big(:,:,jz+1)*(-vort1_big(:,:,jz+1))&
                            +w_big(:,:,jz)*(-vort1_big(:,:,jz)))
end do
! ---
if ( coordz == 0 ) then
    ! There is no convective acceleration of w at wall or at top.
    cross_z_big(:,:,1) = 0._rprec
    jz_min = 2
else
    jz_min = 1
end if

if ( coordz == npz-1 ) then
    cross_z_big(:,:,nz) = 0._rprec
    jz_max = nz - 1
else
    jz_max = nz
end if

do jz = jz_min, jz_max
    cross_z_big(:,:,jz) = 0.5_rprec*(                  &
        (u_big(:,:,jz)+u_big(:,:,jz-1)+ust(jz)+ust(jz-1))*(-vort2_big(:,:,jz))+    &
        (v_big(:,:,jz)+v_big(:,:,jz-1)+vst(jz)+vst(jz-1))*(vort1_big(:,:,jz)) )
end do
! ---
end subroutine cross_product
!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
