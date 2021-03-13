!
!    ################################################################################
!    ################################################################################
!    ##                                                                            ##
!    ##                            dynamic_sc_pasi.f90                             ##
!    ##                                                                            ##
!    ##                               Developed by                                 ##
!    ##                             Marcelo Chamecki                               ##
!    ##                                                                            ##
!    ##                        The Johns Hopkins University                        ##
!    ##                                                                            ##
!    ##                                05/25/2007                                  ##
!    ##                                                                            ##
!    ################################################################################
!    ################################################################################
!
!    ################################################################################
!
!    PURPOSE:  This routine is the plane averaged scale invariant dynamic model
!              for pollen concentration.
!
!    ################################################################################
!
! Modified to fit namelist feature, use allocatable attribution and
! >> remove the parameter attribution (Bicheng Chen 06/12/2016)

SUBROUTINE dynamic_sc_pasi(scalar,scalar_x_tmp,scalar_y_tmp,scalar_z_tmp,u_int_tmp,v_int_tmp,w_int_tmp)

! Modules
  USE types,only:rprec
  USE param,only:ldx,nxt,nynpy,nz,dx, dy,dz, PCon_scheme, periodicbcx, periodicbcy
  USE scalars_module, only: Cs2Sc
  USE sgsmodule,only:magS

! No implicit variables
  IMPLICIT NONE

! Main variables
  REAL(kind=rprec),DIMENSION(ldx,nynpy,nz)       :: scalar                      ! Scalar field
  REAL(kind=rprec),DIMENSION(ldx,nynpy+1,nz)     :: u_int_tmp,v_int_tmp,w_int_tmp      ! Interpolated velocity field
  REAL(kind=rprec),DIMENSION(ldx,nynpy+1,nz)     :: magSx_int,magSy_int,magSz_int! Interpolated |S|
  REAL(kind=rprec),DIMENSION(nxt+1,nynpy+1,nz)   :: scalar_x_tmp,scalar_y_tmp,scalar_z_tmp   ! Interpolated scalar
  REAL(kind=rprec),DIMENSION(nxt+1,nynpy+1,nz)   :: Kx,Ky,Kz                     ! "Leonard" fluxes
  REAL(kind=rprec),DIMENSION(nxt+1,nynpy+1,nz)   :: Xx,Xy,Xz                     ! "Leonard" modeled fluxes
  REAL(kind=rprec),DIMENSION(nxt,nynpy,nz)       :: KiXi,XiXi                    ! Vector dot products
  REAL(kind=rprec),DIMENSION(nz)             :: KX_avg,XX_avg       ! Plane averaged products

! Filtered fields
  REAL(kind=rprec),DIMENSION(nxt,nynpy,nz)       :: scalar_f                     ! Scalar field
  REAL(kind=rprec),DIMENSION(nxt+1,nynpy+1,nz)   :: scalar_x_f,scalar_y_f,scalar_z_f ! Interpolated scalar field
  REAL(kind=rprec),DIMENSION(ldx,nynpy+1,nz)     :: u_int_f,v_int_f,w_int_f      ! Interpolated velocity field
  REAL(kind=rprec),DIMENSION(nxt+1,nynpy+1,nz)   :: dscalar_dx,dscalar_dy,dscalar_dz ! Scalar derivatives

! Other variables
  REAL(kind=rprec),DIMENSION(ldx,nynpy+1,nz)     :: temp,temp_f                  ! Temporary to filter
  INTEGER                                    :: jz                           ! Counter


!
! BEGINNING CODE
!


  ! Check advection scheme
  IF (PCon_scheme/=3) THEN
    PRINT*,'Dynamic Cs2/Sc only available for SMART scheme...'
    STOP
  END IF


  ! Box filter velocity fields at test filter scale
  ! Note: these are the velocity fields already interpolated to surface locations
  temp=u_int_tmp
  CALL box_filter_3D(temp,temp_f,nxt+1,nynpy,nz-1)
  u_int_f=temp_f

  temp=v_int_tmp
  CALL box_filter_3D(temp,temp_f,nxt,nynpy+1,nz-1)
  v_int_f=temp_f

  temp=w_int_tmp
  CALL box_filter_3D(temp,temp_f,nxt,nynpy,nz-1)
  w_int_f=temp_f


  ! Box filter scalar fields at test filter scale
  ! Note: this is the scalar field already interpolated to surface locations
  temp=0._rprec
  temp(1:nxt+1,:,:)=scalar_x_tmp
  CALL box_filter_3D(temp,temp_f,nxt+1,nynpy,nz-1)
  scalar_x_f=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=scalar_y_tmp
  CALL box_filter_3D(temp,temp_f,nxt,nynpy+1,nz-1)
  scalar_y_f=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=scalar_z_tmp
  CALL box_filter_3D(temp,temp_f,nxt,nynpy,nz-1)
  scalar_z_f=temp_f(1:nxt+1,:,:)


  ! First use the vectors Kx, Ky and Kz to calculate the resolved fluxes
  Kx=u_int_tmp(1:nxt+1,1:nynpy+1,1:nz)*scalar_x_tmp
  Ky=v_int_tmp(1:nxt+1,1:nynpy+1,1:nz)*scalar_y_tmp
  Kz=w_int_tmp(1:nxt+1,1:nynpy+1,1:nz)*scalar_z_tmp


  ! Now filter the resolved fluxes
  temp=0._rprec
  temp(1:nxt+1,:,:)=Kx
  CALL box_filter_3D(temp,temp_f,nxt+1,nynpy,nz-1)
  Kx=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=Ky
  CALL box_filter_3D(temp,temp_f,nxt,nynpy+1,nz-1)
  Ky=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=Kz
  CALL box_filter_3D(temp,temp_f,nxt,nynpy,nz-1)
  Kz=temp_f(1:nxt+1,:,:)


  ! Subtract product of filtered fields to get "Leonard" fluxes
  Kx=Kx-u_int_f(1:nxt+1,1:nynpy+1,1:nz)*scalar_x_f
  Ky=Ky-v_int_f(1:nxt+1,1:nynpy+1,1:nz)*scalar_y_f
  Kz=Kz-w_int_f(1:nxt+1,1:nynpy+1,1:nz)*scalar_z_f


  ! Calculate derivatives
  dscalar_dx(2:nxt,1:nynpy,:)=(scalar(2:nxt,1:nynpy,:)-scalar(1:nxt-1,1:nynpy,:))/dx
  dscalar_dy(1:nxt,2:nynpy,:)=(scalar(1:nxt,2:nynpy,:)-scalar(1:nxt,1:nynpy-1,:))/dy
  dscalar_dz(1:nxt,1:nynpy,2:nz-1)=(scalar(1:nxt,1:nynpy,2:nz-1)-scalar(1:nxt,1:nynpy,1:nz-2))/dz


  ! |S| has to be interpolated (same linear interpolation as magS in the pollen_RHS_calc2 routine
  ! For x-diffusion interpolate in x and z directions
  magSx_int(2:nxt,1:nynpy,1)=(magS(1:nxt-1,1:nynpy,1)+magS(2:nxt,1:nynpy,1))/2._rprec
  magSx_int(1,1:nynpy,1)=(magS(nxt,1:nynpy,1)+magS(1,1:nynpy,1))/2._rprec

  magSx_int(2:nxt,1:nynpy,2:nz-1)=(magS(1:nxt-1,1:nynpy,2:nz-1)+magS(2:nxt,1:nynpy,2:nz-1)+ &
     magS(1:nxt-1,1:nynpy,3:nz)+magS(2:nxt,1:nynpy,3:nz))/4._rprec
  magSx_int(1,1:nynpy,2:nz-1)=(magS(nxt,1:nynpy,2:nz-1)+magS(1,1:nynpy,2:nz-1)+ &
     magS(nxt,1:nynpy,3:nz)+magS(1,1:nynpy,3:nz))/4._rprec
  magSx_int(nxt+1,:,:)=magSx_int(1,:,:)


  ! For y-diffusion interpolate in y and z directions
  magSy_int(1:nxt,2:nynpy,1)=(magS(1:nxt,1:nynpy-1,1)+magS(1:nxt,2:nynpy,1))/2._rprec
  magSy_int(1:nxt,1,1)=(magS(1:nxt,nynpy,1)+magS(1:nxt,1,1))/2._rprec

  magSy_int(1:nxt,2:nynpy,2:nz-1)=(magS(1:nxt,1:nynpy-1,2:nz-1)+magS(1:nxt,2:nynpy,2:nz-1)+ &
     magS(1:nxt,1:nynpy-1,3:nz)+magS(1:nxt,2:nynpy,3:nz))/4._rprec
  magSy_int(1:nxt,1,2:nz-1)=(magS(1:nxt,nynpy,2:nz-1)+magS(1:nxt,1,2:nz-1)+ &
     magS(1:nxt,nynpy,3:nz)+magS(1:nxt,1,3:nz))/4._rprec
  magSy_int(:,nynpy+1,:)=magSy_int(:,1,:)


  ! For z-diffusion no interpolion is required
  magSz_int(1:nxt,1:nynpy,1:nz)=magS(1:nxt,1:nynpy,1:nz)



  ! Modeled fluxes at scale delta (temporary storage in Xi vectors)
  Xx=dscalar_dx*magSx_int(1:nxt+1,1:nynpy+1,1:nz)
  Xy=dscalar_dy*magSy_int(1:nxt+1,1:nynpy+1,1:nz)
  Xz=dscalar_dz*magSz_int(1:nxt+1,1:nynpy+1,1:nz)


  ! Now filter the modeled fluxes
  temp=0._rprec
  temp(1:nxt+1,:,:)=Xx
  CALL box_filter_3D(temp,temp_f,nxt+1,nynpy,nz-1)
  Xx=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=Xy
  CALL box_filter_3D(temp,temp_f,nxt,nynpy+1,nz-1)
  Xy=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=Xz
  CALL box_filter_3D(temp,temp_f,nxt,nynpy,nz-1)
  Xz=temp_f(1:nxt+1,:,:)


  ! Box filter scalar field at test filter scale
  ! Note: this is the original scalar field
  temp=0._rprec
  temp(1:nxt,1:nynpy,:)=scalar(1:nxt,1:nynpy,:)
  CALL box_filter_3D(temp,temp_f,nxt,nynpy,nz-1)
  scalar_f(1:nxt,1:nynpy,:)=temp_f(1:nxt,1:nynpy,:)

  ! Calculate derivatives of filtered scalar
  dscalar_dx(2:nxt,1:nynpy,:)=(scalar_f(2:nxt,1:nynpy,:)-scalar_f(1:nxt-1,1:nynpy,:))/dx
  dscalar_dy(1:nxt,2:nynpy,:)=(scalar_f(1:nxt,2:nynpy,:)-scalar_f(1:nxt,1:nynpy-1,:))/dy
  dscalar_dz(1:nxt,1:nynpy,2:nz)=(scalar_f(1:nxt,1:nynpy,2:nz)-scalar_f(1:nxt,1:nynpy,1:nz-1))/dz


  ! Box filter |S| fields at test filter scale (replace original by filtered)
  ! Note: these are the |S| fields already interpolated to surface locations
  temp=magSx_int
  CALL box_filter_3D(temp,temp_f,nxt+1,nynpy,nz-1)
  magSx_int=temp_f

  temp=magSy_int
  CALL box_filter_3D(temp,temp_f,nxt,nynpy+1,nz-1)
  magSy_int=temp_f

  temp=magSz_int
  CALL box_filter_3D(temp,temp_f,nxt,nynpy,nz-1)
  magSz_int=temp_f

  ! Calculate final vector
  DO jz=1,nz
    Xx(:,:,jz)=Xx(:,:,jz)-4._rprec*magSx_int(1:nxt+1,1:nynpy+1,jz)*dscalar_dx(:,:,jz)
    Xy(:,:,jz)=Xy(:,:,jz)-4._rprec*magSy_int(1:nxt+1,1:nynpy+1,jz)*dscalar_dy(:,:,jz)
    Xz(:,:,jz)=Xz(:,:,jz)-4._rprec*magSz_int(1:nxt+1,1:nynpy+1,jz)*dscalar_dz(:,:,jz)
  END DO
  Xx=((dx*dy*dz)**(2._rprec/3._rprec))*Xx
  Xy=((dx*dy*dz)**(2._rprec/3._rprec))*Xy
  Xz=((dx*dy*dz)**(2._rprec/3._rprec))*Xz


  ! Take care of horizontal boundary problems
  ! Couldx modify this to take advantage of periodic BC when (periodicbc==.TRUE.)
  ! x-component
  Xx(1,:,:)=Xx(4,:,:); Xx(2,:,:)=Xx(4,:,:); Xx(3,:,:)=Xx(4,:,:);
  Xx(:,1,:)=Xx(:,4,:); Xx(:,2,:)=Xx(:,4,:); Xx(:,3,:)=Xx(:,4,:);
  Kx(1,:,:)=Kx(4,:,:); Kx(2,:,:)=Kx(4,:,:); Kx(3,:,:)=Kx(4,:,:);
  Kx(:,1,:)=Kx(:,4,:); Kx(:,2,:)=Kx(:,4,:); Kx(:,3,:)=Kx(:,4,:);
  Xx(nxt,:,:)=Xx(nxt-2,:,:);   Xx(nxt+1,:,:)=Xx(nxt-2,:,:); Xx(nxt-1,:,:)=Xx(nxt-2,:,:);
  Xx(:,nynpy-1,:)=Xx(:,nynpy-3,:); Xx(:,nynpy,:)=Xx(:,nynpy-3,:);   Xx(:,nynpy-2,:)=Xx(:,nynpy-3,:);
  Kx(nxt,:,:)=Kx(nxt-2,:,:);   Kx(nxt+1,:,:)=Kx(nxt-2,:,:); Kx(nxt-1,:,:)=Kx(nxt-2,:,:);
  Kx(:,nynpy-1,:)=Kx(:,nynpy-3,:); Kx(:,nynpy,:)=Kx(:,nynpy-3,:);   Kx(:,nynpy-2,:)=Kx(:,nynpy-3,:);

  ! y-component
  Xy(1,:,:)=Xy(4,:,:); Xy(2,:,:)=Xy(4,:,:); Xy(3,:,:)=Xy(4,:,:);
  Xy(:,1,:)=Xy(:,4,:); Xy(:,2,:)=Xy(:,4,:); Xy(:,3,:)=Xy(:,4,:);
  Ky(1,:,:)=Ky(4,:,:); Ky(2,:,:)=Ky(4,:,:); Ky(3,:,:)=Ky(4,:,:);
  Ky(:,1,:)=Ky(:,4,:); Ky(:,2,:)=Ky(:,4,:); Ky(:,3,:)=Ky(:,4,:);
  Xy(nxt-1,:,:)=Xy(nxt-3,:,:); Xy(nxt,:,:)=Xy(nxt-3,:,:);   Xy(nxt-2,:,:)=Xy(nxt-3,:,:);
  Xy(:,nynpy,:)=Xy(:,nynpy-2,:);   Xy(:,nynpy+1,:)=Xy(:,nynpy-2,:); Xy(:,nynpy-1,:)=Xy(:,nynpy-2,:);
  Ky(nxt-1,:,:)=Ky(nxt-3,:,:); Ky(nxt,:,:)=Ky(nxt-3,:,:);   Ky(nxt-2,:,:)=Ky(nxt-3,:,:);
  Ky(:,nynpy,:)=Ky(:,nynpy-2,:);   Ky(:,nynpy+1,:)=Ky(:,nynpy-2,:); Ky(:,nynpy-1,:)=Ky(:,nynpy-2,:);

  ! z-component
  Xz(1,:,:)=Xz(4,:,:); Xz(2,:,:)=Xz(4,:,:); Xz(3,:,:)=Xz(4,:,:);
  Xz(:,1,:)=Xz(:,4,:); Xz(:,2,:)=Xz(:,4,:); Xz(:,3,:)=Xz(:,4,:);
  Kz(1,:,:)=Kz(4,:,:); Kz(2,:,:)=Kz(4,:,:); Kz(3,:,:)=Kz(4,:,:);
  Kz(:,1,:)=Kz(:,4,:); Kz(:,2,:)=Kz(:,4,:); Kz(:,3,:)=Kz(:,4,:);
  Xz(nxt-1,:,:)=Xz(nxt-3,:,:); Xz(nxt,:,:)=Xz(nxt-3,:,:); Xz(nxt-2,:,:)=Xz(nxt-3,:,:);
  Xz(:,nynpy-1,:)=Xz(:,nynpy-3,:); Xz(:,nynpy,:)=Xz(:,nynpy-3,:); Xz(:,nynpy-2,:)=Xz(:,nynpy-3,:);
  Kz(nxt-1,:,:)=Kz(nxt-3,:,:); Kz(nxt,:,:)=Kz(nxt-3,:,:); Kz(nxt-2,:,:)=Kz(nxt-3,:,:);
  Kz(:,nynpy-1,:)=Kz(:,nynpy-3,:); Kz(:,nynpy,:)=Kz(:,nynpy-3,:); Kz(:,nynpy-2,:)=Kz(:,nynpy-3,:);


  ! Calculate contractions (these are for each element - averaging both faces in each direction)
  KiXi(1:nxt,1:nynpy,2:nz-1)=0.5_rprec*(Kx(1:nxt,1:nynpy,2:nz-1)*Xx(1:nxt,1:nynpy,2:nz-1)+Kx(2:nxt+1,1:nynpy,2:nz-1)*Xx(2:nxt+1,1:nynpy,2:nz-1) &
     +Ky(1:nxt,1:nynpy,2:nz-1)*Xy(1:nxt,1:nynpy,2:nz-1)+Ky(1:nxt,2:nynpy+1,2:nz-1)*Xy(1:nxt,2:nynpy+1,2:nz-1) &
     +Kz(1:nxt,1:nynpy,2:nz-1)*Xz(1:nxt,1:nynpy,2:nz-1)+Kz(1:nxt,1:nynpy,3:nz)*Xz(1:nxt,1:nynpy,3:nz))
  XiXi(1:nxt,1:nynpy,2:nz-1)=0.5_rprec*(Xx(1:nxt,1:nynpy,2:nz-1)*Xx(1:nxt,1:nynpy,2:nz-1)+Xx(2:nxt+1,1:nynpy,2:nz-1)*Xx(2:nxt+1,1:nynpy,2:nz-1) &
     +Xy(1:nxt,1:nynpy,2:nz-1)*Xy(1:nxt,1:nynpy,2:nz-1)+Xy(1:nxt,2:nynpy+1,2:nz-1)*Xy(1:nxt,2:nynpy+1,2:nz-1) &
     +Xz(1:nxt,1:nynpy,2:nz-1)*Xz(1:nxt,1:nynpy,2:nz-1)+Xz(1:nxt,1:nynpy,3:nz)*Xz(1:nxt,1:nynpy,3:nz))



  ! Calculate plane averaged terms
  DO jz=1,nz
    KX_avg(jz)=SUM(KiXi(:,:,jz))/(nxt*nynpy)
    XX_avg(jz)=SUM(XiXi(:,:,jz))/(nxt*nynpy)
  END DO

  ! Linear interpolation for the lower boundary nodes
  ! KX and KX2 are iterpolated between their values at jz=4 and 0 at the surface
  KX_avg(3)=MAX((5._rprec/7._rprec)*KX_avg(4),0._rprec)
  KX_avg(2)=MAX((3._rprec/7._rprec)*KX_avg(4),0._rprec)
  KX_avg(1)=MAX((1._rprec/7._rprec)*KX_avg(4),0._rprec)
  ! XX and XX2 are assumed to be constant
  XX_avg(3)=XX_avg(4)
  XX_avg(2)=XX_avg(4)
  XX_avg(1)=XX_avg(4)


  ! Same for upper boundary
  KX_avg(nz-2)=MAX((5._rprec/7._rprec)*KX_avg(nz-3),0._rprec)
  KX_avg(nz-1)=MAX((3._rprec/7._rprec)*KX_avg(nz-3),0._rprec)
  KX_avg(nz)=MAX((1._rprec/7._rprec)*KX_avg(nz-3),0._rprec)
  XX_avg(nz-2)=XX_avg(nz-3)
  XX_avg(nz-1)=XX_avg(nz-3)
  XX_avg(nz)=XX_avg(nz-3)


  !
  ! Determine the values of Cs2/Sc
  !

  ! Calculate average Cs2/Sc
  DO jz=1,nz

    IF (XX_avg(jz)>0._rprec) THEN
      Cs2Sc(:,:,jz)=MAX(0._rprec,(KX_avg(jz)/XX_avg(jz)))
    ELSE
      Cs2Sc(:,:,jz)=0._rprec
    END IF

  END DO


END SUBROUTINE dynamic_sc_pasi






!
!    ################################################################################
!    ################################################################################
!    ##                                                                            ##
!    ##                            dynamic_sc_pasd.f90                             ##
!    ##                                                                            ##
!    ##                               Developed by                                 ##
!    ##                             Marcelo Chamecki                               ##
!    ##                                                                            ##
!    ##                        The Johns Hopkins University                        ##
!    ##                                                                            ##
!    ##                                05/27/2007                                  ##
!    ##                                                                            ##
!    ################################################################################
!    ################################################################################
!
!    ################################################################################
!
!    PURPOSE:  This routine is the plane averaged scale dependent dynamic model
!              for pollen concentration. It uses the procedure proposed by
!              Bou-Zeid et al. (PoF - 2005)
!
!    ################################################################################
!

SUBROUTINE dynamic_sc_pasd(scalar,scalar_x_tmp,scalar_y_tmp,scalar_z_tmp,u_int_tmp,v_int_tmp,w_int_tmp)

! Modules
  USE types,only:rprec
  USE param,only:ldx,nxt,nynpy,nz,dx,p_count,jt_total, PCon_scheme,&
     periodicbcx, periodicbcy
  USE scalars_module
  USE sgsmodule,only:magS

! No implicit variables
  IMPLICIT NONE

! Main variables
  REAL(kind=rprec),DIMENSION(ldx,nynpy,nz)       :: scalar                 	     ! Scalar field
  REAL(kind=rprec),DIMENSION(ldx,nynpy+1,nz)     :: u_int_tmp,v_int_tmp,w_int_tmp 	     ! Interpolated velocity field
  REAL(kind=rprec),DIMENSION(ldx,nynpy+1,nz)     :: magSx_int,magSy_int,magSz_int! Interpolated |S|
  REAL(kind=rprec),DIMENSION(nxt+1,nynpy+1,nz)   :: scalar_x_tmp,scalar_y_tmp,scalar_z_tmp   ! Interpolated scalar
  REAL(kind=rprec),DIMENSION(nxt+1,nynpy+1,nz)   :: Kx,Ky,Kz                     ! "Leonard" fluxes
  REAL(kind=rprec),DIMENSION(nxt+1,nynpy+1,nz)   :: Xx,Xy,Xz                     ! "Leonard" modeled fluxes
  REAL(kind=rprec),DIMENSION(nxt,nynpy,nz)       :: KiXi,XiXi                    ! Vector dot products (2*delta)
  REAL(kind=rprec),DIMENSION(nxt,nynpy,nz)       :: KiXi2,XiXi2                  ! Vector dot products (4*delta)
  REAL(kind=rprec),DIMENSION(nz)             :: KX_avg,XX_avg		     ! Plane averaged products (2*delta)
  REAL(kind=rprec),DIMENSION(nz)             :: KX2_avg,XX2_avg		     ! Plane averaged products (4*delta)
  REAL(kind=rprec),DIMENSION(nz)             :: Cs2Sc_2d,Cs2Sc_4d	     ! Coefficients at 2*delta and 4*delta
  REAL(kind=rprec),DIMENSION(nz)             :: beta_c		             ! Beta_c

! Filtered fields
  REAL(kind=rprec),DIMENSION(nxt,nynpy,nz)       :: scalar_f                     ! Scalar field
  REAL(kind=rprec),DIMENSION(nxt+1,nynpy+1,nz)   :: scalar_x_f,scalar_y_f,scalar_z_f ! Interpolated scalar field
  REAL(kind=rprec),DIMENSION(ldx,nynpy+1,nz)     :: u_int_f,v_int_f,w_int_f      ! Interpolated velocity field
  REAL(kind=rprec),DIMENSION(nxt+1,nynpy+1,nz)   :: dscalar_dx,dscalar_dy,dscalar_dz ! Scalar derivatives

! Other variables
  REAL(kind=rprec),DIMENSION(ldx,nynpy+1,nz)     :: temp,temp_f                  ! Temporary to filter
  INTEGER                                    :: jz                           ! Counter
  CHARACTER(len=50)                          :: filename

!
! BEGINNING CODE
!


  ! Check advection scheme
  IF (PCon_scheme/=3) THEN
    PRINT*,'Dynamic Sc only available for SMART scheme...'
    STOP
  END IF


  !
  ! First test filter scale - 2*delta
  !

  ! Box filter velocity fields at test filter scale
  ! Note: these are the velocity fields already interpolated to surface locations
  temp=u_int_tmp
  CALL box_filter_3D(temp,temp_f,nxt+1,nynpy,nz-1)
  u_int_f=temp_f

  temp=v_int_tmp
  CALL box_filter_3D(temp,temp_f,nxt,nynpy+1,nz-1)
  v_int_f=temp_f

  temp=w_int_tmp
  CALL box_filter_3D(temp,temp_f,nxt,nynpy,nz-1)
  w_int_f=temp_f


  ! Box filter scalar fields at test filter scale
  ! Note: this is the scalar field already interpolated to surface locations
  temp=0._rprec
  temp(1:nxt+1,:,:)=scalar_x_tmp
  CALL box_filter_3D(temp,temp_f,nxt+1,nynpy,nz-1)
  scalar_x_f=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=scalar_y_tmp
  CALL box_filter_3D(temp,temp_f,nxt,nynpy+1,nz-1)
  scalar_y_f=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=scalar_z_tmp
  CALL box_filter_3D(temp,temp_f,nxt,nynpy,nz-1)
  scalar_z_f=temp_f(1:nxt+1,:,:)


  ! First use the vectors Kx, Ky and Kz to calculate the resolved fluxes
  Kx=u_int_tmp(1:nxt+1,1:nynpy+1,1:nz)*scalar_x_tmp
  Ky=v_int_tmp(1:nxt+1,1:nynpy+1,1:nz)*scalar_y_tmp
  Kz=w_int_tmp(1:nxt+1,1:nynpy+1,1:nz)*scalar_z_tmp


  ! Now filter the resolved fluxes
  temp=0._rprec
  temp(1:nxt+1,:,:)=Kx
  CALL box_filter_3D(temp,temp_f,nxt+1,nynpy,nz-1)
  Kx=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=Ky
  CALL box_filter_3D(temp,temp_f,nxt,nynpy+1,nz-1)
  Ky=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=Kz
  CALL box_filter_3D(temp,temp_f,nxt,nynpy,nz-1)
  Kz=temp_f(1:nxt+1,:,:)


  ! Subtract product of filtered fields to get "Leonard" fluxes
  Kx=Kx-u_int_f(1:nxt+1,1:nynpy+1,1:nz)*scalar_x_f
  Ky=Ky-v_int_f(1:nxt+1,1:nynpy+1,1:nz)*scalar_y_f
  Kz=Kz-w_int_f(1:nxt+1,1:nynpy+1,1:nz)*scalar_z_f


  ! Calculate derivatives
  dscalar_dx(2:nxt,1:nynpy,:)=(scalar(2:nxt,1:nynpy,:)-scalar(1:nxt-1,1:nynpy,:))/dx
  dscalar_dy(1:nxt,2:nynpy,:)=(scalar(1:nxt,2:nynpy,:)-scalar(1:nxt,1:nynpy-1,:))/dy
  dscalar_dz(1:nxt,1:nynpy,2:nz-1)=(scalar(1:nxt,1:nynpy,2:nz-1)-scalar(1:nxt,1:nynpy,1:nz-2))/dz


  ! |S| has to be interpolated (same linear interpolation as magS in the pollen_RHS_calc2 routine
  ! For x-diffusion interpolate in x and z directions
  magSx_int(2:nxt,1:nynpy,1)=(magS(1:nxt-1,1:nynpy,1)+magS(2:nxt,1:nynpy,1))/2._rprec
  magSx_int(1,1:nynpy,1)=(magS(nxt,1:nynpy,1)+magS(1,1:nynpy,1))/2._rprec

  magSx_int(2:nxt,1:nynpy,2:nz-1)=(magS(1:nxt-1,1:nynpy,2:nz-1)+magS(2:nxt,1:nynpy,2:nz-1)+ &
     magS(1:nxt-1,1:nynpy,3:nz)+magS(2:nxt,1:nynpy,3:nz))/4._rprec
  magSx_int(1,1:nynpy,2:nz-1)=(magS(nxt,1:nynpy,2:nz-1)+magS(1,1:nynpy,2:nz-1)+ &
     magS(nxt,1:nynpy,3:nz)+magS(1,1:nynpy,3:nz))/4._rprec
  magSx_int(nxt+1,:,:)=magSx_int(1,:,:)


  ! For y-diffusion interpolate in y and z directions
  magSy_int(1:nxt,2:nynpy,1)=(magS(1:nxt,1:nynpy-1,1)+magS(1:nxt,2:nynpy,1))/2._rprec
  magSy_int(1:nxt,1,1)=(magS(1:nxt,nynpy,1)+magS(1:nxt,1,1))/2._rprec

  magSy_int(1:nxt,2:nynpy,2:nz-1)=(magS(1:nxt,1:nynpy-1,2:nz-1)+magS(1:nxt,2:nynpy,2:nz-1)+ &
     magS(1:nxt,1:nynpy-1,3:nz)+magS(1:nxt,2:nynpy,3:nz))/4._rprec
  magSy_int(1:nxt,1,2:nz-1)=(magS(1:nxt,nynpy,2:nz-1)+magS(1:nxt,1,2:nz-1)+ &
     magS(1:nxt,nynpy,3:nz)+magS(1:nxt,1,3:nz))/4._rprec
  magSy_int(:,nynpy+1,:)=magSy_int(:,1,:)


  ! For z-diffusion no interpolion is required
  magSz_int(1:nxt,1:nynpy,1:nz)=magS(1:nxt,1:nynpy,1:nz)



  ! Modeled fluxes at scale delta (temporary storage in Xi vectors)
  Xx=dscalar_dx*magSx_int(1:nxt+1,1:nynpy+1,1:nz)
  Xy=dscalar_dy*magSy_int(1:nxt+1,1:nynpy+1,1:nz)
  Xz=dscalar_dz*magSz_int(1:nxt+1,1:nynpy+1,1:nz)


  ! Now filter the modeled fluxes
  temp=0._rprec
  temp(1:nxt+1,:,:)=Xx
  CALL box_filter_3D(temp,temp_f,nxt+1,nynpy,nz-1)
  Xx=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=Xy
  CALL box_filter_3D(temp,temp_f,nxt,nynpy+1,nz-1)
  Xy=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=Xz
  CALL box_filter_3D(temp,temp_f,nxt,nynpy,nz-1)
  Xz=temp_f(1:nxt+1,:,:)


  ! Box filter scalar field at test filter scale
  ! Note: this is the original scalar field
  temp=0._rprec
  temp(1:nxt,1:nynpy,:)=scalar(1:nxt,1:nynpy,:)
  CALL box_filter_3D(temp,temp_f,nxt,nynpy,nz-1)
  scalar_f(1:nxt,1:nynpy,:)=temp_f(1:nxt,1:nynpy,:)

  ! Calculate derivatives of filtered scalar
  dscalar_dx(2:nxt,1:nynpy,:)=(scalar_f(2:nxt,1:nynpy,:)-scalar_f(1:nxt-1,1:nynpy,:))/dx
  dscalar_dy(1:nxt,2:nynpy,:)=(scalar_f(1:nxt,2:nynpy,:)-scalar_f(1:nxt,1:nynpy-1,:))/dy
  dscalar_dz(1:nxt,1:nynpy,2:nz)=(scalar_f(1:nxt,1:nynpy,2:nz)-scalar_f(1:nxt,1:nynpy,1:nz-1))/dz


  ! Box filter |S| fields at test filter scale (replace original by filtered)
  ! Note: these are the |S| fields already interpolated to surface locations
  temp=magSx_int
  CALL box_filter_3D(temp,temp_f,nxt+1,nynpy,nz-1)
  magSx_int=temp_f

  temp=magSy_int
  CALL box_filter_3D(temp,temp_f,nxt,nynpy+1,nz-1)
  magSy_int=temp_f

  temp=magSz_int
  CALL box_filter_3D(temp,temp_f,nxt,nynpy,nz-1)
  magSz_int=temp_f

  ! Calculate final vector
  DO jz=1,nz
    Xx(:,:,jz)=Xx(:,:,jz)-4._rprec*magSx_int(1:nxt+1,1:nynpy+1,jz)*dscalar_dx(:,:,jz)
    Xy(:,:,jz)=Xy(:,:,jz)-4._rprec*magSy_int(1:nxt+1,1:nynpy+1,jz)*dscalar_dy(:,:,jz)
    Xz(:,:,jz)=Xz(:,:,jz)-4._rprec*magSz_int(1:nxt+1,1:nynpy+1,jz)*dscalar_dz(:,:,jz)
  END DO
  Xx=((dx*dy*dz)**(2._rprec/3._rprec))*Xx
  Xy=((dx*dy*dz)**(2._rprec/3._rprec))*Xy
  Xz=((dx*dy*dz)**(2._rprec/3._rprec))*Xz


  ! Take care of horizontal boundary problems
  ! x-component
  Xx(1,:,:)=Xx(4,:,:); Xx(2,:,:)=Xx(4,:,:); Xx(3,:,:)=Xx(4,:,:);
  Xx(:,1,:)=Xx(:,4,:); Xx(:,2,:)=Xx(:,4,:); Xx(:,3,:)=Xx(:,4,:);
  Kx(1,:,:)=Kx(4,:,:); Kx(2,:,:)=Kx(4,:,:); Kx(3,:,:)=Kx(4,:,:);
  Kx(:,1,:)=Kx(:,4,:); Kx(:,2,:)=Kx(:,4,:); Kx(:,3,:)=Kx(:,4,:);
  Xx(nxt,:,:)=Xx(nxt-2,:,:);   Xx(nxt+1,:,:)=Xx(nxt-2,:,:); Xx(nxt-1,:,:)=Xx(nxt-2,:,:);
  Xx(:,nynpy-1,:)=Xx(:,nynpy-3,:); Xx(:,nynpy,:)=Xx(:,nynpy-3,:);   Xx(:,nynpy-2,:)=Xx(:,nynpy-3,:);
  Kx(nxt,:,:)=Kx(nxt-2,:,:);   Kx(nxt+1,:,:)=Kx(nxt-2,:,:); Kx(nxt-1,:,:)=Kx(nxt-2,:,:);
  Kx(:,nynpy-1,:)=Kx(:,nynpy-3,:); Kx(:,nynpy,:)=Kx(:,nynpy-3,:);   Kx(:,nynpy-2,:)=Kx(:,nynpy-3,:);

  ! y-component
  Xy(1,:,:)=Xy(4,:,:); Xy(2,:,:)=Xy(4,:,:); Xy(3,:,:)=Xy(4,:,:);
  Xy(:,1,:)=Xy(:,4,:); Xy(:,2,:)=Xy(:,4,:); Xy(:,3,:)=Xy(:,4,:);
  Ky(1,:,:)=Ky(4,:,:); Ky(2,:,:)=Ky(4,:,:); Ky(3,:,:)=Ky(4,:,:);
  Ky(:,1,:)=Ky(:,4,:); Ky(:,2,:)=Ky(:,4,:); Ky(:,3,:)=Ky(:,4,:);
  Xy(nxt-1,:,:)=Xy(nxt-3,:,:); Xy(nxt,:,:)=Xy(nxt-3,:,:);   Xy(nxt-2,:,:)=Xy(nxt-3,:,:);
  Xy(:,nynpy,:)=Xy(:,nynpy-2,:);   Xy(:,nynpy+1,:)=Xy(:,nynpy-2,:); Xy(:,nynpy-1,:)=Xy(:,nynpy-2,:);
  Ky(nxt-1,:,:)=Ky(nxt-3,:,:); Ky(nxt,:,:)=Ky(nxt-3,:,:);   Ky(nxt-2,:,:)=Ky(nxt-3,:,:);
  Ky(:,nynpy,:)=Ky(:,nynpy-2,:);   Ky(:,nynpy+1,:)=Ky(:,nynpy-2,:); Ky(:,nynpy-1,:)=Ky(:,nynpy-2,:);

  ! z-component
  Xz(1,:,:)=Xz(4,:,:); Xz(2,:,:)=Xz(4,:,:); Xz(3,:,:)=Xz(4,:,:);
  Xz(:,1,:)=Xz(:,4,:); Xz(:,2,:)=Xz(:,4,:); Xz(:,3,:)=Xz(:,4,:);
  Kz(1,:,:)=Kz(4,:,:); Kz(2,:,:)=Kz(4,:,:); Kz(3,:,:)=Kz(4,:,:);
  Kz(:,1,:)=Kz(:,4,:); Kz(:,2,:)=Kz(:,4,:); Kz(:,3,:)=Kz(:,4,:);
  Xz(nxt-1,:,:)=Xz(nxt-3,:,:); Xz(nxt,:,:)=Xz(nxt-3,:,:); Xz(nxt-2,:,:)=Xz(nxt-3,:,:);
  Xz(:,nynpy-1,:)=Xz(:,nynpy-3,:); Xz(:,nynpy,:)=Xz(:,nynpy-3,:); Xz(:,nynpy-2,:)=Xz(:,nynpy-3,:);
  Kz(nxt-1,:,:)=Kz(nxt-3,:,:); Kz(nxt,:,:)=Kz(nxt-3,:,:); Kz(nxt-2,:,:)=Kz(nxt-3,:,:);
  Kz(:,nynpy-1,:)=Kz(:,nynpy-3,:); Kz(:,nynpy,:)=Kz(:,nynpy-3,:); Kz(:,nynpy-2,:)=Kz(:,nynpy-3,:);


  ! Calculate contractions (these are for each element - averaging both faces in each direction)
  KiXi(1:nxt,1:nynpy,1:nz-1)=0.5_rprec*(Kx(1:nxt,1:nynpy,1:nz-1)*Xx(1:nxt,1:nynpy,1:nz-1)+Kx(2:nxt+1,1:nynpy,1:nz-1)*Xx(2:nxt+1,1:nynpy,1:nz-1) &
     +Ky(1:nxt,1:nynpy,1:nz-1)*Xy(1:nxt,1:nynpy,1:nz-1)+Ky(1:nxt,2:nynpy+1,1:nz-1)*Xy(1:nxt,2:nynpy+1,1:nz-1) &
     +Kz(1:nxt,1:nynpy,1:nz-1)*Xz(1:nxt,1:nynpy,1:nz-1)+Kz(1:nxt,1:nynpy,2:nz)*Xz(1:nxt,1:nynpy,2:nz))
  XiXi(1:nxt,1:nynpy,1:nz-1)=0.5_rprec*(Xx(1:nxt,1:nynpy,1:nz-1)*Xx(1:nxt,1:nynpy,1:nz-1)+Xx(2:nxt+1,1:nynpy,1:nz-1)*Xx(2:nxt+1,1:nynpy,1:nz-1) &
     +Xy(1:nxt,1:nynpy,1:nz-1)*Xy(1:nxt,1:nynpy,1:nz-1)+Xy(1:nxt,2:nynpy+1,1:nz-1)*Xy(1:nxt,2:nynpy+1,1:nz-1) &
     +Xz(1:nxt,1:nynpy,1:nz-1)*Xz(1:nxt,1:nynpy,1:nz-1)+Xz(1:nxt,1:nynpy,2:nz)*Xz(1:nxt,1:nynpy,2:nz))



  !
  ! Second test filter scale - 4*delta
  !


  ! Box filter velocity fields at second test filter scale
  ! Note: these are the velocity fields already interpolated to surface locations
  temp=u_int_tmp
  CALL box_filter_3D_4delta(temp,temp_f,nxt+1,nynpy,nz-1)
  u_int_f=temp_f

  temp=v_int_tmp
  CALL box_filter_3D_4delta(temp,temp_f,nxt,nynpy+1,nz-1)
  v_int_f=temp_f

  temp=w_int_tmp
  CALL box_filter_3D_4delta(temp,temp_f,nxt,nynpy,nz-1)
  w_int_f=temp_f


  ! Box filter scalar fields at second test filter scale
  ! Note: this is the scalar field already interpolated to surface locations
  temp=0._rprec
  temp(1:nxt+1,:,:)=scalar_x_tmp
  CALL box_filter_3D_4delta(temp,temp_f,nxt+1,nynpy,nz-1)
  scalar_x_f=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=scalar_y_tmp
  CALL box_filter_3D_4delta(temp,temp_f,nxt,nynpy+1,nz-1)
  scalar_y_f=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=scalar_z_tmp
  CALL box_filter_3D_4delta(temp,temp_f,nxt,nynpy,nz-1)
  scalar_z_f=temp_f(1:nxt+1,:,:)


  ! First use the vectors Kx, Ky and Kz to calculate the resolved fluxes
  Kx=u_int_tmp(1:nxt+1,1:nynpy+1,1:nz)*scalar_x_tmp
  Ky=v_int_tmp(1:nxt+1,1:nynpy+1,1:nz)*scalar_y_tmp
  Kz=w_int_tmp(1:nxt+1,1:nynpy+1,1:nz)*scalar_z_tmp


  ! Now filter the resolved fluxes
  temp=0._rprec
  temp(1:nxt+1,:,:)=Kx
  CALL box_filter_3D_4delta(temp,temp_f,nxt+1,nynpy,nz-1)
  Kx=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=Ky
  CALL box_filter_3D_4delta(temp,temp_f,nxt,nynpy+1,nz-1)
  Ky=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=Kz
  CALL box_filter_3D_4delta(temp,temp_f,nxt,nynpy,nz-1)
  Kz=temp_f(1:nxt+1,:,:)


  ! Subtract product of filtered fields to get "Leonard" fluxes
  Kx=Kx-u_int_f(1:nxt+1,1:nynpy+1,1:nz)*scalar_x_f
  Ky=Ky-v_int_f(1:nxt+1,1:nynpy+1,1:nz)*scalar_y_f
  Kz=Kz-w_int_f(1:nxt+1,1:nynpy+1,1:nz)*scalar_z_f


  ! Calculate derivatives
  dscalar_dx(2:nxt,1:nynpy,:)=(scalar(2:nxt,1:nynpy,:)-scalar(1:nxt-1,1:nynpy,:))/dx
  dscalar_dy(1:nxt,2:nynpy,:)=(scalar(1:nxt,2:nynpy,:)-scalar(1:nxt,1:nynpy-1,:))/dy
  dscalar_dz(1:nxt,1:nynpy,2:nz-1)=(scalar(1:nxt,1:nynpy,2:nz-1)-scalar(1:nxt,1:nynpy,1:nz-2))/dz


  ! |S| has to be interpolated (same linear interpolation as magS in the pollen_RHS_calc2 routine
  ! For x-diffusion interpolate in x and z directions
  magSx_int(2:nxt,1:nynpy,1)=(magS(1:nxt-1,1:nynpy,1)+magS(2:nxt,1:nynpy,1))/2._rprec
  magSx_int(1,1:nynpy,1)=(magS(nxt,1:nynpy,1)+magS(1,1:nynpy,1))/2._rprec

  magSx_int(2:nxt,1:nynpy,2:nz-1)=(magS(1:nxt-1,1:nynpy,2:nz-1)+magS(2:nxt,1:nynpy,2:nz-1)+ &
     magS(1:nxt-1,1:nynpy,3:nz)+magS(2:nxt,1:nynpy,3:nz))/4._rprec
  magSx_int(1,1:nynpy,2:nz-1)=(magS(nxt,1:nynpy,2:nz-1)+magS(1,1:nynpy,2:nz-1)+ &
     magS(nxt,1:nynpy,3:nz)+magS(1,1:nynpy,3:nz))/4._rprec
  magSx_int(nxt+1,:,:)=magSx_int(1,:,:)


  ! For y-diffusion interpolate in y and z directions
  magSy_int(1:nxt,2:nynpy,1)=(magS(1:nxt,1:nynpy-1,1)+magS(1:nxt,2:nynpy,1))/2._rprec
  magSy_int(1:nxt,1,1)=(magS(1:nxt,nynpy,1)+magS(1:nxt,1,1))/2._rprec

  magSy_int(1:nxt,2:nynpy,2:nz-1)=(magS(1:nxt,1:nynpy-1,2:nz-1)+magS(1:nxt,2:nynpy,2:nz-1)+ &
     magS(1:nxt,1:nynpy-1,3:nz)+magS(1:nxt,2:nynpy,3:nz))/4._rprec
  magSy_int(1:nxt,1,2:nz-1)=(magS(1:nxt,nynpy,2:nz-1)+magS(1:nxt,1,2:nz-1)+ &
     magS(1:nxt,nynpy,3:nz)+magS(1:nxt,1,3:nz))/4._rprec
  magSy_int(:,nynpy+1,:)=magSy_int(:,1,:)


  ! For z-diffusion no interpolion is required
  magSz_int(1:nxt,1:nynpy,1:nz)=magS(1:nxt,1:nynpy,1:nz)



  ! Modeled fluxes at scale delta (temporary storage in Xi vectors)
  Xx=dscalar_dx*magSx_int(1:nxt+1,1:nynpy+1,1:nz)
  Xy=dscalar_dy*magSy_int(1:nxt+1,1:nynpy+1,1:nz)
  Xz=dscalar_dz*magSz_int(1:nxt+1,1:nynpy+1,1:nz)


  ! Now filter the modeled fluxes
  temp=0._rprec
  temp(1:nxt+1,:,:)=Xx
  CALL box_filter_3D_4delta(temp,temp_f,nxt+1,nynpy,nz-1)
  Xx=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=Xy
  CALL box_filter_3D_4delta(temp,temp_f,nxt,nynpy+1,nz-1)
  Xy=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=Xz
  CALL box_filter_3D_4delta(temp,temp_f,nxt,nynpy,nz-1)
  Xz=temp_f(1:nxt+1,:,:)


  ! Box filter scalar field at test filter scale
  ! Note: this is the original scalar field
  temp=0._rprec
  temp(1:nxt,1:nynpy,:)=scalar(1:nxt,1:nynpy,:)
  CALL box_filter_3D_4delta(temp,temp_f,nxt,nynpy,nz-1)
  scalar_f(1:nxt,1:nynpy,:)=temp_f(1:nxt,1:nynpy,:)

  ! Calculate derivatives of filtered scalar
  dscalar_dx(2:nxt,1:nynpy,:)=(scalar_f(2:nxt,1:nynpy,:)-scalar_f(1:nxt-1,1:nynpy,:))/dx
  dscalar_dy(1:nxt,2:nynpy,:)=(scalar_f(1:nxt,2:nynpy,:)-scalar_f(1:nxt,1:nynpy-1,:))/dy
  dscalar_dz(1:nxt,1:nynpy,2:nz)=(scalar_f(1:nxt,1:nynpy,2:nz)-scalar_f(1:nxt,1:nynpy,1:nz-1))/dz


  ! Box filter |S| fields at test filter scale (replace original by filtered)
  ! Note: these are the |S| fields already interpolated to surface locations
  temp=magSx_int
  CALL box_filter_3D_4delta(temp,temp_f,nxt+1,nynpy,nz-1)
  magSx_int=temp_f

  temp=magSy_int
  CALL box_filter_3D_4delta(temp,temp_f,nxt,nynpy+1,nz-1)
  magSy_int=temp_f

  temp=magSz_int
  CALL box_filter_3D_4delta(temp,temp_f,nxt,nynpy,nz-1)
  magSz_int=temp_f

  ! Calculate final vector
  DO jz=1,nz
    Xx(:,:,jz)=Xx(:,:,jz)-16._rprec*magSx_int(1:nxt+1,1:nynpy+1,jz)*dscalar_dx(:,:,jz)
    Xy(:,:,jz)=Xy(:,:,jz)-16._rprec*magSy_int(1:nxt+1,1:nynpy+1,jz)*dscalar_dy(:,:,jz)
    Xz(:,:,jz)=Xz(:,:,jz)-16._rprec*magSz_int(1:nxt+1,1:nynpy+1,jz)*dscalar_dz(:,:,jz)
  END DO
  Xx=((dx*dy*dz)**(2._rprec/3._rprec))*Xx
  Xy=((dx*dy*dz)**(2._rprec/3._rprec))*Xy
  Xz=((dx*dy*dz)**(2._rprec/3._rprec))*Xz


  ! Take care of horizontal boundary problems
  ! x-component
  Xx(1,:,:)=Xx(4,:,:); Xx(2,:,:)=Xx(4,:,:); Xx(3,:,:)=Xx(4,:,:);
  Xx(:,1,:)=Xx(:,4,:); Xx(:,2,:)=Xx(:,4,:); Xx(:,3,:)=Xx(:,4,:);
  Kx(1,:,:)=Kx(4,:,:); Kx(2,:,:)=Kx(4,:,:); Kx(3,:,:)=Kx(4,:,:);
  Kx(:,1,:)=Kx(:,4,:); Kx(:,2,:)=Kx(:,4,:); Kx(:,3,:)=Kx(:,4,:);
  Xx(nxt,:,:)=Xx(nxt-2,:,:);   Xx(nxt+1,:,:)=Xx(nxt-2,:,:); Xx(nxt-1,:,:)=Xx(nxt-2,:,:);
  Xx(:,nynpy-1,:)=Xx(:,nynpy-3,:); Xx(:,nynpy,:)=Xx(:,nynpy-3,:);   Xx(:,nynpy-2,:)=Xx(:,nynpy-3,:);
  Kx(nxt,:,:)=Kx(nxt-2,:,:);   Kx(nxt+1,:,:)=Kx(nxt-2,:,:); Kx(nxt-1,:,:)=Kx(nxt-2,:,:);
  Kx(:,nynpy-1,:)=Kx(:,nynpy-3,:); Kx(:,nynpy,:)=Kx(:,nynpy-3,:);   Kx(:,nynpy-2,:)=Kx(:,nynpy-3,:);

  ! y-component
  Xy(1,:,:)=Xy(4,:,:); Xy(2,:,:)=Xy(4,:,:); Xy(3,:,:)=Xy(4,:,:);
  Xy(:,1,:)=Xy(:,4,:); Xy(:,2,:)=Xy(:,4,:); Xy(:,3,:)=Xy(:,4,:);
  Ky(1,:,:)=Ky(4,:,:); Ky(2,:,:)=Ky(4,:,:); Ky(3,:,:)=Ky(4,:,:);
  Ky(:,1,:)=Ky(:,4,:); Ky(:,2,:)=Ky(:,4,:); Ky(:,3,:)=Ky(:,4,:);
  Xy(nxt-1,:,:)=Xy(nxt-3,:,:); Xy(nxt,:,:)=Xy(nxt-3,:,:);   Xy(nxt-2,:,:)=Xy(nxt-3,:,:);
  Xy(:,nynpy,:)=Xy(:,nynpy-2,:);   Xy(:,nynpy+1,:)=Xy(:,nynpy-2,:); Xy(:,nynpy-1,:)=Xy(:,nynpy-2,:);
  Ky(nxt-1,:,:)=Ky(nxt-3,:,:); Ky(nxt,:,:)=Ky(nxt-3,:,:);   Ky(nxt-2,:,:)=Ky(nxt-3,:,:);
  Ky(:,nynpy,:)=Ky(:,nynpy-2,:);   Ky(:,nynpy+1,:)=Ky(:,nynpy-2,:); Ky(:,nynpy-1,:)=Ky(:,nynpy-2,:);

  ! z-component
  Xz(1,:,:)=Xz(4,:,:); Xz(2,:,:)=Xz(4,:,:); Xz(3,:,:)=Xz(4,:,:);
  Xz(:,1,:)=Xz(:,4,:); Xz(:,2,:)=Xz(:,4,:); Xz(:,3,:)=Xz(:,4,:);
  Kz(1,:,:)=Kz(4,:,:); Kz(2,:,:)=Kz(4,:,:); Kz(3,:,:)=Kz(4,:,:);
  Kz(:,1,:)=Kz(:,4,:); Kz(:,2,:)=Kz(:,4,:); Kz(:,3,:)=Kz(:,4,:);
  Xz(nxt-1,:,:)=Xz(nxt-3,:,:); Xz(nxt,:,:)=Xz(nxt-3,:,:); Xz(nxt-2,:,:)=Xz(nxt-3,:,:);
  Xz(:,nynpy-1,:)=Xz(:,nynpy-3,:); Xz(:,nynpy,:)=Xz(:,nynpy-3,:); Xz(:,nynpy-2,:)=Xz(:,nynpy-3,:);
  Kz(nxt-1,:,:)=Kz(nxt-3,:,:); Kz(nxt,:,:)=Kz(nxt-3,:,:); Kz(nxt-2,:,:)=Kz(nxt-3,:,:);
  Kz(:,nynpy-1,:)=Kz(:,nynpy-3,:); Kz(:,nynpy,:)=Kz(:,nynpy-3,:); Kz(:,nynpy-2,:)=Kz(:,nynpy-3,:);


  ! Calculate contractions (these are for each element - averaging both faces in each direction)
  KiXi2(1:nxt,1:nynpy,1:nz-1)=0.5_rprec*(Kx(1:nxt,1:nynpy,1:nz-1)*Xx(1:nxt,1:nynpy,1:nz-1)+Kx(2:nxt+1,1:nynpy,1:nz-1)*Xx(2:nxt+1,1:nynpy,1:nz-1) &
     +Ky(1:nxt,1:nynpy,1:nz-1)*Xy(1:nxt,1:nynpy,1:nz-1)+Ky(1:nxt,2:nynpy+1,1:nz-1)*Xy(1:nxt,2:nynpy+1,1:nz-1) &
     +Kz(1:nxt,1:nynpy,1:nz-1)*Xz(1:nxt,1:nynpy,1:nz-1)+Kz(1:nxt,1:nynpy,2:nz)*Xz(1:nxt,1:nynpy,2:nz))
  XiXi2(1:nxt,1:nynpy,1:nz-1)=0.5_rprec*(Xx(1:nxt,1:nynpy,1:nz-1)*Xx(1:nxt,1:nynpy,1:nz-1)+Xx(2:nxt+1,1:nynpy,1:nz-1)*Xx(2:nxt+1,1:nynpy,1:nz-1) &
     +Xy(1:nxt,1:nynpy,1:nz-1)*Xy(1:nxt,1:nynpy,1:nz-1)+Xy(1:nxt,2:nynpy+1,1:nz-1)*Xy(1:nxt,2:nynpy+1,1:nz-1) &
     +Xz(1:nxt,1:nynpy,1:nz-1)*Xz(1:nxt,1:nynpy,1:nz-1)+Xz(1:nxt,1:nynpy,2:nz)*Xz(1:nxt,1:nynpy,2:nz))




  ! Calculate plane averaged terms
  DO jz=1,nz
    KX_avg(jz)=SUM(KiXi(:,:,jz))/(nxt*nynpy)
    XX_avg(jz)=SUM(XiXi(:,:,jz))/(nxt*nynpy)
    KX2_avg(jz)=SUM(KiXi2(:,:,jz))/(nxt*nynpy)
    XX2_avg(jz)=SUM(XiXi2(:,:,jz))/(nxt*nynpy)
  END DO

  ! Linear interpolation for the lower boundary nodes
  ! KX and KX2 are iterpolated between their values at jz=4 and 0 at the surface
  KX_avg(3)=MAX((5._rprec/7._rprec)*KX_avg(4),0._rprec)
  KX_avg(2)=MAX((3._rprec/7._rprec)*KX_avg(4),0._rprec)
  KX_avg(1)=MAX((1._rprec/7._rprec)*KX_avg(4),0._rprec)
  KX2_avg(3)=MAX((5._rprec/7._rprec)*KX2_avg(4),0._rprec)
  KX2_avg(2)=MAX((3._rprec/7._rprec)*KX2_avg(4),0._rprec)
  KX2_avg(1)=MAX((1._rprec/7._rprec)*KX2_avg(4),0._rprec)
  ! XX and XX2 are assumed to be constant
  XX_avg(3)=XX_avg(4)
  XX_avg(2)=XX_avg(4)
  XX_avg(1)=XX_avg(4)
  XX2_avg(3)=XX2_avg(4)
  XX2_avg(2)=XX2_avg(4)
  XX2_avg(1)=XX2_avg(4)


  ! Same for upper boundary
  KX_avg(nz-2)=MAX((5._rprec/7._rprec)*KX_avg(nz-3),0._rprec)
  KX_avg(nz-1)=MAX((3._rprec/7._rprec)*KX_avg(nz-3),0._rprec)
  KX_avg(nz)=MAX((1._rprec/7._rprec)*KX_avg(nz-3),0._rprec)
  KX2_avg(nz-2)=MAX((5._rprec/7._rprec)*KX2_avg(nz-3),0._rprec)
  KX2_avg(nz-1)=MAX((3._rprec/7._rprec)*KX2_avg(nz-3),0._rprec)
  KX2_avg(nz)=MAX((1._rprec/7._rprec)*KX2_avg(nz-3),0._rprec)
  XX_avg(nz-2)=XX_avg(nz-3)
  XX_avg(nz-1)=XX_avg(nz-3)
  XX_avg(nz)=XX_avg(nz-3)
  XX2_avg(nz-2)=XX2_avg(nz-3)
  XX2_avg(nz-1)=XX2_avg(nz-3)
  XX2_avg(nz)=XX2_avg(nz-3)


  !
  ! Determine the values of beta_c and Cs2/Sc
  !


  ! Calculate average beta_c
  DO jz=1,nz

    IF (XX2_avg(jz)>0._rprec) THEN
      Cs2Sc_4d(jz)=MAX(0._rprec,KX2_avg(jz)/XX2_avg(jz))
    ELSE
      Cs2Sc_4d(jz)=0._rprec
    END IF

    IF (XX_avg(jz)>0._rprec) THEN
      Cs2Sc_2d(jz)=MAX(0._rprec,KX_avg(jz)/XX_avg(jz))
    ELSE
      Cs2Sc_2d(jz)=0._rprec
    END IF


    IF (XX2_avg(jz)*KX_avg(jz)>0._rprec) THEN
      beta_c(jz)=MAX(0.125_rprec,(KX2_avg(jz)*XX_avg(jz))/(KX_avg(jz)*XX2_avg(jz)))
    ELSE
      beta_c(jz)=0.125_rprec
    END IF
  END DO


  ! Calculate average Cs2/Sc
  DO jz=1,nz

    IF (XX_avg(jz)>0._rprec) THEN
      Cs2Sc(:,:,jz)=MAX(0._rprec,(KX_avg(jz)/XX_avg(jz))/beta_c(jz))
    ELSE
      Cs2Sc(:,:,jz)=0._rprec
    END IF

  END DO


  ! Output for analysis
  IF (mod(jt_total+1,p_count)==0) THEN

    WRITE(filename,'(A,I0.8,A)')'result/dyn_',(jt_total+1),'.txt'
    OPEN(1,FILE=filename,STATUS="unknown")
    DO jz=1,nz
      WRITE(1,'(I8,4F16.8)')jz,Cs2Sc(1,1,jz),beta_c(jz),Cs2Sc_2d(jz),Cs2Sc_4d(jz)
    END DO
    CLOSE(1)

  END IF


END SUBROUTINE dynamic_sc_pasd







!
!    ################################################################################
!    ################################################################################
!    ##                                                                            ##
!    ##                            dynamic_sc_lasi.f90                             ##
!    ##                                                                            ##
!    ##                               Developed by                                 ##
!    ##                             Marcelo Chamecki                               ##
!    ##                                                                            ##
!    ##                        The Johns Hopkins University                        ##
!    ##                                                                            ##
!    ##                                05/27/2007                                  ##
!    ##                                                                            ##
!    ################################################################################
!    ################################################################################
!
!    ################################################################################
!
!    PURPOSE:  This routine is the Lagrangian averaged scale invariant dynamic model
!              for pollen concentration.
!
!    ################################################################################
!

SUBROUTINE dynamic_sc_lasi(scalar,scalar_x_tmp,scalar_y_tmp,scalar_z_tmp,u_int_tmp,v_int_tmp,w_int_tmp)

! Modules
  USE types,only:rprec
  USE param,only:ldx,nxt,nynpy,nz,dx,PCon_scheme, periodicbcx, periodicbcy
  USE scalars_module
  USE sgsmodule,only:magS,epsilon_lag,F_KX,F_XX

! No implicit variables
  IMPLICIT NONE

! Main variables
  REAL(kind=rprec),DIMENSION(ldx,nynpy,nz)       :: scalar                      ! Scalar field
  REAL(kind=rprec),DIMENSION(ldx,nynpy+1,nz)     :: u_int_tmp,v_int_tmp,w_int_tmp      ! Interpolated velocity field
  REAL(kind=rprec),DIMENSION(ldx,nynpy+1,nz)     :: magSx_int,magSy_int,magSz_int! Interpolated |S|
  REAL(kind=rprec),DIMENSION(nxt+1,nynpy+1,nz)   :: scalar_x_tmp,scalar_y_tmp,scalar_z_tmp   ! Interpolated scalar
  REAL(kind=rprec),DIMENSION(nxt+1,nynpy+1,nz)   :: Kx,Ky,Kz                     ! "Leonard" fluxes
  REAL(kind=rprec),DIMENSION(nxt+1,nynpy+1,nz)   :: Xx,Xy,Xz                     ! "Leonard" modeled fluxes
  REAL(kind=rprec),DIMENSION(nxt,nynpy,nz)       :: KiXi,XiXi                    ! Vector dot products (2*delta)

! Filtered fields
  REAL(kind=rprec),DIMENSION(nxt,nynpy,nz)       :: scalar_f                     ! Scalar field
  REAL(kind=rprec),DIMENSION(nxt+1,nynpy+1,nz)   :: scalar_x_f,scalar_y_f,scalar_z_f ! Interpolated scalar field
  REAL(kind=rprec),DIMENSION(ldx,nynpy+1,nz)     :: u_int_f,v_int_f,w_int_f      ! Interpolated velocity field
  REAL(kind=rprec),DIMENSION(nxt+1,nynpy+1,nz)   :: dscalar_dx,dscalar_dy,dscalar_dz ! Scalar derivatives

! Other variables
  REAL(kind=rprec),DIMENSION(ldx,nynpy+1,nz)     :: temp,temp_f                  ! Temporary to filter
  INTEGER                                    :: jx,jy,jz                     ! Counter


!
! BEGINNING CODE
!


  ! Check advection scheme
  IF (PCon_scheme/=3) THEN
    PRINT*,'Dynamic Cs2/Sc only available for SMART scheme...'
    STOP
  END IF


  ! Box filter velocity fields at test filter scale
  ! Note: these are the velocity fields already interpolated to surface locations
  temp=u_int_tmp
  CALL box_filter_3D(temp,temp_f,nxt+1,nynpy,nz-1)
  u_int_f=temp_f

  temp=v_int_tmp
  CALL box_filter_3D(temp,temp_f,nxt,nynpy+1,nz-1)
  v_int_f=temp_f

  temp=w_int_tmp
  CALL box_filter_3D(temp,temp_f,nxt,nynpy,nz-1)
  w_int_f=temp_f


  ! Box filter scalar fields at test filter scale
  ! Note: this is the scalar field already interpolated to surface locations
  temp=0._rprec
  temp(1:nxt+1,:,:)=scalar_x_tmp
  CALL box_filter_3D(temp,temp_f,nxt+1,nynpy,nz-1)
  scalar_x_f=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=scalar_y_tmp
  CALL box_filter_3D(temp,temp_f,nxt,nynpy+1,nz-1)
  scalar_y_f=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=scalar_z_tmp
  CALL box_filter_3D(temp,temp_f,nxt,nynpy,nz-1)
  scalar_z_f=temp_f(1:nxt+1,:,:)


  ! First use the vectors Kx, Ky and Kz to calculate the resolved fluxes
  Kx=u_int_tmp(1:nxt+1,1:nynpy+1,1:nz)*scalar_x_tmp
  Ky=v_int_tmp(1:nxt+1,1:nynpy+1,1:nz)*scalar_y_tmp
  Kz=w_int_tmp(1:nxt+1,1:nynpy+1,1:nz)*scalar_z_tmp


  ! Now filter the resolved fluxes
  temp=0._rprec
  temp(1:nxt+1,:,:)=Kx
  CALL box_filter_3D(temp,temp_f,nxt+1,nynpy,nz-1)
  Kx=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=Ky
  CALL box_filter_3D(temp,temp_f,nxt,nynpy+1,nz-1)
  Ky=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=Kz
  CALL box_filter_3D(temp,temp_f,nxt,nynpy,nz-1)
  Kz=temp_f(1:nxt+1,:,:)


  ! Subtract product of filtered fields to get "Leonard" fluxes
  Kx=Kx-u_int_f(1:nxt+1,1:nynpy+1,1:nz)*scalar_x_f
  Ky=Ky-v_int_f(1:nxt+1,1:nynpy+1,1:nz)*scalar_y_f
  Kz=Kz-w_int_f(1:nxt+1,1:nynpy+1,1:nz)*scalar_z_f


  ! Calculate derivatives
  dscalar_dx(2:nxt,1:nynpy,:)=(scalar(2:nxt,1:nynpy,:)-scalar(1:nxt-1,1:nynpy,:))/dx
  dscalar_dy(1:nxt,2:nynpy,:)=(scalar(1:nxt,2:nynpy,:)-scalar(1:nxt,1:nynpy-1,:))/dy
  dscalar_dz(1:nxt,1:nynpy,2:nz-1)=(scalar(1:nxt,1:nynpy,2:nz-1)-scalar(1:nxt,1:nynpy,1:nz-2))/dz


  ! |S| has to be interpolated (smae linear interpolation as magS in the pollen_RHS_calc2 routine
  ! For x-diffusion interpolate in y and z directions
  magSx_int(1:nxt,1:nynpy-1,1)=(magS(1:nxt,1:nynpy-1,1)+magS(1:nxt,2:nynpy,1))/2._rprec
  magSx_int(1:nxt,1:nynpy-1,2:nz-1)=(magS(1:nxt,1:nynpy-1,2:nz-1)+magS(1:nxt,2:nynpy,2:nz-1)+ &
     magS(1:nxt,1:nynpy-1,3:nz)+magS(1:nxt,2:nynpy,3:nz))/4._rprec
  magSx_int(1:nxt,nynpy,1)=(magS(1:nxt,nynpy,1)+magS(1:nxt,1,1))/2._rprec
  magSx_int(1:nxt,nynpy,2:nz-1)=(magS(1:nxt,nynpy,2:nz-1)+magS(1:nxt,1,2:nz-1)+ &
     magS(1:nxt,nynpy,3:nz)+magS(1:nxt,1,3:nz))/4._rprec
  magSx_int(nxt+1,1:nynpy,1:nz-1)=magSx_int(1,1:nynpy,1:nz-1)

  ! For y-diffusion interpolate in x and z directions
  magSy_int(1:nxt-1,1:nynpy,1)=(magS(1:nxt-1,1:nynpy,1)+magS(2:nxt,1:nynpy,1))/2._rprec
  magSy_int(1:nxt-1,1:nynpy,2:nz-1)=(magS(1:nxt-1,1:nynpy,2:nz-1)+magS(2:nxt,1:nynpy,2:nz-1)+ &
     magS(1:nxt-1,1:nynpy,3:nz)+magS(2:nxt,1:nynpy,3:nz))/4._rprec
  magSy_int(nxt,1:nynpy,1)=(magS(nxt,1:nynpy,1)+magS(1,1:nynpy,1))/2._rprec
  magSy_int(nxt,1:nynpy,2:nz-1)=(magS(nxt,1:nynpy,2:nz-1)+magS(1,1:nynpy,2:nz-1)+ &
     magS(nxt,1:nynpy,3:nz)+magS(1,1:nynpy,3:nz))/4._rprec
  magSy_int(1:nxt,nynpy+1,1:nz-1)=magSy_int(1:nxt,1,1:nz-1)

  ! For z-diffusion interpolate in x and y directions
  magSz_int(1:nxt-1,1:nynpy-1,1:nz-1)=(magS(1:nxt-1,1:nynpy-1,1:nz-1)+magS(2:nxt,1:nynpy-1,1:nz-1)+ &
     magS(1:nxt-1,2:nynpy,1:nz-1)+magS(2:nxt,2:nynpy,1:nz-1))/4._rprec
  magSz_int(nxt,1:nynpy-1,1:nz-1)=(magS(nxt,1:nynpy-1,1:nz-1)+magS(1,1:nynpy-1,1:nz-1)+ &
     magS(nxt,2:nynpy,1:nz-1)+magS(1,2:nynpy,1:nz-1))/4._rprec
  magSz_int(1:nxt-1,nynpy,1:nz-1)=(magS(1:nxt-1,nynpy,1:nz-1)+magS(2:nxt,nynpy,1:nz-1)+ &
     magS(1:nxt-1,1,1:nz-1)+magS(2:nxt,1,1:nz-1))/4._rprec
  magSz_int(nxt,nynpy,1:nz-1)=(magS(nxt,nynpy,1:nz-1)+magS(1,nynpy,1:nz-1)+ &
     magS(nxt,1,1:nz-1)+magS(1,1,1:nz-1))/4._rprec



  ! Modeled fluxes at scale delta (temporary storage in Xi vectors)
  Xx=dscalar_dx*magSx_int(1:nxt+1,1:nynpy+1,1:nz)
  Xy=dscalar_dy*magSy_int(1:nxt+1,1:nynpy+1,1:nz)
  Xz=dscalar_dz*magSz_int(1:nxt+1,1:nynpy+1,1:nz)


  ! Now filter the modeled fluxes
  temp=0._rprec
  temp(1:nxt+1,:,:)=Xx
  CALL box_filter_3D(temp,temp_f,nxt+1,nynpy,nz-1)
  Xx=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=Xy
  CALL box_filter_3D(temp,temp_f,nxt,nynpy+1,nz-1)
  Xy=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=Xz
  CALL box_filter_3D(temp,temp_f,nxt,nynpy,nz-1)
  Xz=temp_f(1:nxt+1,:,:)


  ! Box filter scalar field at test filter scale
  ! Note: this is the original scalar field
  temp=0._rprec
  temp(1:nxt,1:nynpy,:)=scalar(1:nxt,1:nynpy,:)
  CALL box_filter_3D(temp,temp_f,nxt,nynpy,nz-1)
  scalar_f(1:nxt,1:nynpy,:)=temp_f(1:nxt,1:nynpy,:)

  ! Calculate derivatives of filtered scalar
  dscalar_dx(2:nxt,1:nynpy,:)=(scalar_f(2:nxt,1:nynpy,:)-scalar_f(1:nxt-1,1:nynpy,:))/dx
  dscalar_dy(1:nxt,2:nynpy,:)=(scalar_f(1:nxt,2:nynpy,:)-scalar_f(1:nxt,1:nynpy-1,:))/dy
  dscalar_dz(1:nxt,1:nynpy,2:nz)=(scalar_f(1:nxt,1:nynpy,2:nz)-scalar_f(1:nxt,1:nynpy,1:nz-1))/dz


  ! Box filter |S| fields at test filter scale (replace original by filtered)
  ! Note: these are the |S| fields already interpolated to surface locations
  temp=magSx_int
  CALL box_filter_3D(temp,temp_f,nxt+1,nynpy,nz-1)
  magSx_int=temp_f

  temp=magSy_int
  CALL box_filter_3D(temp,temp_f,nxt,nynpy+1,nz-1)
  magSy_int=temp_f

  temp=magSz_int
  CALL box_filter_3D(temp,temp_f,nxt,nynpy,nz-1)
  magSz_int=temp_f

  ! Calculate final vector
  DO jz=1,nz
    Xx(:,:,jz)=Xx(:,:,jz)-4._rprec*magSx_int(1:nxt+1,1:nynpy+1,jz)*dscalar_dx(:,:,jz)
    Xy(:,:,jz)=Xy(:,:,jz)-4._rprec*magSy_int(1:nxt+1,1:nynpy+1,jz)*dscalar_dy(:,:,jz)
    Xz(:,:,jz)=Xz(:,:,jz)-4._rprec*magSz_int(1:nxt+1,1:nynpy+1,jz)*dscalar_dz(:,:,jz)
  END DO
  Xx=((dx*dy*dz)**(2._rprec/3._rprec))*Xx
  Xy=((dx*dy*dz)**(2._rprec/3._rprec))*Xy
  Xz=((dx*dy*dz)**(2._rprec/3._rprec))*Xz


  ! Take care of horizontal boundary problems
  ! Couldx modify this to take advantage of periodic BC when (periodicbc==.TRUE.)
  ! x-component
  Xx(1,:,:)=Xx(4,:,:); Xx(2,:,:)=Xx(4,:,:); Xx(3,:,:)=Xx(4,:,:);
  Xx(:,1,:)=Xx(:,4,:); Xx(:,2,:)=Xx(:,4,:); Xx(:,3,:)=Xx(:,4,:);
  Kx(1,:,:)=Kx(4,:,:); Kx(2,:,:)=Kx(4,:,:); Kx(3,:,:)=Kx(4,:,:);
  Kx(:,1,:)=Kx(:,4,:); Kx(:,2,:)=Kx(:,4,:); Kx(:,3,:)=Kx(:,4,:);
  Xx(nxt,:,:)=Xx(nxt-2,:,:);   Xx(nxt+1,:,:)=Xx(nxt-2,:,:); Xx(nxt-1,:,:)=Xx(nxt-2,:,:);
  Xx(:,nynpy-1,:)=Xx(:,nynpy-3,:); Xx(:,nynpy,:)=Xx(:,nynpy-3,:);   Xx(:,nynpy-2,:)=Xx(:,nynpy-3,:);
  Kx(nxt,:,:)=Kx(nxt-2,:,:);   Kx(nxt+1,:,:)=Kx(nxt-2,:,:); Kx(nxt-1,:,:)=Kx(nxt-2,:,:);
  Kx(:,nynpy-1,:)=Kx(:,nynpy-3,:); Kx(:,nynpy,:)=Kx(:,nynpy-3,:);   Kx(:,nynpy-2,:)=Kx(:,nynpy-3,:);

  ! y-component
  Xy(1,:,:)=Xy(4,:,:); Xy(2,:,:)=Xy(4,:,:); Xy(3,:,:)=Xy(4,:,:);
  Xy(:,1,:)=Xy(:,4,:); Xy(:,2,:)=Xy(:,4,:); Xy(:,3,:)=Xy(:,4,:);
  Ky(1,:,:)=Ky(4,:,:); Ky(2,:,:)=Ky(4,:,:); Ky(3,:,:)=Ky(4,:,:);
  Ky(:,1,:)=Ky(:,4,:); Ky(:,2,:)=Ky(:,4,:); Ky(:,3,:)=Ky(:,4,:);
  Xy(nxt-1,:,:)=Xy(nxt-3,:,:); Xy(nxt,:,:)=Xy(nxt-3,:,:);   Xy(nxt-2,:,:)=Xy(nxt-3,:,:);
  Xy(:,nynpy,:)=Xy(:,nynpy-2,:);   Xy(:,nynpy+1,:)=Xy(:,nynpy-2,:); Xy(:,nynpy-1,:)=Xy(:,nynpy-2,:);
  Ky(nxt-1,:,:)=Ky(nxt-3,:,:); Ky(nxt,:,:)=Ky(nxt-3,:,:);   Ky(nxt-2,:,:)=Ky(nxt-3,:,:);
  Ky(:,nynpy,:)=Ky(:,nynpy-2,:);   Ky(:,nynpy+1,:)=Ky(:,nynpy-2,:); Ky(:,nynpy-1,:)=Ky(:,nynpy-2,:);

  ! z-component
  Xz(1,:,:)=Xz(4,:,:); Xz(2,:,:)=Xz(4,:,:); Xz(3,:,:)=Xz(4,:,:);
  Xz(:,1,:)=Xz(:,4,:); Xz(:,2,:)=Xz(:,4,:); Xz(:,3,:)=Xz(:,4,:);
  Kz(1,:,:)=Kz(4,:,:); Kz(2,:,:)=Kz(4,:,:); Kz(3,:,:)=Kz(4,:,:);
  Kz(:,1,:)=Kz(:,4,:); Kz(:,2,:)=Kz(:,4,:); Kz(:,3,:)=Kz(:,4,:);
  Xz(nxt-1,:,:)=Xz(nxt-3,:,:); Xz(nxt,:,:)=Xz(nxt-3,:,:); Xz(nxt-2,:,:)=Xz(nxt-3,:,:);
  Xz(:,nynpy-1,:)=Xz(:,nynpy-3,:); Xz(:,nynpy,:)=Xz(:,nynpy-3,:); Xz(:,nynpy-2,:)=Xz(:,nynpy-3,:);
  Kz(nxt-1,:,:)=Kz(nxt-3,:,:); Kz(nxt,:,:)=Kz(nxt-3,:,:); Kz(nxt-2,:,:)=Kz(nxt-3,:,:);
  Kz(:,nynpy-1,:)=Kz(:,nynpy-3,:); Kz(:,nynpy,:)=Kz(:,nynpy-3,:); Kz(:,nynpy-2,:)=Kz(:,nynpy-3,:);


  ! Calculate contractions (these are for each element - averaging both faces in each direction)
  KiXi(1:nxt,1:nynpy,2:nz-1)=0.5_rprec*(Kx(1:nxt,1:nynpy,2:nz-1)*Xx(1:nxt,1:nynpy,2:nz-1)+Kx(2:nxt+1,1:nynpy,2:nz-1)*Xx(2:nxt+1,1:nynpy,2:nz-1) &
     +Ky(1:nxt,1:nynpy,2:nz-1)*Xy(1:nxt,1:nynpy,2:nz-1)+Ky(1:nxt,2:nynpy+1,2:nz-1)*Xy(1:nxt,2:nynpy+1,2:nz-1) &
     +Kz(1:nxt,1:nynpy,2:nz-1)*Xz(1:nxt,1:nynpy,2:nz-1)+Kz(1:nxt,1:nynpy,3:nz)*Xz(1:nxt,1:nynpy,3:nz))
  XiXi(1:nxt,1:nynpy,2:nz-1)=0.5_rprec*(Xx(1:nxt,1:nynpy,2:nz-1)*Xx(1:nxt,1:nynpy,2:nz-1)+Xx(2:nxt+1,1:nynpy,2:nz-1)*Xx(2:nxt+1,1:nynpy,2:nz-1) &
     +Xy(1:nxt,1:nynpy,2:nz-1)*Xy(1:nxt,1:nynpy,2:nz-1)+Xy(1:nxt,2:nynpy+1,2:nz-1)*Xy(1:nxt,2:nynpy+1,2:nz-1) &
     +Xz(1:nxt,1:nynpy,2:nz-1)*Xz(1:nxt,1:nynpy,2:nz-1)+Xz(1:nxt,1:nynpy,3:nz)*Xz(1:nxt,1:nynpy,3:nz))


  !
  ! Up to here, no changes from dynamic_sc
  ! Now the lagrangian stuff...
  !

  ! Initialization
  IF (jt==cs_count) THEN
    PRINT *,'F_XX and F_KX initialized'
    F_XX(:,:,:)=XiXi(:,:,:)
    F_KX(:,:,:)=0.03_rprec*XiXi(:,:,:)/Sc
  END IF

  ! Interpolate F_XX and F_KX to position xo=x-u*dt
  CALL interpolag_scalar()

  ! This is the discretization for the evolution equations of the two vector contractions
  ! Note: the factor epsilon (which is related to the time scale) is the same form the
  !       velocity field. So epsilon_lag is actually calculated in lagrange_Sdep.f90
  F_KX(1:nxt,:,:)=MAX((epsilon_lag(1:nxt,:,:)*KiXi(1:nxt,:,:)+(1._rprec-epsilon_lag(1:nxt,:,:))*F_KX(1:nxt,:,:)),0._rprec)
  F_XX(1:nxt,:,:)=(epsilon_lag(1:nxt,:,:)*XiXi(1:nxt,:,:)+(1._rprec-epsilon_lag(1:nxt,:,:))*F_XX(1:nxt,:,:))


  ! Linear interpolation for the lower boundary nodes
  ! KX and KX2 are iterpolated between their values at jz=4 and 0 at the surface
  F_KX(:,:,3)=MAX((5._rprec/7._rprec)*F_KX(:,:,4),0._rprec)
  F_KX(:,:,2)=MAX((3._rprec/7._rprec)*F_KX(:,:,4),0._rprec)
  F_KX(:,:,1)=MAX((1._rprec/7._rprec)*F_KX(:,:,4),0._rprec)
  ! XX and XX2 are assumed to be constant
  F_XX(:,:,3)=F_XX(:,:,4)
  F_XX(:,:,2)=F_XX(:,:,4)
  F_XX(:,:,1)=F_XX(:,:,4)


  ! Same for upper boundary
  F_KX(:,:,nz-2)=MAX((5._rprec/7._rprec)*F_KX(:,:,nz-3),0._rprec)
  F_KX(:,:,nz-1)=MAX((3._rprec/7._rprec)*F_KX(:,:,nz-3),0._rprec)
  F_KX(:,:,nz)=MAX((1._rprec/7._rprec)*F_KX(:,:,nz-3),0._rprec)
  F_XX(:,:,nz-2)=F_XX(:,:,nz-3)
  F_XX(:,:,nz-1)=F_XX(:,:,nz-3)
  F_XX(:,:,nz)=F_XX(:,:,nz-3)


  ! Calculate average Cs2/Sc
  DO jx=1,nxt
    DO jy=1,nynpy
      DO jz=1,nz

        IF (F_XX(jx,jy,jz)>0._rprec) THEN
          Cs2Sc(jx,jy,jz)=MAX(0._rprec,F_KX(jx,jy,jz)/F_XX(jx,jy,jz))
        ELSE
          Cs2Sc(jx,jy,jz)=0._rprec
        END IF

      END DO
    END DO
  END DO


END SUBROUTINE dynamic_sc_lasi





!
!    ################################################################################
!    ################################################################################
!    ##                                                                            ##
!    ##                            dynamic_sc_lasd.f90                             ##
!    ##                                                                            ##
!    ##                               Developed by                                 ##
!    ##                             Marcelo Chamecki                               ##
!    ##                                                                            ##
!    ##                        The Johns Hopkins University                        ##
!    ##                                                                            ##
!    ##                                05/27/2007                                  ##
!    ##                                                                            ##
!    ################################################################################
!    ################################################################################
!
!    ################################################################################
!
!    PURPOSE:  This routine is the plane averaged scale dependent dynamic model
!              for pollen concentration. It uses the procedure proposed by
!              Bou-Zeid et al. (PoF - 2005)
!
!    ################################################################################
!

SUBROUTINE dynamic_sc_lasd(scalar,scalar_x_tmp,scalar_y_tmp,scalar_z_tmp,u_int_tmp,v_int_tmp,w_int_tmp)

! Modules
  USE types,only:rprec
  USE param,only:ldx,nxt,nynpy,nz,dx,PCon_scheme, periodicbcx, periodicbcy
  USE scalars_module
  USE sgsmodule,only:magS,epsilon_lag,epsilon_lag2,F_KX,F_XX,F_KX2,F_XX2

! No implicit variables
  IMPLICIT NONE

! Main variables
  REAL(kind=rprec),DIMENSION(ldx,nynpy,nz)       :: scalar                      ! Scalar field
  REAL(kind=rprec),DIMENSION(ldx,nynpy+1,nz)     :: u_int_tmp,v_int_tmp,w_int_tmp      ! Interpolated velocity field
  REAL(kind=rprec),DIMENSION(ldx,nynpy+1,nz)     :: magSx_int,magSy_int,magSz_int! Interpolated |S|
  REAL(kind=rprec),DIMENSION(nxt+1,nynpy+1,nz)   :: scalar_x_tmp,scalar_y_tmp,scalar_z_tmp   ! Interpolated scalar
  REAL(kind=rprec),DIMENSION(nxt+1,nynpy+1,nz)   :: Kx,Ky,Kz     ! "Leonard" fluxes
  REAL(kind=rprec),DIMENSION(nxt+1,nynpy+1,nz)   :: Xx,Xy,Xz     ! "Leonard" modeled fluxes
  REAL(kind=rprec),DIMENSION(nxt,nynpy,nz)       :: KiXi,XiXi     ! Vector dot products (2*delta)
  REAL(kind=rprec),DIMENSION(nxt,nynpy,nz)       :: KiXi2,XiXi2     ! Vector dot products (4*delta)
  REAL(kind=rprec),DIMENSION(ldx,nynpy,nz)       :: beta_c       ! Beta_c

! Filtered fields
  REAL(kind=rprec),DIMENSION(nxt,nynpy,nz)       :: scalar_f                     ! Scalar field
  REAL(kind=rprec),DIMENSION(nxt+1,nynpy+1,nz)   :: scalar_x_f,scalar_y_f,scalar_z_f ! Interpolated scalar field
  REAL(kind=rprec),DIMENSION(ldx,nynpy+1,nz)     :: u_int_f,v_int_f,w_int_f      ! Interpolated velocity field
  REAL(kind=rprec),DIMENSION(nxt+1,nynpy+1,nz)   :: dscalar_dx,dscalar_dy,dscalar_dz ! Scalar derivatives

! Other variables
  REAL(kind=rprec),DIMENSION(ldx,nynpy+1,nz)     :: temp,temp_f                  ! Temporary to filter
  INTEGER                                    :: jx,jy,jz                     ! Counter


!
! BEGINNING CODE
!


  ! Check advection scheme
  IF (PCon_scheme/=3) THEN
    PRINT*,'Dynamic Sc only available for SMART scheme...'
    STOP
  END IF


  !
  ! First test filter scale - 2*delta
  !

  ! Box filter velocity fields at test filter scale
  ! Note: these are the velocity fields already interpolated to surface locations
  temp=u_int_tmp
  CALL box_filter_3D(temp,temp_f,nxt+1,nynpy,nz-1)
  u_int_f=temp_f

  temp=v_int_tmp
  CALL box_filter_3D(temp,temp_f,nxt,nynpy+1,nz-1)
  v_int_f=temp_f

  temp=w_int_tmp
  CALL box_filter_3D(temp,temp_f,nxt,nynpy,nz-1)
  w_int_f=temp_f


  ! Box filter scalar fields at test filter scale
  ! Note: this is the scalar field already interpolated to surface locations
  temp=0._rprec
  temp(1:nxt+1,:,:)=scalar_x_tmp
  CALL box_filter_3D(temp,temp_f,nxt+1,nynpy,nz-1)
  scalar_x_f=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=scalar_y_tmp
  CALL box_filter_3D(temp,temp_f,nxt,nynpy+1,nz-1)
  scalar_y_f=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=scalar_z_tmp
  CALL box_filter_3D(temp,temp_f,nxt,nynpy,nz-1)
  scalar_z_f=temp_f(1:nxt+1,:,:)


  ! First use the vectors Kx, Ky and Kz to calculate the resolved fluxes
  Kx=u_int_tmp(1:nxt+1,1:nynpy+1,1:nz)*scalar_x_tmp
  Ky=v_int_tmp(1:nxt+1,1:nynpy+1,1:nz)*scalar_y_tmp
  Kz=w_int_tmp(1:nxt+1,1:nynpy+1,1:nz)*scalar_z_tmp


  ! Now filter the resolved fluxes
  temp=0._rprec
  temp(1:nxt+1,:,:)=Kx
  CALL box_filter_3D(temp,temp_f,nxt+1,nynpy,nz-1)
  Kx=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=Ky
  CALL box_filter_3D(temp,temp_f,nxt,nynpy+1,nz-1)
  Ky=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=Kz
  CALL box_filter_3D(temp,temp_f,nxt,nynpy,nz-1)
  Kz=temp_f(1:nxt+1,:,:)


  ! Subtract product of filtered fields to get "Leonard" fluxes
  Kx=Kx-u_int_f(1:nxt+1,1:nynpy+1,1:nz)*scalar_x_f
  Ky=Ky-v_int_f(1:nxt+1,1:nynpy+1,1:nz)*scalar_y_f
  Kz=Kz-w_int_f(1:nxt+1,1:nynpy+1,1:nz)*scalar_z_f


  ! Calculate derivatives
  dscalar_dx(2:nxt,1:nynpy,:)=(scalar(2:nxt,1:nynpy,:)-scalar(1:nxt-1,1:nynpy,:))/dx
  dscalar_dy(1:nxt,2:nynpy,:)=(scalar(1:nxt,2:nynpy,:)-scalar(1:nxt,1:nynpy-1,:))/dy
  dscalar_dz(1:nxt,1:nynpy,2:nz-1)=(scalar(1:nxt,1:nynpy,2:nz-1)-scalar(1:nxt,1:nynpy,1:nz-2))/dz


  ! |S| has to be interpolated (same linear interpolation as magS in the pollen_RHS_calc2 routine)
  ! For x-diffusion interpolate in y and z directions
  magSx_int(1:nxt,1:nynpy-1,1)=(magS(1:nxt,1:nynpy-1,1)+magS(1:nxt,2:nynpy,1))/2._rprec
  magSx_int(1:nxt,1:nynpy-1,2:nz-1)=(magS(1:nxt,1:nynpy-1,2:nz-1)+magS(1:nxt,2:nynpy,2:nz-1)+ &
     magS(1:nxt,1:nynpy-1,3:nz)+magS(1:nxt,2:nynpy,3:nz))/4._rprec
  magSx_int(1:nxt,nynpy,1)=(magS(1:nxt,nynpy,1)+magS(1:nxt,1,1))/2._rprec
  magSx_int(1:nxt,nynpy,2:nz-1)=(magS(1:nxt,nynpy,2:nz-1)+magS(1:nxt,1,2:nz-1)+ &
     magS(1:nxt,nynpy,3:nz)+magS(1:nxt,1,3:nz))/4._rprec
  magSx_int(nxt+1,1:nynpy,1:nz-1)=magSx_int(1,1:nynpy,1:nz-1)

  ! For y-diffusion interpolate in x and z directions
  magSy_int(1:nxt-1,1:nynpy,1)=(magS(1:nxt-1,1:nynpy,1)+magS(2:nxt,1:nynpy,1))/2._rprec
  magSy_int(1:nxt-1,1:nynpy,2:nz-1)=(magS(1:nxt-1,1:nynpy,2:nz-1)+magS(2:nxt,1:nynpy,2:nz-1)+ &
     magS(1:nxt-1,1:nynpy,3:nz)+magS(2:nxt,1:nynpy,3:nz))/4._rprec
  magSy_int(nxt,1:nynpy,1)=(magS(nxt,1:nynpy,1)+magS(1,1:nynpy,1))/2._rprec
  magSy_int(nxt,1:nynpy,2:nz-1)=(magS(nxt,1:nynpy,2:nz-1)+magS(1,1:nynpy,2:nz-1)+ &
     magS(nxt,1:nynpy,3:nz)+magS(1,1:nynpy,3:nz))/4._rprec
  magSy_int(1:nxt,nynpy+1,1:nz-1)=magSy_int(1:nxt,1,1:nz-1)

  ! For z-diffusion interpolate in x and y directions
  magSz_int(1:nxt-1,1:nynpy-1,1:nz-1)=(magS(1:nxt-1,1:nynpy-1,1:nz-1)+magS(2:nxt,1:nynpy-1,1:nz-1)+ &
     magS(1:nxt-1,2:nynpy,1:nz-1)+magS(2:nxt,2:nynpy,1:nz-1))/4._rprec
  magSz_int(nxt,1:nynpy-1,1:nz-1)=(magS(nxt,1:nynpy-1,1:nz-1)+magS(1,1:nynpy-1,1:nz-1)+ &
     magS(nxt,2:nynpy,1:nz-1)+magS(1,2:nynpy,1:nz-1))/4._rprec
  magSz_int(1:nxt-1,nynpy,1:nz-1)=(magS(1:nxt-1,nynpy,1:nz-1)+magS(2:nxt,nynpy,1:nz-1)+ &
     magS(1:nxt-1,1,1:nz-1)+magS(2:nxt,1,1:nz-1))/4._rprec
  magSz_int(nxt,nynpy,1:nz-1)=(magS(nxt,nynpy,1:nz-1)+magS(1,nynpy,1:nz-1)+ &
     magS(nxt,1,1:nz-1)+magS(1,1,1:nz-1))/4._rprec



  ! Modeled fluxes at scale delta (temporary storage in Xi vectors)
  Xx=dscalar_dx*magSx_int(1:nxt+1,1:nynpy+1,1:nz)
  Xy=dscalar_dy*magSy_int(1:nxt+1,1:nynpy+1,1:nz)
  Xz=dscalar_dz*magSz_int(1:nxt+1,1:nynpy+1,1:nz)


  ! Now filter the modeled fluxes
  temp=0._rprec
  temp(1:nxt+1,:,:)=Xx
  CALL box_filter_3D(temp,temp_f,nxt+1,nynpy,nz-1)
  Xx=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=Xy
  CALL box_filter_3D(temp,temp_f,nxt,nynpy+1,nz-1)
  Xy=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=Xz
  CALL box_filter_3D(temp,temp_f,nxt,nynpy,nz-1)
  Xz=temp_f(1:nxt+1,:,:)


  ! Box filter scalar field at test filter scale
  ! Note: this is the original scalar field
  temp=0._rprec
  temp(1:nxt,1:nynpy,:)=scalar(1:nxt,1:nynpy,:)
  CALL box_filter_3D(temp,temp_f,nxt,nynpy,nz-1)
  scalar_f(1:nxt,1:nynpy,:)=temp_f(1:nxt,1:nynpy,:)

  ! Calculate derivatives of filtered scalar
  dscalar_dx(2:nxt,1:nynpy,:)=(scalar_f(2:nxt,1:nynpy,:)-scalar_f(1:nxt-1,1:nynpy,:))/dx
  dscalar_dy(1:nxt,2:nynpy,:)=(scalar_f(1:nxt,2:nynpy,:)-scalar_f(1:nxt,1:nynpy-1,:))/dy
  dscalar_dz(1:nxt,1:nynpy,2:nz)=(scalar_f(1:nxt,1:nynpy,2:nz)-scalar_f(1:nxt,1:nynpy,1:nz-1))/dz


  ! Box filter |S| fields at test filter scale (replace original by filtered)
  ! Note: these are the |S| fields already interpolated to surface locations
  temp=magSx_int
  CALL box_filter_3D(temp,temp_f,nxt+1,nynpy,nz-1)
  magSx_int=temp_f

  temp=magSy_int
  CALL box_filter_3D(temp,temp_f,nxt,nynpy+1,nz-1)
  magSy_int=temp_f

  temp=magSz_int
  CALL box_filter_3D(temp,temp_f,nxt,nynpy,nz-1)
  magSz_int=temp_f

  ! Calculate final vector
  DO jz=1,nz
    Xx(:,:,jz)=Xx(:,:,jz)-4._rprec*magSx_int(1:nxt+1,1:nynpy+1,jz)*dscalar_dx(:,:,jz)
    Xy(:,:,jz)=Xy(:,:,jz)-4._rprec*magSy_int(1:nxt+1,1:nynpy+1,jz)*dscalar_dy(:,:,jz)
    Xz(:,:,jz)=Xz(:,:,jz)-4._rprec*magSz_int(1:nxt+1,1:nynpy+1,jz)*dscalar_dz(:,:,jz)
  END DO
  Xx=((dx*dy*dz)**(2._rprec/3._rprec))*Xx
  Xy=((dx*dy*dz)**(2._rprec/3._rprec))*Xy
  Xz=((dx*dy*dz)**(2._rprec/3._rprec))*Xz


  ! Take care of horizontal boundary problems
  ! x-component
  Xx(1,:,:)=Xx(4,:,:); Xx(2,:,:)=Xx(4,:,:); Xx(3,:,:)=Xx(4,:,:);
  Xx(:,1,:)=Xx(:,4,:); Xx(:,2,:)=Xx(:,4,:); Xx(:,3,:)=Xx(:,4,:);
  Kx(1,:,:)=Kx(4,:,:); Kx(2,:,:)=Kx(4,:,:); Kx(3,:,:)=Kx(4,:,:);
  Kx(:,1,:)=Kx(:,4,:); Kx(:,2,:)=Kx(:,4,:); Kx(:,3,:)=Kx(:,4,:);
  Xx(nxt,:,:)=Xx(nxt-2,:,:);   Xx(nxt+1,:,:)=Xx(nxt-2,:,:); Xx(nxt-1,:,:)=Xx(nxt-2,:,:);
  Xx(:,nynpy-1,:)=Xx(:,nynpy-3,:); Xx(:,nynpy,:)=Xx(:,nynpy-3,:);   Xx(:,nynpy-2,:)=Xx(:,nynpy-3,:);
  Kx(nxt,:,:)=Kx(nxt-2,:,:);   Kx(nxt+1,:,:)=Kx(nxt-2,:,:); Kx(nxt-1,:,:)=Kx(nxt-2,:,:);
  Kx(:,nynpy-1,:)=Kx(:,nynpy-3,:); Kx(:,nynpy,:)=Kx(:,nynpy-3,:);   Kx(:,nynpy-2,:)=Kx(:,nynpy-3,:);

  ! y-component
  Xy(1,:,:)=Xy(4,:,:); Xy(2,:,:)=Xy(4,:,:); Xy(3,:,:)=Xy(4,:,:);
  Xy(:,1,:)=Xy(:,4,:); Xy(:,2,:)=Xy(:,4,:); Xy(:,3,:)=Xy(:,4,:);
  Ky(1,:,:)=Ky(4,:,:); Ky(2,:,:)=Ky(4,:,:); Ky(3,:,:)=Ky(4,:,:);
  Ky(:,1,:)=Ky(:,4,:); Ky(:,2,:)=Ky(:,4,:); Ky(:,3,:)=Ky(:,4,:);
  Xy(nxt-1,:,:)=Xy(nxt-3,:,:); Xy(nxt,:,:)=Xy(nxt-3,:,:);   Xy(nxt-2,:,:)=Xy(nxt-3,:,:);
  Xy(:,nynpy,:)=Xy(:,nynpy-2,:);   Xy(:,nynpy+1,:)=Xy(:,nynpy-2,:); Xy(:,nynpy-1,:)=Xy(:,nynpy-2,:);
  Ky(nxt-1,:,:)=Ky(nxt-3,:,:); Ky(nxt,:,:)=Ky(nxt-3,:,:);   Ky(nxt-2,:,:)=Ky(nxt-3,:,:);
  Ky(:,nynpy,:)=Ky(:,nynpy-2,:);   Ky(:,nynpy+1,:)=Ky(:,nynpy-2,:); Ky(:,nynpy-1,:)=Ky(:,nynpy-2,:);

  ! z-component
  Xz(1,:,:)=Xz(4,:,:); Xz(2,:,:)=Xz(4,:,:); Xz(3,:,:)=Xz(4,:,:);
  Xz(:,1,:)=Xz(:,4,:); Xz(:,2,:)=Xz(:,4,:); Xz(:,3,:)=Xz(:,4,:);
  Kz(1,:,:)=Kz(4,:,:); Kz(2,:,:)=Kz(4,:,:); Kz(3,:,:)=Kz(4,:,:);
  Kz(:,1,:)=Kz(:,4,:); Kz(:,2,:)=Kz(:,4,:); Kz(:,3,:)=Kz(:,4,:);
  Xz(nxt-1,:,:)=Xz(nxt-3,:,:); Xz(nxt,:,:)=Xz(nxt-3,:,:); Xz(nxt-2,:,:)=Xz(nxt-3,:,:);
  Xz(:,nynpy-1,:)=Xz(:,nynpy-3,:); Xz(:,nynpy,:)=Xz(:,nynpy-3,:); Xz(:,nynpy-2,:)=Xz(:,nynpy-3,:);
  Kz(nxt-1,:,:)=Kz(nxt-3,:,:); Kz(nxt,:,:)=Kz(nxt-3,:,:); Kz(nxt-2,:,:)=Kz(nxt-3,:,:);
  Kz(:,nynpy-1,:)=Kz(:,nynpy-3,:); Kz(:,nynpy,:)=Kz(:,nynpy-3,:); Kz(:,nynpy-2,:)=Kz(:,nynpy-3,:);


  ! Calculate contractions (these are for each element - averaging both faces in each direction)
  KiXi(1:nxt,1:nynpy,1:nz-1)=0.5_rprec*(Kx(1:nxt,1:nynpy,1:nz-1)*Xx(1:nxt,1:nynpy,1:nz-1)+Kx(2:nxt+1,1:nynpy,1:nz-1)*Xx(2:nxt+1,1:nynpy,1:nz-1) &
     +Ky(1:nxt,1:nynpy,1:nz-1)*Xy(1:nxt,1:nynpy,1:nz-1)+Ky(1:nxt,2:nynpy+1,1:nz-1)*Xy(1:nxt,2:nynpy+1,1:nz-1) &
     +Kz(1:nxt,1:nynpy,1:nz-1)*Xz(1:nxt,1:nynpy,1:nz-1)+Kz(1:nxt,1:nynpy,2:nz)*Xz(1:nxt,1:nynpy,2:nz))
  XiXi(1:nxt,1:nynpy,1:nz-1)=0.5_rprec*(Xx(1:nxt,1:nynpy,1:nz-1)*Xx(1:nxt,1:nynpy,1:nz-1)+Xx(2:nxt+1,1:nynpy,1:nz-1)*Xx(2:nxt+1,1:nynpy,1:nz-1) &
     +Xy(1:nxt,1:nynpy,1:nz-1)*Xy(1:nxt,1:nynpy,1:nz-1)+Xy(1:nxt,2:nynpy+1,1:nz-1)*Xy(1:nxt,2:nynpy+1,1:nz-1) &
     +Xz(1:nxt,1:nynpy,1:nz-1)*Xz(1:nxt,1:nynpy,1:nz-1)+Xz(1:nxt,1:nynpy,2:nz)*Xz(1:nxt,1:nynpy,2:nz))




  !
  ! Second test filter scale - 4*delta
  !


  ! Box filter velocity fields at second test filter scale
  ! Note: these are the velocity fields already interpolated to surface locations
  temp=u_int_tmp
  CALL box_filter_3D_4delta(temp,temp_f,nxt+1,nynpy,nz-1)
  u_int_f=temp_f

  temp=v_int_tmp
  CALL box_filter_3D_4delta(temp,temp_f,nxt,nynpy+1,nz-1)
  v_int_f=temp_f

  temp=w_int_tmp
  CALL box_filter_3D_4delta(temp,temp_f,nxt,nynpy,nz-1)
  w_int_f=temp_f


  ! Box filter scalar fields at second test filter scale
  ! Note: this is the scalar field already interpolated to surface locations
  temp=0._rprec
  temp(1:nxt+1,:,:)=scalar_x_tmp
  CALL box_filter_3D_4delta(temp,temp_f,nxt+1,nynpy,nz-1)
  scalar_x_f=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=scalar_y_tmp
  CALL box_filter_3D_4delta(temp,temp_f,nxt,nynpy+1,nz-1)
  scalar_y_f=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=scalar_z_tmp
  CALL box_filter_3D_4delta(temp,temp_f,nxt,nynpy,nz-1)
  scalar_z_f=temp_f(1:nxt+1,:,:)


  ! First use the vectors Kx, Ky and Kz to calculate the resolved fluxes
  Kx=u_int_tmp(1:nxt+1,1:nynpy+1,1:nz)*scalar_x_tmp
  Ky=v_int_tmp(1:nxt+1,1:nynpy+1,1:nz)*scalar_y_tmp
  Kz=w_int_tmp(1:nxt+1,1:nynpy+1,1:nz)*scalar_z_tmp


  ! Now filter the resolved fluxes
  temp=0._rprec
  temp(1:nxt+1,:,:)=Kx
  CALL box_filter_3D_4delta(temp,temp_f,nxt+1,nynpy,nz-1)
  Kx=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=Ky
  CALL box_filter_3D_4delta(temp,temp_f,nxt,nynpy+1,nz-1)
  Ky=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=Kz
  CALL box_filter_3D_4delta(temp,temp_f,nxt,nynpy,nz-1)
  Kz=temp_f(1:nxt+1,:,:)


  ! Subtract product of filtered fields to get "Leonard" fluxes
  Kx=Kx-u_int_f(1:nxt+1,1:nynpy+1,1:nz)*scalar_x_f
  Ky=Ky-v_int_f(1:nxt+1,1:nynpy+1,1:nz)*scalar_y_f
  Kz=Kz-w_int_f(1:nxt+1,1:nynpy+1,1:nz)*scalar_z_f


  ! Calculate derivatives
  dscalar_dx(2:nxt,1:nynpy,:)=(scalar(2:nxt,1:nynpy,:)-scalar(1:nxt-1,1:nynpy,:))/dx
  dscalar_dy(1:nxt,2:nynpy,:)=(scalar(1:nxt,2:nynpy,:)-scalar(1:nxt,1:nynpy-1,:))/dy
  dscalar_dz(1:nxt,1:nynpy,2:nz-1)=(scalar(1:nxt,1:nynpy,2:nz-1)-scalar(1:nxt,1:nynpy,1:nz-2))/dz


  ! |S| has to be interpolated (smae linear interpolation as magS in the pollen_RHS_calc2 routine
  ! For x-diffusion interpolate in y and z directions
  magSx_int(1:nxt,1:nynpy-1,1)=(magS(1:nxt,1:nynpy-1,1)+magS(1:nxt,2:nynpy,1))/2._rprec
  magSx_int(1:nxt,1:nynpy-1,2:nz-1)=(magS(1:nxt,1:nynpy-1,2:nz-1)+magS(1:nxt,2:nynpy,2:nz-1)+ &
     magS(1:nxt,1:nynpy-1,3:nz)+magS(1:nxt,2:nynpy,3:nz))/4._rprec
  magSx_int(1:nxt,nynpy,1)=(magS(1:nxt,nynpy,1)+magS(1:nxt,1,1))/2._rprec
  magSx_int(1:nxt,nynpy,2:nz-1)=(magS(1:nxt,nynpy,2:nz-1)+magS(1:nxt,1,2:nz-1)+ &
     magS(1:nxt,nynpy,3:nz)+magS(1:nxt,1,3:nz))/4._rprec
  magSx_int(nxt+1,1:nynpy,1:nz-1)=magSx_int(1,1:nynpy,1:nz-1)

  ! For y-diffusion interpolate in x and z directions
  magSy_int(1:nxt-1,1:nynpy,1)=(magS(1:nxt-1,1:nynpy,1)+magS(2:nxt,1:nynpy,1))/2._rprec
  magSy_int(1:nxt-1,1:nynpy,2:nz-1)=(magS(1:nxt-1,1:nynpy,2:nz-1)+magS(2:nxt,1:nynpy,2:nz-1)+ &
     magS(1:nxt-1,1:nynpy,3:nz)+magS(2:nxt,1:nynpy,3:nz))/4._rprec
  magSy_int(nxt,1:nynpy,1)=(magS(nxt,1:nynpy,1)+magS(1,1:nynpy,1))/2._rprec
  magSy_int(nxt,1:nynpy,2:nz-1)=(magS(nxt,1:nynpy,2:nz-1)+magS(1,1:nynpy,2:nz-1)+ &
     magS(nxt,1:nynpy,3:nz)+magS(1,1:nynpy,3:nz))/4._rprec
  magSy_int(1:nxt,nynpy+1,1:nz-1)=magSy_int(1:nxt,1,1:nz-1)

  ! For z-diffusion interpolate in x and y directions
  magSz_int(1:nxt-1,1:nynpy-1,1:nz-1)=(magS(1:nxt-1,1:nynpy-1,1:nz-1)+magS(2:nxt,1:nynpy-1,1:nz-1)+ &
     magS(1:nxt-1,2:nynpy,1:nz-1)+magS(2:nxt,2:nynpy,1:nz-1))/4._rprec
  magSz_int(nxt,1:nynpy-1,1:nz-1)=(magS(nxt,1:nynpy-1,1:nz-1)+magS(1,1:nynpy-1,1:nz-1)+ &
     magS(nxt,2:nynpy,1:nz-1)+magS(1,2:nynpy,1:nz-1))/4._rprec
  magSz_int(1:nxt-1,nynpy,1:nz-1)=(magS(1:nxt-1,nynpy,1:nz-1)+magS(2:nxt,nynpy,1:nz-1)+ &
     magS(1:nxt-1,1,1:nz-1)+magS(2:nxt,1,1:nz-1))/4._rprec
  magSz_int(nxt,nynpy,1:nz-1)=(magS(nxt,nynpy,1:nz-1)+magS(1,nynpy,1:nz-1)+ &
     magS(nxt,1,1:nz-1)+magS(1,1,1:nz-1))/4._rprec



  ! Modeled fluxes at scale delta (temporary storage in Xi vectors)
  Xx=dscalar_dx*magSx_int(1:nxt+1,1:nynpy+1,1:nz)
  Xy=dscalar_dy*magSy_int(1:nxt+1,1:nynpy+1,1:nz)
  Xz=dscalar_dz*magSz_int(1:nxt+1,1:nynpy+1,1:nz)


  ! Now filter the modeled fluxes
  temp=0._rprec
  temp(1:nxt+1,:,:)=Xx
  CALL box_filter_3D_4delta(temp,temp_f,nxt+1,nynpy,nz-1)
  Xx=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=Xy
  CALL box_filter_3D_4delta(temp,temp_f,nxt,nynpy+1,nz-1)
  Xy=temp_f(1:nxt+1,:,:)

  temp=0._rprec
  temp(1:nxt+1,:,:)=Xz
  CALL box_filter_3D_4delta(temp,temp_f,nxt,nynpy,nz-1)
  Xz=temp_f(1:nxt+1,:,:)


  ! Box filter scalar field at test filter scale
  ! Note: this is the original scalar field
  temp=0._rprec
  temp(1:nxt,1:nynpy,:)=scalar(1:nxt,1:nynpy,:)
  CALL box_filter_3D_4delta(temp,temp_f,nxt,nynpy,nz-1)
  scalar_f(1:nxt,1:nynpy,:)=temp_f(1:nxt,1:nynpy,:)

  ! Calculate derivatives of filtered scalar
  dscalar_dx(2:nxt,1:nynpy,:)=(scalar_f(2:nxt,1:nynpy,:)-scalar_f(1:nxt-1,1:nynpy,:))/dx
  dscalar_dy(1:nxt,2:nynpy,:)=(scalar_f(1:nxt,2:nynpy,:)-scalar_f(1:nxt,1:nynpy-1,:))/dy
  dscalar_dz(1:nxt,1:nynpy,2:nz)=(scalar_f(1:nxt,1:nynpy,2:nz)-scalar_f(1:nxt,1:nynpy,1:nz-1))/dz


  ! Box filter |S| fields at test filter scale (replace original by filtered)
  ! Note: these are the |S| fields already interpolated to surface locations
  temp=magSx_int
  CALL box_filter_3D_4delta(temp,temp_f,nxt+1,nynpy,nz-1)
  magSx_int=temp_f

  temp=magSy_int
  CALL box_filter_3D_4delta(temp,temp_f,nxt,nynpy+1,nz-1)
  magSy_int=temp_f

  temp=magSz_int
  CALL box_filter_3D_4delta(temp,temp_f,nxt,nynpy,nz-1)
  magSz_int=temp_f

  ! Calculate final vector
  DO jz=1,nz
    Xx(:,:,jz)=Xx(:,:,jz)-16._rprec*magSx_int(1:nxt+1,1:nynpy+1,jz)*dscalar_dx(:,:,jz)
    Xy(:,:,jz)=Xy(:,:,jz)-16._rprec*magSy_int(1:nxt+1,1:nynpy+1,jz)*dscalar_dy(:,:,jz)
    Xz(:,:,jz)=Xz(:,:,jz)-16._rprec*magSz_int(1:nxt+1,1:nynpy+1,jz)*dscalar_dz(:,:,jz)
  END DO
  Xx=((dx*dy*dz)**(2._rprec/3._rprec))*Xx
  Xy=((dx*dy*dz)**(2._rprec/3._rprec))*Xy
  Xz=((dx*dy*dz)**(2._rprec/3._rprec))*Xz


  ! Take care of boundary problems
  ! x-component
  Xx(1,:,:)=Xx(4,:,:); Xx(2,:,:)=Xx(4,:,:); Xx(3,:,:)=Xx(4,:,:);
  Xx(:,1,:)=Xx(:,4,:); Xx(:,2,:)=Xx(:,4,:); Xx(:,3,:)=Xx(:,4,:);
  Kx(1,:,:)=Kx(4,:,:); Kx(2,:,:)=Kx(4,:,:); Kx(3,:,:)=Kx(4,:,:);
  Kx(:,1,:)=Kx(:,4,:); Kx(:,2,:)=Kx(:,4,:); Kx(:,3,:)=Kx(:,4,:);
  Xx(nxt,:,:)=Xx(nxt-2,:,:);   Xx(nxt+1,:,:)=Xx(nxt-2,:,:); Xx(nxt-1,:,:)=Xx(nxt-2,:,:);
  Xx(:,nynpy-1,:)=Xx(:,nynpy-3,:); Xx(:,nynpy,:)=Xx(:,nynpy-3,:);   Xx(:,nynpy-2,:)=Xx(:,nynpy-3,:);
  Kx(nxt,:,:)=Kx(nxt-2,:,:);   Kx(nxt+1,:,:)=Kx(nxt-2,:,:); Kx(nxt-1,:,:)=Kx(nxt-2,:,:);
  Kx(:,nynpy-1,:)=Kx(:,nynpy-3,:); Kx(:,nynpy,:)=Kx(:,nynpy-3,:);   Kx(:,nynpy-2,:)=Kx(:,nynpy-3,:);

  ! y-component
  Xy(1,:,:)=Xy(4,:,:); Xy(2,:,:)=Xy(4,:,:); Xy(3,:,:)=Xy(4,:,:);
  Xy(:,1,:)=Xy(:,4,:); Xy(:,2,:)=Xy(:,4,:); Xy(:,3,:)=Xy(:,4,:);
  Ky(1,:,:)=Ky(4,:,:); Ky(2,:,:)=Ky(4,:,:); Ky(3,:,:)=Ky(4,:,:);
  Ky(:,1,:)=Ky(:,4,:); Ky(:,2,:)=Ky(:,4,:); Ky(:,3,:)=Ky(:,4,:);
  Xy(nxt-1,:,:)=Xy(nxt-3,:,:); Xy(nxt,:,:)=Xy(nxt-3,:,:);   Xy(nxt-2,:,:)=Xy(nxt-3,:,:);
  Xy(:,nynpy,:)=Xy(:,nynpy-2,:);   Xy(:,nynpy+1,:)=Xy(:,nynpy-2,:); Xy(:,nynpy-1,:)=Xy(:,nynpy-2,:);
  Ky(nxt-1,:,:)=Ky(nxt-3,:,:); Ky(nxt,:,:)=Ky(nxt-3,:,:);   Ky(nxt-2,:,:)=Ky(nxt-3,:,:);
  Ky(:,nynpy,:)=Ky(:,nynpy-2,:);   Ky(:,nynpy+1,:)=Ky(:,nynpy-2,:); Ky(:,nynpy-1,:)=Ky(:,nynpy-2,:);

  ! z-component
  Xz(1,:,:)=Xz(4,:,:); Xz(2,:,:)=Xz(4,:,:); Xz(3,:,:)=Xz(4,:,:);
  Xz(:,1,:)=Xz(:,4,:); Xz(:,2,:)=Xz(:,4,:); Xz(:,3,:)=Xz(:,4,:);
  Kz(1,:,:)=Kz(4,:,:); Kz(2,:,:)=Kz(4,:,:); Kz(3,:,:)=Kz(4,:,:);
  Kz(:,1,:)=Kz(:,4,:); Kz(:,2,:)=Kz(:,4,:); Kz(:,3,:)=Kz(:,4,:);
  Xz(nxt-1,:,:)=Xz(nxt-3,:,:); Xz(nxt,:,:)=Xz(nxt-3,:,:); Xz(nxt-2,:,:)=Xz(nxt-3,:,:);
  Xz(:,nynpy-1,:)=Xz(:,nynpy-3,:); Xz(:,nynpy,:)=Xz(:,nynpy-3,:); Xz(:,nynpy-2,:)=Xz(:,nynpy-3,:);
  Kz(nxt-1,:,:)=Kz(nxt-3,:,:); Kz(nxt,:,:)=Kz(nxt-3,:,:); Kz(nxt-2,:,:)=Kz(nxt-3,:,:);
  Kz(:,nynpy-1,:)=Kz(:,nynpy-3,:); Kz(:,nynpy,:)=Kz(:,nynpy-3,:); Kz(:,nynpy-2,:)=Kz(:,nynpy-3,:);


  ! Calculate contractions (these are for each element - averaging both faces in each direction)
  KiXi2(1:nxt,1:nynpy,1:nz-1)=0.5_rprec*(Kx(1:nxt,1:nynpy,1:nz-1)*Xx(1:nxt,1:nynpy,1:nz-1)+Kx(2:nxt+1,1:nynpy,1:nz-1)*Xx(2:nxt+1,1:nynpy,1:nz-1) &
     +Ky(1:nxt,1:nynpy,1:nz-1)*Xy(1:nxt,1:nynpy,1:nz-1)+Ky(1:nxt,2:nynpy+1,1:nz-1)*Xy(1:nxt,2:nynpy+1,1:nz-1) &
     +Kz(1:nxt,1:nynpy,1:nz-1)*Xz(1:nxt,1:nynpy,1:nz-1)+Kz(1:nxt,1:nynpy,2:nz)*Xz(1:nxt,1:nynpy,2:nz))
  XiXi2(1:nxt,1:nynpy,1:nz-1)=0.5_rprec*(Xx(1:nxt,1:nynpy,1:nz-1)*Xx(1:nxt,1:nynpy,1:nz-1)+Xx(2:nxt+1,1:nynpy,1:nz-1)*Xx(2:nxt+1,1:nynpy,1:nz-1) &
     +Xy(1:nxt,1:nynpy,1:nz-1)*Xy(1:nxt,1:nynpy,1:nz-1)+Xy(1:nxt,2:nynpy+1,1:nz-1)*Xy(1:nxt,2:nynpy+1,1:nz-1) &
     +Xz(1:nxt,1:nynpy,1:nz-1)*Xz(1:nxt,1:nynpy,1:nz-1)+Xz(1:nxt,1:nynpy,2:nz)*Xz(1:nxt,1:nynpy,2:nz))





  !
  ! Up to here, no changes from dynamic_sc
  ! Now the lagrangian stuff...
  !

  ! Initialization
  IF (jt==cs_count) THEN
    PRINT *,'F_XX and F_KX initialized'
    F_XX(:,:,:)=XiXi(:,:,:)
    F_KX(:,:,:)=0.03_rprec*XiXi(:,:,:)/Sc
    F_XX2(:,:,:)=XiXi2(:,:,:)
    F_KX2(:,:,:)=0.03_rprec*XiXi2(:,:,:)/Sc
  END IF

  ! Interpolate F_XX, F_KX, F_XX2 and F_KX2 to position xo=x-u*dt
  CALL interpolag_scalar_Sdep()


  ! This is the discretization for the evolution equations of the two vector contractions
  ! Note: the factor epsilon (which is related to the time scale) is the same form the
  !       velocity field. So epsilon_lag is actually calculated in lagrange_Sdep.f90
  F_KX(1:nxt,:,:)=MAX((epsilon_lag(1:nxt,:,:)*KiXi(1:nxt,:,:)+(1._rprec-epsilon_lag(1:nxt,:,:))*F_KX(1:nxt,:,:)),0._rprec)
  F_XX(1:nxt,:,:)=(epsilon_lag(1:nxt,:,:)*XiXi(1:nxt,:,:)+(1._rprec-epsilon_lag(1:nxt,:,:))*F_XX(1:nxt,:,:))

  F_KX2(1:nxt,:,:)=MAX((epsilon_lag2(1:nxt,:,:)*KiXi2(1:nxt,:,:)+(1._rprec-epsilon_lag2(1:nxt,:,:))*F_KX2(1:nxt,:,:)),0._rprec)
  F_XX2(1:nxt,:,:)=(epsilon_lag2(1:nxt,:,:)*XiXi2(1:nxt,:,:)+(1._rprec-epsilon_lag2(1:nxt,:,:))*F_XX2(1:nxt,:,:))



  ! Linear interpolation for the lower boundary nodes
  ! KX and KX2 are iterpolated between their values at jz=4 and 0 at the surface
  F_KX(:,:,3)=MAX((5._rprec/7._rprec)*F_KX(:,:,4),0._rprec)
  F_KX(:,:,2)=MAX((3._rprec/7._rprec)*F_KX(:,:,4),0._rprec)
  F_KX(:,:,1)=MAX((1._rprec/7._rprec)*F_KX(:,:,4),0._rprec)
  F_KX2(:,:,3)=MAX((5._rprec/7._rprec)*F_KX2(:,:,4),0._rprec)
  F_KX2(:,:,2)=MAX((3._rprec/7._rprec)*F_KX2(:,:,4),0._rprec)
  F_KX2(:,:,1)=MAX((1._rprec/7._rprec)*F_KX2(:,:,4),0._rprec)
  ! XX and XX2 are assumed to be constant
  F_XX(:,:,3)=F_XX(:,:,4)
  F_XX(:,:,2)=F_XX(:,:,4)
  F_XX(:,:,1)=F_XX(:,:,4)
  F_XX2(:,:,3)=F_XX2(:,:,4)
  F_XX2(:,:,2)=F_XX2(:,:,4)
  F_XX2(:,:,1)=F_XX2(:,:,4)


  ! Same for upper boundary
  F_KX(:,:,nz-2)=MAX((5._rprec/7._rprec)*F_KX(:,:,nz-3),0._rprec)
  F_KX(:,:,nz-1)=MAX((3._rprec/7._rprec)*F_KX(:,:,nz-3),0._rprec)
  F_KX(:,:,nz)=MAX((1._rprec/7._rprec)*F_KX(:,:,nz-3),0._rprec)
  F_KX2(:,:,nz-2)=MAX((5._rprec/7._rprec)*F_KX2(:,:,nz-3),0._rprec)
  F_KX2(:,:,nz-1)=MAX((3._rprec/7._rprec)*F_KX2(:,:,nz-3),0._rprec)
  F_KX2(:,:,nz)=MAX((1._rprec/7._rprec)*F_KX2(:,:,nz-3),0._rprec)
  F_XX(:,:,nz-2)=F_XX(:,:,nz-3)
  F_XX(:,:,nz-1)=F_XX(:,:,nz-3)
  F_XX(:,:,nz)=F_XX(:,:,nz-3)
  F_XX2(:,:,nz-2)=F_XX2(:,:,nz-3)
  F_XX2(:,:,nz-1)=F_XX2(:,:,nz-3)
  F_XX2(:,:,nz)=F_XX2(:,:,nz-3)



  !
  ! Determine the values of beta_c and Cs2/Sc
  !


  DO jx=1,nxt
    DO jy=1,nynpy
      DO jz=1,nz

        ! Calculate beta_c
        IF (F_XX2(jx,jy,jz)>0._rprec .AND. F_KX(jx,jy,jz)>0._rprec) THEN
          beta_c(jx,jy,jz)=MAX(0.125_rprec,(F_KX2(jx,jy,jz)*F_XX(jx,jy,jz))/(F_KX(jx,jy,jz)*F_XX2(jx,jy,jz)))
        ELSE
          beta_c(jx,jy,jz)=0.125_rprec
        END IF

        ! Calculate Cs2/Sc
        IF (F_XX(jx,jy,jz)>0._rprec) THEN
          Cs2Sc(jx,jy,jz)=MAX(0._rprec,(F_KX(jx,jy,jz)/F_XX(jx,jy,jz))/beta_c(jx,jy,jz))
        ELSE
          Cs2Sc(jx,jy,jz)=0._rprec
        END IF

      END DO
    END DO
  END DO


!  ! Calculate Cs2/Sc(4*delta)
!  WHERE (F_XX2>0._rprec)
!    Cs2Sc_4d(1:nxt,:,:)=MAX(0._rprec,F_KX2/F_XX2)
!  ELSEWHERE
!    Cs2Sc_4d(1:nxt,:,:)=0._rprec
!  END WHERE

!  ! Calculate Cs2/Sc(2*delta)
!  WHERE (F_XX>0._rprec)
!    Cs2Sc_2d(1:nxt,:,:)=MAX(0._rprec,F_KX/F_XX)
!  ELSEWHERE
!    Cs2Sc_2d(1:nxt,:,:)=0._rprec
!  END WHERE

!  ! Calculate beta_c
!  WHERE (F_XX2>0._rprec .AND. F_KX>0._rprec)
!    beta_c(1:nxt,:,:)=MAX(0.125_rprec,(F_KX2*F_XX)/(F_KX*F_XX2))
!  ELSEWHERE
!    beta_c(1:nxt,:,:)=0.125_rprec
!  END WHERE


!  ! Calculate Cs2/Sc
!  WHERE (F_XX>0._rprec)
!    Cs2Sc(1:nxt,:,:)=MAX(0._rprec,(F_KX/F_XX)/beta_c(1:nxt,:,:))
!  ELSEWHERE
!    Cs2Sc(1:nxt,:,:)=0._rprec
!  END WHERE



END SUBROUTINE dynamic_sc_lasd






! This is a modified version of interpolag_Sdep.f90

! This subroutine computes the values of F_KX and F_XX
! at positions (x-u*dt) for use in the subroutine lagrangian

SUBROUTINE interpolag_scalar()

  USE types,only:rprec
  USE param
  USE sgsmodule
  use intermediate, only: FF_KX, FF_XX

  IMPLICIT NONE
  INTEGER :: jx, jy, jz
  INTEGER :: jjx, jjy, jjz

  !--saves here: force heap storage
  REAL(kind=rprec) :: frac_x,frac_y,frac_z
  REAL(kind=rprec) :: comp_x,comp_y,comp_z
  !REAL(kind=rprec), SAVE, DIMENSION(nxt+2,nynpy+2,nz+2) :: FF_KX,FF_XX

  INTEGER :: addx, addy, addz
  INTEGER :: jxaddx, jyaddy, jzaddz



  ! creates dummy arrays FF_KX and FF_XX to use in the subroutine

  DO jz=1,nz
    DO jy=1,nynpy
      DO jx=1,nxt
        FF_KX(jx+1,jy+1,jz+1) = F_KX(jx,jy,jz)
        FF_XX(jx+1,jy+1,jz+1) = F_XX(jx,jy,jz)
      END DO
    END DO
  END DO

  !This is a bit like witch craft but the following lines do take care
  !of all the edges including the corners
  FF_KX(1,:,:) = FF_KX(nxt+1,:,:)
  FF_KX(nxt+2,:,:) = FF_KX(2,:,:)

  FF_KX(:,1,:) = FF_KX(:,nynpy+1,:)
  FF_KX(:,nynpy+2,:) = FF_KX(:,2,:)

  FF_KX(:,:,1) = FF_KX(:,:,2)
  FF_KX(:,:,nz+2) = FF_KX(:,:,nz+1)

  FF_XX(1,:,:) = FF_XX(nxt+1,:,:)
  FF_XX(nxt+2,:,:) = FF_XX(2,:,:)

  FF_XX(:,1,:) = FF_XX(:,nynpy+1,:)
  FF_XX(:,nynpy+2,:) = FF_XX(:,2,:)

  FF_XX(:,:,1) = FF_XX(:,:,2)
  FF_XX(:,:,nz+2) = FF_XX(:,:,nz+1)

  ! end of witch craft


  DO jz=1,nz
    jjz = jz+1
    DO jy=1,nynpy
      jjy = jy+1
      DO jx=1,nxt
        jjx = jx+1
        !     the are the values to add to the indices jx, jy, and jz
        !     the are +1 or -1 depending on what cube shouldx be used for interpolation
        addx = int(sign(1._rprec,xlag(jx,jy,jz)))
        addy = int(sign(1._rprec,ylag(jx,jy,jz)))
        addz = int(sign(1._rprec,zlag(jx,jy,jz)))
        jxaddx = jjx + addx
        jyaddy = jjy + addy
        jzaddz = jjz + addz
        !     computes the relative weights given to F_** in the cube depending on point location
        comp_x = abs(xlag(jx,jy,jz))
        comp_y = abs(ylag(jx,jy,jz))
        comp_z = abs(zlag(jx,jy,jz))
        frac_x = 1._rprec - comp_x
        frac_y = 1._rprec - comp_y
        frac_z = 1._rprec - comp_z

        ! Computes interpolated F_KX
        F_KX(jx,jy,jz)=frac_x*frac_y*(FF_KX(jjx,jjy,jjz)*frac_z+FF_KX(jjx,jjy,jzaddz)*comp_z) &
           + frac_x*comp_y*(FF_KX(jjx,jyaddy,jjz)*frac_z+FF_KX(jjx,jyaddy,jzaddz)*comp_z) &
           + comp_x*frac_y*(FF_KX(jxaddx,jjy,jjz)*frac_z+FF_KX(jxaddx,jjy,jzaddz)*comp_z) &
           + comp_x*comp_y*(FF_KX(jxaddx,jyaddy,jjz)*frac_z+FF_KX(jxaddx,jyaddy,jzaddz)*comp_z)

        ! Computes interpolated F_XX
        F_XX(jx,jy,jz)=frac_x*frac_y*(FF_XX(jjx,jjy,jjz)*frac_z+FF_XX(jjx,jjy,jzaddz)*comp_z) &
           + frac_x*comp_y*(FF_XX(jjx,jyaddy,jjz)*frac_z+FF_XX(jjx,jyaddy,jzaddz)*comp_z) &
           + comp_x*frac_y*(FF_XX(jxaddx,jjy,jjz)*frac_z+FF_XX(jxaddx,jjy,jzaddz)*comp_z) &
           + comp_x*comp_y*(FF_XX(jxaddx,jyaddy,jjz)*frac_z+FF_XX(jxaddx,jyaddy,jzaddz)*comp_z)

      END DO
    END DO
  END DO


END SUBROUTINE interpolag_scalar










! This is a modified version of interpolag_Sdep.f90

! This subroutine computes the values of F_KX and F_XX
! at positions (x-u*dt) for use in the subroutine lagrangian

SUBROUTINE interpolag_scalar_Sdep()

  USE types,only:rprec
  USE param
  USE sgsmodule
  use intermediate, only: FF_KX, FF_XX, FF_KX2, FF_XX2

  IMPLICIT NONE
  INTEGER :: jx, jy, jz
  INTEGER :: jjx, jjy, jjz

  !--saves here: force heap storage
  REAL(kind=rprec) :: frac_x,frac_y,frac_z
  REAL(kind=rprec) :: comp_x,comp_y,comp_z
  !REAL(kind=rprec), SAVE, DIMENSION(nxt+2,nynpy+2,nz+2) :: FF_KX,FF_XX
  !REAL(kind=rprec), SAVE, DIMENSION(nxt+2,nynpy+2,nz+2) :: FF_KX2,FF_XX2

  INTEGER :: addx, addy, addz
  INTEGER :: jxaddx, jyaddy, jzaddz



  ! creates dummy arrays FF_KX and FF_XX to use in the subroutine

  DO jz=1,nz
    DO jy=1,nynpy
      DO jx=1,nxt
        FF_KX(jx+1,jy+1,jz+1) = F_KX(jx,jy,jz)
        FF_XX(jx+1,jy+1,jz+1) = F_XX(jx,jy,jz)
        FF_KX2(jx+1,jy+1,jz+1) = F_KX2(jx,jy,jz)
        FF_XX2(jx+1,jy+1,jz+1) = F_XX2(jx,jy,jz)
      END DO
    END DO
  END DO

  !This is a bit like witch craft but the following lines do take care
  !of all the edges including the corners
  FF_KX(1,:,:) = FF_KX(nxt+1,:,:)
  FF_KX(nxt+2,:,:) = FF_KX(2,:,:)

  FF_KX(:,1,:) = FF_KX(:,nynpy+1,:)
  FF_KX(:,nynpy+2,:) = FF_KX(:,2,:)

  FF_KX(:,:,1) = FF_KX(:,:,2)
  FF_KX(:,:,nz+2) = FF_KX(:,:,nz+1)

  FF_XX(1,:,:) = FF_XX(nxt+1,:,:)
  FF_XX(nxt+2,:,:) = FF_XX(2,:,:)

  FF_XX(:,1,:) = FF_XX(:,nynpy+1,:)
  FF_XX(:,nynpy+2,:) = FF_XX(:,2,:)

  FF_XX(:,:,1) = FF_XX(:,:,2)
  FF_XX(:,:,nz+2) = FF_XX(:,:,nz+1)


  FF_KX2(1,:,:) = FF_KX2(nxt+1,:,:)
  FF_KX2(nxt+2,:,:) = FF_KX2(2,:,:)

  FF_KX2(:,1,:) = FF_KX2(:,nynpy+1,:)
  FF_KX2(:,nynpy+2,:) = FF_KX2(:,2,:)

  FF_KX2(:,:,1) = FF_KX2(:,:,2)
  FF_KX2(:,:,nz+2) = FF_KX2(:,:,nz+1)

  FF_XX2(1,:,:) = FF_XX2(nxt+1,:,:)
  FF_XX2(nxt+2,:,:) = FF_XX2(2,:,:)

  FF_XX2(:,1,:) = FF_XX2(:,nynpy+1,:)
  FF_XX2(:,nynpy+2,:) = FF_XX2(:,2,:)

  FF_XX2(:,:,1) = FF_XX2(:,:,2)
  FF_XX2(:,:,nz+2) = FF_XX2(:,:,nz+1)

  ! end of witch craft



  DO jz=1,nz
    jjz = jz+1
    DO jy=1,nynpy
      jjy = jy+1
      DO jx=1,nxt
        jjx = jx+1
        !     the are the values to add to the indices jx, jy, and jz
        !     the are +1 or -1 depending on what cube shouldx be used for interpolation
        addx = int(sign(1._rprec,xlag(jx,jy,jz)))
        addy = int(sign(1._rprec,ylag(jx,jy,jz)))
        addz = int(sign(1._rprec,zlag(jx,jy,jz)))
        jxaddx = jjx + addx
        jyaddy = jjy + addy
        jzaddz = jjz + addz
        !     computes the relative weights given to F_** in the cube depending on point location
        comp_x = abs(xlag(jx,jy,jz))
        comp_y = abs(ylag(jx,jy,jz))
        comp_z = abs(zlag(jx,jy,jz))
        frac_x = 1._rprec - comp_x
        frac_y = 1._rprec - comp_y
        frac_z = 1._rprec - comp_z

        ! Computes interpolated F_KX
        F_KX(jx,jy,jz)=frac_x*frac_y*(FF_KX(jjx,jjy,jjz)*frac_z+FF_KX(jjx,jjy,jzaddz)*comp_z) &
           + frac_x*comp_y*(FF_KX(jjx,jyaddy,jjz)*frac_z+FF_KX(jjx,jyaddy,jzaddz)*comp_z) &
           + comp_x*frac_y*(FF_KX(jxaddx,jjy,jjz)*frac_z+FF_KX(jxaddx,jjy,jzaddz)*comp_z) &
           + comp_x*comp_y*(FF_KX(jxaddx,jyaddy,jjz)*frac_z+FF_KX(jxaddx,jyaddy,jzaddz)*comp_z)

        ! Computes interpolated F_XX
        F_XX(jx,jy,jz)=frac_x*frac_y*(FF_XX(jjx,jjy,jjz)*frac_z+FF_XX(jjx,jjy,jzaddz)*comp_z) &
           + frac_x*comp_y*(FF_XX(jjx,jyaddy,jjz)*frac_z+FF_XX(jjx,jyaddy,jzaddz)*comp_z) &
           + comp_x*frac_y*(FF_XX(jxaddx,jjy,jjz)*frac_z+FF_XX(jxaddx,jjy,jzaddz)*comp_z) &
           + comp_x*comp_y*(FF_XX(jxaddx,jyaddy,jjz)*frac_z+FF_XX(jxaddx,jyaddy,jzaddz)*comp_z)

        ! Computes interpolated F_KX2
        F_KX2(jx,jy,jz)=frac_x*frac_y*(FF_KX2(jjx,jjy,jjz)*frac_z+FF_KX2(jjx,jjy,jzaddz)*comp_z) &
           + frac_x*comp_y*(FF_KX2(jjx,jyaddy,jjz)*frac_z+FF_KX2(jjx,jyaddy,jzaddz)*comp_z) &
           + comp_x*frac_y*(FF_KX2(jxaddx,jjy,jjz)*frac_z+FF_KX2(jxaddx,jjy,jzaddz)*comp_z) &
           + comp_x*comp_y*(FF_KX2(jxaddx,jyaddy,jjz)*frac_z+FF_KX2(jxaddx,jyaddy,jzaddz)*comp_z)

        ! Computes interpolated F_XX2
        F_XX2(jx,jy,jz)=frac_x*frac_y*(FF_XX2(jjx,jjy,jjz)*frac_z+FF_XX2(jjx,jjy,jzaddz)*comp_z) &
           + frac_x*comp_y*(FF_XX2(jjx,jyaddy,jjz)*frac_z+FF_XX2(jjx,jyaddy,jzaddz)*comp_z) &
           + comp_x*frac_y*(FF_XX2(jxaddx,jjy,jjz)*frac_z+FF_XX2(jxaddx,jjy,jzaddz)*comp_z) &
           + comp_x*comp_y*(FF_XX2(jxaddx,jyaddy,jjz)*frac_z+FF_XX2(jxaddx,jyaddy,jzaddz)*comp_z)

      END DO
    END DO
  END DO


END SUBROUTINE interpolag_scalar_Sdep
