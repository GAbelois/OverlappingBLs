!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------   
module stokes_drift
!-----------------------------------------------------------------------
!### The module for Stokes drift (Bicheng Chen 10/20/2016)
!### Are going to move all related variables in the future
!-----------------------------------------------------------------------     
use types, only: rprec, float_missing
implicit none
save
public
! ---
!real(rprec) :: agl_stokes = 0._rprec
real(rprec), parameter :: c_omega = 0.8548_rprec

real(rprec), dimension(:), allocatable :: ust, vst  ! the profile of Sokes velocity in x, y-direction
real(rprec), dimension(:,:), allocatable :: ust_x, vst_y

real(rprec) :: U_stokes = float_missing            ! the Stokes velocity at the surface
real(rprec) :: wavenm_w = float_missing            ! the wave number of the wave
real(rprec) :: omega_w = float_missing             ! the angular frequency of the wave
real(rprec) :: rad_Ust = float_missing             ! degree of agl_Ust 

contains
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine init_stokes_drift()
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    use param
    implicit none
    ! ---
    real(rprec) :: z
    integer :: jx, jy, jz
    ! ---
    allocate (ust(0:nz), vst(0:nz))
    allocate (ust_x(nynpy, 0:nz), vst_y(nxt, 0:nz))
    ! ---
    if (stokes_flag) then     
        select case (ubc_mom)

        case ('stress free')    
            select case (stokes_action)
            case('donelan.pierson1987')
                call DP_spectrum()      
            case('mono')
                call monochromatic()
            case('dynamic')
                call dynamic_stokes
            case default
                write(*,*) 'Stokes drift action not recognized. [donelan.pierson197, read, mono, dynamic]'
                stop
            end select
        
            ! ---
            do jz = 0, nz
                !z = (coordz*(nz - 1) + jz - 0.5_rprec)*dz
                ust(jz) = U_stokes*exp(-2._rprec*wavenm_w*zuvp(jz))*cos(rad_Ust)
                vst(jz) = U_stokes*exp(-2._rprec*wavenm_w*zuvp(jz))*sin(rad_Ust)
            end do
            
            do jx = 1, nxt
                vst_y(jx, :) = vst(:)
            end do
            do jy = 1, nynpy
                ust_x(jy, :) = ust(:)
            end do

        case ('wall')
            call mono_shallow_water()

            do jz = 0, nz
                !z = (coordz*(nz - 1) + jz - 0.5_rprec)*dz
                ust(jz) = U_stokes*cosh(2*wavenm_w*(-zuvp(jz)+lz_tot))/  &
                        2._rprec/(sinh(wavenm_w*lz_tot))**2*cos(rad_Ust)
                vst(jz) = U_stokes*cosh(2*wavenm_w*(-zuvp(jz)+lz_tot))/  &
                        2._rprec/(sinh(wavenm_w*lz_tot))**2*sin(rad_Ust)
            end do

            do jx = 1, nxt
                vst_y(jx, :) = vst(:)
            end do
            do jy = 1, nynpy
                ust_x(jy, :) = ust(:)
            end do

        case default
            write (*, *) 'invalid ubc_mom'
            stop
        end select
        
        !call collect_data_z(ust(1:nz-1), ust_tot)
        !if ( rank == 0 ) then
        !    open(9, file='../output/stokes_drift_profile.csv')
        !    write(9,*) 'z,ust'
        !    do jz = 1, nzt-1
        !        z = -(real(jz)-0.5_rprec)*dz*z_i
        !        write(9,*) z, ',', ust_tot(jz)*u_scale
        !    end do
        !    close(9)
        !    
        !    write(*,*) "wave number:", wavenm_w/z_i
        !    write(*,*) "La_t = ", u_star, U_stokes, sqrt(u_star/U_stokes) 
        !end if
    else
        ust(:) = 0._rprec
        vst(:) = 0._rprec
        U_stokes = 0._rprec
        wavenm_w = 0._rprec
    end if
    ! ---
    end subroutine init_stokes_drift
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    ! real(rprec) function get_angFreq(wn, flag_dim)
    ! ! get angular frequency from wave number
    ! real(rprec), intent(in) :: wn       ! wave number
    ! logical, intent(in) :: flag_dim     ! flag to indicate if dimensional variable or not

    ! if (flag_dim) then
    !     get_angFreq = sqrt(g*wn)
    ! else
    !     get_angFreq = sqrt(g*z_i/u_scale**2*wn)
    ! end if
    ! ! ---
    ! end function get_angFreq
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine monochromatic()
    ! Sets the stokes drift as the usual monochromatic wave (based on one
    ! amplitude and one period
    use param, only: z_i, u_scale, lambda_w, amp_w, pi, g, agl_Ust
    implicit none
    ! ---
    rad_Ust = deg2rad(agl_Ust)
    amp_w = amp_w / z_i
    lambda_w = lambda_w / z_i
    wavenm_w = 2._rprec * pi / lambda_w
    omega_w = sqrt(g*z_i/u_scale**2*wavenm_w)
    U_stokes = omega_w * wavenm_w * amp_w**2
    ! ---
    end subroutine monochromatic
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine mono_shallow_water()
    ! Sets the stokes drift as the usual monochromatic wave (based on one
    ! amplitude and one period
    use param, only: lz_tot, z_i, u_scale, lambda_w, amp_w, pi, g, agl_Ust
    implicit none
    ! ---
    rad_Ust = deg2rad(agl_Ust)
    amp_w = amp_w / z_i
    lambda_w = lambda_w / z_i
    wavenm_w = 2._rprec * pi / lambda_w
    omega_w = sqrt(g*z_i/u_scale**2*wavenm_w*tanh(wavenm_w*lz_tot))
    U_stokes = omega_w * wavenm_w * amp_w**2 / tanh(wavenm_w*lz_tot)
    ! ---
    end subroutine mono_shallow_water
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine DP_spectrum()
    ! TOMAS CHOR
    ! Sets the stokes drift as the integral of a broadband wave spectrum
    ! introduced in Donelan & Pierson, 1987, Journal of Geophysical Research
    ! Subroutine was first programmed by Di Yang and later adapted to this format
    use param
    implicit none
    ! Most of these parameters are given in Sec 5.4 of donelan.pierson1987
    real(rprec), parameter :: KM=370.0_rprec
    real(rprec), parameter :: N1=5.0_rprec,N2=1.15_rprec
    real(rprec), parameter :: A1=22.0_rprec,A2=4.6_rprec,B=3._rprec
    real(rprec), parameter :: RHO_A=1.225_rprec,RHO_W=1000.0_rprec,NU=1.E-6_rprec
    integer, parameter :: NMAX=100000, NTHETA=360
    real(rprec), parameter :: THETA_MAX=PI
    real(rprec), parameter :: DTHETA=THETA_MAX*2._rprec/NTHETA
    real(rprec), parameter :: KAPPA=0.4
    real(rprec), parameter :: CD = 0.17
    real(rprec) :: DKW
  
    !DY U10 needs to be consistent with the friction velocity u_scale
    real(rprec) :: U10
    real(rprec) :: Cdn
    real(rprec) :: US_WIND
    real(rprec) :: Z0
    real(rprec) :: WAVENM_P
  
    real(kind=rprec), dimension(ntheta)::theta
    real(kind=rprec), dimension(nmax,ntheta)::swk
    real(kind=rprec), dimension(nmax)::swk1d
  
    integer :: i,j,k
    real(rprec) :: h,kw,n,a,uk,Z
  
    call get_U10(u_scale, U10) ! Gotta pass u_scale with dimensions
    !print*,'Stokes drift: Calculated U10:', U10,' for u_scale = ', u_scale
  
    Cdn=0.001_rprec*(0.96_rprec+0.041_rprec*U10) ! Eq (14) of donelan.pierson1987
    US_WIND=SQRT(Cdn*U10**2)
    Z0=10._rprec*EXP(-KAPPA/SQRT(CD))
    WAVENM_P = G/(1.2*U10)**2 ! Eq. (22) of donelan.pierson1987
    DKW=WAVENM_P/100._rprec
    wavenm_w = wavenm_p * z_i       ! must normalize the wavenumber
  
    DO J=1,NTHETA
      THETA(J)=-THETA_MAX+J*DTHETA
    ENDDO
  
    DO I=1,NMAX
      KW=I*DKW
      N=(N1-N2)*(ABS((KW**2-KM**2)/(KW**2+KM**2)))**B+N2
      A=EXP((A1-A2)*(ABS((KW**2-KM**2)/(KW**2+KM**2)))**B+A2)
      UK=US_WIND/KAPPA*LOG(PI/KW/Z0)
      IF(KW.LE.0.31_rprec*WAVENM_P) THEN
        H = 1.24_rprec
      ELSE IF(KW.LE.0.9_rprec*WAVENM_P) THEN
        H = 2.61_rprec*(KW/WAVENM_P)**0.65_rprec
      ELSE
        H = 2.28_rprec*(WAVENM_P/KW)**0.65_rprec
      ENDIF
      SWK1D(I)=0.
      DO J=1,NTHETA
        IF(KW.LE.10._rprec*WAVENM_P) THEN
          SWK(I,J) = (1.35E-3_rprec)/KW**3.5_rprec/WAVENM_P**0.5_rprec &
             *EXP(-WAVENM_P**2/KW**2)*1.7_rprec**EXP(-1.22_rprec &
             *(SQRT(KW/WAVENM_P)-1._rprec)**2)*H/(COSH(H*THETA(J)))**2
        ELSE
          SWK(I,J)=(0.194_rprec/A*(RHO_A/RHO_W)*(UK/SQRT(G/KW)-1._rprec)**2 &
             -4._rprec*NU*KW/A/SQRT(G/KW))**(1._rprec/N)/KW**4 &
             /(COSH(H*THETA(J)))**2
          IF(ISNAN(SWK(I,J))) THEN
            print*,'KW, DKW, WAVENM_p, NU, THETA(J), H'
            print*,KW, DKW, WAVENM_p, NU, THETA(J), H
            PRINT*, "SW=NaN!",I,J,SWK(I,J)
            STOP
          ENDIF
        ENDIF
        SWK1D(I) = SWK1D(I) + SWK(I,J)*KW*DTHETA
      ENDDO
    ENDDO
  
    if(coordz.eq.0) then
      WRITE(50,*) 'VARIABLES=K,SWK1D,SWK0'
      DO I=1,NMAX
        KW=I*DKW
        write(50,*) KW/WAVENM_P,SWK1D(I),SWK(I,NTHETA/2)
      ENDDO
    endif
  
    ! DO K=0,nz
    !   z= - (coordz*(nz-1) + k - 0.5_rprec) * dz *z_i
    !   UST_abs(K)=0.
    !   DO I=1,NMAX
    !     KW=I*DKW
    !     DO J=1,NTHETA
    !       UST_abs(K) = UST_abs(K) + 2.*G**0.5*KW**2.5*SWK(I,J) &
    !          *EXP(2.*KW*Z)*DKW*COS(THETA(J))*DTHETA/u_scale
    !     ENDDO
    !   ENDDO
    ! ENDDO
  
    !-------
    ! Calculates U_stokes at z=0 exactly
    U_stokes=0.0_rprec
    DO I=1,NMAX
      KW=I*DKW
      DO J=1,NTHETA
        U_stokes = U_stokes + 2.*G**0.5*KW**2.5*SWK(I,J)*DKW*COS(THETA(J))*DTHETA/u_scale
      ENDDO
    ENDDO
    rad_Ust = deg2rad(agl_Ust)
    !-------
    ! if (coordz==0) then
    !   ust_abs(0) = U_stokes+2.0d0 ! Just in case (otherwise node zero is trash)
    ! endif
    ! ---
    end subroutine DP_spectrum
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    real(rprec) function get_Ustokes(omega, wn, amp)
    ! get Stokes drift velocity at the surface
    ! >> omega - angular frequency
    ! >> wn - wave number
    ! >> amp - amplitude of wave
    real(rprec), intent(in) :: omega, wn, amp
    ! ---
    get_Ustokes = omega*wn*amp**2
    ! ---
    end function get_Ustokes
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine dynamic_stokes()
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------  
    use param
    implicit none
    ! ---
    inquire (iolength=len_gwave) amp_w, lambda_w, agl_Ust
    open (unit=fid_gwave, file=fn_gwave, access='direct', recl=len_gwave)
    read (unit=fid_gwave, rec=1) amp_w, lambda_w, agl_Ust
    close (fid_gwave)
        
    rad_Ust = deg2rad(agl_Ust)
    amp_w = amp_w/z_i
    lambda_w = lambda_w/z_i
    wavenm_w = 2._rprec*pi/lambda_w
    omega_w = sqrt(g*z_i/u_scale**2*wavenm_w)
    U_stokes = omega_w*wavenm_w*amp_w**2
    ! ---
    end subroutine dynamic_stokes
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    ! subroutine get_PM(ustar, cd, rho_air, rho_ocean, Lat, flag_dim, &
    !                   omega, wn, amp)
    ! !!! get the angular frequency and wave number at peak of
    ! !!! >>Pierson-Moskowitz (PM) spectrum
    ! ! >> ustar - friction velocity
    ! ! >> cd - drag coefficient wind stress for air
    ! ! >> rho_air - air density
    ! ! >> rho_ocean - sea water density
    ! ! >> Lat - equilibrium turbulent Langmuir number
    ! ! >> flag_dim - flag to indicate if dimensional variable or not
    ! ! >> omega - angular frequency of PM spectrum
    ! ! >> wn - wave number of PM spectrum
    ! ! >> amp - amplitude of the wave (at Lat level)
    ! ! >> u10 - wind speed at 10 m
    ! real(rprec), intent(in) :: ustar, cd, rho_air, rho_ocean, Lat
    ! logical, intent(in) :: flag_dim
    ! real(rprec), intent(out) :: omega, wn, amp
    
    ! real(rprec) :: u10
    ! ! ---
    ! u10 = ustar*sqrt(rho_ocean/rho_air/cd)
    ! if (flag_dim) then
    !     omega = c_omega*g/u10
    !     wn = omega**2/g
    ! else
    !     omega = c_omega*(g*z_i/u_scale**2)/u10
    !     wn = omega**2/(g*z_i/u_scale**2)
    ! end if
    ! amp = sqrt(ustar/(omega*wn))/Lat
    ! ! ---
    ! end subroutine get_PM
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine get_U10(u_scale, U10)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------    
    real(rprec), intent(in) :: u_scale
    real(rprec), intent(out) :: U10
    real(rprec), parameter :: rho_w = 1031_rprec
    real(rprec), parameter :: rho_a = 1.21_rprec
    real(rprec), parameter :: A = 0.96_rprec
    real(rprec), parameter :: B = 0.041_rprec
    real(rprec), parameter :: D = 1d-3
    complex(rprec) :: coeff
  
    coeff = -(4d0*A**3d0*D*rho_a-27d0*B**2d0*rho_w*u_scale**2d0)*rho_w
    coeff = sqrt(coeff)
  
    U10 = real((-1d0/27d0*A**3d0/B**3d0 + 1d0/2d0*rho_w*u_scale**2d0/(B*D*rho_a) + &
        1d0/6d0*sqrt(1d0/3d0)* coeff *u_scale/(B**2d0*D*rho_a))**(1d0/3d0) &
       - 1d0/3d0*A/B + 1d0/9d0*A**2d0/(B**2d0*(-1d0/27d0*A**3d0/B**3d0 + 1d0/2d0*rho_w*u_scale**2d0/(B*D*rho_a) + &
       1d0/6d0*sqrt(1d0/3d0)* coeff *u_scale/(B**2d0*D*rho_a))**(1d0/3d0)))
    ! ---
    end subroutine get_U10
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    real(rprec) function deg2rad(deg)
    ! Convert from degree to radians
    use param, only: pi
    implicit none
    ! ---
    real(rprec), intent(in) :: deg

    deg2rad = deg / 180._rprec * pi
    ! ---
    end function deg2rad

end module stokes_drift
