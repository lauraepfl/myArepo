module constants! bind(C, "constants")
    use iso_c_binding, only: c_float128, c_long_double, c_double
    implicit None 
    
    integer,parameter  :: dp           = c_double ! default
    ! Numerical constants
    real(dp),parameter :: twopi        = real(6.2831853d0, dp)
    real(dp),parameter :: pi           = twopi/2d0
  
    real(dp),parameter :: mu_mol       = real(1.2195d0,dp )
  
    ! Physical constants
    ! Source:
    ! * SI - SI Brochure (2018)
    ! * PCAD - http://www.astro.wisc.edu/~dolan/constants.html
    ! * NIST - National Institute of Standards and Technology
    ! * IAU - Internatonal Astronomical Union resolution
    real(dp),parameter ::hplanck      = 6.6260702d-27 ! Planck const. [erg s]; SI
    real(dp),parameter ::eV2erg       = 1.6021766d-12 ! Electronvolt [erg]; SI
    real(dp),parameter ::kB           = 1.3806490d-16 ! Boltzmann const. [erg K-1]; SI
    real(dp),parameter ::c_cgs        = 2.9979246d+10 ! Speed of light [cm s-1]; SI
    real(dp),parameter ::a_r          = 7.5657233d-15 ! Radiation density const. [erg cm-3 K-4]; SI (derived)
    real(dp),parameter ::mH           = 1.6605390d-24 ! H atom mass [g] = amu, i.e. atomic mass unit; NIST
    real(dp),parameter ::factG_in_cgs = 6.6740800d-08 ! Gravitational const. [cm3 kg-1 s-2]; NIST
    real(dp),parameter ::sigma_T      = 6.6524587d-25 ! Thomson scattering cross-section [cm2]; NIST
    real(dp),parameter ::M_sun        = 1.9891000d+33 ! Solar Mass [g]; IAU
    real(dp),parameter ::R_sun        = 6.9570000d+10 !solar radius
    real(dp),parameter ::L_sun        = 3.8280000d+33 ! Solar Lum [erg s-1]; IAU
    real(dp),parameter ::rhoc         = 1.8800000d-29 ! Crit. density [g cm-3]
    real(dp),parameter ::sigma_boltz  = 5.6703740d-05 !stefan-boltzman constant
    ! Conversion factors - distance
    ! IAU 2012 convention:
    ! 1 pc = 648000 AU / pi
    ! 1 AU = 14 959 787 070 000 cm
    real(dp),parameter ::pc2cm        = 3.0856776d+18
    real(dp),parameter ::kpc2cm       = 3.0856776d+21
    real(dp),parameter ::Mpc2cm       = 3.0856776d+24
    real(dp),parameter ::Gpc2cm       = 3.0856776d+27
  
    ! Conversion factors - time
    ! Year definition follows IAU recommendation
    ! https://www.iau.org/publications/proceedings_rules/units/
    real(dp),parameter ::yr2sec       = 3.15576000d+07 ! Year [s]
    real(dp),parameter ::kyr2sec      = 3.15576000d+10 ! Kyr [s]
    real(dp),parameter ::Myr2sec      = 3.15576000d+13 ! Myr [s]
    real(dp),parameter ::Gyr2sec      = 3.15576000d+16 ! Gyr [s]
    real(dp),parameter ::sec2Gyr      = 3.16880878d-17 ! sec [Gyr]
    real(dp),parameter :: kms         = 1.0d5
  
  
  end module constants

  module barros_parameters
    use constants
    implicit none

    ! Declaring other constants
    real(dp), parameter :: a1 = 2.106d10 * M_sun
    real(dp), parameter :: a2 = 3.859d0 * kpc2cm
    real(dp), parameter :: a3 = 2.162d10 * M_sun
    real(dp), parameter :: a4 = 9.052d0 * kpc2cm
    real(dp), parameter :: a5 = -1.704d10 *M_sun
    real(dp), parameter :: a6 = 3.107d0 * kpc2cm
    real(dp), parameter :: a7 = 0.243d0 * kpc2cm

    real(dp), parameter :: b1 = 0.056d10 * M_sun
    real(dp), parameter :: b2 = 0.993d0 * kpc2cm
    real(dp), parameter :: b3 = 3.766d10 * M_sun
    real(dp), parameter :: b4 = 6.555d0 * kpc2cm
    real(dp), parameter :: b5 = -3.250d10 * M_sun
    real(dp), parameter :: b6 = 7.651d0 * kpc2cm
    real(dp), parameter :: b7 = 0.776d0 * kpc2cm

    real(dp), parameter :: c1 = 2.046d10 * M_sun
    real(dp), parameter :: c2 = 9.021d0 * kpc2cm
    real(dp), parameter :: c3 = 2.169d10 * M_sun
    real(dp), parameter :: c4 = 9.143d0 * kpc2cm
    real(dp), parameter :: c5 = -3.049d10 * M_sun
    real(dp), parameter :: c6 = 7.758d0 * kpc2cm
    real(dp), parameter :: c7 = 0.168d0 * kpc2cm

    real(dp), parameter :: d1 = 0.928d10 * M_sun
    real(dp), parameter :: d2 = 6.062d0 * kpc2cm
    real(dp), parameter :: d3 = 0.163d10 * M_sun
    real(dp), parameter :: d4 = 3.141d0 * kpc2cm
    real(dp), parameter :: d5 = -0.837d10 * M_sun
    real(dp), parameter :: d6 = 4.485d0 * kpc2cm
    real(dp), parameter :: d7 = 0.128d0 * kpc2cm

    real(dp), parameter :: e1 = 2.61d10 * M_sun
    real(dp), parameter :: e2 = 0.44d0 * kpc2cm

    real(dp), parameter :: f1 = 5.4d0 * kpc2cm
    real(dp), parameter :: f2 = 166.d0 * kms

end module barros_parameters

module barros
  use constants
  implicit none
    
  real(dp), parameter:: R0=8.0*kpc2cm
  private
  public gravity 
  interface gravity 
    module procedure my_gravity
  endinterface 
  contains
  subroutine my_gravity(z, g, scale_l,zmid) bind(C, name="F90gravity")
    use barros_parameters
    implicit None 
    real(dp), intent(in):: z, scale_l,zmid
    real(dp), intent(out):: g
    !g = sqrt(z + 1.0)
    !g = zeta(z, 2*z)!
    !print*, zmid
    !call stop()
    call totg(z, zmid, scale_l, g)
 endsubroutine

 subroutine totg(z, z_midplane, scale_l, gravacc)
    
    use barros_parameters
    implicit NONE
    real(dp), intent(IN)::z, z_midplane, scale_l
    real(dp):: rz, dphi_d,dphi_tot
    real(dp), intent(out):: gravacc
    rz=(z-z_midplane)*scale_l
          ! Analytical gravitational acceleration for disk + bulge + halo
          dphi_d =  dphi_d_thin( rz, a1, a2, a3, a4, a5, a6, a7) + &
  &                 dphi_d_thick(rz, b1, b2, b3, b4, b5, b6, b7) + &
  &                 dphi_d_HI(   rz, c1, c2, c3, c4, c5, c6, c7) + &
  &                 dphi_d_H2(   rz, d1, d2, d3, d4 ,d5 ,d6, d7)
          dphi_tot = dphi_d + dphi_b(rz, e1, e2) + dphi_h(rz, f1, f2)
          gravacc=-dphi_tot


  end subroutine
  real(dp) function dphi_DMhalo(z, rh, vs)
  implicit NONE
  real(dp)::z, rh,vs
  dphi_DMhalo =  (2 *rh*vs*(log(z/rh + 1.) - z/(rh + z)))/z**2
  end function dphi_DMhalo
   real(dp) function zeta(z,b)
   implicit None
     real(dp)::z,b
     zeta = sqrt(z**2 + b**2)
   end function zeta
   
   real(dp) function phi_MN_1(z,M,a,b)
   implicit None
   ! Eq. (14) of Barros et al. (2016)
     real(dp)::z,M,a,b
     phi_MN_1 = (-factG_in_cgs*M)/(sqrt(R0**2+(a+zeta(z,b))**2))
   end function phi_MN_1
   
   real(dp) function dphi_MN_1(z,M,a,b)
   implicit None
   ! z-derivative of Eq. (14) of Barros et al. (2016)
     real(dp)::z,M,a,b
     dphi_MN_1 = factG_in_cgs*M*z*(a+zeta(z,b))/(zeta(z,b)*(sqrt(R0**2+(a+zeta(z,b))**2))**3)
   end function dphi_MN_1
   
   real(dp) function dphi_MN_2(z,M,a,b)
   implicit None
   ! z-derivative of Eq. (15) of Barros et al. (2016)
     real(dp)::z,M,a,b
     dphi_MN_2 = dphi_MN_1(z,M,a,b)*(1.d0+a*(a+zeta(z,b))/(R0**2+(a+zeta(z,b))**2))+&
 &               phi_MN_1(z,M,a,b)*(a*z/zeta(z,b)*(R0**2-(a+zeta(z,b))**2)/(R0**2+(a+zeta(z,b))**2)**2)
   end function dphi_MN_2
   
   real(dp) function dphi_MN_3(z,M,a,b)
   ! z-derivative of Eq. (16) of Barros et al. (2016)
   implicit None
     real(dp)::z,M,a,b
     dphi_MN_3 = dphi_MN_1(z,M,a,b)*(1.d0+a*(a+zeta(z,b))/(R0**2+(a+zeta(z,b))**2)-&
 &               1.d0/3.d0*a**2*(R0**2-2.d0*(a+zeta(z,b))**2)/(R0**2+(a+zeta(z,b))**2)**2)+&
 &               phi_MN_1(z,M,a,b)*(a*z/zeta(z,b)*(R0**2-(a+zeta(z,b))**2)/(R0**2+(a+zeta(z,b))**2)**2+&
 &               4.d0*a**2*(a+zeta(z,b))*z/(3.d0*zeta(z,b))*&
 &               (2.d0*R0**2-(a+zeta(z,b))**2)/(R0**2+(a+zeta(z,b))**2)**3)
   end function dphi_MN_3
   
   real(dp) function dphi_d_thin(z,M1,a1,M2,a2,M3,a3,b)
     real(dp)::z,M1,a1,M2,a2,M3,a3,b
     dphi_d_thin = dphi_MN_3(z,M1,a1,b) + dphi_MN_3(z,M2,a2,b) + dphi_MN_3(z,M3,a3,b)
   end function dphi_d_thin
   
   real(dp) function dphi_d_thick(z,M1,a1,M2,a2,M3,a3,b)
   implicit None
     real(dp)::z,M1,a1,M2,a2,M3,a3,b
     dphi_d_thick = dphi_MN_1(z,M1,a1,b) + dphi_MN_1(z,M2,a2,b) + dphi_MN_1(z,M3,a3,b)
   end function dphi_d_thick
   
   real(dp) function dphi_d_HI(z,M1,a1,M2,a2,M3,a3,b)
   implicit None
     real(dp)::z,M1,a1,M2,a2,M3,a3,b
     dphi_d_HI = dphi_MN_2(z,M1,a1,b) + dphi_MN_2(z,M2,a2,b) + dphi_MN_2(z,M3,a3,b)
   end function dphi_d_HI
   
   real(dp) function dphi_d_H2(z,M1,a1,M2,a2,M3,a3,b)
   implicit None
     real(dp)::z,M1,a1,M2,a2,M3,a3,b
     dphi_d_H2 = dphi_MN_3(z,M1,a1,b) + dphi_MN_3(z,M2,a2,b) + dphi_MN_3(z,M3,a3,b)
   end function dphi_d_H2
   
   real(dp) function dphi_b(z,M_b,a_b)
   ! z-derivative of Eq. (19) of Barros et al. (2016)
    implicit None 
     real(dp)::z,M_b,a_b
     dphi_b = (factG_in_cgs*M_b*z)/(sqrt(R0**2+z**2)*(sqrt(R0**2+z**2)+a_b)**2)
   end function dphi_b
   
   real(dp) function dphi_h(z,r_h,v_h)
   ! z-derivative of Eq. (20) of Barros et al. (2016)
   implicit None
     real(dp)::z,r_h,v_h
     dphi_h = v_h**2*z/(R0**2+z**2+r_h**2)
   end function dphi_h
   
end module