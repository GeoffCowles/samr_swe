!=======================================================================
! Oscar Global Parameter Class
! Copyright:    2009 (c)
!
! THIS IS A DEMONSTRATION RELEASE. THE AUTHOR(S) MAKE NO REPRESENTATION
! ABOUT THE SUITABILITY OF THIS SOFTWARE FOR ANY OTHER PURPOSE. IT IS
! PROVIDED "AS IS" WITHOUT EXPRESSED OR IMPLIED WARRANTY.
!
! THIS ORIGINAL HEADER MUST BE MAINTAINED IN ALL DISTRIBUTED
! VERSIONS.
!
! Authors:      G. Cowles 
!               School for Marine Science and Technology     
!               University of Massachusetts-Dartmouth
!
! Comments:     Global Parameters
!=======================================================================

module gparms  

implicit none

!single precision kind 
integer, parameter :: sp = selected_real_kind(6 , 37)

!double precision kind
integer, parameter :: dp = selected_real_kind(15,307)

!double precision constants
real(dp), parameter :: an8th = 0.125_dp
real(dp), parameter :: a4th  = 0.250_dp
real(dp), parameter :: a3rd  = 1.0_dp/3.0_dp
real(dp), parameter :: ahalf = 0.500_dp
real(dp), parameter :: zero  = 0.000_dp
real(dp), parameter :: one   = 1.000_dp

!----------------------------------------------------------------
!string length
!    fstr:    filename length
!    vstr:    variable name length
!    sstr:    subroutine name length
!    tstr:    text string
!    cstr:    string variable standard length
!----------------------------------------------------------------
integer, parameter  :: fstr = 120  
integer, parameter  :: tstr = 120 
integer, parameter  :: vstr = 30    
integer, parameter  :: sstr = 30   
integer, parameter  :: cstr = 30 

!----------------------------------------------------------------
!trigonemetric
!      pi:    pi
!     d2r:    convert degrees to radians
!     r2d:    convert radians to degrees
!----------------------------------------------------------------

real(dp), parameter :: pi  = 3.14159265358979312 
real(dp), parameter :: d2r = pi/180.0_dp
real(dp), parameter :: r2d = 180.0_dp/pi

!----------------------------------------------------------------
!oceanic parameters
!    gacc        :  gravitational acceleration   [ms^-2]
!    omega_earth :  earths rotation rate         [s^-1]
!----------------------------------------------------------------
real(dp) :: gravity  = 9.8016_dp !default
real(dp), parameter :: omega_earth = 7.292e-5
real(dp), parameter :: rho_w = 1025.0

!----------------------------------------------------------------
!sediment parameters
!    rho_sed :: sediment density [kg/m^3]
!    por_sed :: sediment parosity [-]
!    shields_crit :: critical shields parameter [-]
!    sed_repose_slope :: slope for angle of repose of sediment
!----------------------------------------------------------------
real(dp), parameter :: rho_sed = 2650
real(dp), parameter :: por_sed = .5
real(dp), parameter :: shields_crit = .047
real(dp), parameter :: sed_repose_slope = tan(33.0_dp*d2r)

!----------------------------------------------------------------
!large and small numbers
!    hugenum = largest float
!    tinynum = smallest float
!----------------------------------------------------------------
real(dp), parameter :: hugenum = huge(1.0_dp)
real(dp), parameter :: tinynum = tiny(1.0_dp)

!----------------------------------------------------------------
!params for Sampson test case (see Liang and Borthwick, 2009),pg231
!----------------------------------------------------------------
real(dp) :: a_samp = 3000_dp
real(dp) :: h0_samp = 10.0_dp
real(dp) :: tau_samp = .001_dp
real(dp) :: p_samp 
real(dp) :: s_samp
real(dp) :: B_samp = 5.0_dp



end module gparms
