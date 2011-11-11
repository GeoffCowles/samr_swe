!=======================================================================
!
! Copyright:    2009(c)
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
! Comments:     Control Parameters
!=======================================================================

module cntrl  

  use gparms
  implicit none

  !dimensions
  integer :: NDIMS    !problem dimension
  integer :: NEQUS    !number of state variable equations 
  integer :: NSCAL    !number of scalar variables
  integer :: NCGST    !number of ghosts for cells 
  integer :: NFGST    !number of ghosts for fluxes

  !test case (used for setting initial conditions)
  integer :: caseid  
  integer, parameter :: ramp = 0
  integer, parameter :: dambreakx = 1
  integer, parameter :: dambreaky = 2
  integer, parameter :: rossby = 3
  integer, parameter :: user_defined = 4
  integer, parameter :: dambreak2D = 5
  integer, parameter :: nomotion = 6
  integer, parameter :: nomotiondryx = 7
  integer, parameter :: threehump = 8
  integer, parameter :: threehump_nomotion = 9
  integer, parameter :: nomotiondryy = 10
  integer, parameter :: dambreakx_dry = 11
  integer, parameter :: dambreakx_berm = 12
  integer, parameter :: threehumpy = 13
  integer, parameter :: dambreaky_dry = 14
  integer, parameter :: dambreaky_berm = 15
  integer, parameter :: sampson = 16
  integer, parameter :: tidetest = 17
  integer, parameter :: conrun = 18
  integer, parameter :: heniche = 19
  integer, parameter :: step = 20
  integer, parameter :: supercrit = 21

  !floodying/drying
  real(dp) :: mindepth = 1e-3 
  logical  :: conserve_volume = .false.
  logical  :: wetdry   = .true.

  !other
  
  !friction / roughness
  real(dp) :: C_manning = 0.002_dp !manning coefficient for conrun
!  real(dp) :: C_manning = .018_dp  !manning for threehump
!  real(dp) :: C_manning = .03_dp   !manning for heniche
!  real(dp) :: C_manning =   0.00_dp   !manning for exact Riemann
  real(dp) :: friction_depth = 1e9 !only apply friction at depth less than friction_depth - currently unused

  !flux parameters
  integer :: flux_order = 2    !=1 first order, =2, second order
  integer :: transverse_prop = 2 !=0 no transverse flux prop, = 1 transv prop, no correction, =2 transv prop with correction
  integer, dimension(3) :: mthlim = (/1,1,1/) !=0, no lim, =1 minmod, =2 superbee, = 3 van leer, = 4 monotonized centered, =5 beam-warming

  
 

end module cntrl

