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
  integer, parameter :: roelvink = 22
  integer, parameter :: roelvinky = 23
  integer, parameter :: slosh_inlet = 24
  integer, parameter :: trench = 25
  integer, parameter :: devriend = 26
  integer, parameter :: trenchy = 27
  integer, parameter :: uniform = 28
  integer, parameter :: hibma = 29
  integer, parameter :: channel = 30
  integer, parameter :: exner = 31

  !friction
  integer  :: frictype !=0, no friction, =1 linear, = 2 Manning

  !floodying/drying
  real(dp) :: mindepth 
  logical  :: conserve_volume = .false.
  logical  :: wetdry   = .true.

  !diffusion
  integer :: diffusivity_type = 0
  integer, parameter :: no_diffusivity = 0
  integer, parameter :: constant_diffusivity = 1
  integer, parameter :: smagorinsky_diffusivity = 2
  integer, parameter :: mixing_length_diffusivity = 3
  real(dp) :: diffusivity_coefficient = 5.
  real(dp) :: smagorinsky_coefficient = .0

  !sediment
  integer  :: sedmodel
  real(dp) :: sedinit
  real(dp) :: d50
  real(dp) :: morphfactor
  real(dp) :: porosity = 0. !bed porosity
  integer  :: load_equation = 4
  integer, parameter :: load_MPM = 1 !index for Meyer-Peter Mueller load calc
  integer, parameter :: load_VR  = 2 !index for van Rijn load calc
  integer, parameter :: load_EH  = 3 !index for Engelund + Hansen load calc
  integer, parameter :: load_linear = 4 !linear load = A*u
  logical  :: sed_slope_effect = .false.

  !real(dp) :: d50 = .0001 !.4 mm  !note, .0001 is good for trench case

  !morphology
  real(dp) :: min_morph_depth = .20  !minimum depth for morpho cases (if avoiding wet/dry)

  !mhke
  integer  :: mhke_model
  real(dp) :: mhke_area
  real(dp) :: mhke_cp
  real(dp) :: mhke_xlo
  real(dp) :: mhke_xhi
  real(dp) :: mhke_ylo
  real(dp) :: mhke_yhi

  !trench params
  real(dp), parameter :: trench_half_width = 2.0
  real(dp), parameter :: trench_slope_width = 0.5
  real(dp), parameter :: bath_trench1 = .4
  real(dp), parameter :: bath_trench2 = .56
  real(dp), parameter :: back_slope = 4e-4
  real(dp), parameter :: trench_mid = 11.

  !exner params
  real(dp), parameter :: exner_a0 = 1.0
  real(dp), parameter :: exner_a1 = 1.0
  real(dp), parameter :: exner_Aqf = 1.0
  real(dp), parameter :: exner_lambda = 20.0
  real(dp), parameter :: exner_zeta = 3.0
  real(dp), parameter :: exner_qf = 1.0
  

  !uniform flow params
  real(dp), parameter :: uniform_h = 10.
  real(dp), parameter :: uniform_u = 1.

  !channel parameters
  real(dp), parameter :: channel_mean_depth = 30.
  real(dp), parameter :: channel_max_depth = 30.
  real(dp), parameter :: channel_min_depth = 10.
  real(dp), parameter :: channel_width = 2000.
  real(dp), parameter :: channel_length = 5000.
  real(dp), parameter :: channel_domain = 15000.
  
  !friction / roughness
  real(dp) :: C_manning !manning coefficient, set in input
!  real(dp) :: C_manning = 0.002_dp !manning coefficient for conrun
!  real(dp) :: C_manning = .018_dp  !manning for threehump
!  real(dp) :: C_manning = .03_dp   !manning for heniche
!  real(dp) :: C_manning =   0.00_dp   !manning for exact Riemann
  real(dp) :: friction_depth = 1e9 !only apply friction at depth less than friction_depth - currently unused

  !flux parameters
  integer :: flux_order     !=1 first order, =2, second order
  integer :: transverse_prop  !=0 no transverse flux prop, = 1 transv prop, no correction, =2 transv prop with correction
  integer, dimension(3) :: mthlim = (/1,1,1/) !=0, no lim, =1 minmod, =2 superbee, = 3 van leer, = 4 monotonized centered, =5 beam-warming

  
 

end module cntrl

