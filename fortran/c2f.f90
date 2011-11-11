!==============================================================================
! 
!                                  c2f
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
! Comments:     Set simulation variables in fortran space
!
!==============================================================================

subroutine c2f(ndims_in,nequs_in,nscal_in,ncgst_in,nfgst_in)

  use gparms
  use cntrl
  implicit none
  integer, intent(in) :: ndims_in,nequs_in,nscal_in
  integer, intent(in) :: ncgst_in,nfgst_in

  !set dimensions and ghost cell buffers
  NDIMS = ndims_in    !spatial dimensions of the problem
  NEQUS = nequs_in    !number of dynamic state variables (used for sizing flux vec)
  NSCAL = nscal_in    !number of scalar variables
  NCGST = ncgst_in    !number of ghost cells for cell-centered vars
  NFGST = nfgst_in    !number of ghost cells for face based (flux) vars
 
  return

end subroutine c2f
