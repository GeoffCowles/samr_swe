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

subroutine c2f(ndims_in,nequs_in,nscal_in,ncgst_in,nfgst_in, & 
      caseid_in,cmanning_in,mindepth_in,fluxorder_in,transverse_in, & 
      sedmodel_in, sedinit_in, d50_in, morphfactor_in, frictype_in, &
      mhke_model_in,mhke_area_in,mhke_cp_in, &
      mhke_xlo_in,mhke_xhi_in,mhke_ylo_in,mhke_yhi_in)

  use gparms
  use cntrl
  implicit none
  integer,  intent(in) :: ndims_in,nequs_in,nscal_in
  integer,  intent(in) :: ncgst_in,nfgst_in
  real(dp), intent(in) :: cmanning_in,mindepth_in
  integer,  intent(in) :: caseid_in, fluxorder_in, transverse_in
  integer,  intent(in) :: sedmodel_in,frictype_in
  real(dp), intent(in) :: sedinit_in,d50_in,morphfactor_in
  integer,  intent(in) :: mhke_model_in
  real(dp), intent(in) :: mhke_area_in,mhke_cp_in
  real(dp), intent(in) :: mhke_xlo_in,mhke_xhi_in,mhke_ylo_in,mhke_yhi_in

  !set dimensions and ghost cell buffers
  NDIMS = ndims_in    !spatial dimensions of the problem
  NEQUS = nequs_in    !number of dynamic state variables (used for sizing flux vec)
  NSCAL = nscal_in    !number of scalar variables
  NCGST = ncgst_in    !number of ghost cells for cell-centered vars
  NFGST = nfgst_in    !number of ghost cells for face based (flux) vars

  !problem params 
  caseid = caseid_in

  !set physical parameters
  mindepth   = mindepth_in
  C_manning  = cmanning_in
  flux_order = fluxorder_in
  transverse_prop = transverse_in

  !set sediment parameters
  sedmodel    = sedmodel_in
  sedinit     = sedinit_in
  d50         = d50_in
  morphfactor = morphfactor_in
  frictype    = frictype_in

  !set mhke parameters
  mhke_model = mhke_model_in
  mhke_area  = mhke_area_in
  mhke_cp    = mhke_cp_in
  mhke_xlo   = mhke_xlo_in
  mhke_xhi   = mhke_xhi_in
  mhke_ylo   = mhke_ylo_in
  mhke_yhi   = mhke_yhi_in

  write(*,*) 'interfacing parameters to Fortran'
  write(*,*) 'mindepth: ',mindepth
  write(*,*) 'flux_order: ',flux_order
  write(*,*) 'transverse_prop:',transverse_prop
  write(*,*) 'C_manning: ',C_manning
  write(*,*) 'sedmodel: ',sedmodel
  write(*,*) 'sedinit:  ',sedinit
  write(*,*) 'd50:  ',d50
  write(*,*) 'morphfactor: ',morphfactor

  if(mhke_model==1)then
	  write(*,*)'mhke model:  active'
	  write(*,*)'mhke area:', mhke_area
	  write(*,*)'mhke cp  :', mhke_cp
	  write(*,*)'mhke xlo: ',mhke_xlo
	  write(*,*)'mhke xhi: ',mhke_xhi
	  write(*,*)'mhke ylo: ',mhke_ylo
	  write(*,*)'mhke yhi: ',mhke_yhi
	else
		write(*,*)'mhke model: not active'
	endif	
  !pause
  return

end subroutine c2f
