!==============================================================================
! 
!                                  mark wet/dry
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
! Comments:     Tags dry cells (set depth,u,v, to zero)
!               Redistribute volume to account for added mass, Ref:
!               Brufau et al., I.J. for Num Meth in Fluids, 2004,45:1047:1082
!
!       Input:
!          flow variables on the patch
!          
!       Output:
!          =>  zeta:    
!          =>  veldepth:
!          =>  drycell:  =0, wet, =1, cell is dry
!          =>   addvol:   volume to add later (m^3)
!==============================================================================


subroutine tagwd(dx,i1,i2,j1,j2,igst,jgst,zeta,veldepth,bath,addvol,drycell)
  use gparms
  use cntrl
  implicit none
  
  integer,  intent(in   ) :: i1,i2,j1,j2,igst,jgst
  real(dp), intent(in   ) :: dx(2)
  real(dp), intent(inout) :: zeta(i1-igst:i2+igst,j1-jgst:j2+jgst)
  real(dp), intent(inout) :: veldepth(i1-igst:i2+igst,j1-jgst:j2+jgst,1:NDIMS)
  real(dp), intent(in )   :: bath(i1-igst:i2+igst,j1-jgst:j2+jgst)
  real(dp), intent(inout) :: addvol(i1-igst:i2+igst,j1-jgst:j2+jgst)
  real(dp), intent(out  ) :: drycell(i1-igst:i2+igst,j1-jgst:j2+jgst)
  
 

  !----- local ------------
  integer  :: i,j
  real(dp) :: area,depth
 
  area = dx(1)*dx(2)  !cell area
  
  !initialize wet/dry flag
  drycell = zero

  !loop over cells in the patch, tag dry cells
  do i=i1,i2
	do j=j1,j2
	  depth = zeta(i,j)-bath(i,j)
      if(depth <= mindepth)then
	    drycell(i,j)      = one                        !mark dry   
	    veldepth(i,j,1:2) = zero                       !zero velocity
	    addvol(i,j)       = (mindepth-depth)*area !store volume about to be added
	    zeta(i,j)         = mindepth + bath(i,j)  !set depth to mindepth (potentially adding volume)
	  endif
	end do
  end do

  !if not bothering with conservation, exit
  if(.not. conserve_volume)return

  !if we have artifically added mass
  !   take water from neighboring cells if available
  !   adjust ud/vd to maintain velocity in neighboring cells
  !   if not available, try later
  do i=i1,i2
	do j=j1,j2
      if(drycell(i,j)==one .and. addvol(i,j) > zero)then
	    
	  endif
	end do
  end do

  
  return
end subroutine tagwd

  