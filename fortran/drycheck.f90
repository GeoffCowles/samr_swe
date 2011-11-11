!==============================================================================
! 
!                          treat dry cells
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
!
! Comments:     Calculate Friction using a Manning Formulation
!
!       Input:
!          [h,uh,vh]:  flow variables on the patch
!          
!       Output:
!          => [h,uh,vh]:  updated uh,vh with dry fix
!
!       If h<0, set h to 0 and velocities to 0  
!==============================================================================

subroutine drycheck(dx,dt,i1,i2,j1,j2,h,vh)
    
  use gparms
  use cntrl
  implicit none
  integer,  intent(in   ) :: i1,i2,j1,j2
  real(dp), intent(in   ) :: dx(2)
  real(dp), intent(in   ) :: dt
  real(dp), intent(inout) :: h(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
  real(dp), intent(inout) :: vh(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST,NDIMS)
 
  !local
  integer  :: i,j
  real(dp) :: dtol = 1e-6 !1e-30 
  return  !gwc debug airplane
  do i=i1,i2
	do j=j1,j2
	  if(h(i,j) <= dtol)then  !depth threshold
	      h(i,j)    = dtol
		  vh(i,j,1) = zero
		  vh(i,j,2) = zero
	  endif
    enddo
  end do

  return
end subroutine drycheck

