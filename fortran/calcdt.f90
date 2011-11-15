!==============================================================================
! 
!                                  calcdt
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
! Comments:     Calculate Stable Time Step for SWE equations
!
!       Input:
!          w:  flow vector on the patch
!          
!       Output:
!          =>  stabdt, maximum stable time step on the patch in units of
!              dx/u
!
!       Todo:
!          =>  use real wavespeed at wet/dry interface to calculate time step
!==============================================================================

subroutine calcdt(dx,i1,i2,j1,j2,igst,jgst,h,vh,b,stabdt)

  use gparms
  use cntrl
  implicit none
  
  integer,  intent(in )    :: i1,i2,j1,j2,igst,jgst
  real(dp), intent(in )    :: dx(ndims)
  real(dp), intent(inout ) :: h(i1-igst:i2+igst,j1-jgst:j2+jgst)
  real(dp), intent(inout ) :: vh(i1-igst:i2+igst,j1-jgst:j2+jgst,1:NDIMS)
  real(dp), intent(in )    :: b(i1-igst:i2+igst,j1-jgst:j2+jgst)
  real(dp), intent(out)    :: stabdt

  !----- local ------------
  integer  :: il1,il2,jl1,jl2,i,j
  real(dp) :: maxspeed(2) 
  real(dp) :: cwave
  integer  :: imaxi,imaxj,jmaxi,jmaxj

  !set loop bounds
  il1 = i1 
  il2 = i2 
  jl1 = j1 
  jl2 = j2 

  !initialize 
  maxspeed = zero
 
  !loop over cells (including ghost) and calculate max stable dt in each direction
  !stable time step is the minimum of the time step in either direction
  do i=il1,il2
	do j=jl1,jl2
	  
	    !gwc new
	  	  !wetting/drying adjust velocity and depth
	  	  if(h(i,j) <= mindepth)then
	  		vh(i,j,1) = zero
	  		vh(i,j,2) = zero
	  		h(i,j) = zero !gwc why is this zero?
	  	  endif	  

      cwave = max(sqrt(gravity*h(i,j)),0.0) 
      maxspeed(1) = max(maxspeed(1),abs(vh(i,j,1)/max(h(i,j),1e-9))+cwave)
      maxspeed(2) = max(maxspeed(2),abs(vh(i,j,2)/max(h(i,j),1e-9))+cwave)
    
    enddo
  enddo
  
  stabdt = max(0.,min((dx(1)/maxspeed(1)),(dx(2)/maxspeed(2))))
  
  return

end subroutine calcdt
