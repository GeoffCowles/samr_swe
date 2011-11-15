!==============================================================================
! 
!                      friction terms: Manning Formulation
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
!          [h,uh,vh],bathy:  flow variables on the patch
!          
!       Output:
!          =>  uh,vh:  updated uh,vh with inclusion of friction
!
!       Uses implicit (theta=1) formulation:
!         taubx = -D u |u| = {(-gn^2 [h])/(h^(4/3))}[u]*|u| (n is manning coefficient)     
!         parts evaluated at n+1 are in brackets = [h][u]
!
!==============================================================================

subroutine friction(dx,dt,i1,i2,j1,j2,h,vh,b)
    
  use gparms
  use cntrl
  implicit none
  integer,  intent(in   ) :: i1,i2,j1,j2
  real(dp), intent(in   ) :: dx(2)
  real(dp), intent(in   ) :: dt
  real(dp), intent(in   ) :: h(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
  real(dp), intent(inout) :: vh(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST,NDIMS)
  real(dp), intent(in   ) :: b(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
 
  !local
  integer  :: i,j
  real(dp) :: depth,CD,coef,hvmag,alpha  
  if(C_manning == zero)return   !no friction
  
  coef = gravity*C_manning*C_manning
  do i=i1,i2
	do j=j1,j2
	  depth = h(i,j)
	  if(depth <= friction_depth)then
		if(depth < mindepth)then
		  vh(i,j,1) = zero
		  vh(i,j,2) = zero
		else	 
	      hvmag = sqrt(vh(i,j,1)**2 + vh(i,j,2)**2)   != h*(|u|)
	      alpha = one/(1 + dt*coef*hvmag*(depth**(-7/3)))
	      vh(i,j,1) = alpha*vh(i,j,1) 
	      vh(i,j,2) = alpha*vh(i,j,2) 	 
	    endif
	  endif 
    enddo
  end do

  return
end subroutine friction

!==============================================================================
! 
!                       friction terms: linear formulation
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
!          [h,uh,vh],bathy:  flow variables on the patch
!          
!       Output:
!          =>  uh,vh:  updated uh,vh with inclusion of friction
!
!       Uses implicit (theta=1) formulation:
!         taubx = -h*u*alpha     
!   
!==============================================================================

subroutine linfriction(dx,dt,i1,i2,j1,j2,h,vh,b)
    
  use gparms
  use cntrl
  implicit none
  integer,  intent(in   ) :: i1,i2,j1,j2
  real(dp), intent(in   ) :: dx(2)
  real(dp), intent(in   ) :: dt
  real(dp), intent(in   ) :: h(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
  real(dp), intent(inout) :: vh(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST,NDIMS)
  real(dp), intent(in   ) :: b(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
 
  !local
  integer  :: i,j
  real(dp) :: depth,alpha

  !if(C_manning == zero)return   !no friction
  
  do i=i1,i2
	do j=j1,j2
	  depth = h(i,j)
	  if(depth < mindepth)then  !depth threshold
		vh(i,j,1) = zero
		vh(i,j,2) = zero
	  else	 
	    alpha = one/(1 + dt*tau_samp)
	    vh(i,j,1) = alpha*vh(i,j,1) 
	    vh(i,j,2) = alpha*vh(i,j,2) 	 
	  endif
    enddo
  end do

  return
end subroutine linfriction


