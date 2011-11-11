!==============================================================================
! 
!                                  limit
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
! Comments:     Calculate Slope Limiter - see Toro, pg 208
!
!       Input:  Left/Right Delta(w)
!       
!       Ouput:  Limited Delta(w)
!          
!==============================================================================

subroutine limit(left,rght,lval)

  use gparms
  use cntrl
  implicit none
  
  real(dp), intent(in ) :: left,rght
  real(dp), intent(out) :: lval

  if(rght > zero)then
    lval = max(zero,min(beta*left,rght),min(left,beta*rght))
  else
    lval = min(zero,max(beta*left,rght),max(left,beta*rght))
  endif
 
 
  return 

end subroutine limit
