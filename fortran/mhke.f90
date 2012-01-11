!==============================================================================
! 
!                                  mhke
!
! Copyright:    2012(c)
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
! Comments:     Account for tidal turbines influence on momentum
!
!       Input:
!          dt: time step (s)
!          dx: grid spacing [dx,dy] on the patch
!          xlo: min(x),min(y) of the patch boundary
!          xhi: max(x),max(y) of the patch boundary
!          uh: momentum
!          
!       Output:
!          =>  updated momentum after accounting for turbines
!==============================================================================

subroutine mhke(dt,dx,xlo,xhi,i1,i2,j1,j2,h,vh)

  use gparms
  use cntrl
  implicit none
  
  integer,  intent(in )   :: i1,i2,j1,j2
  real(dp), intent(in   ) :: dt
  real(dp), intent(in   ) :: dx(2),xlo(2),xhi(2)
  real(dp), intent(in   ) :: h(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
  real(dp), intent(inout) :: vh(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST,NDIMS)
  

  !----- local ------------
  integer  :: il1,il2,jl1,jl2,i,j
  real(dp) :: epsilon

  !ensure patch overlaps with MHKE zone, otherwise exit
  !if(mhke_xhi > xlo(1) & mhke_xlo < xhi(1) & mhke_yhi > xlo(2) & mhke_ylo < xhi(2))
  

  !set loop bounds
  il1 = i1 
  il2 = i2 
  jl1 = j1 
  jl2 = j2 
 
  !loop over cells and compute intersection of cell with turbine zone and update momentum
 !  do i=il1,il2
 ! 	do j=jl1,jl2
 ! 
 !       uh(i,j) = uh(i,j) - dt*mhke_cp*mhke_A*
 !     enddo
 !   enddo
  
  return

end subroutine mhke
