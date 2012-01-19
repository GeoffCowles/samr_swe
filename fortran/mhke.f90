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
  integer  :: i,j
  real(dp) :: overlap,minx1,maxx1,miny1,maxy1,minx2,maxx2,miny2,maxy2
  real(dp) :: mhke_region,fac,x_overlap,y_overlap
  real(dp) :: u,v,vmag

  !ensure patch overlaps with MHKE zone, otherwise exit
  !if(mhke_xhi > xlo(1) & mhke_xlo < xhi(1) & mhke_yhi > xlo(2) & mhke_ylo < xhi(2))
  
 
  !loop over cells and compute intersection of cell with turbine zone and update momentum
  minx1 = mhke_xlo
  miny1 = mhke_ylo
  maxx1 = mhke_xhi
  maxy1 = mhke_yhi

  mhke_region = abs((maxx1-minx1)*(maxy1-miny1))
  
  
  fac = (dt/(dx(1)*dx(2)))*(mhke_cp*mhke_area)*(1./mhke_region)
  do i=i1,i2
    minx2 = xlo(1)+dx(1)*dble(i-i1)
    maxx2 = xlo(1)+dx(1)*dble(i-i1+1)
    do j=j1,j2
       miny2 = xlo(2)+dx(2)*dble(j-j1)
       maxy2 = xlo(2)+dx(2)*dble(j-j1+1)

       overlap = zero
       if ( minx1 > maxx2 )then
         overlap = 0.0
       else if ( maxx1 < minx2 )then
         overlap = 0.0
       else if ( miny1 > maxy2 )then
         overlap = 0.0
       else if ( maxy1 < miny2 )then
         overlap = 0.0
       else
       ! calculate overlap area
       x_overlap = MIN( maxx2, maxx1 )-MAX( minx2, minx1 )
       y_overlap = MIN( maxy2, maxy1 )-MAX( miny2, miny1 )
       overlap = x_overlap*y_overlap;
       !write(*,*)i,j,overlap

       !augment momentum
       u = vh(i,j,1)/h(i,j)
       v = vh(i,j,2)/h(i,j)
       vmag = sqrt(u**2 + v**2)
       vh(i,j,1) = vh(i,j,1) - .5*fac*overlap*vmag*u 
       vh(i,j,2) = vh(i,j,2) - .5*fac*overlap*vmag*v 
       !write(*,*)'decreasing x/y mom by',fac*overlap*((vh(i,j,1)/h(i,j))**2),fac*overlap*((vh(i,j,2)/h(i,j))**2)
    
       endif
  
       

    enddo
  enddo
  
  return

end subroutine mhke
