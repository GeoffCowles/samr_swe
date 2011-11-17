!==============================================================================
! 
!                            Exner Equation for Sediment
!
! Copyright:    2011(c)
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
! Comments:     Calculate Sediment fluxes using upwinded Exner Equation
!
!       Input:
!                    dt:  global time step on the patch
!                    dx:  patch size (dx,dy)
!       [h,vh,bedlevel]:  flow variables on the patch
!          
!       Output:
!          =>  ifluxm,jfluxm,ifluxp,jfluxp: fluxes on the edges
!==============================================================================

subroutine flux_sed(cid,dt,dx,i1,i2,j1,j2,h,vh,bedlevel,ifluxm,jfluxm,ifluxp,jfluxp)  

 
  use gparms
  use cntrl
  implicit none

  integer,  intent(in)    :: cid
  integer,  intent(in)    :: i1,i2,j1,j2
  real(dp), intent(in )   :: dt
  real(dp), intent(in )   :: dx(2)
  real(dp), intent(in   ) :: h(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
  real(dp), intent(in   ) :: vh(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST,NDIMS)
  real(dp), intent(inout) :: bedlevel(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
  real(dp), intent(inout) :: ifluxm(i1-NFGST:i2+1+NFGST,j1-NFGST:j2+NFGST)
  real(dp), intent(inout) :: jfluxm(j1-NFGST:j2+1+NFGST,i1-NFGST:i2+NFGST)
  real(dp), intent(inout) :: ifluxp(i1-NFGST:i2+1+NFGST,j1-NFGST:j2+NFGST)
  real(dp), intent(inout) :: jfluxp(j1-NFGST:j2+1+NFGST,i1-NFGST:i2+NFGST)

  !local
  real(dp) :: taux(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
  real(dp) :: tauy(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
  real(dp) :: fac1,fac2
  integer i,j

  !zero out fluxes 
  ifluxm = zero
  ifluxp = zero
  jfluxm = zero
  jfluxp = zero
  taux   = zero
  tauy   = zero
  fac1 = gravity*C_manning*C_manning

  do j=j1-NCGST,j2+NCGST
    do i=i1-NCGST,i2+NCGST
	  fac2 = (h(i,j)**(-7/3))*sqrt(vh(i,j,1)**2 + vh(i,j,2)**2)
	  taux(i,j) = fac1*fac2*vh(i,j,1)
	  tauy(i,j) = fac1*fac2*vh(i,j,2)
	  bedlevel(i,j) = rho_w*sqrt(taux(i,j)**2 + tauy(i,j)**2)
	end do
  end do

  return

end subroutine flux_sed

!==============================================================================
! note: dimensioning on jflux (j before i)
!==============================================================================

subroutine consdiff_sed(dx,i1,i2,j1,j2,ifluxm,jfluxm,ifluxp,jfluxp,bedlevel,b)
    
  use gparms
  use cntrl
  implicit none
  integer,  intent(in   ) :: i1,i2,j1,j2
  real(dp), intent(in   ) :: dx(2)
  real(dp), intent(in   ) :: ifluxm(i1-NFGST:i2+1+NFGST,j1-NFGST:j2+NFGST)
  real(dp), intent(in   ) :: jfluxm(j1-NFGST:j2+1+NFGST,i1-NFGST:i2+NFGST)
  real(dp), intent(in   ) :: ifluxp(i1-NFGST:i2+1+NFGST,j1-NFGST:j2+NFGST)
  real(dp), intent(in   ) :: jfluxp(j1-NFGST:j2+1+NFGST,i1-NFGST:i2+NFGST)
  real(dp), intent(inout) :: bedlevel(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
  real(dp), intent(inout) :: b(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
  
  !local
  real(dp) :: bedlevel_last(j1-NFGST:j2+1+NFGST,i1-NFGST:i2+NFGST)
  integer  :: i,j
  real(dp) :: oodx,oody

  !precompute 1/dx,  1/dy
  oodx = one/dx(1);
  oody = one/dx(2);
  return
  !store previous bed thickness for use in morphodynamic adjustment
  bedlevel_last = bedlevel  

  !loop over internal cells, update state variables using convervative diff on fluxes
  do i=i1,i2
	do j=j1,j2
      bedlevel(i,j) = bedlevel(i,j) &
               -oodx*(ifluxm(i+1,j)-ifluxp(i,j)) &
               -oody*(jfluxm(j+1,i)-jfluxp(j,i)) 
    enddo
  end do

  !update bathymetry using change in bed level morphodynamics is active
  if(sedmodel == 2)then
    b = b + (bedlevel - bedlevel_last)
  endif

  return
end subroutine consdiff_sed

