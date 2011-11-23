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
! Comments:     Calculate Sediment fluxes using upwinded Exner Equation and MPM load
!
!       Input:
!                    dt:  global time step on the patch
!                    dx:  patch size (dx,dy)
!       [h,vh,bedlevel]:  flow variables on the patch
!          
!       Output:
!          =>  ifluxm,jfluxm,ifluxp,jfluxp: fluxes on the edges
!==============================================================================

subroutine flux_sed(cid,dt,dx,xlo,xhi,i1,i2,j1,j2,h,vh,bedlevel,iflux,jflux)  

 
  use gparms
  use cntrl
  implicit none

  integer,  intent(in)    :: cid
  integer,  intent(in)    :: i1,i2,j1,j2
  real(dp), intent(in )   :: dt
  real(dp), intent(in )   :: dx(2)
  real(dp), intent(in )   :: xlo,xhi
  real(dp), intent(in   ) :: h(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
  real(dp), intent(in   ) :: vh(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST,NDIMS)
  real(dp), intent(inout) :: bedlevel(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
  real(dp), intent(inout) :: iflux(i1-NFGST:i2+1+NFGST,j1-NFGST:j2+NFGST)
  real(dp), intent(inout) :: jflux(j1-NFGST:j2+1+NFGST,i1-NFGST:i2+NFGST)

  !local
  real(dp) :: taux(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
  real(dp) :: tauy(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
  real(dp) :: qx(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
  real(dp) :: qy(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
  real(dp) :: taub,tau_edge
  real(dp) :: fac1,fac2,shieldfac,gprime,facMPM,facVR,dstar,loadMPM,loadVR,loadfac
  integer i,j

  !zero out fluxes 
  iflux = zero
  jflux = zero

  !zero out local vars
  qx     = zero
  qy     = zero
  taux   = zero
  tauy   = zero

  !set up the constants
  fac1      = gravity*C_manning*C_manning
  gprime    = gravity*(rho_sed/rho_w-1)  !reduced gravity
  shieldfac = 1./(gprime*d50) !converts stress (in units of m^2/s^2) to nondim stress
  dstar     = 2.0313e+04*d50  !nondimensional grain size [assumes nu=1.36e-6]
  loadfac   = sqrt(gprime*(d50**3)) !converts load from (-) into m^2/s
  facMPM    = 8.0 !multiplier for MPM
  facVR     = .053/(dstar**0.3) !multiplier for van Rijn
  
  
 
  !compute solid load at cell centers (units = m^2/s)
  do j=j1-NCGST,j2+NCGST
    do i=i1-NCGST,i2+NCGST

	  !(re)compute the stresses from the Manning formulation
	  fac2 = (h(i,j)**(-7./3.))*sqrt(vh(i,j,1)**2 + vh(i,j,2)**2)
	  taux(i,j) = fac1*fac2*vh(i,j,1)
	  tauy(i,j) = fac1*fac2*vh(i,j,2)       
	  taub = sqrt(taux(i,j)**2 + tauy(i,j)**2)  !magnitude of bottom stress [m^2/s^2]
	
	  !compute the load at cell centers (ignore slope for now) using MPM formula
	  !see eq. 11 de Swart and Zimmerman 2000.
	  ! loadMPM = facMPM*loadfac*(max(taub*shieldfac - shields_crit,0.0)**1.5d0)
	  ! 	  qx(i,j) = loadMPM*taux(i,j)/taub
	  ! 	  qy(i,j) = loadMPM*tauy(i,j)/taub
	  
	
	  !compute the load at cell centers (ignore slope for now) using Van Rijn 1984
	  loadVR = facVR*loadfac*(max(taub*shieldfac/shields_crit-1,zero)**2.1)
	  qx(i,j) = loadVR*taux(i,j)/taub
	  qy(i,j) = loadVR*tauy(i,j)/taub

	 end do
	end do
	
	!compute load fluxes on edges using upwinded formulation
	!ahalf is for correctly computing flux
	!morphfactor is the acceleration factor for morphodynamics
	!porosity is included as a factor 1./(1-porosity) to account for changes in the bed thickness
	!need to account for fluxes into/out of the water column
	fac1 = ahalf*dt*morphfactor*(1./(1-porosity))
	
	do j=j1,j2
	do i=i1,i2+1
	   tau_edge = ahalf*(taux(i,j)+taux(i-1,j))
		iflux(i,j) = fac1*((1+sign(1.,tau_edge))*qx(i-1,j)+(1-sign(1.,tau_edge))*qx(i,j)) !upwind flux
	   !iflux(i,j) = fac1*(qx(i-1,j)+qx(i,j))  !centered flux
	end do
   end do

   !slope effects (TODO)
 !   sed_angle=DTAN(33.0_r8*pi/180.0_r8)
 !
 !   cff1=0.5_r8*(1.0_r8+SIGN(1.0_r8,FX_r(i,j)))
 !            cff2=0.5_r8*(1.0_r8-SIGN(1.0_r8,FX_r(i,j)))
 !             dzdx=MIN(((h(i+1,j)-h(i  ,j))/om_u(i+1,j)*cff1+             &
 !      &                (h(i  ,j)-h(i-1,j))/om_u(i  ,j)*cff2),0.52_r8)*   &
 !      &                SIGN(1.0_r8,FX_r(i,j))
 !             dzdy=MIN(((h(i,j+1)-h(i,j  ))/on_v(i,j+1)*cff1+             &
 !      &                (h(i,j  )-h(i,j-1))/on_v(i  ,j)*cff2),0.52_r8)*   &
 !      &                SIGN(1.0_r8,FE_r(i,j))
 !             cff=DATAN(dzdx)
 !             a_slopex=sed_angle/(COS(cff)*(sed_angle-dzdx))
 !             cff=DATAN(dzdy)
 !             a_slopey=sed_angle/(COS(cff)*(sed_angle-dzdy))


!
! Add contribution of bed slope to bed load transport fluxes.
!
!             FX_r(i,j)=FX_r(i,j)*a_slopex
!             FE_r(i,j)=FE_r(i,j)*a_slopey





   !set flux to zero at left boundary !FUDGE
   if(cid==trench .or. cid==devriend)then
   if(abs(xlo)<1e-5)then
     iflux(0,:) = iflux(1,:)
     !write(*,'(4F10.3)')h(-1,1),h(0,1),vh(-1,1,1),vh(0,1,1)
   endif
  !     if(abs(xhi-22.)<1e-5)then
  !    	     iflux(i2+1,:) = iflux(i2,:)
  !    	   endif
    endif
   	
  

   do i=i1,i2
  	do j=j1,j2+1
  	  tau_edge = .5*(tauy(i,j)+tauy(i,j-1))
  	  jflux(j,i) = fac1*((1+sign(1.,tau_edge))*qy(i,j-1)+(1-sign(1.,tau_edge))*qy(i,j))
  	end do
   end do
   


  
 !  do j=j1-NCGST,j2+NCGST
 ! 	bedlevel(i1-NCGST:i1-1,j) = bedlevel(i1,j)
 ! 	bedlevel(i2+1:i2+NCGST,j) = bedlevel(i2,j)
 !   end do
 !   do i=i1-NCGST,i2+NCGST
 ! 	bedlevel(i,j1-NCGST:j1-1) = bedlevel(i,j1)
 ! 	bedlevel(i,j2+1:j2+NCGST) = bedlevel(i,j2)
 !   end do
 
 
  

  return

end subroutine flux_sed

!==============================================================================
! note: dimensioning on jflux (j before i)
!==============================================================================

subroutine consdiff_sed(dx,i1,i2,j1,j2,iflux,jflux,bedlevel,b)
    
  use gparms
  use cntrl
  implicit none
  integer,  intent(in   ) :: i1,i2,j1,j2
  real(dp), intent(in   ) :: dx(2)
  real(dp), intent(in   ) :: iflux(i1-NFGST:i2+1+NFGST,j1-NFGST:j2+NFGST)
  real(dp), intent(in   ) :: jflux(j1-NFGST:j2+1+NFGST,i1-NFGST:i2+NFGST)
  real(dp), intent(inout) :: bedlevel(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
  real(dp), intent(inout) :: b(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
  
  !local
  real(dp) :: bedlevel_last(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
  integer  :: i,j
  real(dp) :: oodx,oody

  !precompute 1/dx,  1/dy
  oodx = one/dx(1);
  oody = one/dx(2);
  
  !store previous bed thickness for use in morphodynamic adjustment
  bedlevel_last = bedlevel  

  !loop over internal cells, update state variables using convervative diff on fluxes
  do i=i1,i2
	do j=j1,j2
      bedlevel(i,j) = bedlevel(i,j) &
               -oodx*(iflux(i+1,j)-iflux(i,j)) &
               -oody*(jflux(j+1,i)-jflux(j,i)) 
   !  if(j==1)then
   ! 	  write(*,*)i,bedlevel(i,j)
   ! 	endif
    enddo
  end do
 !  do i=i1,i2
 !   do j=j1,j2
 ! 	 bedlevel(i,j) = .5*bedlevel(i,j) + .25*bedlevel(i-1,j)+.25*bedlevel(i+1,j)
 !   enddo
 !   enddo
  ! write(*,*)'check: ',bedlevel(i1,1),iflux(i1+1,1),iflux(i1,1)
     !stop

  !update bathymetry using change in bed level morphodynamics is active
  if(sedmodel == 2)then
    b = b + (bedlevel - bedlevel_last)
  endif

  return
end subroutine consdiff_sed

