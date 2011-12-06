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
! Comments:     Calculate Sediment fluxes using upwinded formula
!               3 Load options implemented - MPM, VR, Engelund+Hansen
!               Options to modify load with bedslope (Lesser formulation)
!
!       Input:
!                    dt:  global time step on the patch
!                    dx:  patch size (dx,dy)
!       [h,vh,bedlevel]:  flow variables on the patch
!          
!       Output:
!          =>  ifluxm,jfluxm,ifluxp,jfluxp: fluxes on the edges
!==============================================================================

subroutine flux_sed(cid,dt,dx,xlo,xhi,i1,i2,j1,j2,h,vh,bedlevel,b,iflux,jflux)  

 
	use gparms
	use cntrl
	implicit none

	integer,  intent(in)    :: cid
	integer,  intent(in)    :: i1,i2,j1,j2
	real(dp), intent(in )   :: dt
	real(dp), intent(in )   :: dx(2)
	real(dp), intent(in )   :: xlo(2),xhi(2)
	real(dp), intent(inout) :: h(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
	real(dp), intent(inout) :: vh(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST,NDIMS)
	real(dp), intent(in   ) :: bedlevel(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
	real(dp), intent(in   ) :: b(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
	real(dp), intent(inout) :: iflux(i1-NFGST:i2+1+NFGST,j1-NFGST:j2+NFGST)
	real(dp), intent(inout) :: jflux(j1-NFGST:j2+1+NFGST,i1-NFGST:i2+NFGST)

	!local
	real(dp) :: taux(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
	real(dp) :: tauy(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
	real(dp) :: qx(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
	real(dp) :: qy(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
	real(dp) :: taub,tau_edge
	real(dp) :: fac1,fac2,shieldfac,gprime,facMPM,facVR,facEH
	real(dp) :: dstar,loadMPM,loadVR,loadEH,loadfac,loadEH2
	real(dp) :: slope_fac,dbdx,dbdy,slmax,oo2dx,oo2dy
	real(dp) :: dzdx
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
  facEH     = .05*sqrt(gravity)*C_manning*d50
  
  
  !----------------------------------------------------------------------
  !compute solid load at cell centers (units = m^2/s) [=> qx,qy]
  !----------------------------------------------------------------------
	select case(load_equation)

  case(load_MPM)

  do j=j1-NCGST,j2+NCGST
		do i=i1-NCGST,i2+NCGST

			!(re)compute the stresses from the Manning formulation
			fac2 = (h(i,j)**(-7./3.))*sqrt(vh(i,j,1)**2 + vh(i,j,2)**2)
			taux(i,j) = fac1*fac2*vh(i,j,1)
			tauy(i,j) = fac1*fac2*vh(i,j,2)   
			 
			taub = sqrt(taux(i,j)**2 + tauy(i,j)**2)  !magnitude of bottom stress [m^2/s^2]

			!see eq. 11 de Swart and Zimmerman 2000.
			loadMPM = facMPM*loadfac*(max(taub*shieldfac - shields_crit,0.0)**1.5d0)
			qx(i,j) = loadMPM*taux(i,j)/taub
			qy(i,j) = loadMPM*tauy(i,j)/taub

		end do
	end do

	case(load_VR)

	do j=j1-NCGST,j2+NCGST
		do i=i1-NCGST,i2+NCGST

		!(re)compute the stresses from the Manning formulation
		fac2 = (h(i,j)**(-7./3.))*sqrt(vh(i,j,1)**2 + vh(i,j,2)**2)
		taux(i,j) = fac1*fac2*vh(i,j,1)
		tauy(i,j) = fac1*fac2*vh(i,j,2)       
		taub = sqrt(taux(i,j)**2 + tauy(i,j)**2)  !magnitude of bottom stress [m^2/s^2]

		!compute sediment load according to van Rijn
		loadVR  = facVR*loadfac*(max(taub*shieldfac/shields_crit-1,zero)**2.1)
		qx(i,j) = loadVR*taux(i,j)/taub
		qy(i,j) = loadVR*tauy(i,j)/taub

		end do
	end do

	case(load_EH)  ! Engelund and Hansen - see van der wegen and roelvink, JGR 2008, eq. 5
  
	do j=j1-NCGST,j2+NCGST
		do i=i1-NCGST,i2+NCGST

		!(re)compute the stresses from the Manning formulation
		fac2 = (h(i,j)**(-7./3.))*sqrt(vh(i,j,1)**2 + vh(i,j,2)**2)
		taux(i,j) = fac1*fac2*vh(i,j,1)
		tauy(i,j) = fac1*fac2*vh(i,j,2)       
		taub = sqrt(taux(i,j)**2 + tauy(i,j)**2)  !magnitude of bottom stress [m^2/s^2]

    !load using Engelund and Hansen
		loadEH = facEH*taub*(shieldfac**2)*((vh(i,j,1)**2 + vh(i,j,2)**2)**1.5)*(h(i,j)**-3)*(h(i,j)**-(1./6.))
		qx(i,j) = loadEH*taux(i,j)/taub
		qy(i,j) = loadEH*tauy(i,j)/taub

		end do
	end do
  
	end select
	
	! set loads and stresses to zero in dry points (this may not be necessary)
	do j=j1-NCGST,j2+NCGST
		do i=i1-NCGST,i2+NCGST 
			if(h(i,j) < mindepth)then
				qx(i,j) = zero
				taux(i,j) = zero
				qy(i,j) = zero
				tauy(i,j) = zero
			endif
		end do
	end do
	
	!-----------------------------------------------------------------------------------------
	!modify loads using bedslope formulations (e.g. lesser et al. in Warner, VDW & ROELVINK)
	!note the minus sign on dbdx/dbdy.  Our stresses here are in the sense of direction of the flow
	!a +u flow will produce a +taux.  We would then want a positive bed slope (rising in x) to 
	!have an adverse effect (< 0) on the load 
	!-----------------------------------------------------------------------------------------
	
  if(sed_slope_effect)then
		oo2dx = 1./(2*dx(1))
		oo2dy = 1./(2*dx(2))
		slmax = .9*sed_repose_slope  !limit the bedslope to 90% of the angle of repose
		do j=j1-1,j2+1
	    do i=i1-1,i2+1
		
      	dbdx = -sign(1.,taux(i,j))*min( (b(i+1,j)-b(i-1,j))*oo2dx , slmax)
        slope_fac = sed_repose_slope/( cos(atan(dbdx))*(sed_repose_slope - dbdx) )
        qx(i,j) = qx(i,j)*slope_fac

        dbdy = -sign(1.,tauy(i,j))*min( (b(i,j+1)-b(i,j-1))*oo2dy , slmax)
        slope_fac = sed_repose_slope/( cos(atan(dbdy))*(sed_repose_slope - dbdy) )
        qy(i,j) = qy(i,j)*slope_fac
       !  if(j==1)then
       ! 	        write(*,*)i,qx(i,j),qy(i,j)
       ! 	      endif
			end do
		end do
	endif
	
	
	!-----------------------------------------------------------------------------------------
	!compute load fluxes on edges using upwinded formulation
	!ahalf is for correctly computing flux
	!morphfactor is the acceleration factor for morphodynamics
	!porosity is included as a factor 1./(1-porosity) to account for changes in the bed thickness
	!need to account for fluxes into/out of the water column
	!-----------------------------------------------------------------------------------------
	
	fac1 = ahalf*dt*morphfactor*(1./(1-porosity))
	
	! idir
	do j=j1,j2
		do i=i1,i2+1
	  	tau_edge = ahalf*(taux(i,j)+taux(i-1,j))
			iflux(i,j) = fac1*((1+sign(1.,tau_edge))*qx(i-1,j)+(1-sign(1.,tau_edge))*qx(i,j)) !upwind flux
		end do
	end do
	
	! jdir  Note that SAMRAI stores jfluxes in arrays with dimensions (jl,il)
	do i=i1,i2
		do j=j1,j2+1
			tau_edge = .5*(tauy(i,j)+tauy(i,j-1))
			jflux(j,i) = fac1*((1+sign(1.,tau_edge))*qy(i,j-1)+(1-sign(1.,tau_edge))*qy(i,j))
	    !if(tauy(i,j)*tauy(i,j-1) <= 0.0) jflux(j,i)=zero
	  end do
	end do

	
	!-----------------------------------------------------------------------------------------
	! set fluxes (zero or extrapolate) on certain physical boundaries for specific cases
	!-----------------------------------------------------------------------------------------
	
	!set j-dir flux to zero (fudge for testing devriend)
	!if(cid == devriend) jflux = zero

	!extrapolate i flux at inflow boundary !FUDGE
	if(cid==trench .or. cid==devriend)then
  if(abs(xlo(1))<1e-5)then
      iflux(0,:) = iflux(1,:)
    endif
   !    if(abs(xhi(1)-22.)<1)then
   !             	     iflux(i2+1,:) = iflux(i2,:)
   !             	   endif
  endif
   
	!set jflux to zero at bottom boundary for y-direction trench !FUDGE
	if(cid==trenchy)then
		if(abs(xlo(2))<1e-5)then
			jflux(0,:) = jflux(1,:)
		endif
	endif

  !  !set flux to zero along line of symmetry in devriend or roelvink case case
  !    if(cid==devriend)then
     if(cid==roelvink)then
      if(abs(xlo(2))<.1)then
        jflux(0,:) = zero
      endif
      endif
      if(cid==devriend) jflux = zero

if(cid==hibma)then
	if(abs(xlo(1))<1)then
      iflux(0,:) = iflux(1,:)
    endif
    if(abs(xhi(1)-80000.)<1)then
              	     iflux(i2+1,:) = zero
              	   endif
endif

  return

end subroutine flux_sed

!==============================================================================
! note: dimensioning on jflux (j before i)
!==============================================================================

subroutine consdiff_sed(dx,i1,i2,j1,j2,xlo,xhi,iflux,jflux,bedlevel,b,h)
    
  use gparms
  use cntrl
  implicit none
  integer,  intent(in   ) :: i1,i2,j1,j2
  real(dp), intent(in   ) :: dx(2),xlo(2),xhi(2)
  real(dp), intent(in   ) :: iflux(i1-NFGST:i2+1+NFGST,j1-NFGST:j2+NFGST)
  real(dp), intent(in   ) :: jflux(j1-NFGST:j2+1+NFGST,i1-NFGST:i2+NFGST)
  real(dp), intent(inout) :: bedlevel(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
  real(dp), intent(inout) :: b(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
  real(dp), intent(in   ) :: h(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
  
  !local
  real(dp) :: delta(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
  integer  :: i,j
  real(dp) :: oodx,oody,xc

  !precompute 1/dx,  1/dy
  oodx = one/dx(1);
  oody = one/dx(2);
  
  !store previous bed thickness for use in morphodynamic adjustment
  delta = zero
  !bedlevel_last = bedlevel  

  !loop over internal cells, update state variables using convervative diff on fluxes
  do i=i1,i2
	do j=j1,j2
      delta(i,j) =  &
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
 

  !ensure bedlevel does not change at inflow in Hibma case
   if(caseid==hibma)then
		if(abs(xlo(1))<1)then
	      delta(1,:) = zero
	    endif
	  endif

   !ensure that after modification we are maintaining some kind of minimum depth
   if(caseid .ne. hibma)then
	 delta = min(delta,h-min_morph_depth)  !kind of a fudge
   endif
!  endif
  ! 	do i=i1,i2
  ! 		do j=j1,j2
  ! 			if(bedlevel(i,j)-bedlevel_last(i,j) > (h(i,j)+.10) )then
  ! 			  bedlevel(i,j) = bedlevel_last(i,j)
  ! 			endif
  ! 		end do
  ! 	end do
	
  !update bedlevel
  bedlevel = bedlevel + delta 

 !  !!fudge, we are getting an accretion singularity at the wall of the inlet leading to bottom penetrating the FS
 !    if(caseid==roelvink)then
 !      	do i=i1,i2
 !    	   xc = xlo(1)+dx(1)*dble(i-i1)+dx(1)/2
 !    	   if(xc >=5000 .and. xc <= 5250)then
 !    		  do j=j1,j2
 !    	       bedlevel(i,j) = min(bedlevel(i,j),.8)
 !    	     end do
 !    	   endif
 !    	 end do
 !    	endif

  !update bathymetry using change in bed level morphodynamics is active
  if(sedmodel == 2)then
    b = b + delta
  endif



  return
end subroutine consdiff_sed

