!==============================================================================
! 
!                                  init
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
! Comments:     Initialize Data on the Patch
!
!       Output:
!          [h,uh,vh]:  flow variables on the patch
!          b        :  bathymetry (positive down)
!          
!==============================================================================

subroutine initflow(cid,dx,xlo,xhi,i1,i2,j1,j2,igst,jgst,h,vh,b,bedlevel)  

  use gparms
  use cntrl
  implicit none
  
  integer,  intent(in   ) :: cid,i1,i2,j1,j2,igst,jgst
  real(dp), intent(in   ) :: dx(2),xlo(2),xhi(2)
  real(dp), intent(inout) :: h(i1-igst:i2+igst,j1-jgst:j2+jgst)
  real(dp), intent(inout) :: vh(i1-igst:i2+igst,j1-jgst:j2+jgst,1:NDIMS)
  real(dp), intent(in   ) :: b(i1-igst:i2+igst,j1-jgst:j2+jgst)
  real(dp), intent(inout) :: bedlevel(i1-igst:i2+igst,j1-jgst:j2+jgst)

  !----- local ------------
  integer  :: i,j
  real(dp) :: xhalf,xc,yc,depth_left,depth_right,x0,xberm,fac,fac2,zeta,u_left
  real(dp) :: theta_scrit,t1,t2,t3,t4,f1,f2
  real(dp) :: devriend_x,devriend_y,devriend_amp,devriend_rad,dist
  

 

  
  caseid = cid

  bedlevel = 0.
  
 
  select case(caseid)
  

  !------------------------------------------------------------
  case(rossby)  !Rossby soliton
  !------------------------------------------------------------

 ! integer ic0,ic1,imid
 !  REAL B,A,xc,yc,expy,phi,dphidx 
 ! 
 !  imid = .5*(ilast0+ifirst0+1)
 !  B = .395
 !  A =  .771*B*B 
 !  do il1,il2
 !    do j=jl1,jl2
 !      xc = xlo(0) + dx(0)*(dble(ic0-ifirst0)+half) 
 !      yc = xlo(1) + dx(1)*(dble(ic1-ifirst1)+half) 
 !      phi = A*(1./cosh(B*xc))*(1./cosh(B*xc))
 !      dphidx = -2*B*tanh(B*xc)*phi
 !      expy = exp(-yc*yc/2.) 
 !      depth(ic0,ic1)      =  phi*((3+6*yc*yc)/4.)*expy+1.
 !        velocity(ic0,ic1,0) =  phi*((-9.+6*yc*yc)/4)*expy
 !        velocity(ic0,ic1,1) =  2.*yc*dphidx*expy 
 !      end do
 !    end do
  
  !------------------------------------------------------------
  case(dambreakx) !Dambreak  (Toro Test 1)
  !------------------------------------------------------------
  depth_left  = 1.0
  u_left = 2.5
  depth_right = 0.1
  x0 = 10.0
  do i=i1,i2
	do j=j1,j2
      xc = xlo(1) + dx(1)*dble(i-i1)+dx(1)/2
      if(xc <= x0)then
	    h(i,j)    = depth_left
	    vh(i,j,1) = u_left*depth_left
	    vh(i,j,2) = zero
	    ! b(i,j)    = zero
	  else
		h(i,j)    = depth_right
	    vh(i,j,1) = zero
	    vh(i,j,2) = zero
	    ! b(i,j)    = zero
	  endif
	end do 
  end do


  !------------------------------------------------------------
  case(dambreakx_dry) !Dambreak into dry region Toro page 120 Case3
  !------------------------------------------------------------
  depth_left  = 1.0
  depth_right = 0.0
  x0 = 20.0
  do i=i1,i2
	do j=j1,j2
      xc = xlo(1) + dx(1)*dble(i-i1)+dx(1)/2
      if(xc <= x0)then
	    h(i,j)    = depth_left
	    vh(i,j,1) = zero
	    vh(i,j,2) = zero
	    ! b(i,j)    = zero
	  else
		h(i,j)    = depth_right
	    vh(i,j,1) = zero
	    vh(i,j,2) = zero
	    ! b(i,j)    = zero
	  endif
	end do 
  end do

  !------------------------------------------------------------
  case(dambreaky_dry) !Dambreak into dry region Toro page 120 Case3
  !------------------------------------------------------------
  depth_left  = 1.0
  depth_right = 0.0
  x0 = 20.0
  do i=i1,i2
	do j=j1,j2
      yc = xlo(2) + dx(2)*dble(j-j1)+dx(2)/2
      if(yc <= x0)then
	    h(i,j)    = depth_left
	    vh(i,j,1) = zero
	    vh(i,j,2) = zero
	    ! b(i,j)    = zero
	  else
		h(i,j)    = depth_right
	    vh(i,j,1) = zero
	    vh(i,j,2) = zero
	    ! b(i,j)    = zero
	  endif
	end do 
  end do

  !------------------------------------------------------------
  case(dambreakx_berm) !Dambreak into dry with berm
  !------------------------------------------------------------
  depth_left  = 1.0
  depth_right = 0.0
  x0    = 10.
  xberm = 25.
  do i=i1,i2
	do j=j1,j2
      xc = xlo(1) + dx(1)*dble(i-i1)+dx(1)/2
	  ! b(i,j)   = max(zero, 0.25-an8th*abs(xc-xberm))
      if(xc <= x0)then
	    h(i,j)    = depth_left
	    vh(i,j,1) = zero
	    vh(i,j,2) = zero
	  else
		h(i,j)    = depth_right
	    vh(i,j,1) = zero
	    vh(i,j,2) = zero
	  endif
	end do 
  end do

  !------------------------------------------------------------
  case(dambreaky_berm) !Dambreak into dry with berm
  !------------------------------------------------------------
  depth_left  = 1.0
  depth_right = 0.0
  x0    = 10.
  xberm = 25.
  do i=i1,i2
	do j=j1,j2
      yc = xlo(2) + dx(2)*dble(j-j1)+dx(2)/2
	  ! b(i,j)   = max(zero, 0.25-an8th*abs(yc-xberm))
      if(yc <= x0)then
	    h(i,j)    = depth_left
	    vh(i,j,1) = zero
	    vh(i,j,2) = zero
	  else
		h(i,j)    = depth_right
	    vh(i,j,1) = zero
	    vh(i,j,2) = zero
	  endif
	end do 
  end do
  

  !------------------------------------------------------------
  case(dambreaky) !Dambreak  (Sod Problem for SWE)
  !------------------------------------------------------------
  depth_left  = 1.0
  depth_right = 0.5
  xhalf = 0.0
  do i=i1,i2
  	do j=j1,j2
      xc = xlo(2) + dx(2)*dble(j-j1)+dx(2)/2 
      if(xc <= xhalf)then
   	    h(i,j)    = depth_left
	    vh(i,j,1) = zero
	    vh(i,j,2) = zero
	    ! b(i,j)    = zero
	  else
		h(i,j)    = depth_right
	    vh(i,j,1) = zero
	    vh(i,j,2) = zero
	    ! b(i,j)    = zero
	  endif
	end do 
  end do


!------------------------------------------------------------
 case(dambreak2D) !Dambreak  (2D Sod Problem for SWE)
!------------------------------------------------------------
  depth_left  = 1.0
  depth_right = 0.5
  xhalf = 0.0
  do i=i1,i2
    do j=j1,j2
      xc = xlo(1) + dx(1)*dble(i-i1)+dx(1)/2 !x direction
      yc = xlo(2) + dx(2)*dble(j-j1)+dx(2)/2 !y direction
      if(sqrt(xc**2 + yc**2) < 10.0)then
        h(i,j)    = depth_left
        vh(i,j,1) = zero
        vh(i,j,2) = zero
        ! b(i,j)    = zero
      else
        h(i,j)    = depth_right
        vh(i,j,1) = zero
        vh(i,j,2) = zero
        ! b(i,j)    = zero
      endif
    end do 
  end do

  case(ramp) !Supercritical flow past Ramp
    write(*,*)'init for ramp not yet setup'
    stop

  case(user_defined) !User defined
  depth_left = 1.0
  depth_right = 0.5
  do i=i1,i2
	do j=j1,j2
      xc = xlo(1) + dx(1)*dble(i-i1)+dx(1)/2
      yc = xlo(2) + dx(2)*dble(j-j1)+dx(1)/2
      write(*,*)'init for ramp not yet setup'
      stop
	end do 
  end do

  !---------------------------------------------------------------------
  case(nomotion) !quiescent flow with bathy gradients, all submerged
  !---------------------------------------------------------------------
  do i=i1,i2
	do j=j1,j2
      xc = xlo(1) + dx(1)*dble(i-i1)+dx(1)/2
      ! xc = xc**2
      !       yc = xlo(2) + dx(2)*(dble(j-j1)+dx(2)/2)
      !       yc = yc**2
      !       h(i,j)      = 0.5*(yc+xc)
      ! b(i,j)      = xc
      h(i,j)       = 2.0-b(i,j)
      vh(i,j,1) = zero
      vh(i,j,2) = zero
	end do 
  end do


  !---------------------------------------------------------------------
  case(nomotiondryx) !quiescent flow with bathy gradients and dry region
  !---------------------------------------------------------------------
  do i=i1,i2
	do j=j1,j2
      xc = xlo(1) + dx(1)*dble(i-i1)+dx(1)/2
      ! b(i,j)    = xc
      h(i,j)    = ahalf-b(i,j)
      if(h(i,j) < 0) h(i,j)=zero
      vh(i,j,1) = zero
      vh(i,j,2) = zero
	end do 
  end do

  !---------------------------------------------------------------------
  case(nomotiondryy) !quiescent flow with bathy gradients and dry region
  !---------------------------------------------------------------------
  do i=i1,i2
	do j=j1,j2
      yc = xlo(2) + dx(2)*dble(j-j1)+dx(2)/2
      ! b(i,j)    = yc
      h(i,j)    = ahalf-b(i,j)
      if(h(i,j) < 0) h(i,j)=zero
      vh(i,j,1) = zero
      vh(i,j,2) = zero
	end do 
  end do

  !---------------------------------------------------------------------
  case(threehump) !Liang and Borthwick
  !---------------------------------------------------------------------
  do i=i1,i2
   	do j=j1,j2
      xc = xlo(1)+dx(1)*dble(i-i1)+dx(1)/2
      yc = xlo(2)+dx(2)*dble(j-j1)+dx(2)/2
      ! b(i,j)   = max(zero, 1-an8th*sqrt((xc-30)**2 + (yc-6 )**2), & 
      !                                1-an8th*sqrt((xc-30)**2 + (yc-24)**2), & 
      !                                3-.3*sqrt((xc-47.5)**2 + (yc-15)**2))
      h(i,j)      = 1.875-b(i,j)
      if(xc > 16.)then
	    h(i,j) = zero
	  endif
      vh(i,j,1) = zero
      vh(i,j,2) = zero
 	end do 
  end do

  !---------------------------------------------------------------------
  case(threehumpy) !Liang and Borthwick in y-direction
  !---------------------------------------------------------------------
  do i=i1,i2
   	do j=j1,j2
      xc = xlo(1)+dx(1)*dble(i-i1)+dx(1)/2
      yc = xlo(2)+dx(2)*dble(j-j1)+dx(2)/2
      ! b(i,j)   = max(zero, 1-an8th*sqrt((yc-30)**2 + (xc-6 )**2), & 
      !                                1-an8th*sqrt((yc-30)**2 + (xc-24)**2), & 
      !                                3-.3*sqrt((yc-47.5)**2 + (xc-15)**2))
      h(i,j)      = 1.875-b(i,j)
      if(yc > 16.)then
	    h(i,j) = zero
	  endif
      vh(i,j,1) = zero
      vh(i,j,2) = zero
 	end do 
  end do

  !---------------------------------------------------------------------
  case(threehump_nomotion) !Liang and Borthwick bathymetry no motion test
  !---------------------------------------------------------------------
  do i=i1,i2
    do j=j1,j2
      xc = xlo(1)+dx(1)*dble(i-i1)+dx(1)/2
      yc = xlo(2)+dx(2)*dble(j-j1)+dx(2)/2
      ! b(i,j)   = max(zero, 1-an8th*sqrt((xc-30)**2 + (yc-6 )**2), & 
      !                            1-an8th*sqrt((xc-30)**2 + (yc-24)**2), & 
      !                            3-.3*sqrt((xc-47.5)**2 + (yc-15)**2))
      !h(i,j)       = 4.0 - b(i,j) !all wet
      h(i,j)      = max(1.0-b(i,j),zero) !some dry
      vh(i,j,1) = zero
      vh(i,j,2) = zero
 	end do 
  end do

  !---------------------------------------------------------------------
  case(sampson) !Sampson parabolic bottom topography analytical
  !see Liang & Borthwick, Computers and Fluids, 38, 2009 page 231
  !---------------------------------------------------------------------
  
   p_samp = sqrt(8.0_dp*gravity*h0_samp)/a_samp
   s_samp = ahalf*sqrt(p_samp**2 - tau_samp**2)   
   fac = (a_samp**2*B_samp**2)/(8*gravity**2*h0_samp) 
   fac2 = (a4th*tau_samp**2 -s_samp**2)
   
   write(*,*)fac,fac2,tau_samp,s_samp
   do i=i1,i2
      do j=j1,j2
        xc = xlo(1)+dx(1)*dble(i-i1)+dx(1)/2
        yc = xlo(2)+dx(2)*dble(j-j1)+dx(2)/2
        zeta = h0_samp + fac*fac2 - B_samp**2/(4*gravity) - (xc/gravity)*B_samp*s_samp
        h(i,j) = max(zeta-b(i,j),zero)
        vh(i,j,1) = zero
        vh(i,j,2) = zero
   	end do 
    end do

  !---------------------------------------------------------------------
  case(tidetest) !time-dependent bcs
  !---------------------------------------------------------------------
  

   do i=i1,i2
      do j=j1,j2
	    h(i,j) = 5.
        vh(i,j,1) = zero
        vh(i,j,2) = zero
      end do 
    end do

  !---------------------------------------------------------------------
  case(conrun) !briggs et al. conical runup test - see begnudelli/valiani
  !---------------------------------------------------------------------


   do i=i1,i2
      do j=j1,j2
	    h(i,j) = max(.32-b(i,j),zero)
        vh(i,j,1) = zero
        vh(i,j,2) = zero
      end do 
    end do

  !---------------------------------------------------------------------
  case(heniche) !heniche non-constant slope test, see Brufau etal, IJNMF v39
  !---------------------------------------------------------------------

   do i=i1,i2
      do j=j1,j2
	    h(i,j) = max(1.75-b(i,j),zero)
        vh(i,j,1) = zero
        vh(i,j,2) = zero
      end do 
    end do

  !---------------------------------------------------------------------
  case(step) !from SAMRAI Euler step-2d example, but for SWE
  !---------------------------------------------------------------------

  do i=i1,i2
    do j=j1,j2
	  h(i,j) = 0.1
      vh(i,j,1) = h(i,j)*3.0
      vh(i,j,2) = 0.0
    end do 
  end do

  !---------------------------------------------------------------------
  case(supercrit) !supercritical flow over ramp
  !---------------------------------------------------------------------
  theta_scrit = 8.95*3.14159/180.

  do i=i1,i2
    do j=j1,j2
	  h(i,j) = 1
      vh(i,j,1) = 8.57
      vh(i,j,2) = 0.0
      xc = xlo(1)+dx(1)*dble(i-i1)+dx(1)/2
	  yc = xlo(2)+dx(2)*dble(j-j1)+dx(2)/2
	  if(xc > 10. .and. yc < (xc-10)*sin(theta_scrit))then
        h(i,j) = zero
        vh(i,j,1) = zero
      endif
    end do 
  end do

  !---------------------------------------------------------------------
  case(roelvink) !inlet roelvink test (Coastal Engineering 53, 2006)
  !---------------------------------------------------------------------
   do i=i1,i2
      do j=j1,j2
        xc = xlo(1)+dx(1)*dble(i-i1)+dx(1)/2
		if(xc < 5000)then
		  fac    = xc/5000.
		  h(i,j) = 10*(1-fac) + 2*(fac)
		else
	      h(i,j) = 2
		end if
        vh(i,j,1) = zero
        vh(i,j,2) = zero
      end do
    end do

	!---------------------------------------------------------------------
	case(roelvinky) !inlet roelvink test (Coastal Engineering 53, 2006)
	!---------------------------------------------------------------------
	 do i=i1,i2
	    do j=j1,j2
	      yc = xlo(2)+dx(2)*dble(j-j1)+dx(2)/2
			if(yc < 5000)then
			  fac    = yc/5000.
			  h(i,j) = 10*(1-fac) + 2*(fac)
			else
		      h(i,j) = 2
			end if
	      vh(i,j,1) = zero
	      vh(i,j,2) = zero
	    end do
	  end do
	
	!---------------------------------------------------------------------
	case(slosh_inlet) !roelvink test with no open boundary and flat bottom
	!---------------------------------------------------------------------
	 do i=i1,i2
	    do j=j1,j2
			xc = xlo(1)+dx(1)*dble(i-i1)+dx(1)/2
			if(xc < 2500)then
			  h(i,j) = 1.2
			else
		      h(i,j) = 1.
			end if
			!h(i,j) = 1.  !uncomment this for zero motion slosh_inlet
	      vh(i,j,1) = zero
	      vh(i,j,2) = zero
	    end do
	  end do
	
	!---------------------------------------------------------------------
	case(trench) !warner migrating trench case
	!---------------------------------------------------------------------
	t1 = trench_mid-(trench_half_width+trench_slope_width);
	t2 = trench_mid-trench_half_width;
	t3 = trench_mid+trench_half_width;
	t4 = trench_mid+(trench_half_width+trench_slope_width);
	
	 do i=i1,i2
	    do j=j1,j2
			xc = xlo(1)+dx(1)*dble(i-i1)+dx(1)/2
			
			if(xc > t2 .and. xc < t3)then
			  h(i,j) = bath_trench2
			else if (xc > t1 .and. xc < t2)then
				f1 = (xc-t1)/trench_slope_width
				f2 = 1-f1;
				h(i,j) = f2*bath_trench1 + f1*bath_trench2
			else if (xc > t3 .and. xc < t4)then
				f1 = (xc-t3)/trench_slope_width;
				f2 = 1-f1;
				h(i,j) = f1*bath_trench1 + f2*bath_trench2
			else
				h(i,j) = bath_trench1
			end if
	      vh(i,j,1) = 0.2
	      vh(i,j,2) = 0.0
	    end do
	   
	  end do

	!---------------------------------------------------------------------
	case(devriend) !de Vriend hump morphodynamic case
	!---------------------------------------------------------------------
	
	devriend_amp = 5.    !amplitude of hump
   devriend_x   = 5000. !x location of hump center
   devriend_y   = 0000. !y location of hump center
   devriend_rad = 1000. !radius of hump
   do i=i1-igst,i2+igst
	 do j=j1-jgst,j2+jgst
      xc = xlo(1)+dx(1)*dble(i-i1)+dx(1)/2
      yc = xlo(2)+dx(2)*dble(j-j1)+dx(2)/2
      dist = (xc-devriend_x)**2 +(yc-devriend_y)**2
      h(i,j) = 10.-devriend_amp*exp(-dist/(2*devriend_rad*devriend_rad))
      vh(i,j,1) = 10.
      vh(i,j,2) = 0.0
 	 end do 
   end do

	!---------------------------------------------------------------------
	case(trenchy) !warner migrating trench case in y-direction
	!---------------------------------------------------------------------
	t1 = trench_mid-(trench_half_width+trench_slope_width);
	t2 = trench_mid-trench_half_width;
	t3 = trench_mid+trench_half_width;
	t4 = trench_mid+(trench_half_width+trench_slope_width);

	 do i=i1,i2
	    do j=j1,j2
		   yc = xlo(2)+dx(2)*dble(j-j1)+dx(2)/2
			if(yc > t2 .and. yc < t3)then
			  h(i,j) = bath_trench2
			else if (yc > t1 .and. yc < t2)then
				f1 = (yc-t1)/trench_slope_width
				f2 = 1-f1;
				h(i,j) = f2*bath_trench1 + f1*bath_trench2
			else if (yc > t3 .and. yc < t4)then
				f1 = (yc-t3)/trench_slope_width;
				f2 = 1-f1;
				h(i,j) = f1*bath_trench1 + f2*bath_trench2
			else
				h(i,j) = bath_trench1
			end if
	      vh(i,j,1) = 0.0
	      vh(i,j,2) = 0.2
	    end do
   
  end do

	!---------------------------------------------------------------------
	case(uniform) !simple uniform flow
	!---------------------------------------------------------------------
	

	 do i=i1,i2
	    do j=j1,j2
		   xc = xlo(1)+dx(1)*dble(i-i1)+dx(1)/2
		   yc = xlo(2)+dx(2)*dble(j-j1)+dx(2)/2
		   vh(i,j,1) = uniform_h*uniform_u
		   vh(i,j,2) = 0.0
		   h(i,j) = uniform_h
		end do
  end do
	
 
  end select

  return

end subroutine initflow
