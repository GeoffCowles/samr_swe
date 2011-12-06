!==============================================================================
! 
!                             set bathymetry
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

subroutine setbathy(cid,dx,xlo,xhi,i1,i2,j1,j2,igst,jgst,b)  

  use gparms
  use cntrl
  implicit none
  
  integer,  intent(in   ) :: cid,i1,i2,j1,j2,igst,jgst
  real(dp), intent(in   ) :: dx(2),xlo(2),xhi(2)
  real(dp), intent(inout) :: b(i1-igst:i2+igst,j1-jgst:j2+jgst)

  !----- local ------------
  integer  :: i,j
  real(dp) :: xhalf,xc,yc,depth_left,depth_right,x0,xberm
  real(dp) :: inner_rad,outer_rad,con_height,fac,rad,theta_scrit
  real(dp) :: t1,t2,t3,t4,f1,f2
  real(dp) :: devriend_x,devriend_y,devriend_amp,devriend_rad,dist
  real(dp) :: slope
  integer*4 timeArray(3)    ! Holds the hour, minute, and second
  real rand
  caseid = cid
 
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
  case(dambreakx) !Dambreak  (Sod Problem for SWE)
  !------------------------------------------------------------
  
  b = 0.0
 


  !------------------------------------------------------------
  case(dambreakx_dry) !Dambreak into dry region Toro page 120 Case3
  !------------------------------------------------------------
  
  b = 0.0
 

  !------------------------------------------------------------
  case(dambreaky_dry) !Dambreak into dry region Toro page 120 Case3
  !------------------------------------------------------------
  
  b = 0.0
 

  !------------------------------------------------------------
  case(dambreakx_berm) !Dambreak into dry with berm
  !------------------------------------------------------------
  depth_left  = 1.0
  depth_right = 0.0
  x0    = 10.
  xberm = 25.
  do i=i1-igst,i2+igst
	do j=j1-jgst,j2+jgst
      xc = xlo(1) + dx(1)*dble(i-i1)+dx(1)/2
	  b(i,j)   = max(zero, 0.35-an8th*abs(xc-xberm))
	end do 
  end do

  !------------------------------------------------------------
  case(dambreaky_berm) !Dambreak into dry with berm
  !------------------------------------------------------------
  depth_left  = 1.0
  depth_right = 0.0
  x0    = 10.
  xberm = 25.
  do i=i1-igst,i2+igst
	do j=j1-jgst,j2+jgst
      yc = xlo(2) + dx(2)*dble(j-j1)+dx(2)/2
	  b(i,j)   = max(zero, 0.35-an8th*abs(yc-xberm))
	end do 
  end do
  

  !------------------------------------------------------------
  case(dambreaky) !Dambreak  (Sod Problem for SWE)
  !------------------------------------------------------------
  
  b = zero
  
!------------------------------------------------------------
 case(dambreak2D) !Dambreak  (2D Sod Problem for SWE)
!------------------------------------------------------------
  b = zero

  case(ramp) !Supercritical flow past Ramp
    write(*,*)'init for ramp not yet setup'
    stop

  case(user_defined) !User defined
  depth_left = 1.0
  depth_right = 0.5
  do i=i1-igst,i2+igst
	do j=j1-jgst,j2+jgst
      xc = xlo(1) + dx(1)*dble(i-i1)+dx(1)/2
      yc = xlo(2) + dx(2)*dble(j-j1)+dx(1)/2
      write(*,*)'init for ramp not yet setup'
      stop
	end do 
  end do

  !---------------------------------------------------------------------
  case(nomotion) !quiescent flow with bathy gradients, all submerged
  !---------------------------------------------------------------------
  
 do i=i1-igst,i2+igst
	do j=j1-jgst,j2+jgst
      xc = xlo(1) + dx(1)*dble(i-i1)+dx(1)/2
      b(i,j)      = xc
	end do 
  end do


  !---------------------------------------------------------------------
  case(nomotiondryx) !quiescent flow with bathy gradients and dry region
  !---------------------------------------------------------------------
  do i=i1-igst,i2+igst
	do j=j1-jgst,j2+jgst
      xc = xlo(1) + dx(1)*dble(i-i1)+dx(1)/2
      b(i,j)    = xc
	end do 
  end do

  !---------------------------------------------------------------------
  case(nomotiondryy) !quiescent flow with bathy gradients and dry region
  !---------------------------------------------------------------------
  do i=i1-igst,i2+igst
	do j=j1-jgst,j2+jgst
      yc = xlo(2) + dx(2)*dble(j-j1)+dx(2)/2
      b(i,j) = yc
	end do 
  end do

  !---------------------------------------------------------------------
  case(threehump) !Liang and Borthwick
  !---------------------------------------------------------------------
  do i=i1-igst,i2+igst
	do j=j1-jgst,j2+jgst
      xc = xlo(1)+dx(1)*dble(i-i1)+dx(1)/2
      yc = xlo(2)+dx(2)*dble(j-j1)+dx(2)/2
      b(i,j)   = max(zero, 1-an8th*sqrt((xc-30)**2 + (yc-6 )**2), & 
                               1-an8th*sqrt((xc-30)**2 + (yc-24)**2), & 
                               3-.3*sqrt((xc-47.5)**2 + (yc-15)**2))
 	end do 
  end do

  !---------------------------------------------------------------------
  case(threehumpy) !Liang and Borthwick in y-direction
  !---------------------------------------------------------------------
  do i=i1-igst,i2+igst
	do j=j1-jgst,j2+jgst
      xc = xlo(1)+dx(1)*dble(i-i1)+dx(1)/2
      yc = xlo(2)+dx(2)*dble(j-j1)+dx(2)/2
      b(i,j)   = max(zero, 1-an8th*sqrt((yc-30)**2 + (xc-6 )**2), & 
                               1-an8th*sqrt((yc-30)**2 + (xc-24)**2), & 
                               3-.3*sqrt((yc-47.5)**2 + (xc-15)**2))
 	end do 
  end do

  !---------------------------------------------------------------------
  case(threehump_nomotion) !Liang and Borthwick bathymetry no motion test
  !---------------------------------------------------------------------
  do i=i1-igst,i2+igst
	do j=j1-jgst,j2+jgst
      xc = xlo(1)+dx(1)*dble(i-i1)+dx(1)/2
      yc = xlo(2)+dx(2)*dble(j-j1)+dx(2)/2
      b(i,j)   = max(zero, 1-an8th*sqrt((xc-30)**2 + (yc-6 )**2), & 
                             1-an8th*sqrt((xc-30)**2 + (yc-24)**2), & 
                             3-.3*sqrt((xc-47.5)**2 + (yc-15)**2))
 	end do 
  end do

  !---------------------------------------------------------------------
  case(sampson) !Liang and Borthwick sampson test
  !---------------------------------------------------------------------
  do i=i1-igst,i2+igst
	do j=j1-jgst,j2+jgst
      xc = xlo(1)+dx(1)*dble(i-i1)+dx(1)/2
      b(i,j) = h0_samp*((xc/a_samp)**2)
 	end do 
  end do

  !---------------------------------------------------------------------
  case(tidetest) 
  !---------------------------------------------------------------------
  do i=i1-igst,i2+igst
	do j=j1-jgst,j2+jgst
      b(i,j) = 0.
 	end do 
  end do

  !---------------------------------------------------------------------
  case(conrun) 
  !---------------------------------------------------------------------
  inner_rad = 1.1
  outer_rad = 3.6
  con_height = .625
 
  do i=i1-igst,i2+igst
	do j=j1-jgst,j2+jgst
	  xc = xlo(1)+dx(1)*dble(i-i1)+dx(1)/2
	  yc = xlo(2)+dx(2)*dble(j-j1)+dx(2)/2
	  rad = sqrt(xc**2 + yc**2)
	  if(rad > outer_rad)then
	    b(i,j) = zero
	  elseif(rad < inner_rad)then
	    b(i,j) = con_height
	  else
	    fac  = (rad-inner_rad)/(outer_rad-inner_rad)
	    b(i,j) = (1-fac)*(con_height)
	  endif
 	end do 
  end do
  
  !---------------------------------------------------------------------
  case(heniche) 
  !---------------------------------------------------------------------
  inner_rad = 1.1
  outer_rad = 3.6
  con_height = .625
 
  do i=i1-igst,i2+igst
	do j=j1-jgst,j2+jgst
	  xc = xlo(1)+dx(1)*dble(i-i1)+dx(1)/2
	  yc = xlo(2)+dx(2)*dble(j-j1)+dx(2)/2
	  if(xc < 300)then
		b(i,j) = .001*(xc-zero)
	  else if(xc >=300 .and. xc < 400 )then
		b(i,j) = (.001*300) + (xc-300)*.01
	  else
		b(i,j) = (.001*300) + (.01*100) + (xc-400)*.001
	  endif
 	end do 
  end do

  !---------------------------------------------------------------------
  case(step) 
  !---------------------------------------------------------------------
  do i=i1-igst,i2+igst
	do j=j1-jgst,j2+jgst
      b(i,j) = 0.
 	end do 
  end do

  !---------------------------------------------------------------------
  case(supercrit) 
  !---------------------------------------------------------------------
  theta_scrit = 8.95*3.14159/180.
  do i=i1-igst,i2+igst
	do j=j1-jgst,j2+jgst
	  xc = xlo(1)+dx(1)*dble(i-i1)+dx(1)/2
	  yc = xlo(2)+dx(2)*dble(j-j1)+dx(2)/2
	  b(i,j) = 0.0
	  if(xc > 10. .and. yc < (xc-10)*sin(theta_scrit))then
        b(i,j) = 20.
      endif
 	end do 
  end do

  !---------------------------------------------------------------------
  case(roelvink)
  !---------------------------------------------------------------------
  do i=i1-igst,i2+igst
     do j=j1-jgst,j2+jgst
	  xc = xlo(1)+dx(1)*dble(i-i1)+dx(1)/2
	  if(xc < 5000)then
	    fac    = xc/5000.
	    b(i,j) = -10*(1-fac) + -2*(fac)
	  else
		b(i,j) = -2
	  end if
     end do
  end do

	!---------------------------------------------------------------------
	case(roelvinky)
	!---------------------------------------------------------------------
	do i=i1-igst,i2+igst
	   do j=j1-jgst,j2+jgst
		  yc = xlo(2)+dx(2)*dble(j-j1)+dx(2)/2
		  if(yc < 5000)then
		    fac    = yc/5000.
		    b(i,j) = -10*(1-fac) + -2*(fac)
		  else
			b(i,j) = -2
		  end if
	   end do
	end do

	!---------------------------------------------------------------------
	case(slosh_inlet)
	!---------------------------------------------------------------------
	do i=i1-igst,i2+igst
	   do j=j1-jgst,j2+jgst
		  b(i,j) = 0.0
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
			  b(i,j) = bath_trench2
			else if (xc > t1 .and. xc < t2)then
				f1 = (xc-t1)/trench_slope_width
				f2 = 1-f1;
				b(i,j) = f2*bath_trench1 + f1*bath_trench2
			else if (xc > t3 .and. xc < t4)then
				f1 = (xc-t3)/trench_slope_width;
				f2 = 1-f1;
				b(i,j) = f1*bath_trench1 + f2*bath_trench2
			else
				b(i,j) = bath_trench1
			endif
	      b(i,j) = -b(i,j)+bath_trench1
	    end do
	  end do
	
  !---------------------------------------------------------------------
  case(devriend) !morphodynamics of gaussian hump
  !x = 0:10:20000;
  ! y = 0:10:10000;
  ! [X,Y] = meshgrid(x,y);
  ! dist = (X-5000).^2 + (Y-5000).^2;
  ! Z = 5.*exp(-dist/(2*1000*1000));
  ! pcolor(X,Y,Z)
  ! shading interp
  ! axis equal
  !---------------------------------------------------------------------
  devriend_amp = 5.    !amplitude of hump
  devriend_x   = 5000. !x location of hump center
  devriend_y   = 5000. !y location of hump center
  devriend_rad = 1000. !radius of hump
  if(frictype==0)then
	 fac = 0.
  else
	 fac = 1.
  endif
  do i=i1-igst,i2+igst
	do j=j1-jgst,j2+jgst
      xc = xlo(1)+dx(1)*dble(i-i1)+dx(1)/2
      yc = xlo(2)+dx(2)*dble(j-j1)+dx(2)/2
      dist = (xc-devriend_x)**2 + (yc-devriend_y)**2
      b(i,j) = devriend_amp*exp(-dist/(2*devriend_rad*devriend_rad)) 
      b(i,j) = b(i,j) - fac*xc*4.64e-5 !add slope to counteract friction? (tau = .4667 at C_manning=.01)
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
			  b(i,j) = bath_trench2
			else if (yc > t1 .and. yc < t2)then
				f1 = (yc-t1)/trench_slope_width
				f2 = 1-f1;
				b(i,j) = f2*bath_trench1 + f1*bath_trench2
			else if (yc > t3 .and. yc < t4)then
				f1 = (yc-t3)/trench_slope_width;
				f2 = 1-f1;
				b(i,j) = f1*bath_trench1 + f2*bath_trench2
			else
				b(i,j) = bath_trench1
			endif
	      b(i,j) = -b(i,j)+bath_trench1
	    end do
	  end do
	
	!---------------------------------------------------------------------
	case(uniform) !simple uniform, subcritical flow
	!---------------------------------------------------------------------
    
    slope = -gravity*C_manning*C_manning*(uniform_h**(-7./3.))*((uniform_u*uniform_h)**2)
    slope =  slope/(uniform_h*gravity)
    
	 do i=i1,i2
	    do j=j1,j2  
		   xc = xlo(1)+dx(1)*dble(i-i1)+dx(1)/2
		   yc = xlo(2)+dx(2)*dble(j-j1)+dx(2)/2
	      b(i,j) = slope*xc
	      if(mod(i,2)==0)then
		       fac = .25
		    else
			     fac = -.25
			  endif
	      b(i,j) = b(i,j) + fac
	    end do
	  end do
  

  !---------------------------------------------------------------------
	case(hibma) !hibma channel shoal formation case in 80x2.5 km domain
	!---------------------------------------------------------------------
   fac = .75
   call itime(timeArray)     ! Get the current time
	 i = rand ( timeArray(1)+timeArray(2)+timeArray(3) )
	 do i=i1,i2
	    do j=j1,j2  
		   xc = xlo(1)+dx(1)*dble(i-i1)+dx(1)/2
		   yc = xlo(2)+dx(2)*dble(j-j1)+dx(2)/2
	      b(i,j) = -15 + 15*(xc/80000.) + 2*fac*(rand(0)-.5)
	      !set odd/even bathymetry in i-direction only
	      ! if(mod(i,2)==0)then
	      ! 		       fac = .25
	      ! 		    else
	      ! 			     fac = -.25
	      ! 			  endif
	      ! 	      b(i,j) = -15 + 15*(xc/80000.) + fac   !odd/even bathymetry
	    end do
	  end do
  
  end select

  return

end subroutine setbathy
