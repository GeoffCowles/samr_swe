!==============================================================================
! 
!                                  calc_diffusion
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
! Comments:     Add diffusion term to equations
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

subroutine calc_diffusion(dt,dx,i1,i2,j1,j2,h,vh)

  use gparms
  use cntrl
  implicit none
  
  integer,  intent(in )    :: i1,i2,j1,j2
  real(dp), intent(in )    :: dt
  real(dp), intent(in )    :: dx(2)
  real(dp), intent(in   )  :: h(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
	real(dp), intent(inout)  :: vh(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST,NDIMS)
  

  !----- local ------------
  integer  :: i,j,diff_equation
  real(dp) :: ediff(i1-1:i2+1,j1-1:j2+1)  ! eddy diffusivity in m^2/s 
  real(dp) :: taux(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
	real(dp) :: tauy(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
  
  real(dp) :: fac1,fac2,oodx2,oody2,d2udx2,d2udy2,d2vdx2,d2vdy2,taub,maxdiff
  real(dp) :: oo2dx,oo2dy,dbar,dudx,dvdx,dudy,dvdy,oodx,oody


  if(diffusivity_type == 0) return

  !initialize 
  ediff = zero
  fac1  = gravity*C_manning*C_manning
  
 
  !set the eddy viscosity
  select case(diffusivity_type)

  case(constant_diffusivity) !use a constant user-prescribed eddy viscosity
    ediff = diffusivity_coefficient
  case(smagorinsky_diffusivity) !coef * L^2 * sqrt( (du/dx-dv/dy)^2 + (du/dy+dv/dx)^2 )
    write(*,*)'this needs more debugging'
    stop
    fac2 = dx(1)*dx(2)  !L^2
    oo2dx = 1./(2*dx(1))
    oo2dy = 1./(2*dx(2))
    oodx  = 1./(dx(1))
    oody  = 1./(dx(2))
    do i=i1-1,i2+1
		do j=j1-1,j2+1
			dudx = oo2dx*(vh(i+1,j,1)/h(i+1,j) - vh(i-1,j,1)/h(i-1,j))
			dvdy = oo2dy*(vh(i,j+1,2)/h(i,j+1) - vh(i,j-1,2)/h(i,j-1))
			dudy = oody*max(abs(vh(i,j+1,1)/h(i,j+1)-vh(i,j,1)/h(i,j)),abs(vh(i,j,1)/h(i,j)-vh(i,j-1,1)/h(i,j-1)))
			dvdx = oodx*max(abs(vh(i+1,j,2)/h(i+1,j)-vh(i,j,2)/h(i,j)),abs(vh(i,j,2)/h(i,j)-vh(i-1,j,2)/h(i-1,j)))
			!dbar = sqrt( (dudx-dvdy)**2 + (dudy+dvdx)**2)  !http://mitgcm.org/public/r2_manual/latest/online_documents/node86.html
			dbar = sqrt( (dudy+dvdx)*(dudy+dvdx)) !Kantha and Clayson, page 125, shear only
			ediff(i,j) = smagorinsky_coefficient*fac2*dbar + diffusivity_coefficient
	  end do
    end do
  case(mixing_length_diffusivity) !see Castanedo paper referenced below
  	do i=i1-1,i2+1
		do j=j1-1,j2+1
		  fac2 = (h(i,j)**(-7./3.))*sqrt(vh(i,j,1)**2 + vh(i,j,2)**2)
			taux(i,j) = fac1*fac2*vh(i,j,1)
			tauy(i,j) = fac1*fac2*vh(i,j,2)   
			taub = sqrt(taux(i,j)**2 + tauy(i,j)**2)  !magnitude of bottom stress [m^2/s^2]
			ediff(i,j) = sqrt(taub)*h(i,j)*diffusivity_coefficient   !ustar * H * coefficient
	  end do
    end do
  end select

  !compute the diffusion.  There are three ways this is done depending on the treatment of H
  ! SEE: Castanedo et al., Models for the Turbulent Diffusion Terms of Shallow Water Equations, JHE 131 2005
  !    In their analysis, the second formulation most closely matches experiment + analytical theory
  !  D = d/dx(A*H*(du/dx))
  !  D = H(d/dx(A(du/dx))
  !  D = d/dx(A(duH/dx))
 
  oodx2 = 1./(dx(1)*dx(1))
  oody2 = 1./(dx(2)*dx(2))
  diff_equation = 2
  select case(diff_equation)
  case(1)
 	do j=j1,j2
 		do i=i1,i2
 			d2udx2 = oodx2*(ediff(i-1,j)*vh(i-1,j,1)/h(i-1,j) - 2*ediff(i,j)*vh(i,j,1)/h(i,j) + ediff(i+1,j)*vh(i+1,j,1)/h(i+1,j))
 			d2udy2 = oody2*(ediff(i,j-1)*vh(i,j-1,1)/h(i,j-1) - 2*ediff(i,j)*vh(i,j,1)/h(i,j) + ediff(i,j+1)*vh(i,j+1,1)/h(i,j+1))
 			d2vdx2 = oodx2*(ediff(i-1,j)*vh(i-1,j,2)/h(i-1,j) - 2*ediff(i,j)*vh(i,j,2)/h(i,j) + ediff(i+1,j)*vh(i+1,j,2)/h(i+1,j))
 			d2vdy2 = oody2*(ediff(i,j-1)*vh(i,j-1,2)/h(i,j-1) - 2*ediff(i,j)*vh(i,j,2)/h(i,j) + ediff(i,j+1)*vh(i,j+1,2)/h(i,j+1))
 			vh(i,j,1) = vh(i,j,1) + dt*h(i,j)*(d2udx2 + d2udy2)
 			vh(i,j,2) = vh(i,j,2) + dt*h(i,j)*(d2vdx2 + d2vdy2)
 		end do
 	end do
  case(2)
  do j=j1,j2
 		do i=i1,i2
 			d2udx2 = oodx2*(ediff(i-1,j)*vh(i-1,j,1) - 2*ediff(i,j)*vh(i,j,1) + ediff(i+1,j)*vh(i+1,j,1))
 			d2udy2 = oody2*(ediff(i,j-1)*vh(i,j-1,1) - 2*ediff(i,j)*vh(i,j,1) + ediff(i,j+1)*vh(i,j+1,1))
 			d2vdx2 = oodx2*(ediff(i-1,j)*vh(i-1,j,2) - 2*ediff(i,j)*vh(i,j,2) + ediff(i+1,j)*vh(i+1,j,2))
 			d2vdy2 = oody2*(ediff(i,j-1)*vh(i,j-1,2) - 2*ediff(i,j)*vh(i,j,2) + ediff(i,j+1)*vh(i,j+1,2))
 			vh(i,j,1) = vh(i,j,1) + dt*(d2udx2 + d2udy2)
 			vh(i,j,2) = vh(i,j,2) + dt*(d2vdx2 + d2vdy2)
 		end do
 	end do
  end select
  

  return

end subroutine calc_diffusion
