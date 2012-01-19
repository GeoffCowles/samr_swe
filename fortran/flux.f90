!==============================================================================
! 
!                                  hlle flux
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
! References:
!   "Finite Volume Methods and Adaptive Refinement for Tsunami Propagation
!	and Inundation." University of Washington, July 2006. 
!	David L George
!
! Comments:     Calculate flux for SWE using HLLE method
!
!       Input:
!                    dt:  global time step on the patch
!                    dx:  patch size (dx,dy)
!          [zeta,uh,vh]:  flow variables on the patch
!                    zb:  bathymetry measured from reference
!          
!       Output:
!          =>  ifluxm,jfluxm,ifluxp,jfluxp: fluctuations on the edges 
!==============================================================================

subroutine george_flux(cid,dt,dx,i1,i2,j1,j2,h,vh,b,ifluxm,jfluxm,ifluxp,jfluxp)  

 
  use gparms
  use cntrl
  implicit none

  integer,  intent(in)    :: cid
  integer,  intent(in)    :: i1,i2,j1,j2
  real(dp), intent(in )   :: dt
  real(dp), intent(in )   :: dx(2)
  real(dp), intent(in   ) :: h(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
  real(dp), intent(in   ) :: vh(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST,NDIMS)
  real(dp), intent(in )   :: b(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
  real(dp), intent(inout) :: ifluxm(i1-NFGST:i2+1+NFGST,j1-NFGST:j2+NFGST,NEQUS)
  real(dp), intent(inout) :: jfluxm(j1-NFGST:j2+1+NFGST,i1-NFGST:i2+NFGST,NEQUS)
  real(dp), intent(inout) :: ifluxp(i1-NFGST:i2+1+NFGST,j1-NFGST:j2+NFGST,NEQUS)
  real(dp), intent(inout) :: jfluxp(j1-NFGST:j2+1+NFGST,i1-NFGST:i2+NFGST,NEQUS)

  !----- local ------------
  integer  :: mx,my,maxm,mbc,meqn,mwaves,maux
!  integer, allocatable :: mthlim(:)
  integer  :: i,j,m
  integer  :: iclaw,jclaw,s2ci,s2cj
  real(dp) :: dtdx,dtdy
  real(dp), allocatable :: faddm(:,:)
  real(dp), allocatable :: faddp(:,:)
  real(dp), allocatable :: gaddm(:,:,:) 
  real(dp), allocatable :: gaddp(:,:,:)
  real(dp), allocatable :: q1d(:,:)
  real(dp), allocatable :: amdq(:,:)
  real(dp), allocatable :: apdq(:,:)
  real(dp), allocatable :: bmasdq(:,:)
  real(dp), allocatable :: bpasdq(:,:)
  real(dp), allocatable :: cqxx(:,:)
  real(dp), allocatable :: aux(:,:)
  real(dp), allocatable :: s(:,:)
  real(dp), allocatable :: fwave(:,:,:)
  
  !zero out fluxes (this might not be a good idea)
  ifluxm = zero
  ifluxp = zero
  jfluxm = zero
  jfluxp = zero

  !bind samrai environment to clawpack
  mx = i2-i1+1  !internal cells in x-dir
  my = j2-j1+1  !internal cells in y-dir
  maxm = max(mx,my) !max number of internal cells in any dir
  mbc = NCGST   !# of ghost layers
  meqn = 3      !# equations
  mwaves = 3    !# waves
  maux = 1      !# auxiliary vars
  dtdx = dt/dx(1)
  dtdy = dt/dx(2)
!  allocate(mthlim(mwaves)); mthlim = 4  !#limiter type on each wave

  !mapping of indices
  !clawpack internal x-dir cell index runs 1:mx whereas samrai runs i1:i2
  !clawpack internal y-dir cell index runs 1:my whereas samrai runs j1:j2
  !loop index is over samrai count, but we need to access claw-dimensioned vars
  s2ci = (1-i1)
  s2cj = (1-j1)
  
  !allocate workspace data and initialize
  allocate(faddm(1-mbc:maxm+mbc, meqn))
  allocate(faddp(1-mbc:maxm+mbc, meqn))
  allocate(gaddm(1-mbc:maxm+mbc, meqn, 2))
  allocate(gaddp(1-mbc:maxm+mbc, meqn, 2))
 
  allocate(q1d(1-mbc:maxm+mbc, meqn))
  allocate(amdq(1-mbc:maxm+mbc, meqn))
  allocate(apdq(1-mbc:maxm+mbc, meqn))
  allocate(bmasdq(1-mbc:maxm+mbc, meqn))
  allocate(bpasdq(1-mbc:maxm+mbc, meqn))
  allocate(cqxx(1-mbc:maxm+mbc, meqn))

  allocate(aux(1-mbc:maxm+mbc, maux))
  allocate(s(1-mbc:maxm+mbc, mwaves))
  allocate(fwave(1-mbc:maxm+mbc, meqn, mwaves))

!--------------------------------------------------------------------------
!  perform x-sweeps
!--------------------------------------------------------------------------
   
  do j=j1-1,j2+1  !j=0,my+1
	jclaw = j+s2cj
 
    !copy i-dir pencil at fixed j into 1d array
    do i=i1-NCGST,i2+NCGST
	  iclaw    = i+s2ci
	  q1d(iclaw,1) = h(i,j)     !h
	  q1d(iclaw,2) = vh(i,j,1)  !uh
	  q1d(iclaw,3) = vh(i,j,2)  !vh
	  aux(iclaw,1) = b(i,j)     !b
	end do 
!	write(*,*)'flux',j,i1,i2,j1-1,j2+1
! 	if(j==1)then
! 	write(*,'(4F10.3)')q1d(0,1:2),q1d(1,1:2)	
! endif
! 	 if(j==j2/2)then
! 	    write(*,'(9F12.4)')q1d(-1,1:3),q1d(0,1:3),q1d(1,1:3)
! 	 endif
	 
	!not this is dangerous because the following is being set on all patches
	!regardless if the left boundary is the inflow.  In debugging the roelvink test
	!case we discovered that we can reset the scalar_bcond to DIRICHLET and vector_bcond
	!to FlOW_BC and this zeroth order extrapolation we enforce here below happens automatically
	!but only to the correct patch.  At some point we need to modify HENICHE in the same manner
	!(this is in swe::setPhysicalBoundaryConditions) and get rid of this fudge below.
	if(cid==heniche)then !.or. cid==roelvink)then
    !set hv and hu in the ghost cells of an open boundary
! 	!geoff - heniche fudge - gwc FUDGE
   write(*,*)'fudging'
  	q1d(-1,2) = q1d(1,2)
   	q1d( 0,2) = q1d(1,2)
    endif

 !    if(j==j2/2)then
!	    write(*,*)'dirichlet values at 0',q1d
 ! 		do i=i1-NCGST,i2+NCGST
 ! 		  iclaw    = i+s2ci
 !     write(*,*)i,iclaw,q1d(iclaw,1:3),aux(iclaw,1)
 ! end do
 !     stop
 !     end if


    

    !compute modifications fadd and gadd to fluxes along this slice:
    call flux2(1,maxm,meqn,maux,mbc,mx,mwaves,mthlim, &
              dtdx,q1d,aux, &
              faddm,faddp,gaddm,gaddp, &
              fwave,s, &
              amdq,apdq,cqxx,bmasdq,bpasdq, &
              flux_order, transverse_prop)
	 
    !update fluctuations - seems samrai fluxes use odd index ordering - careful
    do m=1,meqn
	do i=i1,i2+1
	  iclaw = i+s2ci
	  ifluxm(i,j,m)   = ifluxm(i,j,m)   + dt*faddm(iclaw,m)
      ifluxp(i,j,m)   = ifluxp(i,j,m)   + dt*faddp(iclaw,m)
      jfluxm(j,i,m)   = jfluxm(j,i,m)   + dt*gaddm(iclaw,m,1)
      jfluxp(j,i,m)   = jfluxp(j,i,m)   + dt*gaddp(iclaw,m,1)
      jfluxm(j+1,i,m) = jfluxm(j+1,i,m) + dt*gaddm(iclaw,m,2)
      jfluxp(j+1,i,m) = jfluxp(j+1,i,m) + dt*gaddp(iclaw,m,2)
    end do
    end do

  end do

!--------------------------------------------------------------------------
!  perform y-sweeps
!--------------------------------------------------------------------------
   
  do i=i1-1,i2+1  !i=0,mx+1
	iclaw = i+s2ci
 
    !copy j-dir pencil at fixed i into 1d array
    do j=j1-NCGST,j2+NCGST
	  jclaw    = j+s2cj
	  q1d(jclaw,1) = h(i,j)     !h
	  q1d(jclaw,2) = vh(i,j,1)  !uh
	  q1d(jclaw,3) = vh(i,j,2)  !vh
	  aux(jclaw,1) = b(i,j)     !b
	end do 
	
! 	if(i==0)then
! 	do j=j1-NCGST,j2+NCGST
! 	  jclaw    = j+s2cj
! 	  write(*,'(I4,4F10.5)')jclaw,q1d(jclaw,1:3),aux(jclaw,1)
! 	end do
! endif
	 
	
	!    if(j==j2/2)then
	!	    write(*,*)'dirichlet values at 0',q1d
	 ! 		do i=i1-NCGST,i2+NCGST
	 ! 		  iclaw    = i+s2ci
	 !     write(*,*)i,iclaw,q1d(iclaw,1:3),aux(iclaw,1)
	 ! end do
	 !     stop
	 !     end if
	

    !compute modifications fadd and gadd to fluxes along this slice:
    call flux2(2,maxm,meqn,maux,mbc,my,mwaves,mthlim, &
              dtdy,q1d,aux, &
              faddm,faddp,gaddm,gaddp, &
              fwave,s, &
              amdq,apdq,cqxx,bmasdq,bpasdq, &
              flux_order, transverse_prop)

	 
    !update fluctuations - seems samrai fluxes use odd index ordering - careful
    do m=1,meqn
	do j=j1,j2+1
	  jclaw = j+s2cj
      jfluxm(j,i,m)   = jfluxm(j,i,m)   + dt*faddm(jclaw,m)
      jfluxp(j,i,m)   = jfluxp(j,i,m)   + dt*faddp(jclaw,m)
      ifluxm(i,j,m)   = ifluxm(i,j,m)   + dt*gaddm(jclaw,m,1)
      ifluxp(i,j,m)   = ifluxp(i,j,m)   + dt*gaddp(jclaw,m,1)
      ifluxm(i+1,j,m) = ifluxm(i+1,j,m) + dt*gaddm(jclaw,m,2)
      ifluxp(i+1,j,m) = ifluxp(i+1,j,m) + dt*gaddp(jclaw,m,2)
    end do
    end do

  end do !outer loop over i

  !deallocate temporary workspace
  deallocate(faddm)
  deallocate(faddp)
  deallocate(gaddm)
  deallocate(gaddp)
  deallocate(q1d)
  deallocate(amdq)
  deallocate(apdq)
  deallocate(bmasdq)
  deallocate(bpasdq)
  deallocate(aux)
  deallocate(cqxx)
  deallocate(s)
  deallocate(fwave)
!  deallocate(mthlim)

  return

end subroutine george_flux

!==============================================================================
! note: dimensioning on jflux (j before i)
!==============================================================================

subroutine consdiff(dx,i1,i2,j1,j2,ifluxm,jfluxm,ifluxp,jfluxp,h,vh,b)
    
  use gparms
  use cntrl
  implicit none
  integer,  intent(in   ) :: i1,i2,j1,j2
  real(dp), intent(in   ) :: dx(2)
  real(dp), intent(in   ) :: ifluxm(i1-NFGST:i2+1+NFGST,j1-NFGST:j2+NFGST,NEQUS)
  real(dp), intent(in   ) :: jfluxm(j1-NFGST:j2+1+NFGST,i1-NFGST:i2+NFGST,NEQUS)
  real(dp), intent(in   ) :: ifluxp(i1-NFGST:i2+1+NFGST,j1-NFGST:j2+NFGST,NEQUS)
  real(dp), intent(in   ) :: jfluxp(j1-NFGST:j2+1+NFGST,i1-NFGST:i2+NFGST,NEQUS)
  real(dp), intent(inout) :: h(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
  real(dp), intent(inout) :: vh(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST,NDIMS)
  !real(dp)  :: shit(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST,NDIMS)
  real(dp), intent(in )   :: b(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
     
  !local
  integer  :: i,j,k
  real(dp) :: oodx,oody
  logical  :: ifound = .false.

  ifound = .false.
  !precompute 1/dx,  1/dy
  oodx = one/dx(1);
  oody = one/dx(2);
  !shit = vh;

  !loop over internal cells, update state variables using convervative diff on fluxes
  do i=i1,i2
	do j=j1,j2
      h(i,j) = h(i,j) &
                                   -oodx*(ifluxm(i+1,j,1)-ifluxp(i,j,1)) &
                                   -oody*(jfluxm(j+1,i,1)-jfluxp(j,i,1)) 
      do k=1,NDIMS
       vh(i,j,k) = vh(i,j,k) &
               -oodx*(ifluxm(i+1,j,k+1)-ifluxp(i,j,k+1)) &
               -oody*(jfluxm(j+1,i,k+1)-jfluxp(j,i,k+1))
      enddo
      !gwc new - 1e-3 seems to stabilize sampson - depth threshold
          if(h(i,j)< mindepth)then
              	h(i,j) = zero
              	vh(i,j,1:2) = zero
              endif
!       if(j==j2/2)then
! 	    write(*,'(I5,3F10.6)')i,h(i,j),vh(i,j,1),vh(i,j,2)
! 	  endif
!        if(i==0)then
!  	    write(*,'(I5,3F10.6)')j,h(i,j),vh(i,j,1),vh(i,j,2),shit(i,j,2)
!         write(*,'(I5,4F12.2)')j,ifluxm(i+1,j,2),ifluxp(i,j,2),jfluxm(j+1,i,2),jfluxp(j,i,2)  
!  	  endif

!       !dry bed adjustment
!       depth = zeta(i,j)-zb(i,j)
!       if(depth <=mindepth)then  
!         zeta(i,j) = zb(i,j)+mindepth
!         veldepth(i,j,1:2) = zero
!       endif
      
    enddo
  end do
 


 
  return
end subroutine consdiff

