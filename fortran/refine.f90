!==============================================================================
! 
!                                 refinement
!                conservative linear interpolation of h,uh,vh
!             preserved quiescent flow for refinement near shoreline
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
! Comments:     Refine conserved variables on  patch
!
!       h:      zf is conservatively refined using limited gradients, hf => max(zf-bf,0)
!      hu:      hu is conservatively refined using limited gradients 
!      hv:      hv is conservatively refined using limited gradients
!                   
!       Output:
!          =>  hf,huf,hvf: conserved variables on fine mesh 
!==============================================================================

subroutine refine(ic1,ic2,jc1,jc2,if1,if2,jf1,jf2, &
                  ratio,dxc,dxf,  &
                  hc,vhc,bc, &
                  hf,vhf,bf)
  
   use gparms
   use cntrl
   implicit none
   integer,  intent(in)  :: ic1,ic2,jc1,jc2,if1,if2,jf1,jf2
   integer,  intent(in)  :: ratio(2)
   real(dp), intent(in)  :: dxc(2),dxf(2)
   real(dp), intent(in ) :: hc(ic1-NCGST:ic2+NCGST,jc1-NCGST:jc2+NCGST)
   real(dp), intent(in ) :: vhc(ic1-NCGST:ic2+NCGST,jc1-NCGST:jc2+NCGST,NDIMS)
   real(dp), intent(in ) :: bc(ic1-NCGST:ic2+NCGST,jc1-NCGST:jc2+NCGST)
   real(dp), intent(in ) :: hf(if1-NCGST:if2+NCGST,jf1-NCGST:jf2+NCGST)
   real(dp), intent(out) :: vhf(if1-NCGST:if2+NCGST,jf1-NCGST:jf2+NCGST,NDIMS)
   real(dp), intent(out) :: bf(if1-NCGST:if2+NCGST,jf1-NCGST:jf2+NCGST)

   !local
   real(dp), allocatable :: idiff(:)
   real(dp), allocatable :: jdiff(:)
   real(dp), allocatable :: islope(:,:)
   real(dp), allocatable :: jslope(:,:)
   real(dp), allocatable :: deltax(:)
   real(dp), allocatable :: deltay(:)
   integer,  allocatable :: if2c(:)
   integer,  allocatable :: jf2c(:)
   integer  :: ic,jc,i,j,k,ir,jr
   real(dp) :: bound,coef2

!-----------------------------------------------------------------------
!  allocate local data
!-----------------------------------------------------------------------
   allocate(idiff(ic1:ic2+1))
   allocate(jdiff(jc1:jc2+1))
   allocate(islope(ic1:ic2,jc1:jc2))
   allocate(jslope(ic1:ic2,jc1:jc2))
   allocate(deltax(0:ratio(1)-1))
   allocate(deltay(0:ratio(2)-1))
  
 
!-----------------------------------------------------------------------
! set deltax for reconstruction
!-----------------------------------------------------------------------
    
   do i=0,ratio(1)-1
     deltax(i)=(dble(i)+ahalf)*dxf(1)-dxc(1)*ahalf
   enddo

   do j=0,ratio(2)-1
     deltay(j)=(dble(j)+ahalf)*dxf(1)-dxc(1)*ahalf
   enddo

!-----------------------------------------------------------------------
! build fine->coarse index mapping
!-----------------------------------------------------------------------
   allocate(if2c(if1:if2))
   allocate(jf2c(jf1:jf2))

   do i=if1,if2
	 if(i < 0)then
	   if2c(j) = (i+1)/ratio(1)-1
	 else
	   if2c(i) = i/ratio(1)
	 endif
   end do

   do j=jf1,jf2
	 if(j < 0)then
	   jf2c(j) = (j+1)/ratio(2)-1
	 else  
	   jf2c(j) = j/ratio(2)
	 endif
   end do

!-----------------------------------------------------------------------
! set eta in each coarse cell
!  => b_c + h_c if the coarse cell is wet
!  => set to the value of surrounding wet cells if cell is dry
!-----------------------------------------------------------------------
   do jc=jc1-1,jc2+1
     do ic=ic1-1,ic2+1
   	   if(hc(ic,jc)>drytol)then
	     zc(ic,jc) = bc(ic,jc)+hc(ic,jc)
	   else
		
	   endif
	enddo
  end do

!-----------------------------------------------------------------------
! conservative slope-limited refinement of zeta_c into zeta_f
!-----------------------------------------------------------------------
  !slopes in i
  do jc=jc1,jc2
	do ic=ic1,ic2+1
	  idiff(ic) = zetac(ic,jc,k)-zetac(ic-1,jc,k)
	end do

	do ic=ic1,ic2
       coef2=ahalf*(idiff(ic+1)+idiff(ic))
       bound=2.0*min(abs(idiff(ic+1)),abs(idiff(ic)))
       if (idiff(ic)*idiff(ic+1).gt.zero) then
         islope(ic,jc)=sign(min(abs(coef2),bound),coef2)/dxc(1)                
       else
         islope(ic,jc)=zero
       endif
    end do

  end do  

  !slopes in j
  do ic=ic1,ic2
	do jc=jc1,jc2+1
	  jdiff(jc) = vhc(ic,jc,k)-vhc(ic,jc-1,k)
	end do

	do jc=jc1,jc2
       coef2=ahalf*(jdiff(jc+1)+jdiff(jc))
       bound=2.0*min(abs(jdiff(jc+1)),abs(jdiff(jc)))
       if (jdiff(jc)*jdiff(jc+1).gt.zero) then
         jslope(ic,jc)=sign(min(abs(coef2),bound),coef2)/dxc(2)                
       else
         jslope(ic,jc)=zero
       endif
    end do

  end do

  !interpolate from coarse => fine
  do j=jf1,jf2	
	jc = jf2c(j)
    jr=j-jc*ratio(2)

	do i=if1,if2
	  ic = if2c(i)
	  ir = i-ic*ratio(1)
	  zetaf(i,j,k) = zetac(i,j,k) + islope(ic,jc)*deltax(ir) +  jslope(ic,jc)*deltay(jr)
    enddo
  enddo

!-----------------------------------------------------------------------
! compute depth in fine mesh cells
! gwc => should this be max with drytol?
!-----------------------------------------------------------------------
  do j=jf1,jf2
	do i=if1,if2
	  hf(i,j) = max(zetaf(i,j)-b(i,j),zero)
	end do
  end do

!-----------------------------------------------------------------------
! precompute max velocities in each coarse cell
! used later to adjust fine cell momentum 
!-----------------------------------------------------------------------
   do jc=jc1,jc2
	 do ic=ic1,ic2
	   do k=-1,1
		 do l=-1,1
		  
	  
!-----------------------------------------------------------------------
! conservative slope-limited refinement of (x,y)-momentum (hu,hv)
!-----------------------------------------------------------------------
  do k=1,2

  !slopes in i
  do jc=jc1,jc2
	do ic=ic1,ic2+1
	  idiff(ic) = vhc(ic,jc,k)-vhc(ic-1,jc,k)
	end do
	
	do ic=ic1,ic2
       coef2=ahalf*(idiff(ic+1)+idiff(ic))
       bound=2.0*min(abs(idiff(ic+1)),abs(idiff(ic)))
       if (idiff(ic)*idiff(ic+1).gt.zero) then
         islope(ic,jc)=sign(min(abs(coef2),bound),coef2)/dxc(1)                
       else
         islope(ic,jc)=zero
       endif
    end do

  end do  
  
  !slopes in j
  do ic=ic1,ic2
	do jc=jc1,jc2+1
	  jdiff(jc) = vhc(ic,jc,k)-vhc(ic,jc-1,k)
	end do
	
	do jc=jc1,jc2
       coef2=ahalf*(jdiff(jc+1)+jdiff(jc))
       bound=2.0*min(abs(jdiff(jc+1)),abs(jdiff(jc)))
       if (jdiff(jc)*jdiff(jc+1).gt.zero) then
         jslope(ic,jc)=sign(min(abs(coef2),bound),coef2)/dxc(2)                
       else
         jslope(ic,jc)=zero
       endif
    end do

  end do

  !interpolate from coarse => fine
  do j=jf1,jf2	
	jc = jf2c(j)
    jr=j-jc*ratio(2)
	
	do i=if1,if2
	  ic = if2c(i)
	  ir = i-ic*ratio(1)
	  vhf(i,j,k) = vhc(i,j,k) + islope(ic,jc)*deltax(ir) +  jslope(ic,jc)*deltay(jr)
    enddo
  enddo

  end do !loop over 2 directions

!-----------------------------------------------------------------------
!  check velocity in fine cells and ensure it does not exceed the max
!  value set in the coarse parent cell
!  make sure velocities are zero in dry cells
!-----------------------------------------------------------------------
  do j=1jf1,jf2
	jc = jf2c(j)
	do i=1if1,if2
	  ic = if2c(i)
	  if(hf(i,k) > drytol)then
		uf = vhf(i,j,1)/hf(i,k)
		if(uf > umax(ic,jc)) vhf(i,j,1) = 
		
!-----------------------------------------------------------------------
!  set velocities in dry cells to zero
!  value set in the coarse parent cell
!-----------------------------------------------------------------------		

!-----------------------------------------------------------------------
!  deallocate local data
!-----------------------------------------------------------------------
   deallocate(idiff)
   deallocate(jdiff)
   deallocate(islope)
   deallocate(jslope)
   deallocate(deltax)
   deallocate(deltay)
   deallocate(if2c)
   deallocate(jf2c)

end subroutine refine
!==============================================================================


