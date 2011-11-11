!==============================================================================
! 
!                          dump a junk probe file
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
!
! 
!==============================================================================

subroutine junkprobe(cid,dx,xlo,xhi,time,i1,i2,j1,j2,h,vh,b)
    
  use gparms
  use cntrl
  implicit none
  integer,  intent(in) :: cid
  integer,  intent(in) :: i1,i2,j1,j2
  real(dp), intent(in) :: dx(2),xlo(2),xhi(2)
  real(dp), intent(in) :: time
  real(dp), intent(in) :: h(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
  real(dp), intent(in) :: vh(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST,NDIMS)
  real(dp), intent(in) :: b(i1-NCGST:i2+NCGST,j1-NCGST:j2+NCGST)
 
  !local
  integer  :: i,j,iloc,jloc
  real(dp) :: probe_x,probe_y,rmin,rad,xc,yc,u

  caseid = cid

  select case(caseid)
 
  !---------------------------------------------------------------------
  case(conrun) 
  !---------------------------------------------------------------------
  probe_x = -2.6
    probe_y = 0
  
 
  rmin = 1e9
  do i=i1-1,i2+1
    do j=j1-1,j2+1 
	  xc = xlo(1)+dx(1)*dble(i-i1)+dx(1)/2
	  yc = xlo(2)+dx(2)*dble(j-j1)+dx(2)/2
	  rad = sqrt((xc-probe_x)**2 + (yc-probe_y)**2)
	  if(rad < rmin)then
		iloc = i
		jloc = j
		rmin = rad
	  endif
 	end do 
  end do
  if(time<1e-8)then  !start a new file
	open(unit=33,file='probe.dat',position='rewind')
  else
    open(unit=33,file='probe.dat',position='append')
  endif
  !write(33,*)time,h(-1,jloc)+b(-1,jloc),h(0,jloc)+b(0,jloc),b(-1,jloc),vh(-1,jloc,1)/h(-1,jloc)
  write(33,*)time,b(iloc,jloc)+h(iloc,jloc)
 
  close(33)

  !---------------------------------------------------------------------
  case(heniche) 
  !---------------------------------------------------------------------
 !  probe_x = -2.6
 !   probe_y = 0
  probe_x = 430
  probe_y = 250

 !set to first order on return- GWC DEBUG FUDGE
 if(time > 2000) flux_order = 1

  rmin = 1e9
  do i=i1-1,i2+1
    do j=j1-1,j2+1 
	  xc = xlo(1)+dx(1)*dble(i-i1)+dx(1)/2
	  yc = xlo(2)+dx(2)*dble(j-j1)+dx(2)/2
	  rad = sqrt((xc-probe_x)**2 + (yc-probe_y)**2)
	  if(rad < rmin)then
		iloc = i
		jloc = j
		rmin = rad
	  endif
 	end do 
  end do
  if(time<1e-8)then  !start a new file
	open(unit=33,file='probe.dat',position='rewind')
  else
    open(unit=33,file='probe.dat',position='append')
  endif
  
  if(h(iloc,jloc) > 0)then
    u = vh(iloc,jloc,1)/h(iloc,jloc)
  else
	u = zero
  endif
  write(33,*)time,(gravity*C_manning*C_manning/(h(iloc,jloc)**(a3rd)))*(u**2),u,h(iloc,jloc)
  
  close(33)
  
  end select

  return
end subroutine junkprobe

