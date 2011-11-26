!==============================================================================
! 
!                                  detectgrad
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
! Comments:     Calculate Stable Time Step for SWE equations
!
!       Input:
!          var:  flow variable used for tagging cells
!          
!       Output:
!          temptags: cell markers tagging cell for refinement 
!==============================================================================


subroutine detectgrad(time,i1,i2,j1,j2,ngi,ngj,ngtagi,ngtagj,ngtti,ngttj,dx,xlo,xhi,& 
                      gradtol,yestag,notag,var,tags,temptags)  
  use gparms
  use cntrl
  implicit none
  
  real(dp), intent(in ) :: time
  integer,  intent(in ) :: i1,i2,j1,j2,ngi,ngj,ngtagi,ngtagj,ngtti,ngttj
  real(dp), intent(in ) :: dx(2),xlo(2),xhi(2)
  real(dp), intent(in ) :: gradtol
  integer,  intent(in ) :: yestag,notag
  real(dp), intent(in ) :: var(i1-ngi:i2+ngi,j1-ngj:j2+ngj)
  integer,  intent(out) :: tags(i1-ngtagi:i2+ngtagi,j1-ngtagj:j2+ngtagj)
  integer,  intent(out) :: temptags(i1-ngtti:i2+ngtti,j1-ngttj:j2+ngttj)
 
  !----- local ------------
  integer  :: i,j
  real(dp) :: tol,diag01,loctol,facejump,varm1,varp1,xc,yc
  logical  :: tagcell
  tol = gradtol
  diag01 = sqrt(dx(1)**2+dx(2)**2)
  ! do j=j1,j2
  ! 	 write(*,*)var(-1,j),var(0,j),var(1,j)
  !   end do
  do i=i1,i2
	do j=j1,j2
      
      
      if (tags(i,j) .ne. 0) then
        loctol = 0.125*tol
      else 
        loctol = tol
      endif
   
      tagcell = .false.
      

      !i-direction jump
      if (.not.tagcell) then
        varm1 = var(i-1,j)
        varp1 = var(i+1,j)
        facejump = abs(var(i,j)-varm1)
        facejump = max(facejump,abs(var(i,j)-varp1))
        tagcell = ((facejump).gt.(loctol*dx(1)))
        !write(*,*)'doing i face',i,j,tagcell,facejump,varm1,varp1,var(i,j)
      endif
      !j-direction jump
      if (.not.tagcell) then
        varm1 = var(i,j-1)
        varp1 = var(i,j+1)
        facejump = abs(var(i,j)-varm1)
        facejump = max(facejump,abs(var(i,j)-varp1))
        tagcell = ((facejump).gt.(loctol*dx(2)))
      endif
      !diagonal 1 jump
      if (.not.tagcell) then
        varm1 = var(i-1,j-1)
        varp1 = var(i+1,j+1)
        facejump = abs(var(i,j)-varm1)
        facejump = max(facejump,abs(var(i,j)-varp1))
        tagcell = ((facejump).gt.(loctol*diag01))
      endif
      !diagonal 2 jump
      if (.not.tagcell) then
        varm1 = var(i-1,j+1)
        varp1 = var(i+1,j-1)
        facejump = abs(var(i,j)-varm1)
        facejump = max(facejump,abs(var(i,j)-varp1))
        tagcell = ((facejump).gt.(loctol*diag01))
      endif

      if ( tagcell ) then
        temptags(i,j) = yestag
      endif
      !write(*,*)i,j,temptags(i,j)
      
 
    enddo
  enddo

  if(caseid==roelvink)then
	do i=i1,i2
		xc = xlo(1) + dx(1)*dble(i-i1)+dx(1)/2
		
		do j=j1,j2
	      yc = xlo(2) + dx(2)*dble(j-j1)+dx(2)/2
		  temptags(i,j) = notag
		  if(xc >= 2500 .and. xc <=8000 .and. yc >5000 .and. yc < 10000) temptags(i,j) = yestag
		end do
	end do
  endif
! if(caseid==roelvink)then
! 	do i=i1,i2
! 		xc = xlo(1) + dx(1)*dble(i-i1)+dx(1)/2		
! 		do j=j1,j2
! 	      yc = xlo(2) + dx(2)*dble(j-j1)+dx(2)/2
! 		  if(xc <= 2500 .and. xc >=7500 .and. yc <6000 .and. yc > 10000) temptags(i,j) = notag
! 		end do
! 	end do
! endif

! below is useful if we want to just test AMR.  Advance adaption front x(t), set ICS for slosh_inlet to h=1,b=0,u=0
! if(caseid==slosh_inlet)then
! 	do i=i1,i2
! 		xc = xlo(1) + dx(1)*dble(i-i1)+dx(1)/2		
! 		do j=j1,j2
! 			temptags(i,j) = notag
! 			if(xc < time) temptags(i,j) = yestag
! 	   end do
! 	end do
! endif
  
  return

end subroutine detectgrad

  