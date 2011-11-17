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


subroutine detectgrad(i1,i2,j1,j2,ngi,ngj,ngtagi,ngtagj,ngtti,ngttj,dx,gradtol, &
                      yestag,notag,var,tags,temptags)  
  use gparms
  implicit none
  
  integer,  intent(in ) :: i1,i2,j1,j2,ngi,ngj,ngtagi,ngtagj,ngtti,ngttj
  real(dp), intent(in ) :: dx(2)
  real(dp), intent(in ) :: gradtol
  integer,  intent(in ) :: yestag,notag
  real(dp), intent(in ) :: var(i1-ngi:i2+ngi,j1-ngj:j2+ngj)
  integer,  intent(out) :: tags(i1-ngtagi:i2+ngtagi,j1-ngtagj:j2+ngtagj)
  integer,  intent(out) :: temptags(i1-ngtti:i2+ngtti,j1-ngttj:j2+ngttj)
 
  !----- local ------------
  integer  :: i,j
  real(dp) :: tol,diag01,loctol,facejump,varm1,varp1
  logical  :: tagcell
 
  tol = gradtol
  diag01 = sqrt(dx(1)**2+dx(2)**2)
  
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
  !stop
  return

end subroutine detectgrad

  