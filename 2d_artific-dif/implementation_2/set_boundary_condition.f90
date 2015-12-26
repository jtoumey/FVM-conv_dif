SUBROUTINE SET_BOUNDARY_CONDITION(np,nx,ny,Fx,Fy,dx,dy,ap,an,as,aw,ae,Su,Sp)
!
IMPLICIT NONE
!
!   variables passed in
integer, intent(inout) :: np,nx,ny
real, intent(in) :: Fx,Fy,dx,dy
real, dimension(:), intent(inout) :: an(np),as(np),aw(np),ae(np)
real, dimension(:), intent(inout) :: Su(np),Sp(np)
real, dimension(:), intent(inout) :: ap(np)
!
!   variables used only in this subroutine
integer ii
!

!
!   West boundary, running S to N
do ii = 1,ny
   aw(ii) = 0.
   Su(ii) = Su(ii) + 100.*Fx*dy 

end do
!   North boundary, running W to E 
do ii = ny,np,ny
   an(ii) = 0.
   Su(ii) = Su(ii) + 100.*Fy*dx
end do
!   South boundary, running W to E 
do ii = 1,np,ny
   as(ii) = 0.
   Su(ii) = Su(ii) + 0.*Fy*dx
end do
!   East boundary, running S to N
do ii = np,np-ny,-1
   ae(ii) = 0.
   Su(ii) = Su(ii) + 0.*Fx*dy 

end do
!do ii = 1,np
!   ap(ii) = ae(ii) + aw(ii) + an(ii) + as(ii) - Sp(ii)
!end do
!
END SUBROUTINE SET_BOUNDARY_CONDITION
