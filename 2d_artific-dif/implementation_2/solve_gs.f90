subroutine solve_gs(np,nx,ny,as,aw,ap,ae,an,Su,phi)
!
implicit none
!
! variables passed in
integer, intent(inout) :: np,nx,ny
!real, intent(in) :: Fx,Fy,dx,dy
real, dimension(:), intent(inout) :: as(np),aw(np),ap(np),ae(np),an(np)
real, dimension(:), intent(inout) :: Su(np)
real, dimension(:), intent(inout) :: phi(np)
!
! variables used only in this subroutine
integer ii
!

!----------------------------------------------------------------------!
!   Start at SW corner
!
ii = 1
phi(ii) = (an(ii)*phi(ii+1) + ae(ii)*phi(ii+ny))/ap(ii)
!
end subroutine solve_gs
