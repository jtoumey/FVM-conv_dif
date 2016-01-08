SUBROUTINE CALC_RESIDUAL(np,nx,ny,as,aw,ap,ae,an,Su,phi_prev)
!
implicit none
!
! variables passed in
integer :: nx,ny,np
real resid
real, dimension(np) :: as,aw,ap,ae,an 
real, dimension(np) :: Su
real, dimension(np) :: phi_prev
!
! variables used only in this subroutine
integer ii,jj
!
!----------------------------------------------------------------------!
!
!   We are solving Ax = b; \sum{a_nb \phi_nb} = Su; 
!      * A is equivalent to a_nb
!      * x is equivalent to \phi_nb
!      * b is equivalent to Su
!     (* Note that a_P = \sum{a_nb} - Sp)
!
!   Calculate the residual by substituting \phi_prev (known from the
!   previous iteration) into the discretized equation and comparing this 
!   solution to Su. 
!
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
!
!   Residual contribution from interior cells 
!
!----------------------------------------------------------------------!
!!do jj = 1,34
   !write(*,*)'p'
!end do
!----------------------------------------------------------------------!
!
!   Residual contribution from the four corners 
!
!----------------------------------------------------------------------!
! SW Corner (only count contribution from N and E cells)

ii = 1 ! keep track of index 
resid = resid + abs(ae(ii)*phi_prev(ii+ny) + an(ii)*phi_prev(ii+1) + Su(ii) - ap(ii)*phi_prev(ii))
!
! NW Corner (only count contribution from S and E cells)

ii = ny ! keep track of index
resid = resid + abs(ae(ii)*phi_prev(ii+ny) + as(ii)*phi_prev(ii-1) + Su(ii) - ap(ii)*phi_prev(ii))
!
! SE Corner (only count contribution from N and W cells)

ii = (nx-1)*ny + 1 ! keep track of index
resid = resid + abs(aw(ii)*phi_prev(ii-ny) + an(ii)*phi_prev(ii+1) + Su(ii) - ap(ii)*phi_prev(ii))
!
! NE Corner (only count contribution from S and W cells)

ii = np ! keep track of index
resid = resid + abs(aw(ii)*phi_prev(ii-ny) + as(ii)*phi_prev(ii-1) + Su(ii) - ap(ii)*phi_prev(ii))
!
END SUBROUTINE CALC_RESIDUAL
