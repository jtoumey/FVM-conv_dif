SUBROUTINE CALC_RESIDUAL(np,nx,ny,as,aw,ap,ae,an,Su,phi_prev,resid)
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
integer :: ii,jj
integer :: l_bnd,u_bnd
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
resid = 0.
!----------------------------------------------------------------------!
!
!   Residual contribution from interior cells 
!
!----------------------------------------------------------------------!
do ii = 2,nx-1
   l_bnd = (ii - 1)*ny + 2
   u_bnd = l_bnd + (ny - 3)
   !
   do jj = l_bnd,u_bnd
      resid = resid + abs(aw(jj)*phi_prev(jj-ny) + ae(jj)*phi_prev(jj+ny) + an(jj)*phi_prev(jj+1) + as(jj)*phi_prev(jj-1) + &
      Su(jj) - ap(jj)*phi_prev(jj))
   end do
end do
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
!----------------------------------------------------------------------!
!
!   Residual contribution from the cells adjacent to boundaries,
!   excluding corners 
!
!----------------------------------------------------------------------!
!
! West Boundary, no corners (only count contribution from N, S, and E cells)

do ii = 2,ny-1
   resid = resid + abs(ae(ii)*phi_prev(ii+ny) + an(ii)*phi_prev(ii+1) + as(ii)*phi_prev(ii-1) + Su(ii) - ap(ii)*phi_prev(ii))
end do
!
! South Boundary, no corners (only count contribution from N, E, W cells 

do ii = ny+1,(nx-2)*ny+1,ny
   resid = resid + abs(aw(ii)*phi_prev(ii-ny) + ae(ii)*phi_prev(ii+ny) + an(ii)*phi_prev(ii+1) + Su(ii) - ap(ii)*phi_prev(ii))
end do
!
! East Boundary, no corners (only count contribution from N, S, and W cells)

do ii = (nx-1)*ny+2,np-1
   resid = resid + abs(aw(ii)*phi_prev(ii-ny) + an(ii)*phi_prev(ii+1) + as(ii)*phi_prev(ii-1) + Su(ii) - ap(ii)*phi_prev(ii))
end do
!
! North Boundary, no corners (only count contribution from S, E, W cells 

do ii = 2*ny,np-ny,ny
   resid = resid + abs(aw(ii)*phi_prev(ii-ny) + ae(ii)*phi_prev(ii+ny) + as(ii)*phi_prev(ii-1) + Su(ii) - ap(ii)*phi_prev(ii))
end do
!
END SUBROUTINE CALC_RESIDUAL
