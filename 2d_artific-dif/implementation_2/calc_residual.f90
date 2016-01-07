SUBROUTINE CALC_RESIDUAL(np,nx,ny,as,aw,ap,ae,an,phi_prev)
!
implicit none
!
! variables passed in
integer :: nx,ny,np
real resid
real, dimension(np) :: as,aw,ap,ae,an 
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



END SUBROUTINE CALC_RESIDUAL
