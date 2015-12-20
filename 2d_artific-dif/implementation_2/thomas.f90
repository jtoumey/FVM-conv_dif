!*************************************************************************!
!                                                                         !
!  Module:       THOMAS.F90                                               !
!                                                                         !
!  Programmer:   Julian M. Toumey                                         !
!                Madison, WI                                              !
!                                                                         !
!  Date:         26 January 2015                                          !
!                                                                         !
!  Language:     FORTRAN90                                                !
!                                                                         !
!  Description:  Solves a tri-diagonal system of linear equations using   !
!                the Thomas Algorithm (TDMA). The algorithm in this       !
!                function comes from Versteeg and Malalasekera, 2nd ed.,  !
!                Section 7.2                                              !
!                                                                         !
!                  a     Sub diagonal                                     !
!                  b     Main diagonal                                    !
!                  c     Super diagonal                                   !
!                  d     Right side of Linear System                      !
!                  phi   Solution vector                                  !
!                                                                         !
!*************************************************************************!
SUBROUTINE THOMAS(n,a,b,c,d,phi)
integer i,n
real a(n),b(n),c(n),d(n),phi(n)
!
!...Forward sweep
!
do i = 2,n
   b(i) = b(i) - (a(i)*c(i-1))/b(i-1) 
   d(i) = d(i) - (a(i)*d(i-1))/b(i-1)
end do
!
!...Back substitution
!
phi(n) = d(n)/b(n)
!
do i = n-1,1,-1
   phi(i) = (d(i) - c(i)*phi(i+1))/b(i)
end do
end SUBROUTINE THOMAS