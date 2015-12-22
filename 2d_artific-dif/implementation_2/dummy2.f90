subroutine dummy2(np,ap)
integer ii
integer, intent(inout) :: np
real, dimension(:) :: ap(np)


do ii = 1,np
   ap(ii) = ii
end do
!
end subroutine dummy2
