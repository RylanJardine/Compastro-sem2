module density
  use kernel
  implicit none



contains
  subroutine get_density(m,x,rho,nx,n,h)
    real,intent(in) :: m(nx), x(nx),h(nx)
    real,intent(out) :: rho(nx)
    real :: rhop(nx),w(nx),hin,xin
    integer :: i,j
    integer,intent(in) :: nx,n

    do i=1,n
      call kern(h(i),nx,x,x(i),w,n)
      do j=1,n

        rhop(j)=m(j)*w(j)
        ! rho(i)=rho(i)+rho(j)
        ! print*,rhop(j)
        if (rhop(j) .NE. 0) then
          print*,rhop(j),j
        endif

      enddo
      ! print*,rhop
      rho(i)=sum(rhop)

    !   if (i==5) then
    !   !!!!!!!!!!
    !     stop
    !   !!!!!!!!
    ! endif
    enddo

  end subroutine
end module
