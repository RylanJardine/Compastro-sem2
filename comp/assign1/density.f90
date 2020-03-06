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

    do i=1,nx
      do j=1,nx
        call kern(h(i),nx,x,x(i),w)
        rhop(j+1)=rhop(j)+m(j)*w(j)
        rho(i)=rho(i)+rho(j)


      enddo
      rho(i)=sum(rhop)
    enddo

  end subroutine
end module
