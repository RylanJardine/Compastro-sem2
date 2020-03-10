module ghost
  implicit none

contains
  subroutine set_ghosts(rho,nx,x,v,xmin,dx,m,cs,n,h,xmax)
    real, intent(in) :: xmin,xmax,dx
    real,intent(inout) :: rho(nx),cs(nx)
    real,intent(out) :: x(nx),m(nx),v(nx),h(nx)
    integer,intent(in) :: n,nx
    integer :: i,ng
    real,parameter :: pi=4.*atan(1.)

    ng=(nx-n)/2
    h(n:)=1.2*dx
    cs(n:)=1.

    do i=n+1,n+ng
      x(i)=x(i-1)+dx
      v(i)=v(i-n)
      rho(i)=rho(i-n)
      m(i)=m(i-n)
    enddo

    do i=1,ng
      x(i+n+ng)=-dx*i
      v(i+n+ng)=v(n-i)
      rho(i+n+ng)=rho(n-i)
      m(i+n+ng)=m(n-i)
    enddo

  end subroutine
end module
