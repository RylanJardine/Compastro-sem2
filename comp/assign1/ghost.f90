module ghost
  implicit none

contains
  subroutine set_ghosts(rho,nx,x,v,xmin,dx,m,cs,n,h,xmax)
    real, intent(in) :: xmin,xmax,dx
    real,intent(inout) :: rho(nx),cs(nx)
    ! integer, intent(in) :: nx
    real,intent(out) :: x(nx),m(nx),v(nx),h(nx)
    integer,intent(in) :: n,nx
    integer :: i,ng
    real,parameter :: pi=4.*atan(1.)


    ! x(n+1)=x(1)

    ng=(nx-n)/2
    ! v(n+1)=v(1)

    h(n:)=1.2*dx
    cs(n:)=1.
    !iterate over all i for x to create a grid of positions
    ! write(1,*) x(1),v(1)
    do i=n+1,n+ng

      x(i)=x(i-1)+dx
      v(i)=cs(i)*10.**(-4)*sin(2*pi*x(i))
      !Smoothing length???
      ! h(i)=1.2*(x(i)-x(i-1))
      ! write(1,*) x(i),v(i)
    enddo

    do i=n+ng,nx
      x(i)=

    print*,'pre-rho',rho
    do i=1,ng
      rho(n+i)=rho(i)
    enddo
    print*,'post-rho',rho
    m(n:)=rho(n:)*dx



  end subroutine
end module
