module set
  implicit none


contains
  subroutine setup(rho,nx,x,v,xmin,dx,m,rho0,cs,n,h,xmax)
    real, intent(in) :: rho(nx),xmin,rho0,cs(nx),xmax
    ! integer, intent(in) :: nx
    real,intent(out) :: x(nx),m(nx),v(nx),h(nx),dx
    integer,intent(out) :: n
    integer,intent(in) :: nx
    integer :: i
    real,parameter :: pi=4.*atan(1.)
    n=100
    !set initial x
    open(1,file='results.dat', status='replace',action='write')
    x(1)=xmin


    v(1)=0

    h(1)=1.2*xmin
    dx=(xmax-xmin)/(n)
    !iterate over all i for x to create a grid of positions
    write(1,*) x(1),v(1)
    do i=2,nx

      x(i)=x(i-1)+dx
      ! v(i)=cs(i)*10.**(-4)*sin(x(i))
      v(i)=sin(2*pi*x(i))
      h(i)=1.2*x(i)
      write(1,*) x(i),v(i)
      print*,'x(i)',x(i),'v(i)',v(i)
    enddo

    m(:)=rho(1)*x(:)
    close(1)



  end subroutine
end module
