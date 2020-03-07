module set
  implicit none


contains
  subroutine setup(rho,nx,x,v,xmin,dx,m,cs,n,h,xmax)
    real, intent(in) :: rho(n),xmin,cs(n),xmax
    ! integer, intent(in) :: nx
    real,intent(out) :: x(n),m(n),v(n),h(n),dx
    integer,intent(inout) :: n
    integer,intent(in) :: nx
    integer :: i
    real,parameter :: pi=4.*atan(1.)
    n=100
    !set initial x
    ! open(1,file='results.dat', status='replace',action='write')
    x(1)=xmin


    v(1)=0


    dx=(xmax-xmin)/(n-1)
    h(:n)=1.2*dx
    !iterate over all i for x to create a grid of positions
    ! write(1,*) x(1),v(1)
    do i=2,n

      x(i)=x(i-1)+dx
      v(i)=cs(i)*10.**(-4)*sin(2*pi*x(i))
      !Smoothing length???
      ! h(i)=1.2*(x(i)-x(i-1))
      ! write(1,*) x(i),v(i)
    enddo

    m(:n)=rho(:n)*(xmax-xmin)/(n)
    ! print*,m,sum(m(1:100))
    ! close(1)



  end subroutine
end module
