module set2
  implicit none


contains
  subroutine setup_shock(rho,nx,x,v,xmin,dx,m,cs,n,h,xmax)

        real, intent(out) :: xmax,xmin
        integer,intent(inout) :: n
        integer,intent(in) :: nx
        ! integer, intent(in) :: nx
        real,intent(out) :: x(:),m(:),v(:),h(:),dx,rho(:),cs(:)
        ! real, intent (inout) :: a(n)


        integer :: i
        real,parameter :: pi=4.*atan(1.)
        ! n=100

        xmin=-0.5
        xmax=0.5
        cs=1.
        rho=1.
        ! print*,rho(:), 5
        x(1)=xmin


        do i=1,n
          if (x(i)<0) then
            rho(i)=1.
          else if (x(i)>=0) then
            rho(i)=0.1
          endif
        enddo




        dx=(xmax-xmin)/(n-1)

        h=1.2*dx
        do i=2,n

          x(i)=x(i-1)+dx
          v(i)=cs(i)*10.**(-4)*sin(2*pi*x(i))

        enddo

        m=rho*dx

        ! call set_ghosts(rho,nx,x,v,dx,m,cs,n,h)

      end subroutine



end module
