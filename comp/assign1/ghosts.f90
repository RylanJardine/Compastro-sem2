module set
  implicit none

contains
    subroutine setup(rho,nx,x,v,xmin,dx,m,cs,n,h,xmax,ng)

      real, intent(out) :: xmax,xmin
      integer,intent(inout) :: n,ng
      integer,intent(in) :: nx
      ! integer, intent(in) :: nx
      real,intent(out) :: x(:),m(:),v(:),h(:),dx,rho(:),cs(:)
      ! real, intent (inout) :: a(n)      real,intent


      integer :: i
      real,parameter :: pi=4.*atan(1.)
      n=100

      xmin=0.
      xmax=1.
      cs=1.
      rho=1.
      ! print*,rho(:), 5
      x(1)=xmin


      v(1)=0.


      dx=(xmax-xmin)/(n-1)

      h=1.2*dx
      do i=2,n

        x(i)=x(i-1)+dx
        v(i)=cs(i)*10.**(-4)*sin(2*pi*x(i))

      enddo

      m=rho*dx

      ! call set_ghosts(rho,nx,x,v,dx,m,cs,n,h)
      call set_ghosts(n,ng,x,v,rho,m,h)


    end subroutine


  subroutine set_ghosts(n,ng,x,v,rho,m,h)

    ! call set_ghosts(rho,nx,x,v,dx,m,cs,n,h)

    integer,intent(in) :: n
    integer,intent(inout) :: ng
    real, intent (inout) ::x(:),v(:),rho(:),m(:),h(:)
    integer :: i,j
    real, parameter :: radkern=2.,xmax=1.,xmin=0.
    real :: dx


    dx=xmax-xmin
    if (abs(xmax-xmin) < tiny(0.)) print*,'Error xmax-xmin is 0'

    ! ng=0

    do i=1,n
      if (x(i)+radkern*h(i)>xmax) then
        ng=ng+1
        j=n+ng
        x(j)=x(i)-dx
        v(j)=v(i)
        m(j)=m(i)
      else if (x(i)-radkern*h(i)<xmin) then
        ng=ng+1
        j=n+ng
        x(j)=x(i)+dx
        v(j)=v(i)
        m(j)=m(i)
      endif
    enddo
    print*,ng

  end subroutine

end module
