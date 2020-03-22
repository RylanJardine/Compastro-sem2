module set
  implicit none


contains
  subroutine setup(rho,nx,x,v,xmin,dx,m,cs,n,h,xmax)

    real, intent(out) :: xmax,xmin
    integer,intent(inout) :: n
    integer,intent(in) :: nx
    ! integer, intent(in) :: nx
    real,intent(out) :: x(:),m(:),v(:),h(:),dx,rho(:),cs(:)
    ! real, intent (inout) :: a(n)


    integer :: i
    real,parameter :: pi=4.*atan(1.)
    n=100

    xmin=0.
    xmax=1.
    cs(1:n)=1.
    rho(1:n)=1.
    ! print*,rho(:), 5
    x(1)=xmin


    v(1)=0.


    dx=(xmax-xmin)/(n-1)

    h(1:n)=1.2*dx
    do i=2,n

      x(i)=x(i-1)+dx
      v(i)=cs(i)*10.**(-4)*sin(2*pi*x(i))

    enddo

    m(1:n)=rho(1:n)*dx

    call set_ghosts(rho,nx,x,v,dx,m,cs,n,h)

  end subroutine


  subroutine set_ghosts(rho,nx,x,v,dx,m,cs,n,h)
    real, intent(in) :: dx
    integer,intent(in) :: n,nx
    real,intent(inout) :: rho(nx),cs(nx)
    real,intent(out) :: x(nx),m(nx),v(nx),h(nx)

    integer :: i,ng
    ! real,parameter :: pi=4.*atan(1.)

    ng=(nx-n)/2
    h(n:)=1.2*dx
    cs(n:)=1.


    ! the ghost points which lead immediately after the particles e.g. 101, 102 ...
    do i=n+1,n+ng
      x(i)=x(i-1)+dx
      v(i)=v(i-n+1)
      rho(i)=rho(i-n+1)
      m(i)=m(i-n+1)
      ! cs(i)=cs(i-n+1)
      ! p(i)=p(i-n+1)
      ! a(i)=a(i-n+1)
    enddo


    ! the ghost particles which exist before the first particles e.g. -1, -2...
    do i=1,ng
      ! x(i+n+ng)=-dx*i
      x(i+n+ng)=x(1)-dx*i
      v(i+n+ng)=v(n-i)
      rho(i+n+ng)=rho(n-i)
      m(i+n+ng)=m(n-i)
      ! cs(i+n+ng)=cs(n-i)
      ! p(i+n+ng)=p(n-i)
      ! a(i+n+ng)=a(n-i)
    enddo

  end subroutine

end module
