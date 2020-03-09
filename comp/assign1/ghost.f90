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
      v(i)=v(i-n)
      rho(i)=rho(i-n)
      m(i)=m(i-n)
      !Smoothing length???
      ! h(i)=1.2*(x(i)-x(i-1))
      ! write(1,*) x(i),v(i)
      ! rho(n+i)=rho(i)
    enddo

    do i=1,ng
      x(i+n+ng)=-dx*i
      v(i+n+ng)=v(n-i)
      rho(i+n+ng)=rho(n-i)
      m(i+n+ng)=m(n-i)
      ! rho(i+n+ng)=rho(0-i)
    enddo

    ! print*,'pre-rho',v
    ! do i=1,nx
    !   rho(n+i)=rho(i)
    ! enddo
    print*,'post-rho',v
    ! m(n:ng)=rho(n:ng)*dx
    ! m(ng:nx)=rho(ng:nx)*dx



  end subroutine
end module
