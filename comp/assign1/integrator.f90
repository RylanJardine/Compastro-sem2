module integrator
  use density
  implicit none

contains
  subroutine tim(x,v,a,nx,dt,cs,rho,p,n,m,h,dx,ek)
    integer,intent(in) :: nx,n
    real, intent(inout) :: x(nx),v(nx),a(nx),h(nx),rho(nx),m(nx),p(nx),cs(nx),dx,ek
    real :: vs(nx),dt,a0(nx)
    ! dt=0.001
    dt=0.2*h(1)/cs(1)

    a0=a
    ! x(:)=x(:)+dt*v(:)+0.5*(dt)**2*a0(:)
    !
    ! vs(:)=v(:)+dt*a(:)
    ! call derivs(cs,rho,p,n,a,nx,x,m,h,dx,vs)
    !
    ! v(:)=vs(:)+0.5*dt*(a(:)-a0(:))

    ! a=a0


    x(1:n)=x(1:n)+dt*v(1:n)+0.5*(dt)**2*a0(1:n)

    vs(1:n)=v(1:n)+dt*a(1:n)
    call derivs(cs,rho,p,n,a,nx,x,m,h,dx,v)

    v(1:n)=vs(1:n)+0.5*dt*(a(1:n)-a0(1:n))

    ! print*,p

    ek=sum(0.5*m(1:n)*v(1:n)**2)

  end subroutine
end module
