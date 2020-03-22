module integrator
  use density
  implicit none

contains
  subroutine tim(x,v,a,nx,dt,cs,rho,p,n,m,h,dx)
    integer,intent(in) :: nx,n
    real, intent(inout) :: x(nx),v(nx),a(nx),h(nx),rho(nx),m(nx),p(nx),cs(nx),dx
    real :: vs(nx),dt,a0(nx)
    dt=0.005
    a0=a
    x(:n)=x(:n)+dt*v(:n)+0.5*(dt)**2*a(:n)
    ! print*,a0(1),a(1)
    vs(:n)=v(:n)+dt*a(:n)
    call derivs(cs,rho,p,n,a,nx,x,m,h,dx,v)
    ! a(1)=0.
    ! a(n)=0.
    ! print*,
    ! v(:)=vs(:)+0.5*dt(a(:)-ab(:))
    v(:n)=vs(:n)+0.5*dt*(a(:n)-a0(:n))
    ! v(:n)=v(:n)+dt*a(:n)+0.5*dt*(a(:n)-a0(:n))
    ! print*,(a(1)-a0(1))
    ! print*,v(1),a(1),a0(1),dt
    a=a0
    ! a(1)=0.
    ! a(n)=0.
    !
    v(1)=0.
    v(n)=0.
    ! print*,a
    !
    ! read*,

    ! a(1)=0.
    ! a(n)=0.



  end subroutine
end module
