module physics
  use set2
  use param
  implicit none
  public :: dw


contains
  !find density of particles
  subroutine get_density(m,x,rho,nx,n,h,ng)
    ! define type of variables
    integer,intent(in) :: nx,n,ng
    real,intent(in) :: m(nx), x(nx)
    real,intent(out) :: rho(nx)
    real,intent(inout) ::h(nx)
    real ::rhop(nx),q
    integer :: i,j




    !loop over i,j to find density per particle
    do i=1,n+ng
      rho(i)=0.
      do j=1,n+ng
      ! define q input of kernel function
      q=abs(x(i)-x(j))/h(i)
      !create each iteration of density sum
      rhop(j)=m(j)*w(q)*2./(3.*h(i))
      ! sum over all elements to find density per particle
      rho(i)=rho(i)+rhop(j)
      enddo
    enddo


  end subroutine

  !calculate kernel
  real function w(q)
    real,INTENT(IN) :: q

      ! set kernel hybrid function
      if (0. .LE. q .and. q <= 1.) then
        w=(1.-3./2.*q**2.*(1-q/2.))
      elseif (1. < q .and. q <= 2.) then
        w=1./4.*(2.-q)**3.
      else
        w=0.
      endif

  end function


  !find pressure/sound speed
  subroutine equation_of_state(cs,rho,p,nx,n,ng,u)
    ! setup density, pressure and sound speed arrays
    integer,intent(in) :: nx,n,ng
    real,intent(in) :: rho(nx),u(nx)
    real,intent(inout) :: p(nx), cs(nx)

    if (y==3) then
      !calculate adiabatic EOS (gamma N.E. 1)
      p(1:n+ng)=(gamma-1)*rho(1:n+ng)*u(1:n+ng)
      cs(1:n+ng)=sqrt(gamma*p(1:n+ng)/rho(1:n+ng))
    else
      ! calculate isothermal EOS (gamma=1)
      p(1:n+ng)=rho(1:n+ng)
      cs(1:n+ng)=sqrt(gamma*p(1:n+ng)/rho(1:n+ng))
    endif


  end subroutine

  !get accel, du/dt
  subroutine get_accel(rho,p,n,a,nx,x,m,h,v,cs,ng,du)
    integer,intent(in) :: n,nx,ng
    real,intent(in) :: rho(nx),m(nx),v(nx),cs(nx)
    real,intent(inout) ::  p(nx),x(nx),h(nx)
    real, intent(out) :: a(nx),du(nx)
    integer :: i,j
    real :: qa,qb,ap(nx),vab,dub(nx)

    do i=1,n+ng

      a(i)=0.
      du(i)=0.

        do j=1,n+ng
          !define viscosity terms qa and qb
          vab=(v(i)-v(j))*((x(i)-x(j))/abs(x(i)-x(j)))
          qa=visc(rho(i),vab,cs(i))
          qb=visc(rho(j),vab,cs(j))
          !sum over individual components of accel to find accel per particle
          ap(j)=m(j)*((p(i)+qa)/rho(i)**2*dw(x(i),x(j),h(i))+(p(j)+qb)/rho(j)**2*dw(x(i),x(j),h(j)))
          a(i)=a(i)-ap(j)
          !sum over individual components of du/dt to find du/dt per particle
          dub(j)=m(j)*(p(i)+qa)/rho(i)**2*(v(i)-v(j))*dw(x(i),x(j),h(i))
          du(i)=du(i)+dub(j)
        enddo

      enddo


  end subroutine

  !define grad w
  real function dw(xin,xon, hin)
    real,intent(in) :: xin,xon, hin
    real :: sig3,q,dq

    !define normalisation
    sig3=2./(3.*hin)

      ! define parameter q,dq
      q=abs(xin-xon)/hin
      dq=(xin-xon)/(hin*abs(xin-xon))


      !criteria to remove divide by zero errors
      if (abs(xin-xon) .le. tiny(1.)) then
        dq=0.
      endif

      ! set kernel hybrid function
      if (0. .LE. q .and. q <= 1.) then
        dw=sig3*(-3.*q*dq+9./4.*q**2*dq)
      elseif (1. < q .and. q <= 2.) then
        dw=3.*sig3/4.*(2.-q)**2*(-dq)
      else
        dw=0.
      endif




  end function

  !set ghosts, density, EOS, var smoothing length, and accel
  subroutine derivs(cs,rho,p,n,a,nx,x,m,h,v,ng,du,u)
    integer,intent(in) :: n,nx
    real, intent(inout) :: rho(nx), cs(nx),p(nx), x(nx),h(nx),v(nx),m(nx),du(nx),u(nx)
    real,intent(out) :: a(nx)
    integer :: i
    integer,intent(inout) :: ng

    !set ghost particles
    call set_ghosts(rho,nx,x,v,m,cs,n,h,a,p,ng,u,du)

    !loop over density and ghosts
    do i=1,3

      call get_density(m,x,rho,nx,n,h,ng)
      call set_ghosts(rho,nx,x,v,m,cs,n,h,a,p,ng,u,du)

      !variable smoothing length if chosen
      if (vvhh==1) then
        h(1:n+ng)=1.2*m(1:n+ng)/rho(1:n+ng)
      endif

    enddo

    !set pressure sound speed, accel and du/dt
    call equation_of_state(cs,rho,p,nx,n,ng,u)

    call get_accel(rho,p,n,a,nx,x,m,h,v,cs,ng,du)


  end subroutine

  !define viscous terms
  real function visc(rho,vab,cs)
    real,intent(in) :: vab,rho,cs
    real :: vsig

    !define vsig and viscous parameters
    if (vab<0.) then
      vsig=alpha*cs-beta*(vab)
      visc=-0.5*rho*vsig*vab
    else
      visc=0.
    endif


  end function

end module
