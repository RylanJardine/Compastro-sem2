module density
  ! use set
  use set2
  implicit none
  public :: dw


contains
  subroutine get_density(m,x,rho,nx,n,h,ng)
    ! set arrays as inputs and outputs of subroutine
    integer,intent(in) :: nx,n,ng
    real,intent(in) :: m(nx), x(nx)
    real,intent(out) :: rho(nx)
    real,intent(inout) ::h(nx)
    ! define place holder density, kernel as well as specific smoothing length and particle position
    real :: w(nx)
    ! define integers
    integer :: i,j


    !
    ! call kernel function and to assign density to each particle

    do i=1,n+ng
      ! call kernel
      call kern(h(i),nx,x,x(i),w,n,ng)
      ! sum over all elements to find density per particle
      rho(i)=sum(m(1:n+ng)*w(1:n+ng))

    enddo

    ! h=m/rho
    ! print*,'c',rho

    ! do i=1,n+ng
    !   rho(i)=0.
    !   do j=1,n+ng
    !     rho(i)=rho(i)+m(j)*w(x(i),x(j),h(j))
    !
    !   enddo
    ! enddo


  end subroutine

  real function w(xi,xj,hin)

    ! specific smoothing length and specific and general positions

    real,intent(in) :: hin,xi,xj


    ! define additional parameters
    real :: sig3,q

    ! set normalising factor for 1D
    sig3=2./(3.*hin)

    ! create array of kernel values
      ! define parameter q
      q=abs(xi-xj)/hin


      ! set kernel hybrid function
      if (0. .LE. q .and. q <= 1.) then
        w=sig3*(1.-3./2.*q**2.*(1-q/2.))
      elseif (1. < q .and. q <= 2.) then
        w=sig3/4.*(2.-q)**3.
      else
        w=0.
      endif


  end function

  subroutine kern(hin,nx,x,xin,w,n,ng)
    ! specific smoothing length and specific and general positions
    integer,intent(in) :: nx,n,ng
    real,intent(in) :: hin,x(nx),xin

    ! kernel
    real,intent(out) :: w(nx)
    ! set integers
    integer :: i

    ! define additional parameters
    real :: sig3,q(nx)

    ! set normalising factor for 1D
    sig3=2./(3.*hin)

    ! create array of kernel values
    do i=1,n+ng
      ! define parameter q
      q(i)=abs(xin-x(i))/hin


      ! set kernel hybrid function
      if (0. .LE. q(i) .and. q(i) <= 1.) then
        w(i)=sig3*(1.-3./2.*q(i)**2.*(1-q(i)/2.))
      elseif (1. < q(i) .and. q(i) <= 2.) then
        w(i)=sig3/4.*(2.-q(i))**3.
      else
        w(i)=0.
      endif

    enddo

  end subroutine


  subroutine equation_of_state(cs,rho,p,nx,n,ng,u)
    ! setup density, pressure and sound speed arrays
    integer,intent(in) :: nx,n,ng
    real,intent(in) :: rho(nx),u(nx)
    real,intent(inout) :: p(nx), cs(nx)
    ! real, parameter :: gamma=1.

    ! calculate isothermal pressure
    ! p(1:n+ng)=rho(1:n+ng)

    if (y==3) then
      p(1:n+ng)=(gamma-1)*rho(1:n+ng)*u(1:n+ng)
      cs(1:n+ng)=sqrt(gamma*p(1:n+ng)/rho(1:n+ng))
    else
      p(1:n+ng)=rho(1:n+ng)
      cs(1:n+ng)=sqrt(gamma*p(1:n+ng)/rho(1:n+ng))
    endif

    ! if (y==3) then
    !   p(1:n+ng)=(gamma-1)*rho(1:n+ng)*u(1:n+ng)
    !   cs(1:n+ng)=sqrt(gamma*p(1:n+ng)/rho(1:n+ng))
    ! elseif (y==2) then
      ! p(1:n+ng)=rho(1:n+ng)
      ! cs(1:n+ng)=sqrt(gamma*p(1:n+ng)/rho(1:n+ng))
    ! endif


  end subroutine

  subroutine get_accel(rho,p,n,a,nx,x,m,h,v,cs,ng,du)
    integer,intent(in) :: n,nx,ng
    real,intent(in) :: rho(nx),m(nx),v(nx),cs(nx)
    real,intent(inout) ::  p(nx),x(nx),h(nx)
    real, intent(out) :: a(nx),du(nx)
    integer :: i,j
    real :: qa,qb,ap(nx),vab,dub(nx)

    ! qa=0.
    ! qb=0.

    do i=1,n+ng

      a(i)=0.
      du(i)=0.
        ! compute the indivdual terms of the density sum
        do j=1,n+ng
          vab=(v(i)-v(j))*((x(i)-x(j))/abs(x(i)-x(j)))
          qa=visc(rho(i),vab,cs(i))
          qb=visc(rho(j),vab,cs(j))
          ap(j)=m(j)*((p(i)+qa)/rho(i)**2*dw(x(i),x(j),h(i),j)+(p(j)+qb)/rho(j)**2*dw(x(i),x(j),h(j),j))
          a(i)=a(i)-ap(j)

          dub(j)=m(j)*(p(i)+qa)/rho(i)**2*(v(i)-v(j))*dw(x(i),x(j),h(i),j)
          du(i)=du(i)+dub(j)
        enddo
        ! sum over all elements to find density per particle
        ! a(i)=-sum(ap(1:n+ng))
      enddo
      ! print*,sum(m(1:n)*a(1:n))

  end subroutine


  real function dw(xin,xon, hin,i)
    real,intent(in) :: xin,xon, hin
    real :: sig3,q,dq
    integer, intent(in) :: i
    ! integer :: i

    sig3=2./(3.*hin)

    ! create array of kernel values
    ! do i=1,nx
      ! define parameter q
      q=abs(xin-xon)/hin
      dq=(xin-xon)/(hin*abs(xin-xon))

      if (abs(xin-xon) .le. 10.**(-16)) then
        dq=0.
      endif
      ! print*,xin-xon


      ! print*,q,dq
      ! set kernel hybrid function
      if (0. .LE. q .and. q <= 1.) then
        dw=sig3*(-3.*q*dq+9./4.*q**2*dq)
      elseif (1. < q .and. q <= 2.) then
        dw=3.*sig3/4.*(2.-q)**2*(-dq)
      else
        dw=0.
      endif
      if (dw < 10.**(-18)) then
        ! print*,dw,i
      endif
    ! enddo


  end function


  subroutine derivs(cs,rho,p,n,a,nx,x,m,h,dx,v,ng,du,u)
    integer,intent(in) :: n,nx
    real, intent(in) ::  dx
    real, intent(inout) :: rho(nx), cs(nx),p(nx), x(nx),h(nx),v(nx),m(nx),du(nx),u(nx)
    real,intent(out) :: a(nx)
    integer :: i
    integer,intent(inout) :: ng

    if (y==1) then
      call set_ghosts(rho,nx,x,v,m,cs,n,h,a,p,ng,u,du)
    else
      call set_ghosts2(rho,nx,x,v,m,cs,n,h,a,p,ng,u,du)
    endif


    do i=1,3

      call get_density(m,x,rho,nx,n,h,ng)
      if (y==1) then
        call set_ghosts(rho,nx,x,v,m,cs,n,h,a,p,ng,u,du)
      else
        call set_ghosts2(rho,nx,x,v,m,cs,n,h,a,p,ng,u,du)
      endif

      h(1:n+ng)=1.2*m(1:n+ng)/rho(1:n+ng)

    enddo


    call equation_of_state(cs,rho,p,nx,n,ng,u)

    call get_accel(rho,p,n,a,nx,x,m,h,v,cs,ng,du)


  end subroutine

  real function visc(rho,vab,cs)
    real,intent(in) :: vab,rho,cs
    real,parameter :: alpha=1.,beta=2.
    real :: vsig

    if (vab<0.) then
      vsig=alpha*cs-beta*(vab)
      visc=-0.5*rho*vsig*vab
    else
      visc=0.
    endif


  end function

end module
