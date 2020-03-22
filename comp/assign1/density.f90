module density
  use set
  implicit none
  public :: dw


contains
  subroutine get_density(m,x,rho,nx,n,h)
    ! set arrays as inputs and outputs of subroutine
    integer,intent(in) :: nx,n
    real,intent(in) :: m(nx), x(nx),h(nx)
    real,intent(out) :: rho(nx)
    ! define place holder density, kernel as well as specific smoothing length and particle position
    real :: w(nx)
    ! define integers
    integer :: i



    !call kernel function and to assign density to each particle

    do i=1,nx
      ! call kernel
      call kern(h(i),nx,x,x(i),w,n)
      ! sum over all elements to find density per particle

      rho(i)=sum(m*w)
      ! print*,h(i),nx,x,x(i),w,n,rho(i)
      ! read*,
    enddo

    ! print*,'c',rho


  end subroutine



  subroutine kern(hin,nx,x,xin,w,n)
    ! specific smoothing length and specific and general positions
    integer,intent(in) :: nx,n
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
    do i=1,nx
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


  subroutine equation_of_state(cs,rho,p,nx,n)
    ! setup density, pressure and sound speed arrays
    integer,intent(in) :: nx,n
    real,intent(in) :: rho(nx)
    real,intent(out) :: p(nx), cs(nx)

    ! calculate isothermal pressure

    p=cs**2*rho
    ! print*,rho


  end subroutine

  subroutine get_accel(rho,p,n,a,nx,x,m,h)
    integer,intent(in) :: n,nx
    real,intent(in) :: rho(nx),m(nx)
    real,intent(inout) ::  p(nx),x(nx),h(nx)
    real, intent(out) :: a(nx)
    integer :: i,j
    real :: qa,qb,ap(nx)

    qa=0.
    qb=0.
    ! p=rho
    ! call equation_of_state(cs,rho,p,nx,n)

    do i=1,nx
      ! do j=1,n
        ! call kernel
        ! call kern(h(i),nx,x,x(i),w,n)
        ! compute the indivdual terms of the density sum
        do j=1,nx
          ap(j)=m(j)*((p(i)+qa)/rho(i)**2*dw(x(i),x(j),h(i),j)+(p(j)+qa)/rho(j)**2*dw(x(i),x(j),h(j),j))

        enddo
        ! print*,'bob'
        ! sum over all elements to find density per particle
        a(i)=-sum(ap)
      enddo

      ! read*,

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

      if (abs(xin-xon) .le. 10.**(-8)) then
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


  subroutine derivs(cs,rho,p,n,a,nx,x,m,h,dx,v)
    integer,intent(in) :: n,nx
    real, intent(in) ::  dx
    real, intent(inout) :: rho(nx), cs(nx),p(nx), x(nx),h(nx),v(nx),m(nx)
    real,intent(out) :: a(nx)



    call get_density(m,x,rho,nx,n,h)

    call equation_of_state(cs,rho,p,nx,n)
    ! print*,p
    call get_accel(rho,p,n,a,nx,x,m,h)
    if (p(23)==rho(23)) then
      print*, 'yes'
    else
      print*,'no'
    endif
    ! print*,a
    call set_ghosts(rho,nx,x,v,dx,m,cs,n,h)



  end subroutine

end module
