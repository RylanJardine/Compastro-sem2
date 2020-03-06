module kernel
  implicit none

contains
  subroutine kern(hin,nx,x,xin,w)
    real,intent(in) :: hin,x(nx),xin
    real,intent(out) :: w(nx)
    integer :: i
    real :: sig3,q(nx)
    integer,intent(in) :: nx

    sig3=2./(3.*hin)

    ! do i=1,nx
    !   q(i)=abs(xin-x(i))
    !   if (0 .LE. q(i) .LE. 1) then
    !     w(i)=sig3*(1-3./2.*q(i)**2*(1-q(i)/2.))
    !   elseif (1 < q(i) .LE. 2) then
    !     w(i)=sig3/4.*(2.-q(i))**3
    !   else
    !     w(i)=0.
    !   endif
    ! enddo


    do i=1,nx
      q(i)=abs(xin-x(i))
      if (0 .LE. q(i)/hin .and. q(i)/hin <= 1) then
        w(i)=sig3*(1-3./2.*q(i)**2*(1-q(i)/2.))
      elseif (1 < q(i)/hin .and. q(i)/hin <= 2) then
        w(i)=sig3/4.*(2.-q(i))**3
      else
        w(i)=0.
      endif
    enddo

  end subroutine
end module
