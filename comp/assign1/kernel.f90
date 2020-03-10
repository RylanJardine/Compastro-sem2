module kernel
  implicit none

contains
  subroutine kern(hin,nx,x,xin,w,n)
    real,intent(in) :: hin,x(nx),xin
    real,intent(out) :: w(nx)
    integer :: i
    real :: sig3,q(nx)
    integer,intent(in) :: nx,n

    sig3=2./(3.*hin)

    do i=1,nx
      q(i)=abs(xin-x(i))/hin
      if (0. .LE. q(i) .and. q(i) <= 1.) then
        w(i)=sig3*(1.-3./2.*q(i)**2.*(1-q(i)/2.))
      elseif (1. < q(i) .and. q(i) <= 2.) then
        w(i)=sig3/4.*(2.-q(i))**3.
      else
        w(i)=0.
      endif
    enddo

  end subroutine
end module
