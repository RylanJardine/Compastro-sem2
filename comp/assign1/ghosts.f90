module set_ghosts
  implicit none

contains
  subroutine set_ghosts(n,ng,x,v,rho,m,h,p)
    integer,intent(in) :: n
    integer,intent(out) :: ng
    real, intent (inout) ::x,v,rho,m,h,p
    integer :: i,j
    real, parameter :: radkern=2.

    dx=xmax-xmin
    if (abs(xmax-xmin) < tiny(0.)) print*,'Error xmax-xmin is 0'

    ng=0

    do i=1,n
      if (x(i)+radkern*h(i)>xmax) then
        ng=ng+1
        j=n+ng
        x(j)=x(i)-dx
        v(j)=v(i)
        m(j)=m(i)
      else if (x(i)-radkern*h(i)<xmin)
        ng=ng+1
        j=n+ng
        x(j)=x(i)+dx
        v(j)=v(i)
        m(j)=m(i)
      endif
    enddo

  end subroutine

end module ghosts
