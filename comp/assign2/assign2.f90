program assign
  use integrator
  implicit none
  real :: x(2),v(2),a(2),r,t,l,u,dl
  integer :: i,j
  real,parameter :: dt=0.1




  t=0.
  x(1)=1.-e
  x(2)=0.

  v(1)=0.
  v(2)=sqrt((1+e)/(1-e))


  r=sqrt(x(1)**2+x(2)**2)
  a(1)=-G*M*x(1)/r**3
  a(2)=-G*M*x(2)/r**3
  l=x(1)*v(2)-x(2)*v(1)
  dl=l
  u=0.5*(v(1)**2+v(2)**2)-1./r


  print*,'Please select Leapfrog (1) or Runge Kutta (2)'
  read*,j
  ! if (j==1) then
  !   call tim(x,v,a,dt)
  ! else if (j==2) then
  !   call rung(x,v,a,dt)
  ! else
  !   print*, 'Default = Leapfrog'
  !   call tim(x,v,a,dt)
  ! endif

  open(1,file='param.dat',status='replace',action='write')
  write(1,*)'# x y v_x v_y a_x a_y t l e dl'
  write(1,*) x(1),x(2),v(1),v(2),a(1),a(2),t,l,u,l-dl
  do i=1,5000
    t=t+dt
    if (j==1) then
      call tim(x,v,a,dt)
    else if (j==2) then
      call rung(x,v,a,dt)
    else
      call tim(x,v,a,dt)
    endif
    l=x(1)*v(2)-x(2)*v(1)
    r=sqrt(x(1)**2+x(2)**2)
    u=0.5*(v(1)**2+v(2)**2)-1./r
  write(1,*) x(1),x(2),v(1),v(2),a(1),a(2),t,l,u,l-dl
  enddo
  close(1)

end program
