program assign
  ! Import modules for use
  ! use set
  use set2
  use inout
  use physics
  use integrator
  use param
  ! use ghost
  ! use equation
  implicit none
  ! set parameters and arrays of length nx
  integer, parameter :: nx=1500
  real :: x(nx), v(nx), m(nx), rho(nx),  p(nx), cs(nx),h(nx),a(nx), ek,po,mt,du(nx),u(nx)
  ! define the max and min of your grid and number of particles, n
  real :: t,dt,tprint,dtout,tstop,bet,alp
  integer :: n,ifile,z,ng,vvh,b,vh




  !DO not ever use my global integers z and y


  print*,'Select Standing Wave (1), Isothermal Shock tube (2) or Adiabatic Shock Tube Problem (3)'
  read*,z
  if (z==1) then
    call setup(rho,x,v,m,cs,n,h,u,z)
    tstop=5.
    dtout=0.05
  else if ((z==2) .or. (z==3)) then
    call setup_shock(rho,x,v,m,cs,n,h,a,p,u,z)
    tstop=0.2
    dtout=0.01
  else
    print*, 'Please Select (1),(2) or (3)'
  endif

  print*, 'Add viscosity (yes=1, no=0)?'
  read*,b
  if (b==1) then
    print*, 'Please select Alpha'
    read*,alp
    print*,'Please select Beta'
    read*,bet
  else if (b==0) then
    alp=0.
    bet=0.
  else
    print*, 'Default without viscosity'
    alp=0.
    bet=0.
  endif

  Print*,'Variable smoothing length (yes=1, no=0)?'
  read*,vh
  if (vh==1) then
    vvh=1
  else if (vh==0) then
    vvh=0
  else
    print*, 'Default with variable smoothing length'
    vvh=0
  endif

  call vars(alp,bet,vvh)

  ! call to setup initial conditions at time t=0



  call derivs(cs,rho,p,n,a,nx,x,m,h,v,ng,du,u)

  open(1,file='kin.dat',status='replace',action='write')
  po=sum(m(1:n)*a(1:n))
  ek=sum(0.5*m(1:n)*v(1:n)**2)
  mt=sum(m(1:n))
  t=0.
  write(1,*)'#t ek po mt'
  write(1,*)t,ek,po,mt


  ifile=0
  tprint=(ifile+1)*dtout
  call output(n,x,v,h,nx,rho,m,p,cs,a,t,ifile,ng,u,du)
  call tim(x,v,a,nx,dt,cs,rho,p,n,m,h,ng,du,u)

  do while (t<tstop)
    t=t+dt
    if (t>tprint) then
      ifile=ifile+1
      tprint=(ifile+1)*dtout
      call output(n,x,v,h,nx,rho,m,p,cs,a,t,ifile,ng,u,du)
    end if
    call tim(x,v,a,nx,dt,cs,rho,p,n,m,h,ng,du,u)
    po=sum(m(1:n)*a(1:n))
    mt=sum(m(1:n))
    ek=sum(0.5*m(1:n)*v(1:n)**2)
    write(1,*)t,ek,po,mt
  enddo
  close(1)

  print*,'hello world'


end program
