module param
  implicit none
  !define GLOBAL variables
  real :: alpha,beta
  integer ::vvhh
  real :: dx,dx2
  real :: gamma,xmax,xmin,xmid,xmid2,xmin2
  integer :: y



contains
  subroutine vars(alp,bet,vvh,z)
    real,INTENT(IN) :: alp,bet
    integer,intent(in) ::vvh,z
    !define visc and variable smoothing length options
    alpha=alp
    beta=bet
    vvhh=vvh
    y=z

    !select grid size and gamma value for selected problem

    !Isothermal shock tube
    if (z==2) then
      dx=0.001
      dx2=0.01
      gamma=1.
      xmin=-1.5
      xmid=-1.
      xmid2=0.
      xmax=0.5
      xmin2=-0.5
      !Adiabatic shock tube
    else if (z==3) then
      dx=0.001
      dx2=0.008
      gamma=1.4
      xmin=-1.5
      xmid=-1.
      xmid2=0.
      xmax=0.5
      xmin2=-0.5
      !Standing wave
    else if (z==1) then
      gamma=1.
      xmin=0.
      xmax=1.
    endif


  end subroutine


end module
