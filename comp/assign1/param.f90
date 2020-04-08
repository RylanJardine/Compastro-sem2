module param
  implicit none
  real :: alpha,beta
  integer ::vvhh


  !For Problem 1
  ! real, parameter :: gamma=1.
  ! integer, parameter :: n=100

contains
  subroutine vars(alp,bet,vvh)
    real,INTENT(IN) :: alp,bet
    integer,intent(in) ::vvh
    alpha=alp
    beta=bet
    vvhh=vvh

  end subroutine


end module
