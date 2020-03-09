module equation
  implicit none

contains
  subroutine equation_of_state(cs,rho,p,n)
    integer,intent(in) :: n
    real,intent(in) :: rho(n)
    real,intent(out) ::p(n), cs(n)


    p=cs**2*rho


  end subroutine

end module
