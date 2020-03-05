module inout
  use set
  implicit none


contains
  subroutine output()
  real,intent(in) ::
  real,


  write(lu,*) ’# x y z ... rest of column labels’
  write(lu,*) time
  do i=1,n
    write(lu,*) x(i),y(i),z(i), rest of particle properties
  enddo
