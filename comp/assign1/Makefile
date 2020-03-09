FC=gfortran
# FFLAGS=-O3 -Wall -Wextra -fdefault-real-8 -fdefault-double-8
FFLAGS=-O3 -Wall -Wextra -fdefault-real-8 -fdefault-double-8
# FFLAGS=-O3 -Wall -Wextra
# FFLAGS=-O3 -Wall -Wextra -freal-4-real-16
SRC= equation.f90 kernel.f90 density.f90 inout.f90 set.f90 ghost.f90 assign.f90
OBJ=${SRC:.f90=.o}


%.o: %.f90
		$(FC) $(FFLAGS) -o $@ -c $<

solve: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

clean:
	rm *.o *.mod
