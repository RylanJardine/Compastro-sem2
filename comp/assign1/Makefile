FC=gfortran
FFLAGS=-O3 -Wall -Wextra -fdefault-real-8 -fdefault-double-8
# FFLAGS=-O3 -Wall -Wextra -freal-4-real-16
SRC=output.f90 set.f90 assign.f90
OBJ=${SRC:.f90=.o}


%.o: %.f90
		$(FC) $(FFLAGS) -o $@ -c $<

solve: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

clean:
	rm *.o *.mod
