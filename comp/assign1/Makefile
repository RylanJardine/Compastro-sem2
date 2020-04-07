FC=gfortran

# FFLAGS=-O3 -Wall -Wextra -fdefault-real-8 -fdefault-double-8
# FFLAGS=-O3 -Wall -Wextra
# FFLAGS=-O3 -Wall -Wextra -freal-4-real-16

FFLAGS=-O3 -Wall -Wextra -fdefault-real-8 -fdefault-double-8 -fcheck=all -std=f2018 -finit-real=NaN
SRC= set2.f90 density.f90 inout.f90  integrator.f90 assign.f90
# SRC= set.f90 set2.f90 density.f90 inout.f90  integrator.f90 assign.f90
OBJ=${SRC:.f90=.o}


%.o: %.f90
		$(FC) $(FFLAGS) -o $@ -c $<

solve: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

run: solve
		./solve

clean:
	rm *.o *.mod snap*

git:
	git add --all
	git commit -m 'lazy commit'
	git push
