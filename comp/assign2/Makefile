FC=gfortran

# FFLAGS=-O3 -Wall -Wextra -fdefault-real-8 -fdefault-double-8
# FFLAGS=-O3 -Wall -Wextra -freal-4-real-16

FFLAGS=-O3 -Wall -Wextra -fdefault-real-16 -fdefault-double-8 -fcheck=all -std=f2018 -finit-real=NaN
SRC= integrate.f90 assign2.f90

OBJ=${SRC:.f90=.o}


%.o: %.f90
		$(FC) $(FFLAGS) -o $@ -c $<

solve: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

run: solve
		./solve

clean:
	rm *.o *.mod param.dat

git:
	git add --all
	git commit -m 'lazy commit'
	git push
