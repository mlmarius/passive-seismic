f77=gfortran

# to check  implicit variables use option `-fimplicit-none`

.f.o:
	$(f77) -c $*.f -ffixed-line-length-80 --std=legacy -g -Wall -Wextra -Warray-temporaries -Wconversion -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan -lm

#  path of executables, leave empty for current dir
BIN =

RAYP=rayP
SOLVEP=solveP

all: rayP solveP

rayP: rayP.o
	$(f77) -o $(BIN)$(RAYP) $(RAYP).o

solveP: solveP.o
	$(f77) -o $(BIN)$(SOLVEP) $(SOLVEP).o

clean:
	rm -f *.o *.so rayP
