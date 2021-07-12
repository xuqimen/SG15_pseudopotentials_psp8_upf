# specify the path MKLROOT if necessary, e.g.
# MKLROOT = /opt/intel/compilers_and_libraries_2017.4.196/linux/mkl

#FFLAGS = 

# if need to run with DEBUG mode, add -DDBUG
#CPPFLAGS = -m64 -I${MKLROOT}/include -I include/ -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -ldl -lrt -O3 -Wall -g -DDEBUG

# if run locally, use the following (add -DDBUG to turn on debug mode)
#CPPFLAGS = -m64 -I${MKLROOT}/include -I include/ -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -lrt -O3 -Wall -g -DDEBUG
CPPFLAGS = -lm -ldl -lrt -Wall -g
#FPPFLAGS = 

#SOURCEC = main.c
#OBJSC = main.o initialization.o readfiles.o structRelaxation.o md.o electrostatics.o electronicGroundState.o electronDensInit.o MatVecRoutines.o nonlocal.o occupation.o tools.o linearSolver.o mixing.o eigenSolver.o ExchangeCorrelation.o energy.o forces.o finalization.o
OBJSC = TM2PSP.o


#LIBBASE = ../lib/sparc

override CC=gcc

all: TM2PSP

TM2PSP: $(OBJSC)
	$(CC) $(CPPFLAGS) -o TM2PSP $^

.PHONY: clean
clean:
	rm -f  *.o TM2PSP
