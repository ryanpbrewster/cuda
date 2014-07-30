
SRC = qc_kind.f90 \
      qc_constant.f90 \
      qc_rnd.f90 \
      qc_step.f90 \
      qc_geom.f90 \
      qc_basis.f90 \
      qc_input.f90 \
      nw_vectors.f90 \
      qc_psi.f90 \
      qc_mcmp2.f90 \
      qc_monte.f90 \
      main.f90 


LIBS =
FC       = gfortran
CC       = g++
CFLAGS1 = -O2 -Wall -pedantic -fbounds-check

include $(patsubst %, %/module.mk, $(MODULES))
OBJ :=  $(patsubst %.f90, %.o, $(filter %.f90, $(SRC))) \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(SRC))) 


mcmp2.exe: $(OBJ)
	$(FC) -o mcmp2.exe $(OBJ) $(LIBS)


.SUFFIXES :
.SUFFIXES : .f90 .o .cpp
.f90.o:
	$(FC) $(CFLAGS1) -c $< -o $@
.cpp.o:
	$(CC) $(CFLAGS1) -c $< -o $@

.PHONY: clean
clean:
	rm -f *.o *.mod mcmp2.exe *~
