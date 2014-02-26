
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


PROJECT := mcmp2.exe

RM := rm -f

LIBS :=
FC       = gfortran
CFLAGS1 = -O2    -Wall -pedantic -fbounds-check -ffpe-trap=invalid,zero,overflow -fall-intrinsics

include $(patsubst %, %/module.mk, $(MODULES))
OBJ :=  $(patsubst %.f90, %.o, $(filter %.f90, $(SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(SRC))) 


$(PROJECT) : $(OBJ) 
	$(FC) $(LFLAGS1) -o $@ $(OBJ) $(LIBS)


.SUFFIXES :
.SUFFIXES : .f90 .o .f
.f90.o:
	$(FC) $(CFLAGS1) -c $< -o $@

.f.o:
	$(FC) $(CFLAGS1) -c $< -o $@

.PHONY: veryclean
veryclean:
	$(RM) *.o *.mod mcmp2.exe 

.PHONY: clean
clean:
	$(RM) $(OBJ) *.mod mcmp2.exe

