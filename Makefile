# fortran compiler
.SUFFIXES: .f90 .f .c .o
FC = mpif90
CCC = mpicc


FFLAGS := -autodouble -fpconstant -O0 -132 -DUNDERSCORE -DNONHOOK # -cpp -ipo  
CCFLAGS := -O2 -DUNDERSCORE -DNONHOOK 
#FFLAGS := -O3 -ffixed-line-length-none -fdefault-real-8 

BIG := # -mcmodel=medium
DBG :=# -g -traceback
PROF :=# -pg
OMP := #-openmp

# 2decomp&fft
include libs/2decomp_fft/src/Makefile.inc
INCLUDE = -I libs/2decomp_fft/include 
#
LIB = -L libs/2decomp_fft/lib -l2decomp_fft -L libs/fft/lib -lfft -L libs/nonroot/SPHEREPACK/spherepack3.2_intel_openmpi/lib  -lspherepack -L libs/nonroot/gsl/gsl-1.16/lib -lgsl -lgslcblas


TARGET = inho

SRC = param.f90 Update_Pos.f90 common.f90 init.f90 initmpi.f90 initsolver.f90 bound.f90 boundc.f90 chkdiv.f90 chkdivc.f90 chkdt.f90 mom.f90 momc.f90 crnk.f90 crnkc.f90 fftuvw.f90 fftuvwc.f90 fillps.f90 fillpsc.f90 zredistribute.f90 zredistributec.f90 solver.f90 solverc.f90 correc.f90 correcc.f90 initparticles.f90 kernel.f90 forcing.f90 interp_spread.f90 interp_spreadc.f90 interpolation.f90 output.f90 outputc.f90 post.f90 main.f90 surf_vol.f90 sph1.f90 sph2.f90 sph3.f90 math.f90 sph.f90 capsulesolid.f90 loadd.f90 indicator.f90 collisions.f90 

# SRCF =  ##fortran 77 rule for .f file

SRCC = sphharm.c membmodel.c


OBJ  = $(SRC:.f90=.o)
OBJF = $(SRCF:.f=.o)
OBJC = $(SRCC:.c=.o)

.f.o:
	$(FC) $(FFLAGS) $(BIG) $(DBG) $(PROF) $(OMP)  -c $< $(LIB)	

.f90.o:
	$(FC) $(FFLAGS) $(BIG) $(DBG) $(PROF) $(OMP)  -c $< $(LIB)	

.c.o:
	$(CCC) $(CCFLAGS) -c $< $(LIB) -I /scratch/mehd/nonroot/gsl/gsl-1.16/include/gsl

all: $(TARGET)

$(TARGET): $(OBJC) $(OBJF) $(OBJ)
	$(FC) $(FFLAGS) $(BIG) $(DBG) $(PROF) $(OMP) -DUNDERSCORE -DNONHOOK -o $@ $(OBJC) $(OBJF) $(OBJ) $(LIB)

.PHONY: clean  
clean:
	rm -f *.o *.mod $(TARGET) .depend

veryclean: 
	rm -f *.o *.mod $(TARGET) .depend libs/fft/fft.o /libs/fft/lib/libfft.a; cd libs/2decomp_fft/src/; make clean; cd ../../../

libraries:
	cd libs/fft; $(FC) $(FFLAGS) $(BIG) $(DBG) $(PROF) $(OMP) -c fft.f; ar qc libfft.a fft.o; mv libfft.a lib; cd ../../;
	cd libs/2decomp_fft/src/; make; cd ../../../

%.o: %.f90
	$(FC) $(INCLUDE) -c -o $@ $(FFLAGS) $(BIG) $(DBG) $(PROF) $(OMP) $<




.depend dep:
	./.makedepo $(SRC) > .depend

include .depend
