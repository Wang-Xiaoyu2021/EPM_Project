# Makefile for Unified fitting program

#module load hdf5 libraries/intel-mkl-11.1.064 libraries/intel_fftw-11.1.064

FC=ifort
MFC=mpif90
#FC=pathf90 -Ofast
#FFLAGS= -dp -Ounroll2
FFLAGS=-O3
LDFLAGS=
INCLUDE=
#LIB = /opt/LAPACK/lapack_LINUX.a  /opt/LAPACK/blas_LINUX.a
#LIB=-L/opt/acml4.0.0/pathscale64/lib -lacml
#LIB=-L/copt/atlas/3.6.0/lib64 -llapack -lf77blas -lcblas -latlas -lg2c
#LIB=/copt/atlas/3.6.0/lib64/liblapack.a /copt/atlas/3.6.0/lib64/libf77blas.a /copt/atlas/3.6.0/lib64/libatlas.a
#LIB=-L/copt/2.0/pgi/7.1/linux86-64/7.1-3/lib -llapack -lblas 
#LIB=-llapack -lblas 
LIB=${BLASLIB}


HEAD = epm.inc


OBJ=monitor.o \
funx.o \
funxsqs.o \
pwk.o \
recvec.o \
gaussj.o \
gridg.o \
pwset.o \
fft.o \
vcell.o \
settp.o \
readdata.o \
formf.o \
skip.o \
solve_so.o \
solve_nso.o \
getline.o\
cfft.o\
cfftd.o\
funjwl.o \
funk.o 
 

OBJ1=band_str.o

OBJ2=parabolic.o

OBJ3=fit_dqed.o dqed.o dqedhdjwl.o

OBJ4=fit_sa_dqed.o dqed.o dqedhdjwl.o amebsa.o amotsa.o funklocal.o

OBJ5=fit_pgapack.o evaluate_pgapack.o funk_pgapack.o timestamp.o

.SUFFIXES:
%.o: %.f
	$(FC) ${FFLAGS} -c $<
%.o: %.f90
	$(FC) ${FFLAGS} -c $<

all: band_str.x parabolic.x  \
     multifit_dqed multifit_sa_dqed


band_str.x: $(OBJ) $(OBJ1) 
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(OBJ1) $(LIB)

parabolic.x: $(OBJ) $(OBJ2) 
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(OBJ2) $(LIB)

multifit_dqed: $(OBJ) $(OBJ3)
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(OBJ3) $(LIB)

multifit_sa_dqed: $(OBJ) $(OBJ4) 
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(OBJ4) $(LIB)

#multifit_ga: $(OBJ) $(OBJ5)
#	$(MFC) $(FFLAGS) -o $@ $(OBJ) $(OBJ5) $(LIB) $(LDFLAGS)


clean:
	rm -f *.o *.mod

