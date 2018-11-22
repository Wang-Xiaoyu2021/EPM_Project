# Makefile for Unified fitting program

#FC=pgf90  -fast -pc 64 -Mvect=sse -Mcache_align -Kieee
FC=ifort -O3 -axN
#FFLAGS= -dp -Ounroll2
FFLAGS= 
#LIB =-lapp -lc -lm -lu -lf -lsci -lfi -lsma -lcomm
LIB = /opt/LAPACK/lapack_LINUX.a  /opt/LAPACK/blas_LINUX.a
#LIB=-llapack -lblas /usr/lib/libg2c.so.0
#LIB=-L /usr/local/pgi/linux86/lib -llapack -lblas

HEAD = epm.inc


OBJ=monitor.o \
funx.o \
pwk.o \
recvec.o \
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
funk.o \
dqed.o\
funklocal.o\
dqedhdjwl.o\
amebsa.o\
amotsa.o\

SRC1=band_str.f

SRC2=parabolic.f

SRC3=fit_dqed.f

SRC4=fit_sa_dqed.f

SRC5=fit_pgapack.f evaluate_pgapack.f funk_pgapack.f

.SUFFIXES:
%.o: %.f
	$(FC) ${FFLAGS} -c $<
%.o: %.f90
	$(FC) ${FFLAGS} -c $<

all: band_str.x parabolic.x  \
     multifit_dqed multifit_sa_dqed\
     multifit_ga


band_str.x: $(OBJ) $(SRC1) 
	$(FC) -c band_str.o band_str.f
	$(FC) -o band_str.x band_str.o $(OBJ) $(LIB)


parabolic.x: $(OBJ) $(SRC2) 
	$(FC) -c parabolic.o parabolic.f
	$(FC) -o parabolic.x  parabolic.o $(OBJ) $(LIB)

multifit_dqed: $(OBJ) $(SRC3)
	$(FC) -c fit_dqed.o fit_dqed.f
	$(FC) -o multifit_dqed fit_dqed.o $(OBJ) $(LIB)

multifit_sa_dqed: $(OBJ) $(SRC4) 
	$(FC) -c fit_sa_dqed.o fit_sa_dqed.f
	$(FC) -o multifit_sa_dqed fit_sa_dqed.o $(OBJ) $(LIB)

multifit_ga: $(OBJ) $(SRC5)
	$(FC) -c ga170.o ga170.f
	$(FC) -o multifit_ga ga170.o $(OBJ) $(LIB)


