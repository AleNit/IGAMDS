########################################## DEFINITIONS

EXEC = exe.e

FC  = gfortran

# optimization 
FFLAGS = -O3 -msse2
#FFLAGS = -O0 

# precision
FFLAGS += -fdefault-real-8  -fdefault-double-8

# flags for debug at runtime
#FFLAGS += -g
#FFLAGS += -ffpe-trap=zero,overflow,underflow 
#FFLAGS += -ffpe-summary=all -fcheck=all 
#FFLAGS += -ffpe-trap=invalid -fcheck=all  
#FFLAGS += -fbacktrace -fbounds-check
#FFLAGS += -fcheck=all -Wall

# libreries
LIB =   -llapack 
LIB +=  -L/usr/lib/x86_64-linux-gnu/hdf5/serial/ -lz -lhdf5_fortran -lhdf5 -ldl

# include
INCL = -I/usr/include/hdf5/serial/

OBJS = utilities.o mgmres.o refinement.o metrics.o bcs.o \
       std_mat.o alloc.o postpr.o solgw.o

default: executable 

executable:
	@echo
	@echo '--------------------------------------------------------------'
	@echo '  making executable : ' $(EXEC)
	@echo '--------------------------------------------------------------'
	@echo
	@echo   
	make .e
	@echo
#	make partclean
		

.e: main.F03 $(OBJS)
	$(FC) $(FFLAGS) $(INCL) -o $(EXEC) main.F03 $(OBJS)  $(LIB)

utilities.o: utilities.f03
	$(FC) $(FFLAGS) -c utilities.f03

mgmres.o: mgmres.f90
	$(FC) $(FFLAGS) -c mgmres.f90

refinement.o: refinement.f03 utilities.f03
	$(FC) $(FFLAGS) -c refinement.f03 

metrics.o: metrics.f03 utilities.f03 refinement.f03
	$(FC) $(FFLAGS) -c metrics.f03 

bcs.o: bcs.f03 metrics.f03 refinement.f03 utilities.f03
	$(FC) $(FFLAGS) -c bcs.f03

std_mat.o: std_mat.f03 bcs.f03 metrics.f03 \
                      refinement.f03 utilities.f03
	$(FC) $(FFLAGS) -c std_mat.f03              

alloc.o: refinement.f03 bcs.f03 std_mat.f03
	$(FC) $(FFLAGS) -c alloc.f03

postpr.o: postpr.F03 std_mat.f03 metrics.f03
	$(FC) $(FFLAGS) $(INCL) -c postpr.F03

solgw.o: solgw.f90
	$(FC) $(FFLAGS) -c solgw.f90
   
clean:
	rm -rf *.o *.e *.mod

partclean:
	rm *.mod *.o 


