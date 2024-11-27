.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# PATH options
srcdir = src
objdir = libyelmo/include
bindir = libyelmo/bin
libdir = libs
testdir = tests

# Command-line options at make call
debug  ?= 0
openmp ?= 0 
petsc  ?= 0 

## COMPILER CONFIGURATION ##
# (should be loaded from config directory)

<COMPILER_CONFIGURATION>

# Determine whether to use normal flags or debugging flags
DFLAGS = $(DFLAGS_NODEBUG)
ifeq ($(debug), 1)
	DFLAGS = $(DFLAGS_DEBUG)
endif

# Debugging flags with profiling output enabled
ifeq ($(debug), 2)
	DFLAGS = $(DFLAGS_PROFILE)
endif

###############################################
##							
## List of yelmo rules and source files
##
###############################################

include config/Makefile_yelmo.mk

###############################################
##
## Compilation of complete programs
##
###############################################

check-hash:
	./check_githash > git_yelmo.txt

# Static library compilation
yelmo-static: $(yelmo_libs) $(yelmo_physics) $(yelmo_base) $(yelmo_tests)
	ar rc $(objdir)/libyelmo.a $(yelmo_libs) $(yelmo_physics) $(yelmo_base) $(yelmo_tests)
	ranlib $(objdir)/libyelmo.a
	./check_githash > git_yelmo.txt
	@echo " "
	@echo "    $(objdir)/libyelmo.a is ready."
	@echo " "

# Test programs that use yelmo-static
initmip : yelmo-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_LIS) -o $(bindir)/yelmo_initmip.x tests/yelmo_initmip.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS) $(objdir)/nml.o
		@echo " "
		@echo "    yelmo_initmip.x is ready."
		@echo " "

benchmarks : yelmo-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_LIS) -o $(bindir)/yelmo_benchmarks.x tests/yelmo_benchmarks.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS) $(objdir)/nml.o
		@echo " "
		@echo "    yelmo_benchmarks.x is ready."
		@echo " "

calving : yelmo-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_LIS) -o $(bindir)/yelmo_calving.x tests/yelmo_calving.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS) $(objdir)/nml.o
		@echo " "
		@echo "    yelmo_calving.x is ready."
		@echo " "

mismip : yelmo-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_LIS) -o $(bindir)/yelmo_mismip.x tests/yelmo_mismip.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS) $(objdir)/nml.o
		@echo " "
		@echo "    yelmo_mismip.x is ready."
		@echo " "

ismiphom : yelmo-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_LIS) -o $(bindir)/yelmo_ismiphom.x tests/yelmo_ismiphom.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS) $(objdir)/nml.o
		@echo " "
		@echo "    yelmo_ismiphom.x is ready."
		@echo " "

regridding : yelmo-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_LIS) -o $(bindir)/yelmo_ismiphom.x tests/yelmo_test_regridding.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS) $(objdir)/nml.o
		@echo " "
		@echo "    yelmo_test_regridding.x is ready."
		@echo " "

icetemp : yelmo-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_LIS) -o $(bindir)/test_icetemp.x tests/test_icetemp.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS) $(objdir)/nml.o
		@echo " "
		@echo "    test_icetemp.x is ready."
		@echo " "

roots : $(objdir)/root_finder.o tests/test_roots.f90
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_roots.x tests/test_roots.f90 \
			$(objdir)/root_finder.o $(LFLAGS) $(objdir)/nml.o
		@echo " "
		@echo "    test_roots.x is ready."
		@echo " "

levelset : yelmo-static tests/test_levelset.f90
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_levelset.x tests/test_levelset.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS) $(objdir)/nml.o
		@echo " "
		@echo "    test_levelset.x is ready."
		@echo " "

mismip1 : yelmo-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_LIS) -o $(bindir)/yelmo_mismip.x tests/yelmo_mismip_new.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS) $(objdir)/nml.o
		@echo " "
		@echo "    yelmo_mismip.x is ready."
		@echo " "

opt : yelmo-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_LIS) -o $(bindir)/yelmo_opt.x tests/yelmo_opt.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS) $(objdir)/nml.o
		@echo " "
		@echo "    yelmo_opt.x is ready."
		@echo " "

trough : yelmo-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_LIS) -o $(bindir)/yelmo_trough.x tests/yelmo_trough.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS) $(objdir)/nml.o
		@echo " "
		@echo "    yelmo_trough.x is ready."
		@echo " "

slab : yelmo-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_LIS) -o $(bindir)/yelmo_slab.x tests/yelmo_slab.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS) $(objdir)/nml.o
		@echo " "
		@echo "    yelmo_slab.x is ready."
		@echo " "

check_sim : $(objdir)/ncio.o
		$(FC) $(DFLAGS) $(FFLAGS) -o ./check_sim.x check_sim.f90 \
			$(LFLAGS) $(objdir)/ncio.o
		@echo " "
		@echo "    check_sim.x is ready."
		@echo " "

.PHONY : usage
usage:
	@echo ""
	@echo "    * USAGE * "
	@echo ""
	@echo " make initmip    : compiles yelmo_initmip.x, for running yelmo on a given domain defined in param file."
	@echo " make benchmarks : compiles yelmo_benchmarks.x, for running various EISMINT domain experiments.""
	@echo " make mismip     : compiles yelmo_mismip.x, for running various MISMIP domain experiments.""
	@echo " make opt        : compiles yelmo_opt.x, for running various MISMIP domain experiments.""
	@echo " make trough     : compiles yelmo_trough.x, for running trough domain experiments.""
	@echo " make clean      : cleans object files"
	@echo ""

clean:
	rm -f $(bindir)/*.x
	rm -f  *.x gmon.out $(objdir)/*.o $(objdir)/*.mod $(objdir)/*.a $(objdir)/*.so
	rm -rf *.x.dSYM
	rm -f ./git_yelmo.txt


# Experimental
yelmo-shared: $(yelmo_libs) $(yelmo_physics) $(yelmo_base) $(yelmo_tests)
	$(FC) -shared -o $(objdir)/libyelmo.so -fPIC $(yelmo_libs) $(yelmo_physics) $(yelmo_base) $(yelmo_tests) \
		$(LFLAGS) $(INC_LIS) 
	./check_githash > git_yelmo.txt
	@echo " "
	@echo "    $(objdir)/libyelmo.so is ready."
	@echo " "
