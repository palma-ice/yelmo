.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# PATH options
srcdir = src
objdir = libyelmo/include
bindir = libyelmo/bin
libdir = libs
testdir = tests

# Command-line options at make call
debug ?= 0

## COMPILER CONFIGURATION ##
# (should be loaded from config directory)

<COMPILER_CONFIGURATION>

# Determine whether to use normal flags or debugging flags
DFLAGS   = $(DFLAGS_NODEBUG)
ifeq ($(debug), 1)
	DFLAGS   = $(DFLAGS_DEBUG)
endif

# Debugging flags with profiling output enabled
ifeq ($(debug), 2)
	DFLAGS   = $(DFLAGS_PROFILE)
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

# Static library compilation
yelmo-static: $(yelmo_libs) $(yelmo_physics) $(yelmo_base) $(yelmo_tests)
	ar rc $(objdir)/libyelmo.a $(yelmo_libs) $(yelmo_physics) $(yelmo_base) $(yelmo_tests)
	ranlib $(objdir)/libyelmo.a
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

mismip : yelmo-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_LIS) -o $(bindir)/yelmo_mismip.x tests/yelmo_mismip.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS) $(objdir)/nml.o
		@echo " "
		@echo "    yelmo_mismip.x is ready."
		@echo " "

icetemp : yelmo-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_LIS) -o $(bindir)/test_icetemp.x tests/test_icetemp.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS) $(objdir)/nml.o
		@echo " "
		@echo "    test_icetemp.x is ready."
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
