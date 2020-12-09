#############################################################
##							
## Rules for individual libraries or modules
##
#############################################################

## EXTERNAL LIBRARIES #######################################

$(objdir)/ncio.o: $(libdir)/ncio.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_NC) -c -o $@ $<

$(objdir)/nml.o: $(libdir)/nml.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/gaussian_filter.o: $(libdir)/gaussian_filter.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/climate_adjustments.o: $(libdir)/climate_adjustments.f90 $(objdir)/yelmo_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/ice_enhancement.o: $(libdir)/ice_enhancement.f90 $(objdir)/yelmo_defs.o $(objdir)/nml.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/ice_optimization.o: $(libdir)/ice_optimization.f90 $(objdir)/yelmo_defs.o $(objdir)/gaussian_filter.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/interp1D.o: $(libdir)/interp1D.f90 $(objdir)/yelmo_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

## INTERNAL PHYSICS LIBRARIES ###############################

$(objdir)/basal_dragging.o: $(srcdir)/physics/basal_dragging.f90 $(objdir)/yelmo_defs.o \
							$(objdir)/yelmo_tools.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/grounding_line_flux.o: $(srcdir)/physics/grounding_line_flux.f90 $(objdir)/yelmo_defs.o \
							$(objdir)/yelmo_tools.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/calving.o: $(srcdir)/physics/calving.f90 $(objdir)/yelmo_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/deformation.o: $(srcdir)/physics/deformation.f90 $(objdir)/yelmo_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/thermodynamics.o : $(srcdir)/physics/thermodynamics.f90 $(objdir)/yelmo_defs.o $(objdir)/gaussian_filter.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/ice_tracer.o : $(srcdir)/physics/ice_tracer.f90 $(objdir)/yelmo_defs.o $(objdir)/solver_tridiagonal.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/ice_enthalpy.o : $(srcdir)/physics/ice_enthalpy.f90 $(objdir)/yelmo_defs.o \
					  $(objdir)/solver_tridiagonal.o $(objdir)/thermodynamics.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/mass_conservation.o : $(srcdir)/physics/mass_conservation.f90 $(objdir)/yelmo_defs.o \
								$(objdir)/solver_advection.o $(objdir)/solver_advection_sico.o \
								$(objdir)/solver_advection_new.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/solver_ssa_sico5.o: $(srcdir)/physics/solver_ssa_sico5.F90 $(objdir)/yelmo_defs.o \
							$(objdir)/yelmo_tools.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_LIS) -c -o $@ $<

$(objdir)/solver_tridiagonal.o: $(srcdir)/physics/solver_tridiagonal.f90 $(objdir)/yelmo_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/solver_advection.o: $(srcdir)/physics/solver_advection.f90 $(objdir)/yelmo_defs.o \
								$(objdir)/yelmo_tools.o $(objdir)/solver_advection_sico.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/solver_advection_sico.o : $(srcdir)/physics/solver_advection_sico.F90 $(objdir)/yelmo_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_LIS) -c -o $@ $<

$(objdir)/solver_advection_new.o: $(srcdir)/physics/solver_advection_new.f90 $(objdir)/yelmo_defs.o \
								$(objdir)/yelmo_tools.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/topography.o: $(srcdir)/physics/topography.f90 $(objdir)/yelmo_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/velocity_general.o: $(srcdir)/physics/velocity_general.f90 $(objdir)/yelmo_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/velocity_sia.o: $(srcdir)/physics/velocity_sia.f90 \
						  	$(objdir)/yelmo_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/velocity_hybrid.o: $(srcdir)/physics/velocity_hybrid.f90 \
						  	$(objdir)/yelmo_defs.o $(objdir)/yelmo_tools.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/velocity_diva.o: $(srcdir)/physics/velocity_diva.f90 \
						  	$(objdir)/yelmo_defs.o $(objdir)/yelmo_tools.o $(objdir)/basal_dragging.o \
						  	$(objdir)/solver_ssa_sico5.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/velocity_l1l2.o: $(srcdir)/physics/velocity_l1l2.f90 \
						  	$(objdir)/yelmo_defs.o $(objdir)/yelmo_tools.o $(objdir)/basal_dragging.o \
						  	$(objdir)/solver_ssa_sico5.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/velocity_hybrid_pd12.o: $(srcdir)/physics/velocity_hybrid_pd12.f90 \
						  	$(objdir)/yelmo_defs.o $(objdir)/yelmo_tools.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<



## YELMO BASE ###############################################

$(objdir)/yelmo_defs.o: $(srcdir)/yelmo_defs.f90 $(objdir)/nml.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_LIS) -c -o $@ $<

$(objdir)/yelmo_grid.o: $(srcdir)/yelmo_grid.f90 $(objdir)/yelmo_defs.o $(objdir)/nml.o $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/yelmo_regridding.o : $(srcdir)/yelmo_regridding.f90 $(objdir)/yelmo_defs.o $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_LIS) -c -o $@ $<

$(objdir)/yelmo_tools.o: $(srcdir)/yelmo_tools.f90 $(objdir)/yelmo_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/yelmo_timesteps.o : $(srcdir)/yelmo_timesteps.f90 $(objdir)/yelmo_defs.o $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_LIS) -c -o $@ $<

$(objdir)/yelmo_topography.o: $(srcdir)/yelmo_topography.f90 $(objdir)/yelmo_defs.o \
							  $(objdir)/yelmo_grid.o $(objdir)/yelmo_tools.o  \
 							  $(objdir)/mass_conservation.o $(objdir)/calving.o \
 							  $(objdir)/topography.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/yelmo_dynamics.o: $(srcdir)/yelmo_dynamics.f90 $(objdir)/yelmo_defs.o \
							$(objdir)/velocity_general.o \
							$(objdir)/velocity_sia.o \
							$(objdir)/velocity_hybrid.o \
							$(objdir)/velocity_diva.o \
							$(objdir)/velocity_l1l2.o \
							$(objdir)/velocity_hybrid_pd12.o \
							$(objdir)/solver_ssa_sico5.o \
							$(objdir)/basal_dragging.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/yelmo_material.o: $(srcdir)/yelmo_material.f90 $(objdir)/yelmo_defs.o $(objdir)/deformation.o \
							$(objdir)/ice_tracer.o 
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/yelmo_thermodynamics.o: $(srcdir)/yelmo_thermodynamics.f90 $(objdir)/yelmo_defs.o \
								  $(objdir)/ice_enthalpy.o $(objdir)/thermodynamics.o \
								  $(objdir)/solver_advection.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/yelmo_boundaries.o: $(srcdir)/yelmo_boundaries.f90 $(objdir)/yelmo_defs.o $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/yelmo_data.o: $(srcdir)/yelmo_data.f90 $(objdir)/yelmo_defs.o $(objdir)/nml.o $(objdir)/ncio.o \
						$(objdir)/topography.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/yelmo_regions.o: $(srcdir)/yelmo_regions.f90 $(objdir)/ncio.o $(objdir)/yelmo_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/yelmo_io.o: $(srcdir)/yelmo_io.f90 $(objdir)/ncio.o $(objdir)/yelmo_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/control.o: control.f90 $(objdir)/nml.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/yelmo_ice.o: $(srcdir)/yelmo_ice.f90 $(objdir)/yelmo_defs.o  \
				   	   $(objdir)/yelmo_timesteps.o \
				   	   $(objdir)/yelmo_topography.o \
	                   $(objdir)/yelmo_dynamics.o \
	                   $(objdir)/velocity_hybrid_pd12.o \
	                   $(objdir)/velocity_sia.o \
	                   $(objdir)/yelmo_material.o \
	                   $(objdir)/yelmo_thermodynamics.o \
	                   $(objdir)/yelmo_boundaries.o \
	                   $(objdir)/yelmo_data.o \
	                   $(objdir)/yelmo_regions.o \
	                   $(objdir)/yelmo_io.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/yelmo.o: $(srcdir)/yelmo.f90 $(objdir)/yelmo_ice.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

## YELMO TESTS ###############################################

$(objdir)/ice_benchmarks.o: $(testdir)/ice_benchmarks.f90 $(objdir)/yelmo_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/mismip3D.o: $(testdir)/mismip3D.f90 $(objdir)/ncio.o $(objdir)/yelmo_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

#############################################################
##							
## List of yelmo files
##
#############################################################

yelmo_libs = 		   $(objdir)/gaussian_filter.o \
					   $(objdir)/climate_adjustments.o \
					   $(objdir)/ice_enhancement.o \
					   $(objdir)/ice_optimization.o \
					   $(objdir)/interp1D.o \
					   $(objdir)/nml.o \
			 		   $(objdir)/ncio.o

yelmo_physics =  	   $(objdir)/basal_dragging.o \
					   $(objdir)/grounding_line_flux.o \
					   $(objdir)/calving.o \
					   $(objdir)/deformation.o \
					   $(objdir)/thermodynamics.o \
					   $(objdir)/ice_tracer.o \
					   $(objdir)/ice_enthalpy.o \
					   $(objdir)/mass_conservation.o \
					   $(objdir)/solver_ssa_sico5.o \
					   $(objdir)/solver_tridiagonal.o \
					   $(objdir)/solver_advection.o \
					   $(objdir)/solver_advection_sico.o \
					   $(objdir)/solver_advection_new.o \
					   $(objdir)/topography.o \
					   $(objdir)/velocity_general.o \
					   $(objdir)/velocity_sia.o \
					   $(objdir)/velocity_hybrid.o \
					   $(objdir)/velocity_diva.o \
					   $(objdir)/velocity_l1l2.o \
					   $(objdir)/velocity_hybrid_pd12.o
					   
					   

yelmo_base = 		   $(objdir)/yelmo_defs.o \
					   $(objdir)/yelmo_grid.o \
					   $(objdir)/yelmo_regridding.o \
	                   $(objdir)/yelmo_tools.o \
			 		   $(objdir)/yelmo_timesteps.o \
					   $(objdir)/yelmo_ice.o \
			 	       $(objdir)/yelmo_topography.o \
	         		   $(objdir)/yelmo_dynamics.o \
	         		   $(objdir)/yelmo_material.o \
	         		   $(objdir)/yelmo_thermodynamics.o \
	         		   $(objdir)/yelmo_boundaries.o \
	                   $(objdir)/yelmo_data.o \
	                   $(objdir)/yelmo_regions.o \
	                   $(objdir)/yelmo_io.o \
	         		   $(objdir)/yelmo.o

yelmo_tests = 		   $(objdir)/ice_benchmarks.o \
					   $(objdir)/mismip3D.o

# Extras for testing 
yelmo_thermo =         $(objdir)/yelmo_defs.o \
					   $(objdir)/yelmo_grid.o \
	                   $(objdir)/yelmo_tools.o \
	                   $(objdir)/interp1D.o \
			 		   $(objdir)/thermodynamics.o \
					   $(objdir)/ice_enthalpy.o

