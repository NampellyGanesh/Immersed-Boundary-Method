all:
	 @ifort -r8 -traceback -c global_variables.f90
	 @ifort -r8 -traceback -c geometric_pros.f90                                                                                    	
	 @ifort -r8 -traceback -c IBM.f90                                                                                               	
	 @ifort -r8 -traceback -c flux_calculation_viscous.f90                                                                          	
	 @ifort -r8 -traceback -check all global_variables.o geometric_pros.o IBM.o flux_calculation_viscous.o main_program.f90 -o b.out	

# @gfortran -fdefault-real-8 -fbacktrace -c global_variables.f90
# @gfortran -fdefault-real-8 -fbacktrace -c geometric_pros.f90
# @gfortran -fdefault-real-8 -fbacktrace -c IBM.f90
# @gfortran -fdefault-real-8 -fbacktrace -c flux_calculation_viscous.f90
# @gfortran -fdefault-real-8 -fbacktrace -fbounds-check -fcheck=all global_variables.o geometric_pros.o IBM.o flux_calculation_viscous.o main_program.f90 -o b.out



	
	
clean:
	@rm *.o  *.mod b.out
