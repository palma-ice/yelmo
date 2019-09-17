# Yelmo

This is the code repository for the ice sheet model Yelmo.
While the model has been designed to be easy to use, there
are many parameters that require knowledge of ice-sheet
physics and numerous parameterizations. The (growing) model documentation
is provide to provide help with proper use of the model.

The model is described in the following article:
Robinson et al., in prep.

The model documentation can be found here:
https://palma-ice.github.io/yelmo-docs

## Dependencies

- [NetCDF](https://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html)
- [Library of Iterative Solvers for Linear Systems](http://www.ssisc.org/lis/)

See:
https://palma-ice.github.io/yelmo-docs/dependencies/

## Test cases

The published model description includes several test simulations for validation
of the model's performance. The following section describes how to perform these
tests using the same model version documented in the article. From this point,
it is assumed that the user has already configured the model for their system
(see https://palma-ice.github.io/yelmo-docs) and is ready to compile the mode.

1. EISMINT1 moving margin experiment
To perform the moving margin experiment, compile the benchmarks
executable and call it with the EISMINT parameter file:
```
make benchmarks
python run_yelmo.py -r -e benchmarks output/test-moving par/gmd/yelmo_EISMINT_moving.nml
```

2. EISMINT2 EXPA
To perform Experiment A from the EISMINT2 benchmarks, compile the benchmarks
executable and call it with the EXPA parameter file:
```
make benchmarks
python run_yelmo.py -r -e benchmarks output/test-expa par/gmd/yelmo_EISMINT_expa.nml
```

3. EISMINT2 EXPF
To perform Experiment F from the EISMINT2 benchmarks, compile the benchmarks
executable and call it with the EXPF parameter file:
```
make benchmarks
python run_yelmo.py -r -e benchmarks output/test-expf par/gmd/yelmo_EISMINT_expf.nml
```

4. MISMIP RF
To perform the MISMIP rate factor experiment, compile the mismip executable
and call it with the MISMIP parameter file:
```
make mismip
python run_yelmo.py -r -e mismip output/test-mismip par/gmd/yelmo_MISMIP3D.nml
```

5. Age profile experiments
To perform the age profile experiments, compile the Fortran program `tests/test_icetemp.f90`
and run it. To perform the different permutations, it is necessary to recompile for
single or double precision after changing the precision parameter `prec` in the file
`src/yelmo_defs.f90`. The number of vertical grid points can be specified in the main
program file, as well as the output filename.
```
make icetemp
./libyelmo/bin/test_icetemp.x
```

6. Antarctica present-day and glacial simulations
To perform the Antarctica simulations as presented in the paper, it is necessary
to compile the `initmip` executable and run with the present-day (pd) and
glacial (lgm) parameter files:
```
make initmip
python run_yelmo.py -r -e initmip output/test-ant-pd par/gmd/yelmo_Antarctica_pd.nml
python run_yelmo.py -r -e initmip output/test-ant-lgm par/gmd/yelmo_Antarctica_lgm.nml
```
