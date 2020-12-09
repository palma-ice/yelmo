# Yelmo

This is the code repository for Yelmo, a 3D ice-sheet-shelf model solving
for the coupled dynamics and thermodynamics of the ice sheet system. Yelmo
can be used for idealized simulations, stand-alone ice sheet simulations
and fully coupled ice-sheet and climate simulations.
The physics and design of the model are described in the following article:
Robinson et al., in prep.

The (growing) model documentation is provided help with proper use of the model,
and can be found at [yelmo-docs](https://palma-ice.github.io/yelmo-docs).
While the model has been designed to be easy to use, there
are many parameters that require knowledge of ice-sheet
physics and numerous parameterizations. It is not recommended to use the ice
sheet model as a black box without understanding of the key parameters that
affect its performance.

To get started with comping and running the model, see the quick-start
instructions below in the section "Usage".

The test cases shown by Robinson et al. (in prep) can be reproduced following the
instructions below in the section "Test cases".

## Quickstart with Docker and VS Code

The easiest way to use Yelmo is with VS Code and Docker. It works on any plattform and uses a Linux based container. You don't need to know Docker or VS Code to get started. Just install the following:

1) [Docker](https://docs.docker.com/engine/install/)
2) [VS Code](https://code.visualstudio.com) 
3) [install the remote development extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.vscode-remote-extensionpack)
4) get the code (see below)

Then make sure that Docker is running and start VS Code. 
Open the folder with the Yelmo code. Say Yes, when VS Code asks you if you want to open it in the container.

Now you can directly go to step 3 below, just make sure that you use the terminal in VS Code.

This will give you everything you need to develop and run Yelmo.

## Dependencies

- NetCDF (tested with version 4.0 and higher): [NetCDF library](https://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html)
- LIS: [Library of Iterative Solvers for Linear Systems](http://www.ssisc.org/lis/)

See: [Dependencies](https://palma-ice.github.io/yelmo-docs/dependencies/)

OPTIONAL:
- Python 3.x, which is only needed for automatic configuration of the Makefile
and the use of the script `run_yelmo.py` for job preparation and submission.

## Directory structure
```fortran
    config/
        Configuration files for compilation on different systems.
    input/
        Location of any input data needed by the model.
    libs/
        Auxiliary libraries nesecessary for running the model.
    libyelmo/
        Folder containing all compiled files in a standard way with
        lib/, include/ and bin/ folders.
    output/
        Default location for model output.
    par/
        Default parameter files that manage the model configuration.
    src/
        Source code for Yelmo.
    tests/
        Source code and analysis scripts for specific model benchmarks and tests.
```

## Usage

Follow the steps below to (1) obtain the code, (2) configure the Makefile for your system,
(3) compile the Yelmo static library and an executable program and (4) run a test simulation.

### 1. Get the code.

Clone the repository from [https://github.com/palma-ice/yelmo](https://github.com/palma-ice/yelmo):
```
git clone git@github.com:palma-ice/yelmo.git $YELMOROOT
cd $YELMOROOT
```
where `$YELMOROOT` is the installation directory.

If you plan to make changes to the code, it is wise to check out a new branch:
```
git checkout -b user-dev
```
You should now be working on the branch `user-dev`.

### 2. Create the system-specific Makefile.

To compile Yelmo, you need to generate a Makefile that is appropriate for your system. In the folder `config`, you need to specify a configuration file that defines the compiler and flags, including definition of the paths to the `NetCDF` and `LIS` libraries. You can use another file in the config folder as a template, e.g.,

```
cd config
cp pik_ifort myhost_mycompiler
```
then modify the file `myhost_mycompiler` to match your paths. Back in `$YELMOROOT`, you can then generate your Makefile with the provided python configuration script:

```
cd $YELMOROOT
python config.py config/myhost_mycompiler
```
The result should be a Makefile in `$YELMOROOT` that is ready for use.

### 3. Compile the code.

Now you are ready to compile Yelmo as a static library:
```
make clean    # This step is very important to avoid errors!!
make yelmo-static [debug=1]
```
This will compile all of the Yelmo modules and libraries (as defined in `config/Makefile_yelmo.mk`),
and link them in a static library. All compiled files can be found in the folder `libyelmo/`.

Once the static library has been compiled, it can be used inside of external Fortran programs and modules
via the statement `use yelmo`.
To include/link yelmo-static during compilation of another program, its location must be defined:
```
INC_YELMO = -I${YELMOROOT}/include
LIB_YELMO = -L${YELMOROOT}/include -lyelmo
```

Alternatively, several test programs exist in the folder `tests/` to run Yelmo
as a stand-alone ice sheet.
For example, it's possible to run different EISMINT benchmarks, MISMIP benchmarks and the
ISIMIP6 INITMIP simulation for Greenland, respectively:
```
make benchmarks    # compiles the program `libyelmo/bin/yelmo_benchmarks.x`
make mismip        # compiles the program `libyelmo/bin/yelmo_mismip.x`
make initmip       # compiles the program `libyelmo/bin/yelmo_initmip.x`
```

The Makefile additionally allows you to specify debugging compiler flags with the option `debug=1`, in case you need to debug the code (e.g., `make benchmarks debug=1`). Using this option, the code will run much slower, so this option is not recommended unless necessary.

### 4. Run the model.

Once an executable has been created, you can run the model. This can be
achieved via the included Python job submission script `run_yelmo.py`. The following steps
are carried out via the script:

1. The output directory is created.
2. The executable is copied to the output directory
3. The relevant parameter files are copied to the output directory.
4. Links to the input data paths (`input` and `ice_data`) are created in the output directory. (Note `ice_data` is typically linked to an external data repository, and is not necessary.)
4. The executable is run from the output directory, either as a background process or it is submitted to the queue (the script currently supports `qsubmit` and `sbatch` commands).

To run a benchmark simulation, for example, use the following command:

```
python run_yelmo.py -r -e benchmarks output/test par/yelmo_EISMINT.nml
```
where the option `-r` implies that the model should be run as a background process. If this is omitted, then the output directory will be populated, but no executable will be run, while `-s` instead will submit the simulation to cluster queue system instead of running in the background. The option `-e` lets you specify the executable. For some standard cases, shortcuts have been created:
```
benchmarks = libyelmo/bin/yelmo_benchmarks.x
mismip     = libyelmo/bin/yemo_mismip.x
initmip    = libyelmo/bin/yelmo_initmip.x
```
The last two mandatory arguments are always the output/run directory and the parameter file to be used for this simulation. In the case of the above simulation, the output directory is defined as `output/test`, where all model parameters (loaded from the file `par/yelmo_EISMINT.nml`) and model output can be found.

## Test cases

The published model description includes several test simulations for validation
of the model's performance. The following section describes how to perform these
tests using the same model version documented in the article. From this point,
it is assumed that the user has already configured the model for their system
(see https://palma-ice.github.io/yelmo-docs) and is ready to compile the mode.

### 1. EISMINT1 moving margin experiment
To perform the moving margin experiment, compile the benchmarks
executable and call it with the EISMINT parameter file:
```
make benchmarks
python run_yelmo.py -r -e benchmarks output/test-moving par-gmd/yelmo_EISMINT_moving.nml
```

### 2. EISMINT2 EXPA
To perform Experiment A from the EISMINT2 benchmarks, compile the benchmarks
executable and call it with the EXPA parameter file:
```
make benchmarks
python run_yelmo.py -r -e benchmarks output/test-expa par-gmd/yelmo_EISMINT_expa.nml
```

### 3. EISMINT2 EXPF
To perform Experiment F from the EISMINT2 benchmarks, compile the benchmarks
executable and call it with the EXPF parameter file:
```
make benchmarks
python run_yelmo.py -r -e benchmarks output/test-expf par-gmd/yelmo_EISMINT_expf.nml
```

### 4. MISMIP RF
To perform the MISMIP rate factor experiment, compile the mismip executable
and call it with the MISMIP parameter file:
```
make mismip
python run_yelmo.py -r -e mismip output/test-mismip par-gmd/yelmo_MISMIP3D.nml
```
To perform the different permutations, it is necessary to change the resolution
`&mismip::dx` in the parameter file (reducing the timestep `&mismip::dtt` for higher resolutions).
By default, no subgrid grounding line treatment or scaling is applied. To add subgrid treatment,
set the following parameters:
```
&ydyn
    beta_gl_sep         = -1
    beta_gl_stag        =  3
/
```
Meanwhile, to apply linear scaling to basal friction as it approaches the grounding line,
set `&ydyn::beta_gl_scale = 2`.

### 5. Age profile experiments
To perform the age profile experiments, compile the Fortran program `tests/test_icetemp.f90`
and run it:
```
make icetemp
./libyelmo/bin/test_icetemp.x
```
To perform the different permutations, it is necessary to recompile for
single or double precision after changing the precision parameter `prec` in the file
`src/yelmo_defs.f90`. The number of vertical grid points can be specified in the main
program file, as well as the output filename.

### 6. Antarctica present-day and glacial simulations
To perform the Antarctica simulations as presented in the paper, it is necessary
to compile the `initmip` executable and run with the present-day (pd) and
glacial (lgm) parameter values:
```
make initmip
# In par-gmd/yelmo_Antarctica.nml, set control.clim_nm="clim_pd"
python run_yelmo.py -r -e initmip output/test-ant-pd par-gmd/yelmo_Antarctica.nml
# In par-gmd/yelmo_Antarctica.nml, set control.clim_nm="clim_lgm"
python run_yelmo.py -r -e initmip output/test-ant-lgm par-gmd/yelmo_Antarctica.nml
```
