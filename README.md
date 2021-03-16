# Yelmo

Yelmo is a 3D ice-sheet-shelf model solving
for the coupled dynamics and thermodynamics of the ice sheet system. Yelmo
can be used for idealized simulations, stand-alone ice sheet simulations
and fully coupled ice-sheet and climate simulations.

The physics and design of the model are described in the following article:

> Robinson, A., Alvarez-Solas, J., Montoya, M., Goelzer, H., Greve, R., and Ritz, C.: Description and validation of the ice-sheet model Yelmo (version 1.0), Geosci. Model Dev., 13, 2805–2823, [https://doi.org/10.5194/gmd-13-2805-2020](https://doi.org/10.5194/gmd-13-2805-2020), 2020.

The (growing) model documentation is provided help with proper use of the model,
and can be found at:

 [https://palma-ice.github.io/yelmo-docs](https://palma-ice.github.io/yelmo-docs)
 
While the model has been designed to be easy to use, there
are many parameters that require knowledge of ice-sheet
physics and numerous parameterizations. It is not recommended to use the ice
sheet model as a black box without understanding of the key parameters that
affect its performance.

To get started with comping and running the model, see the quick-start
instructions below in the section "Usage". Or go to the documentation directly: [https://palma-ice.github.io/yelmo-docs/getting-started/](https://palma-ice.github.io/yelmo-docs/getting-started/).

The test cases shown by Robinson et al. (2020) can be reproduced following the
instructions below in the section "Test cases".

The sections below can also be found in the documentation here: [Getting started](https://palma-ice.github.io/yelmo-docs/getting-started/).

# Getting started

Here you can find the basic information and steps needed to get **Yelmo** running.

## Super-quick start

A summary of commands to get started is given below. For more detailed information see subsequent sections.

```
# Clone repository
git clone git@github.com:palma-ice/yelmo.git

# Enter directory and run configuration script
cd yelmo
python config.py config/pik_ifort 

# Compile the benchmarks program
make clean 
make benchmarks 

# Run a test simulation of the EISMINT1-moving experiment
./runylmo -r -e benchmarks -o output/eismint1-moving -n par-gmd/yelmo_EISMINT-moving.nml

# Compile the initmip program and run a simulation of Antarctica
make initmip 
./runylmo -r -e initmip -o output/ant-pd -n par/yelmo_initmip.nml -p ctrl.clim_nm="clim_pd"
```

## Dependencies

- NetCDF library (preferably version 4.0 or higher)
- LIS: [Library of Iterative Solvers for Linear Systems](http://www.ssisc.org/lis/)

See: [Dependencies](https://palma-ice.github.io/yelmo-docs/dependencies/)

OPTIONAL:

- Python 3.x, which is only needed for automatic configuration of the Makefile and the use of the script `runylmo` for job preparation and submission.
- Python library `runner` for changing parameters at the command line using `runylmo`, and for running ensembles. Installation instructions here [https://github.com/alex-robinson/runner](https://github.com/alex-robinson/runner)

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

#### Alternative - quickstart with Docker and VS Code

Instead of a manual install, one way to get up and running quickly with Yelmo is with VS Code and Docker. It works on any plattform and uses a Linux based container. You don't need to know Docker or VS Code to get started. Just install the following:

1) [Docker](https://docs.docker.com/engine/install/)
2) [VS Code](https://code.visualstudio.com) 
3) [install the remote development extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.vscode-remote-extensionpack)
4) get the code (see below)

Then make sure that Docker is running and start VS Code. 
Open the folder with the Yelmo code. Say Yes, when VS Code asks you if you want to open it in the container.

Now you can directly go to step 3 below, just make sure that you use the terminal in VS Code.

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
achieved via the included Python job submission script `runylmo`. The following steps
are carried out via the script:

1. The output directory is created.
2. The executable is copied to the output directory
3. The relevant parameter files are copied to the output directory.
4. Links to the input data paths (`input` and `ice_data`) are created in the output directory. Note that many simulations, such as benchmark experiments, do not depend on these external data sources, but the links are made anyway.
4. The executable is run from the output directory, either as a background process or it is submitted to the queue via `sbatch` (the SLURM workload manager).

To run a benchmark simulation, for example, use the following command:

```
./runylmo -r -e benchmarks -o output/test -n par/yelmo_EISMINT.nml
```

where the option `-r` implies that the model should be run as a background process. If this is omitted, then the output directory will be populated, but no executable will be run, while `-s` instead will submit the simulation to cluster queue system instead of running in the background. The option `-e` lets you specify the executable. For some standard cases, shortcuts have been created:

```
benchmarks = libyelmo/bin/yelmo_benchmarks.x
mismip     = libyelmo/bin/yemo_mismip.x
initmip    = libyelmo/bin/yelmo_initmip.x
```
The last two mandatory arguments `-o OUTDIR` and `-n PAR_PATH` are the output/run directory and the parameter file to be used for this simulation, respectively. In the case of the above simulation, the output directory is defined as `output/test`, where all model parameters (loaded from the file `par/yelmo_EISMINT.nml`) and model output can be found.

It is also possible to modify parameters inline via the option `-p KEY=VAL [KEY=VAL ...]`. The parameter should be specified with its namelist group and its name. E.g., to change the resolution of the EISMINT benchmark experiment to 10km, use:

```
./runylmo -r -e benchmarks -o output/test -n par/yelmo_EISMINT.nml -p ctrl.dx=10
```

See `runylmo -h` for more details on the run script. 

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
./runylmo -r -e benchmarks -o output/eismint-moving -n par-gmd/yelmo_EISMINT_moving.nml
```

### 2. EISMINT2 EXPA
To perform Experiment A from the EISMINT2 benchmarks, compile the benchmarks
executable and call it with the EXPA parameter file:

```
make benchmarks
./runylmo -r -e benchmarks -o output/eismint-expa -n par-gmd/yelmo_EISMINT_expa.nml
```

### 3. EISMINT2 EXPF
To perform Experiment F from the EISMINT2 benchmarks, compile the benchmarks
executable and call it with the EXPF parameter file:

```
make benchmarks
./runylmo -r -e benchmarks -o output/eismint-expf -n par-gmd/yelmo_EISMINT_expf.nml
```

### 4. MISMIP RF
To perform the MISMIP rate factor experiment, compile the mismip executable
and call it with the MISMIP parameter file the three parameter permutations of interest (default, subgrid and subgrid+gl-scaling):

```
make mismip
./runylmo -r -e mismip -o output/mismip-rf-0 -n par-gmd/yelmo_MISMIP3D.nml -p ydyn.beta_gl_stag=0 ydyn.beta_gl_scale=0
./runylmo -r -e mismip -o output/mismip-rf-1 -n par-gmd/yelmo_MISMIP3D.nml -p ydyn.beta_gl_stag=3 ydyn.beta_gl_scale=0
./runylmo -r -e mismip -o output/mismip-rf-2 -n par-gmd/yelmo_MISMIP3D.nml -p ydyn.beta_gl_stag=3 ydyn.beta_gl_scale=2
```
To additionally change the resolution of the simulations change the parameter `mismip.dx`, e.g. for the default simulation with 10km resolution , call:

```
./runylmo -r -e mismip -o output/mismip-rf-0-10km -n par-gmd/yelmo_MISMIP3D.nml -p ydyn.beta_gl_stag=0 ydyn.beta_gl_scale=0 mismip.dx=10
```

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
./runylmo -r -e initmip -o output/ant-pd -n par-gmd/yelmo_Antarctica.nml -p ctrl.clim_nm="clim_pd"
./runylmo -r -e initmip -o output/ant-lgm -n par-gmd/yelmo_Antarctica.nml -p ctrl.clim_nm="clim_lgm"
```
