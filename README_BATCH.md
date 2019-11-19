# Running ensemble simulations

The script run_yelmo.py is designed to prepare and run one
Yelmo simulation with parameter loaded from a specified 
parameter file. To perform an ensemble of Yelmo simulations
with various combinations of parameter values, run_yelmo.py 
can be called via the [Python `runner` module](https://github.com/perrette/runner).
See below for installation instructions.

## Simple ensemble

To run a simple ensemble, call run_yelmo.py via the `job run` command. For example,
to run an ensemble of benchmark simulations with different resolutions of 30.0 and 50.0 km,
you can run the following command:
```
job run --shell -f -o output/run -p eismint.dx=30.0,50.0 -- python run_yelmo.py -x -r -e benchmarks {} par/gmd/yelmo_HALFAR.nml
```
Everything after the `--` corresponds to the normal call to run_yelmo.py, with two key differences. It 
is necessary to add the option `-x` to signal to run_yelmo.py that it should use capabilities for loading
Namelist files from the `runner` module. And the first argument, which is normally the output directory,
is replaced by the placeholder `{}` - this will be populated with the actual output directory of each
simulation by the `job run` command. 

Everything before the `--` corresponds to the options of the `job run` command. `--shell` ensures that the output is printed to the screen as well as to a log file to see which simulations are submitted and how, `-f` means that any previous simulations will be overwritten in the same ensemble output folder, 
`-o output/run` sets the ensemble output directory location
 (individual simulations will appear as numbered folder `0,1,2,...` inside of `output/run`).
`-p` is the option to specify parameter options to modify, in this case setting the parameter values of `eismint.dx` (group.parameter) to `30.0` and `50.0`. 

So, the above command would make an ensemble of two simulations with different values of `eismint.dx` in the output location `output/run`. 

It's also possible to acheive the above in two steps, first sampling the parameter values and storing them
in a file, then loading that file to generate the ensemble:
```
job product eismint.dx=30.0,50.0 > ensemble1.txt 
job run --shell -f -o output/run -i ensemble1.txt -- python run_yelmo.py -x -r -e benchmarks {} par/gmd/yelmo_HALFAR.nml
```

## Parameter sampling 

Using the command `job sample`, it is possible to generate an ensemble drawing from different
parameter distributions instead:
```
job sample eismint.dx=U?30.0,50.0 --size 10 > ensemble2.txt
job run --shell -f -o output/run -i ensemble2.txt -- python run_yelmo.py -x -r -e benchmarks {} par/gmd/yelmo_HALFAR.nml
```
The above command samples 10 values of `eismint.dx` from a uniform distribution ranging from `30.0` to `50.0` following a Latin-Hypercube sampling method. It is also possible to sample from a normal distribution, e.g., 
with mean and standard deviation of (40.0,10.0):
```
job sample eismint.dx=N?40.0,10.0 --size 10 > ensemble3.txt
```

## Installing runner 

  1. Download or clone the runner package from the repository:
[https://github.com/perrette/runner](https://github.com/perrette/runner)

  2. Run `python setup.py install` to install the package. 

Once completed successfully, a system command `job` should be available
at the command line. 