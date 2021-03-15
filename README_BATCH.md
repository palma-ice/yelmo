# Running ensemble simulations

The script `runylmo` is designed to prepare and run one
Yelmo simulation with parameters loaded from a specified 
parameter file, and potentially updated on the command-line via `-p`. To perform an ensemble of Yelmo simulations
with various combinations of parameter values, `runylmo` 
can be called via the [Python `runner` module](https://github.com/alex-robinson/runner).
See below for installation instructions.

## Simple ensemble

To run a simple ensemble, call `runylmo` via the `jobrun` command. For example,
to run an ensemble of benchmark simulations with different resolutions of 30.0 and 50.0 km,
you can run the following command:
```
jobrun ./runylmo -r -e benchmarks -n par-gmd/yelmo_HALFAR.nml -- -o output/run -p eismint.dx=30.0,50.0
```
Everything between `jobrun` and `--` corresponds to the normal call to `runylmo`, with one key difference. The output directory argument `-o output/run` is specified after the `--` as this is now the ensemble output directory. Individual simulation directories will appear inside this main ensemble directory.

Everything after the `--` corresponds to the options of the `jobrun` command, namely you need to specify the ensemble output directory via `-o` and the parameter combinations to run via `-p KEY=VAL [KEY=VAL...]`. In this case,
`-o output/run` sets the ensemble output directory location
 (individual simulations will appear as numbered folder `0,1,2,...` inside of `output/run`).
`-p` is the option to specify parameter options to modify, in this case setting the parameter values of `eismint.dx` (group.parameter) to `30.0` and `50.0`. 

So, the above command would make an ensemble of two simulations with different values of `eismint.dx` in the output location `output/run`. 

It's also possible to acheive the above in two steps, first sampling the parameter values and storing them
in a file, then loading that file to generate the ensemble:
```
job product eismint.dx=30.0,50.0 > ensemble1.txt 
jobrun ./runylmo -r -e benchmarks -n par-gmd/yelmo_HALFAR.nml -- -o output/run -i ensemble1.txt
```

## Parameter sampling 

Using the command `job sample`, it is possible to generate an ensemble drawing from different
parameter distributions instead:
```
job sample eismint.dx=U?30.0,50.0 --size 10 > ensemble2.txt
jobrun ./runylmo -r -e benchmarks -n par-gmd/yelmo_HALFAR.nml -- -o output/run -i ensemble2.txt
```
The above command samples 10 values of `eismint.dx` from a uniform distribution ranging from `30.0` to `50.0` following a Latin-Hypercube sampling method. It is also possible to sample from a normal distribution, e.g., 
with mean and standard deviation of (40.0,10.0):
```
job sample eismint.dx=N?40.0,10.0 --size 10 > ensemble3.txt
```

## Installing runner 

  1. Download or clone the runner package from the repository:
[https://github.com/alex-robinson/runner](https://github.com/alex-robinson/runner)

  2. From inside the main directory, install the package via `pip`:
```
cd runner
pip install ./
```

Once completed successfully, the system commands `job`, `job sample`, `job run` and `jobrun` should be available at the command line. Check via `job -h` for example.