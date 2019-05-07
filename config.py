'''
Script to generate Makefile with the correct compiler
configuration for a given machine.
'''

from subprocess import *
import argparse
import sys 

# Manage command-line arguments to load compiler option
parser = argparse.ArgumentParser()
parser.add_argument(
    'config', metavar='CONFIG', type=str,
     help='Name of config file to use (e.g., config/aire_gfortran)')
args = parser.parse_args()


# Determine locations
target_dir  = "./"
config_dir  = "config/"
config_path = args.config 

# Load template Makefile and insert compiler configuration
makefile0    = open(config_dir+"Makefile").read()
compile_info = open(config_path).read()
makefile1 = makefile0.replace("<COMPILER_CONFIGURATION>",compile_info)

# Write the new Makefile to the main directory
open(target_dir+"Makefile","w").write(makefile1)

print( "".join(["\n\nMakefile configuration complete for configuration file: ",config_path,"\n"]) )

instructions = '''==== How to run yelmo ====\n
# Make a link to the ice_data path, if available.
ln -s path/to/ice_data ice_data 

# Compile a test program
make clean 
make yelmo_benchmarks      # makes executable libyelmo/bin/yelmo_benchmarks.x 

# Run the program using run_yelmo.py (see run_yelmo.py -h for details) 
python run_yelmo.py -r -e benchmarks output/test par/yelmo_EISMINT.nml 

# Check the output 
cd output/test 
ncview yelmo2D.nc 
'''

print(instructions)



### Old script commands to automatically determine the config_path from host + compiler:

# Determine the current hostname
#proc = Popen("hostname",shell=True,stdout=PIPE,stderr=PIPE)
#host = proc.communicate()[0].strip().decode("utf-8")

#if host == "manager.cei.ucm.local": host = "eolo"
#if host in ["login01","login02"]:   host = "pik"

#print("{}{}_{}".format(config_dir,host,compiler))

