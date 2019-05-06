#!/bin/bash
#$ -V                             # Ensure user enivronment variables are available
#$ -cwd                           # To use the current directory
#$ -m ae                          # Send mail when job is aborted (a), begins (b) and ends (e)
#$ -N  yelmo                      # (nombre del trabajo)
#$ -o ./out.out                   # (fichero de salida)   $JOB_ID
#$ -e ./out.err                   # (fichero de error)

# Run the job
./libyelmo/bin/yelmo_benchmarks.x 

