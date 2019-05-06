#! /bin/bash
#SBATCH --qos=short
#SBATCH --job-name=yelmo
#SBATCH --account=anthroia
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-user=robinson@pik-potsdam.de
#SBATCH --output=out.out
#SBATCH --error=out.err
#SBATCH --time=20:00:00
#SBATCH --mem=50000 

# Run the job
./libyelmo/bin/yelmo_benchmarks.x 

##SBATCH --partition=ram_gpu
##SBATCH --mem=50000 
