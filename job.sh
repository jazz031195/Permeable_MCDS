#!/bin/bash
#SBATCH --mail-type=END
#SBATCH --mail-user=jasmine.nguyen-duc@chuv.ch
#SBATCH --job-name sim_job
#SBATCH --account rad
#SBATCH --partition cluster
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=7G
#SBATCH --time=47:59:00
#SBATCH -o /OUT/output/%N.%j.%a.out
#SBATCH -e /OUT/err/%N.%j.%a.err

./run.sh