#!/bin/bash
#SBATCH -A CSC113                                # charge from account CSC113
#SBATCH --job-name="Dbg"                         # set up job name
#SBATCH --output="bash/output/Dbg.%j.%N.out"     # set up the name of output
#SBATCH --error="bash/output/Dbg.%j.%N.err"      # set up the name of output
#SBATCH --partition=debug                        # assign jobs to the share node
#SBATCH --nodes=1                                # number of nodes 
#SBATCH --ntasks-per-node=24                     # number of tasks
#SBATCH --export=ALL                             # set up export
#SBATCH -t 00:30:00                              # set up maximum run time
#SBATCH --mail-user=michael.arciniaga@gmail.com  # set up email address
#SBATCH --mail-type=begin                        # email me when the job starts
#SBATCH --mail-type=end                          # email me when the job finishes

#This job runs with 2 nodes, 24 cores per node for a total of 48 cores.
#ibrun in verbose mode will give binding detail


# declare -a pmatrix

# [d Nk nd tau tp da eta J1 J2 mup u0]

pmatrix[0]="10 16 0.85 0.038 -0.2 0.0 0.17 -0.78 0.69 7 30 0.0001"


./params.out ${pmatrix[0]}
