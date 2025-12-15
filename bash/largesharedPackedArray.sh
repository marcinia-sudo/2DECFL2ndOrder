#!/bin/bash


#SBATCH -A CSC113                                   # charge from account CSC113
#SBATCH --job-name="LargeShared"                    # set up job name
#SBATCH --output="bash/output/LSARR_%A_%a_%N.out"   # set up the name of std output
#SBATCH --error="bash/output/LSARR_%A_%a_%N.err"    # set up the name of error output
#SBATCH --array=0-3%4                               # n%m an array of n jobs where m jobs are run in parallel
#SBATCH --mem=128G                                  # RAM requested
#SBATCH --partition=large-shared                    # assign jobs to the share node
#SBATCH --nodes=1                                   # number of nodes 
#SBATCH --ntasks-per-node=128                       # number of tasks
#SBATCH --cpus-per-task=1                           # number cpu per task
#SBATCH --export=ALL                                # set up export
#SBATCH -t 48:00:00                                 # set up maximum run time
#SBATCH --mail-user=michael.arciniaga@gmail.com     # set up email address
#SBATCH --mail-type=begin                           # email me when the job starts
#SBATCH --mail-type=end                             # email me when the job finishes
#


pmatrix[0]="14 128 0.85 0.012 -0.2 0.0 0.17 -0.78 0.69 7 30 0.0001"
pmatrix[1]="14 128 0.80 0.012 -0.2 0.0 0.17 -0.78 0.69 7 30 0.0001"
pmatrix[2]="16 64 0.85 0.012 -0.2 0.0 0.17 -0.78 0.69 7 30 0.0001"
pmatrix[3]="16 64 0.80 0.012 -0.2 0.0 0.17 -0.78 0.69 7 30 0.0001"


#Set the number of runs that each SLURM task should do

PER_TASK=1

# Ex. Let say we want to run a.out 20 times with various parameters
# If we set PER_TASK = 4 this script will block run 0-3 in one task, 
# 4-8 in the next task, and so on for the job array. 


# Calculate the starting and ending values for this task based 
# on the SLURM task and the number of runs per tasks.
START_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))
END_NUM=$(( ($SLURM_ARRAY_TASK_ID + 1) * $PER_TASK - 1))

# Print the task and run range 

echo "MY SLURM_ARRAY_JOB_ID is $SLURM_ARRAY_JOB_ID"
echo This is task $SLURM_ARRAY_TASK_ID, which will do runs $START_NUM to $END_NUM

# RUN The loop of runs for this task.
for (( run=$START_NUM; run<=$END_NUM; run++ ))
do
    echo This is SLURM task $SLURM_ARRAY_TASK_ID, run number $run
    ./params.out ${pmatrix[$run]}
    echo -e "\n"
done

