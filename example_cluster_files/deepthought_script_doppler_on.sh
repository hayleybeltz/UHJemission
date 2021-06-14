#!/bin/bash

# main file used to submit scripts to deepthought2 for running Eliza's 3D transmission code.

# The -a flag is the job array specifier. 1,000 individual files are submitted to the 
# cluster, each with a number between 1 and 1000 assigned to SLURM_ARRAY_TASK_ID.

# The -t flag is the time specifier. I usually give it ~15% longer than I expect it to
# take on the given CPU, just in case. The only penalty for making your time too long is
# that it takes longer to get off the queue, whereas the penalty for making your time
# too short is that you lose all your output!

# The --mem-per-cpu flag species how much memory (in megabytes) you want assigned
# to your job. I arrived at this number by trial and error!

# The ntasks flag refers to the number of processes that will be spawned. We aren't doing
# any fancy parallelization, so this can be kept to 1.

# The --share flag lets the supercomputer scheduler know that we can share our resources with
# other people. This is fine so long as we're not doing anything top-secret, and it makes
# scheduling our job a bit easier.

# The --constraint line is specific to deepthought; there's a particular architecture that 
# works well with my code.

# Useful reference with a bunch of sbatch commands: https://slurm.schedmd.com/sbatch.html
# author: @arjunsavel

#SBATCH -a 1-1000
#SBATCH -t 63:00
#SBATCH --mem-per-cpu=5120
#SBATCH --ntasks=1
#SBATCH --share
#SBATCH --constraint="rhel8"

# each job runs the below bash commmands.

# go to my home directory
cd /lustre/asavel/wasp_76

# activate my environment
source wasp_76_env/bin/activate

# run my script! the args are the job array, doppler on/off, and the number of 
# jobs that have already been run.
python3 run_RT_deepthought.py $SLURM_ARRAY_TASK_ID 1 56000

