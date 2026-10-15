#!/bin/bash
#SBATCH --account=torch_pr_921_general
#SBATCH --mem=16G
#SBATCH --cpus-per-task=32
#SBATCH --job-name=generate-life_spd-mu
#SBATCH --output=generate-life_spd-mu_%A_%a.out
#SBATCH --error=generate-life_spd-mu_%A_%a.err
#SBATCH --array=1-441%21

# Give Julia the same number of threads as allocated CPUs
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Print some info for debugging
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"
echo "Julia threads: $JULIA_NUM_THREADS"

# Convert task ID -> (i,j) in range 1:21 × 1:21
task_id=$((SLURM_ARRAY_TASK_ID - 1))

i=$(( task_id / 21 + 1 ))
j=$(( task_id % 21 + 1 ))

echo "Running pair: ($i, $j)"

# Run Julia script
julia generate-life.jl $i $j
