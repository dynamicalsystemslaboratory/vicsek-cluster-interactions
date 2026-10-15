#!/bin/bash
#SBATCH --account=torch_pr_921_general
#SBATCH --mem=16G
#SBATCH --cpus-per-task=32
#SBATCH --job-name=size-generate-order_spd-rho-k
#SBATCH --output=size-generate-order_spd-rho-k_%A_%a.out
#SBATCH --error=size-generate-order_spd-rho-k_%A_%a.err
#SBATCH --array=1-484%24    # = 121 * N (here N = 4)

set -euo pipefail

Nr=11
Ns=11

# -----------------------------
# Julia runtime settings
# -----------------------------
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export JULIA_PKG_PRECOMPILE_AUTO=0

# IMPORTANT: shared depot (prevents recompilation every run)
export JULIA_DEPOT_PATH="$HOME/.julia"

# -----------------------------
# Debug info
# -----------------------------
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"
echo "Julia threads: $JULIA_NUM_THREADS"
echo "Julia depot: $JULIA_DEPOT_PATH"

# -----------------------------
# Index mapping
# -----------------------------
task_id=$((SLURM_ARRAY_TASK_ID - 1))

k=$((task_id / (Nr*Ns) + 1))
grid=$((task_id % (Nr*Ns)))

i=$((grid / Nr + 1))
j=$((grid % Ns + 1))

echo "Running indices: i=$i, j=$j, k=$k"

# Run Julia script
sleep $((RANDOM % 30))
julia --project=env generate-order.jl "$i" "$j" "$k"
