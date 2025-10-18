# Working Setup Singularity

**Last Updated**: October 15, 2025  
**Status**: ✅ Training successfully running  
**Location**: `/home/md.rasheduzzaman/data/molestus_Meta/NT_LM/env_cuda/nt_container`

---
## `bashrc` file content goes first

```bash
# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
    . /etc/bashrc
fi

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/data/users/md.rasheduzzaman/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/data/users/md.rasheduzzaman/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/data/users/md.rasheduzzaman/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/data/users/md.rasheduzzaman/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup

if [ -f "/data/users/md.rasheduzzaman/miniconda3/etc/profile.d/mamba.sh" ]; then
    . "/data/users/md.rasheduzzaman/miniconda3/etc/profile.d/mamba.sh"
fi
# <<< conda initialize <<<

# --- MAMBA INTEGRATION ---
# (Mamba is a fast drop-in replacement for Conda)
if command -v mamba &> /dev/null; then
    alias conda='mamba'
fi

# Redirect pip cache
export PIP_CACHE_DIR=/data/users/md.rasheduzzaman/.cache/pip
export HF_HOME=/data/users/md.rasheduzzaman/.cache/huggingface
export NUCLEOTIDE_TRANSFORMER_CACHE=/data/users/md.rasheduzzaman/.cache/nucleotide_transformer

# Temporary directories under /data (fast and quota-friendly)
export TMPDIR=/data/users/md.rasheduzzaman/tmp

# Ensure cache directories exist
mkdir -p /data/users/md.rasheduzzaman/.cache/{pip,huggingface,nucleotide_transformer}
mkdir -p /data/users/md.rasheduzzaman/tmp

# Add local bin last
export PATH=$HOME/.local/bin:$PATH
```

## Problem Summary

- **Original Issue**: glibc 2.17, cuda11.1 on HPC cluster, needed glibc 2.27+ and cuda111.8+ for PyTorch/transformers
- **Solution**: Use Singularity container with NVIDIA PyTorch base image
- **System**: Tesla V100S-PCIE-32GB (32GB VRAM), sbatch scheduler

## Container Setup (Completed)

### 1. Pull Pre-built Container
```bash
module load singularity

# Pull NVIDIA PyTorch container (Python 3.12, PyTorch 2.6.0)
singularity pull docker://nvcr.io/nvidia/pytorch:24.12-py3

# Rename
mv pytorch_24.12-py3.sif nt_v2_finetune.sif
```

### 2. Install Python Packages
```bash
# Install to user directory (persists across container runs)
singularity exec --nv nt_v2_finetune.sif pip install --user \
    transformers==4.44.2 \
    peft>=0.11.1 \
    datasets>=3.0.1 \
    accelerate>=1.0.0 \
    evaluate \
    scikit-learn \
    numpy \
    pandas \
    tqdm \
    matplotlib \
    seaborn
```
**Note**: Packages are installed to `~/.local/` and persist across container runs.

**Note**: Used `transformers==4.44.2` (not 4.45.0) due to PyTorch 2.6.0a compatibility issues.

### 3. Verify Installation
```bash
singularity exec --nv \
    --bind /data/users/md.rasheduzzaman:/data/users/md.rasheduzzaman \
    --bind /home/md.rasheduzzaman:/home/md.rasheduzzaman \
    --env PYTHONUSERBASE=/home/md.rasheduzzaman/.local \
    nt_v2_finetune.sif python -c "
import torch
import transformers, peft, datasets, accelerate
print(f'PyTorch: {torch.__version__}')
print(f'CUDA: {torch.cuda.is_available()}')
print(f'GPU: {torch.cuda.get_device_name(0)}')
print(f'Transformers: {transformers.__version__}')
print(f'PEFT: {peft.__version__}')
print('✓ All packages loaded successfully')
"
```

**Expected Output**:
```
PyTorch: 2.6.0a0+df5bbc09d1.nv24.12
CUDA: True
GPU: Tesla V100S-PCIE-32GB
Transformers: 4.44.2
PEFT: 0.17.1
✓ All packages loaded successfully
```

---

## Issues Encountered & Fixes

### Issue 1: No Mapping Entry for Fakeroot
**Error**: `could not use fakeroot: no mapping entry found in /etc/subuid`  
**Solution**: Used pre-built container instead of building from scratch

### Issue 2: Module Not Found (numpy)
**Error**: `ModuleNotFoundError: No module named 'numpy'`  
**Cause**: Conda base environment interfering with container Python  
**Solution**: Use `--cleanenv` flag to isolate container environment


## Working Training Command

### Interactive Session
```bash
srun --pty -p gpu --gres=gpu:1 -c 32 --mem=280g -t 24:00:00 /bin/bash
#wait for resource allocation

# First: deactivate conda or use --cleanenv
conda deactivate
cd /data/users/md.rasheduzzaman/molestus_Meta/NT_LM/env_cuda/nt_container

# Run training
singularity exec --nv --cleanenv \
    --bind /data/users/md.rasheduzzaman:/data/users/md.rasheduzzaman \
    --bind /home/md.rasheduzzaman:/home/md.rasheduzzaman \
    --env PATH=/home/md.rasheduzzaman/.local/bin:/usr/local/bin:/usr/bin:/bin \
    --env PYTHONUSERBASE=/home/md.rasheduzzaman/.local \
    --env HF_HOME=/data/users/md.rasheduzzaman/.cache/huggingface \
    --env TMPDIR=/data/users/md.rasheduzzaman/tmp \
    --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
    nt_v2_finetune.sif python3 -u scripts/train_lora_ntv2_50m.py \
        --train-csv processed_data/train_chunks.csv \
        --valid-csv processed_data/valid_chunks.csv \
        --output-dir results/ntv2_50m_lora \
        --num-train-epochs 5 \
        --per-device-train-batch-size 2 \
        --per-device-eval-batch-size 4 \
        --fp16 \
        --gradient-accumulation-steps 6
```
## slurm (sbatch) scheduler script: script name `submit_train.sbatch`

```bash
#!/bin/bash
#SBATCH --job-name=nt_lora_optimized
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:1
#SBATCH --mem=320G
#SBATCH --time=120:00:00
#SBATCH --output=logs/train_optimized_%j.log
#SBATCH --error=logs/train_optimized_%j.err.txt
#SBATCH --exclusive

set -e  # Exit on error
set -x  # Print commands (for debugging)

# ==============================================================================
# Optimized Training Script - Maximum Resource Utilization
# ==============================================================================

cd /home/md.rasheduzzaman/data/molestus_Meta/NT_LM/env_cuda/nt_container

module load singularity
conda deactivate 2>/dev/null || true

mkdir -p logs

echo "=========================================="
echo "Job Information"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Job started: $(date)"
echo "Node: $(hostname)"
echo "Working directory: $PWD"

# Resource allocation
echo ""
echo "=========================================="
echo "Resource Allocation"
echo "=========================================="
echo "CPUs requested: 32"
echo "Memory requested: 320 GB"
echo "GPUs requested: 1"

# GPU information
echo ""
echo "GPU Information:"
nvidia-smi --query-gpu=name,memory.total,driver_version,compute_cap --format=csv,noheader

echo "=========================================="

# Set CPU affinity and threading
export OMP_NUM_THREADS=32
export MKL_NUM_THREADS=32
export OPENBLAS_NUM_THREADS=32
export NUMEXPR_NUM_THREADS=32

# Optimal NCCL settings for single-GPU
export NCCL_DEBUG=WARN
export NCCL_P2P_DISABLE=1

# PyTorch optimizations
export CUDA_LAUNCH_BLOCKING=0
export TORCH_CUDNN_V8_API_ENABLED=1

# Auto-resume from latest checkpoint
CHECKPOINT=""
if [ -d "results/ntv2_50m_lora" ]; then
    LATEST=$(ls -td results/ntv2_50m_lora/checkpoint-* 2>/dev/null | head -1)
    if [ ! -z "$LATEST" ]; then
        CHECKPOINT="--resume-from-checkpoint $LATEST"
        echo ""
        echo "Resuming from checkpoint: $LATEST"
        echo ""
    fi
fi

echo "=========================================="
echo "Starting Training"
echo "=========================================="

# Run training with maximum resource utilization
singularity exec --nv --cleanenv \
    --bind /data/users/md.rasheduzzaman:/data/users/md.rasheduzzaman \
    --bind /home/md.rasheduzzaman:/home/md.rasheduzzaman \
    --env PATH=/home/md.rasheduzzaman/.local/bin:/usr/local/bin:/usr/bin:/bin \
    --env PYTHONUSERBASE=/home/md.rasheduzzaman/.local \
    --env HF_HOME=/data/users/md.rasheduzzaman/.cache/huggingface \
    --env TMPDIR=/data/users/md.rasheduzzaman/tmp \
    --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
    --env OMP_NUM_THREADS=32 \
    --env MKL_NUM_THREADS=32 \
    --env TORCH_CUDNN_V8_API_ENABLED=1 \
    --env DATALOADER_NUM_WORKERS=16 \
    nt_v2_finetune.sif python3 -u scripts/train_lora_ntv2_50m.py \
        --train-csv processed_data/train_chunks.csv \
        --valid-csv processed_data/valid_chunks.csv \
        --output-dir results/ntv2_50m_lora \
        --num-train-epochs 5 \
        --per-device-train-batch-size 2 \
        --per-device-eval-batch-size 8 \
        --fp16 \
        --gradient-accumulation-steps 6 \
        --save-steps 10000 \
        --eval-steps 2000 \
        --dataloader-num-workers 16 \
        --dataloader-prefetch-factor 4 \
        $CHECKPOINT

EXIT_CODE=$?

echo ""
echo "=========================================="
echo "Job Summary"
echo "=========================================="
echo "Job finished: $(date)"
echo "Exit code: $EXIT_CODE"

# Final GPU stats
echo ""
echo "Final GPU State:"
nvidia-smi

# Checkpoint summary
if [ -d "results/ntv2_50m_lora" ]; then
    echo ""
    echo "Saved Checkpoints:"
    ls -lth results/ntv2_50m_lora/checkpoint-* 2>/dev/null | head -5
fi

echo "=========================================="

exit $EXIT_CODE
```
Now we need to make it executable and submit it:
```bash
chmod +x submit_train.sbatch
sbatch submit_train.sbatch
```