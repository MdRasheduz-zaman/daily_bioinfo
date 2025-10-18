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