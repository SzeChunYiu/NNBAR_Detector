#!/bin/bash
# Step 1: Install CUDA Toolkit
set -e

echo "Installing CUDA Toolkit 12.4..."

# Add NVIDIA repository
cd /tmp
wget -q https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-keyring_1.1-1_all.deb
sudo dpkg -i cuda-keyring_1.1-1_all.deb
rm cuda-keyring_1.1-1_all.deb

# Update and install CUDA toolkit
sudo apt-get update
sudo apt-get install -y cuda-toolkit-12-4

# Add to PATH
if ! grep -q "cuda" ~/.bashrc; then
    echo 'export PATH=/usr/local/cuda/bin:$PATH' >> ~/.bashrc
    echo 'export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH' >> ~/.bashrc
fi

export PATH=/usr/local/cuda/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH

echo ""
echo "CUDA installed successfully!"
nvcc --version
