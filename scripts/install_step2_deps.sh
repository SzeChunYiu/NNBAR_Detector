#!/bin/bash
# Step 2: Install system dependencies
set -e

echo "Installing system dependencies for Geant4 and GPU packages..."

sudo apt-get update
sudo apt-get install -y \
    build-essential \
    cmake \
    git \
    wget \
    libxerces-c-dev \
    libexpat1-dev \
    zlib1g-dev \
    libgl1-mesa-dev \
    libglu1-mesa-dev \
    libxt-dev \
    libxmu-dev \
    libxi-dev \
    qtbase5-dev \
    libqt5opengl5-dev \
    libgsl-dev \
    libboost-dev

echo "System dependencies installed!"
