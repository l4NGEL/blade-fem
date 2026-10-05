# Build/test environment: identical package set to the CI runner (Ubuntu 24.04).
FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential cmake ninja-build git ca-certificates \
        libeigen3-dev libboost-dev libboost-program-options-dev \
        libhdf5-dev catch2 \
        gmsh \
        python3 python3-numpy python3-matplotlib python3-h5py \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
