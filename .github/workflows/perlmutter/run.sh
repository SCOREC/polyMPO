#!/bin/bash

name=polyMPO
cd $SCRATCH/globus-compute/$name-test

export root=$PWD
module load cmake
module load cray-hdf5
module load cray-netcdf
export MPICH_CXX=$root/kokkos/bin/nvcc_wrapper

cd build-$name
salloc --time 00:10:00 --constrain=gpu --qos=interactive --nodes=1 --ntasks-per-node=40 --cpus-per-task=1 --gpus=1 --account=m4564 ctest
cat $PWD/Testing/Temporary/LastTest.log