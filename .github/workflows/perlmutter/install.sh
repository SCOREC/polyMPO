#!/bin/bash

branch=$1

cd $SCRATCH/globus-compute/polyMPO-test

export root=$PWD
module load cmake/3.24.3
module load cray-hdf5
module load cray-netcdf

export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:\
$root/build-kokkos/install:\
$root/build-omegah/install:\
$root/build-cabana/install:\
$root/build-engpar/install:\
$root/build-pumipic/install

export MPICH_CXX=$root/kokkos/bin/nvcc_wrapper

# #kokkos
# git clone -b 4.5.00 https://github.com/kokkos/kokkos.git
# cmake -S kokkos -B build-kokkos \
#   -DCMAKE_INSTALL_PREFIX=build-kokkos/install \
#   -DCMAKE_CXX_COMPILER=$root/kokkos/bin/nvcc_wrapper \
#   -DKokkos_ARCH_AMPERE80=ON \
#   -DKokkos_ENABLE_SERIAL=ON \
#   -DKokkos_ENABLE_OPENMP=off \
#   -DKokkos_ENABLE_CUDA=ON \
#   -DKokkos_ENABLE_CUDA_LAMBDA=ON
# cmake --build build-kokkos -j 24 --target install

# #engpar
# unset MPICH_CXX #don't want nvcc_wrapper for engpar
# git clone https://github.com/SCOREC/EnGPar.git
# cmake -S EnGPar -B build-engpar \
#   -DCMAKE_INSTALL_PREFIX=build-engpar/install \
#   -DCMAKE_BUILD_TYPE=Release \
#   -DCMAKE_C_COMPILER=cc \
#   -DCMAKE_CXX_COMPILER=CC \
#   -DCMAKE_CXX_FLAGS="-std=c++11" \
#   -DENABLE_PARMETIS=OFF \
#   -DENABLE_PUMI=OFF \
#   -DIS_TESTING=OFF
# cmake --build build-engpar -j 24 --target install
# export MPICH_CXX=$root/kokkos/bin/nvcc_wrapper #restore use of nvcc_wrapper

# #omegah
# git clone -b scorec-v10.8.4 https://github.com/SCOREC/omega_h.git
# cmake -S omega_h -B build-omegah \
#   -DCMAKE_INSTALL_PREFIX=build-omegah/install \
#   -DCMAKE_BUILD_TYPE=Release \
#   -DBUILD_SHARED_LIBS=OFF \
#   -DOmega_h_USE_Kokkos=ON \
#   -DOmega_h_USE_CUDA=on \
#   -DOmega_h_CUDA_ARCH=80 \
#   -DOmega_h_USE_MPI=on  \
#   -DBUILD_TESTING=off  \
#   -DCMAKE_C_COMPILER=cc \
#   -DCMAKE_CXX_COMPILER=CC
# cmake --build build-omegah -j 24 --target install

# #cabana
# git clone -b 0.6.1 https://github.com/ECP-copa/Cabana.git cabana
# cmake -S cabana -B build-cabana \
#   -DCMAKE_INSTALL_PREFIX=build-cabana/install \
#   -DCMAKE_BUILD_TYPE=Release \
#   -DCMAKE_CXX_COMPILER=$root/kokkos/bin/nvcc_wrapper \
#   -DCabana_ENABLE_TESTING=OFF \
#   -DCabana_ENABLE_EXAMPLES=OFF
# cmake --build build-cabana -j 24 --target install

# #pumipic
# git clone --recursive https://github.com/SCOREC/pumi-pic.git
# cmake -S pumi-pic -B build-pumipic \
#   -DCMAKE_INSTALL_PREFIX=build-pumipic/install \
#   -DCMAKE_BUILD_TYPE=Release \
#   -DCMAKE_CXX_COMPILER=CC \
#   -DENABLE_CABANA=ON \
#   -DCMAKE_CXX_STANDARD=20 \
#   -DTEST_DATA_DIR=$root/pumi-pic/pumipic-data \
#   -DIS_TESTING=ON \
#   -DPS_IS_TESTING=ON
# cmake --build build-pumipic -j 24 --target install

# polyMPO
rm polyMPO -rf
git clone -b cws/pumipicDps https://github.com/SCOREC/polyMPO.git
cd polyMPO && git checkout $branch && cd -
rm build-polyMPO -rf
cmake -S polyMPO -B build-polyMPO \
  -DCMAKE_INSTALL_PREFIX=build-polyMPO/install \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=$root/kokkos/bin/nvcc_wrapper \
  -DCMAKE_Fortran_COMPILER=ftn \
  -DIS_TESTING=on
cmake --build build-polyMPO --target install -j4