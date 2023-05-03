## code layout

- src - source files
- test - examples and tests that exercise APIs

## read this first

The instructions in the 'setup' section should be run once to get the stack of software built.  After that, the instructions in the 'develop and rebuild' section should be followed to rebuild polyMPO after making source code changes.

CUDA flags for the NVIDIA Turing GPUs are included below.  The scorec workstations `cranium` and `blockade` have Turing GPUs.

The following assumes that a valid C and C++ compiler, and `cmake`, are in your PATH.  On SCOREC systems these are provided by `module` commands.  If you are not on a SCOREC system these must be edited accordingly.

`CMAKE_INSTALL_PREFIX` is the path where the library, headers, and test binary
are installed.

## install dependencies

Create a directory to work from.  It will contain all source code and build directories.

```
mkdir polyMpoDev #this can be any name - just be consistent
cd polyMpoDev
```

Create an environment script `setupEnvironment.sh` with the following contents.  **It contains SCOREC specific `module` commands that will have to be modified if you are building on a non-SCOREC system.**

```
export root=$PWD 
module unuse /opt/scorec/spack/lmod/linux-rhel7-x86_64/Core 
module use /opt/scorec/spack/v0154_2/lmod/linux-rhel7-x86_64/Core 
module load gcc/10.1.0 cuda/11.4
module load mpich/3.3.2
module load cmake
module load netcdf-c/4.7.3

function getname() {
  name=$1
  machine=`hostname -s`
  buildSuffix=${machine}-cuda
  echo "build-${name}-${buildSuffix}"
}
export engpar=$root/`getname engpar`/install # This is where engpar will be (or is) installed
export kk=$root/`getname kokkos`/install   # This is where kokkos will be (or is) installed
export oh=$root/`getname omegah`/install  # This is where omega_h will be (or is) installed
export cab=$root/`getname cabana`/install # This is where cabana will be (or is) installed
export pumipic=$root/`getname pumipic`/install # This is where PumiPIC will be (or is) installed
export CMAKE_PREFIX_PATH=$engpar:$kk:$kk/lib64/cmake:$oh:$cab:$pumipic:$CMAKE_PREFIX_PATH
export MPICH_CXX=$root/kokkos/bin/nvcc_wrapper
```

Create a file named `buildAll.sh` with the following contents. **It contains compiler flags specific to NVIDIA GPUs with the Turing architecture (i.e., `-DKokkos_ARCH_TURING75=ON` and `-DOmega_h_CUDA_ARCH=75`) that need to be modified to match the architecture of the GPU in your system.**

```
#!/bin/bash -e

cd $root

#kokkos
git clone -b 3.1.00 https://github.com/kokkos/kokkos.git
mkdir -p $kk
cmake -S kokkos -B ${kk%%install} \
  -DCMAKE_INSTALL_PREFIX=$kk \
  -DCMAKE_CXX_COMPILER=$root/kokkos/bin/nvcc_wrapper \
  -DKokkos_ARCH_TURING75=ON \
  -DKokkos_ENABLE_SERIAL=ON \
  -DKokkos_ENABLE_OPENMP=off \
  -DKokkos_ENABLE_CUDA=on \
  -DKokkos_ENABLE_CUDA_LAMBDA=on \
  -DKokkos_ENABLE_DEBUG=on
cmake --build ${kk%%install} -j 24 --target install

#engpar
mkdir -p $engpar
git clone https://github.com/SCOREC/EnGPar.git
cmake -S EnGPar -B ${engpar%%install} \
  -DCMAKE_INSTALL_PREFIX=$engpar \
  -DCMAKE_BUILD_TYPE="Release" \
  -DCMAKE_C_COMPILER="mpicc" \
  -DCMAKE_CXX_COMPILER="mpicxx" \
  -DCMAKE_CXX_FLAGS="-std=c++11" \
  -DENABLE_PARMETIS=OFF \
  -DENABLE_PUMI=OFF \
  -DIS_TESTING=OFF
cmake --build ${engpar%%install} -j 24 --target install

#omegah
mkdir -p $oh
git clone https://github.com/SCOREC/omega_h.git
cmake -S omega_h -B ${oh%%install} \
  -DCMAKE_INSTALL_PREFIX=$oh \
  -DCMAKE_BUILD_TYPE="Release" \
  -DBUILD_SHARED_LIBS=OFF \
  -DOmega_h_USE_Kokkos=ON \
  -DOmega_h_USE_CUDA=on \
  -DOmega_h_CUDA_ARCH=75 \
  -DOmega_h_USE_MPI=on  \
  -DBUILD_TESTING=off  \
  -DCMAKE_C_COMPILER=`which mpicc` \
  -DCMAKE_CXX_COMPILER=`which mpicxx` \
  -DKokkos_PREFIX=$kk/lib64/cmake
cmake --build ${oh%%install} -j 24 --target install

#cabana
mkdir -p $cab
git clone -b 0.3.0 https://github.com/ECP-copa/Cabana.git cabana
cmake -S cabana -B ${cab%%install} \
  -DCMAKE_INSTALL_PREFIX=$cab \
  -DCMAKE_BUILD_TYPE="Release" \
  -DCMAKE_CXX_COMPILER=$root/kokkos/bin/nvcc_wrapper \
  -DCabana_ENABLE_Cuda=ON \
  -DCabana_ENABLE_TESTING=OFF \
  -DCabana_ENABLE_EXAMPLES=OFF
cmake --build ${cab%%install} -j 24 --target install

#pumipic
mkdir -p $pumipic
git clone --recursive https://github.com/SCOREC/pumi-pic.git
cmake -S pumi-pic -B ${pumipic%%install} \
  -DCMAKE_INSTALL_PREFIX=$pumipic \
  -DCMAKE_BUILD_TYPE="Debug" \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DENABLE_CABANA=ON \
  -DTEST_DATA_DIR=$root/pumi-pic/pumipic-data \
  -DOmega_h_PREFIX=$oh \
  -DEnGPar_PREFIX=$engpar \
  -DIS_TESTING=OFF \
  -DPS_IS_TESTING=OFF
cmake --build ${pumipic%%install} -j 24 --target install
```

Make the script executable:

```
chmod +x buildAll.sh
```


Source the environment script from this work directory:

```
source setupEnvironment.sh
```

Run the build script:

```
./buildAll.sh
```

### Build polyMPO


Create a file named `doConfigPolyMpo-GPU.sh` with the following contents:

```
bdir=$PWD/buildPolyMPO-GPU
cmake -S polyMPO -B $bdir \
-DKokkos_DIR=$kk/lib64/cmake/Kokkos \
-DCMAKE_CXX_COMPILER=$kk/bin/nvcc_wrapper \
-DCMAKE_INSTALL_PREFIX=$bdir/install
```

Create a file named `buildPolyMpo-GPU.sh` with the following contents:

```
bdir=$PWD/buildPolyMPO-GPU
cmake --build $bdir --target install -j4
```

Make them executable:

```
chmod +x doConfigPolyMpo-GPU.sh buildPolyMpo-GPU.sh
```

Clone the repo

```
git clone -b cws/pumipicDps https://github.com/SCOREC/polyMPO.git
```

Run the configure script then run the build script:

```
./doConfigPolyMpo-GPU.sh
./buildPolyMpo-GPU.sh
```

## Run polyMPO tests

```
ctest --test-dir $root/buildPolyMPO-GPU
```

Show stdout and stderr when the tests run:

```
ctest --test-dir $root/buildPolyMPO-GPU -V
```

Show stdout and stderr *and* run only the `unitTest`:

```
ctest --test-dir $root/buildPolyMPO-GPU -V -R unit
```

## polyMPO: develop and rebuild

To resume work on polyMPO (i.e., when starting a new shell session/terminal) run the following commands to setup your environment:

```
cd /path/to/polyMpoDev
source setupEnvironment.sh
```

Note, `setupEnvironment.sh` **MUST** be sourced from the top-level work directory (`polyMpoDev`).

### C++ source changes only

Assuming changes existing polyMPO C++ source/header files were made you can just run make as follows:

```
./buildPolyMpo-GPU.sh
```

### CMake changes

If CMake files were changed (i.e., to add new C++ source files) then you should delete the contents of the build directory, then rerun the configure and build scripts:

```
cd $root
rm -rf buildPolyMPO-GPU
./doConfigPolyMpo-GPU.sh
./buildPolyMpo-GPU.sh
```

