
## layout

- src - source files
- test - examples and tests that exercise APIs

## setup

```
git clone https://github.com/SCOREC/polyMpmTest.git
```

## build

The following assumes that a valid C and C++ compiler, and `cmake`, are in your PATH.

`CMAKE_INSTALL_PREFIX` is the path where the library, headers, and test binary
are installed.

`kk` is the path where kokkos is installed (either GPU or OpenMP based installation).

Set path as `export kk=/path/to/kokkos/install`

`kk_compiler` is the path to kokkos compiler
Set path as `export kk_compiler=/path/to/nvcc_wrapper` for GPU build.

Ignore this for OpenMP build

Load necessary modules:
```
module unuse /opt/scorec/spack/lmod/linux-rhel7-x86_64/Core
module use /opt/scorec/spack/v0154_2/lmod/linux-rhel7-x86_64/Core
module load \
gcc/10.1.0 \
openmpi \
cmake/3.20.0  \
cuda/11.4
```

### Install Kokkos

Clone the repo

```
git clone https://github.com/Kokkos/kokkos.git
```

Create a file named `installKokkos.sh` with the following contents:

```
bkend=$1
d=buildKokkos${bkend}
omp=""
cuda=""
[[ ${bkend} == "OPENMP" ]] && omp="-DKokkos_ENABLE_OPENMP=on"
[[ ${bkend} == "CUDA" ]] && cuda="-DKokkos_ENABLE_CUDA=on -DKokkos_ENABLE_CUDA_LAMBDA=on"
[[ ${bkend} == "OPENMP_CUDA" ]] && ompcuda="-DKokkos_ENABLE_OPENMP=on -DKokkos_ENABLE_CUDA=on -DKokkos_ENABLE_CUDA_LAMBDA=on"
cmake -S kokkos -B $d \
  -DCMAKE_CXX_COMPILER=g++ \
  -DKokkos_ARCH_TURING75=ON \
  -DBUILD_SHARED_LIBS=ON \
  $omp \
  $cuda \
  $ompcuda \
  -DKokkos_ENABLE_SERIAL=ON \
  -DKokkos_ENABLE_DEBUG=on \
  -DKokkos_ENABLE_TESTS=off \
  -DCMAKE_INSTALL_PREFIX=$d/install
cmake --build $d -j 8 --target install
```

Make it executable:

```
chmod +x installKokkos.sh
```

Install with CUDA backend

```
./installKokkos.sh CUDA
```

Install with OpenMP backend

```
./installKokkos.sh OPENMP
```

### Building for GPU execution

Create a file named `doConfigPolyMpmSheathGpu.sh` with the following contents:

```
bdir=$PWD/build-GPU
cmake -S polyMpmTest -B $bdir \
-DKokkos_DIR=$PWD/buildKokkosCUDA/install/lib64/cmake/Kokkos \
-DCMAKE_CXX_COMPILER=$PWD/buildKokkosCUDA/install/bin/nvcc_wrapper \
-DCMAKE_INSTALL_PREFIX=$bdir/install
```

Create a file named `buildPolyMpmSheathGpu.sh` with the following contents:

```
bdir=$PWD/build-GPU
cmake --build $bdir --target install
```

Make them executable:

```
chmod +x doConfigPolyMpmSheathGpu.sh buildPolyMpmSheathGpu.sh
```

Run the configure script then run the build script:

```
./doConfigPolyMpmSheathGpu.sh
./buildPolyMpmSheathGpu.sh
```

To rebuild the code after making some changes:

```
./buildPolyMpmSheathGpu.sh
```


### Building for CPU execution using OpenMP

Create a file named `doConfigPolyMpmSheathOmp.sh` with the following contents:

```
bdir=$PWD/build-omp
cmake -S polyMpmTest -B $bdir \
-DKokkos_DIR=$PWD/buildKokkosOPENMP/install/lib64/cmake/Kokkos \
-DCMAKE_INSTALL_PREFIX=$bdir/install
```

Create a file named `buildPolyMpmSheathOmp.sh` with the following contents:

```
bdir=$PWD/build-omp
cmake --build $bdir --target install
```

Make them executable:

```
chmod +x doConfigPolyMpmSheathOmp.sh buildPolyMpmSheathOmp.sh
```

Run the configure script then run the build script:

```
./doConfigPolyMpmSheathOmp.sh
./buildPolyMpmSheathOmp.sh
```

To rebuild the code after making some changes:

```
./buildPolyMpmSheathOmp.sh
```

## Run test

Run with OpenMP on CPU

```
./build-omp/test/PolyMpmSheathWachspress_Demo
```

Run with CUDA on GPU

```
./build-GPU/test/PolyMpmSheathWachspress_Demo
```

