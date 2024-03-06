
#Instructions
# 1. Build polyMPO
# 2. Modify polyMPO/test/testFortranInit.f90
# 3. cd to buildPolyMPO and run "make VERBOSE=1"
# 4. Copy libs and flags from the output of step 3
# 5. Replace varibles below
# 6. cd to polyMPO/test and run "./use_at_your_own_risk/generateMakefile.sh"
# 7. run "make -f use_at_your_own_risk/MakefileGenerated.testFortranInit"
# 8. run "./testFortranInit", there should be no output

#Modify these

libs="/opt/scorec/spack/v0154_2/install/linux-rhel7-x86_64/gcc-6.5.0/gcc-10.1.0-tf5jjaditemasrbsl7tz6pnqa6duqwkg/bin/gfortran -Wl,-rpath -Wl,/opt/scorec/spack/v0154_2/install/linux-rhel7-x86_64/gcc-10.1.0/mpich-3.3.2-gi3wrjquyo564rk27x6r2c6ilr7ndmpl/lib -L/opt/scorec/spack/v0154_2/install/linux-rhel7-x86_64/gcc-10.1.0/hwloc-2.2.0-yuypimarf5tpp2wn4ljlgtnogpoxkzge/lib -g CMakeFiles/testFortranInit.dir/testFortranInit.f90.o -o testFortranInit   -L/opt/scorec/spack/v0154_2/install/linux-rhel7-x86_64/gcc-10.1.0/hwloc-2.2.0-yuypimarf5tpp2wn4ljlgtnogpoxkzge/lib  libpmpoUtils.a libpmpoUtils.a ../src/libpolympo_fortran.a libreadMPAS.a libpmpoUtils.a libreadMPAS.a ../src/libpolyMPO-core.a /usr/lib64/libm.so /lore/castia5/cranium/cpu-pumi-pic/build-pumipic-blockade-host/install/lib/libpumipic-core.a /lore/castia5/cranium/cpu-pumi-pic/build-pumipic-blockade-host/install/lib/libparticleStructs.a /lore/castia5/cranium/cpu-pumi-pic/build-pumipic-blockade-host/install/lib/libsupport.a /lore/castia5/cranium/cpu-pumi-pic/build-omegah-blockade-host/install/lib64/libomega_h.a /lore/castia5/cranium/cpu-pumi-pic/build-kokkos-blockade-host/install/lib64/libkokkoscontainers.a /lore/castia5/cranium/cpu-pumi-pic/build-kokkos-blockade-host/install/lib64/libkokkoscore.a -ldl /lore/castia5/cranium/cpu-pumi-pic/build-kokkos-blockade-host/install/lib64/libkokkossimd.a /usr/lib64/libz.so /lore/castia5/cranium/cpu-pumi-pic/build-engpar-blockade-host/install/lib/libengpar.a /lore/castia5/cranium/cpu-pumi-pic/build-engpar-blockade-host/install/lib/libdiffusive.a /lore/castia5/cranium/cpu-pumi-pic/build-engpar-blockade-host/install/lib/libmultilevel.a /lore/castia5/cranium/cpu-pumi-pic/build-engpar-blockade-host/install/lib/libengpar_metrics.a /lore/castia5/cranium/cpu-pumi-pic/build-engpar-blockade-host/install/lib/libcoloring.a /lore/castia5/cranium/cpu-pumi-pic/build-engpar-blockade-host/install/lib/libagi.a /lore/castia5/cranium/cpu-pumi-pic/build-engpar-blockade-host/install/lib/libengpar_support.a /lore/castia5/cranium/cpu-pumi-pic/build-engpar-blockade-host/install/lib/libpcu.a ../src/libpolympo_fortran.a /opt/scorec/spack/v0154_2/install/linux-rhel7-x86_64/gcc-10.1.0/mpich-3.3.2-gi3wrjquyo564rk27x6r2c6ilr7ndmpl/lib/libmpifort.so /opt/scorec/spack/v0154_2/install/linux-rhel7-x86_64/gcc-10.1.0/mpich-3.3.2-gi3wrjquyo564rk27x6r2c6ilr7ndmpl/lib/libmpi.so /opt/scorec/spack/v0154_2/install/linux-rhel7-x86_64/gcc-10.1.0/netcdf-fortran-4.5.2-eaneg3ac2jmh57rusifdxw65vitm6e65/lib/libnetcdff.so /opt/scorec/spack/v0154_2/install/linux-rhel7-x86_64/gcc-10.1.0/netcdf-c-4.7.3-jrc35lqeialwtfqnfxitlqesrgkajws7/lib/libnetcdf.so -lmpicxx -lmpi -lstdc++"

flags="/opt/scorec/spack/v0154_2/install/linux-rhel7-x86_64/gcc-6.5.0/gcc-10.1.0-tf5jjaditemasrbsl7tz6pnqa6duqwkg/bin/gfortran -DFP64 -DKOKKOS_ENABLED -DLOCAL_ID_FTYPE=C_INT32_T -DLOCAL_ID_TYPE=int32_t -DPOLYMPO_HAS_NETCDF -DPP_ENABLE_CAB -I/opt/scorec/spack/v0154_2/install/linux-rhel7-x86_64/gcc-10.1.0/netcdf-c-4.7.3-jrc35lqeialwtfqnfxitlqesrgkajws7/include -I/opt/scorec/spack/v0154_2/install/linux-rhel7-x86_64/gcc-10.1.0/netcdf-fortran-4.5.2-eaneg3ac2jmh57rusifdxw65vitm6e65/include -I/lore/castia5/cranium/cpu-pumi-pic/polyMPO/src -I/lore/castia5/cranium/cpu-pumi-pic/polyMPO/test -I/lore/castia5/cranium/cpu-pumi-pic/buildPolyMPO-CPU/src -I/lore/castia5/cranium/cpu-pumi-pic/build-pumipic-blockade-host/install/include -I/lore/castia5/cranium/cpu-pumi-pic/build-kokkos-blockade-host/install/include -I/usr/include -I/lore/castia5/cranium/cpu-pumi-pic/build-cabana-blockade-host/install/include -I/lore/castia5/cranium/cpu-pumi-pic/build-omegah-blockade-host/install/include -I/lore/castia5/cranium/cpu-pumi-pic/build-engpar-blockade-host/install/include -I/opt/scorec/spack/v0154_2/install/linux-rhel7-x86_64/gcc-10.1.0/mpich-3.3.2-gi3wrjquyo564rk27x6r2c6ilr7ndmpl/include -g -J../mod_files -c /lore/castia5/cranium/cpu-pumi-pic/polyMPO/test/testFortranInit.f90 -o CMakeFiles/testFortranInit.dir/testFortranInit.f90.o"

path=/lore/castia5/cranium/cpu-pumi-pic/buildPolyMPO-CPU/test

#Leave these alone

file=use_at_your_own_risk/MakefileGenerated.testFortranInit
printf "testFortranInit: testFortranInit.o \n\t $libs \n\ntestFortranInit.o: testFortranInit.f90 \n\t $flags" > $file

# replace ".." with path
sed -i "s@\.\.@$path/\.\.@g" $file

# set to build in current directory
sed -i "s@CMakeFiles/testFortranInit.dir/@@g" $file

# replace .f90.o with .o file
sed -i "s@\.f90\.o@\.o@g" $file

# add path to .a files without a path
sed -i "s@\ \([a-zA-Z]*\)\.a@\ $path/\1.a@g" $file