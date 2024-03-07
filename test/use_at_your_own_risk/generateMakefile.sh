#Instructions
# 1. Build polyMPO
# 2. (Optional) Below replace path to build dir
# 3. cd to polyMPO/test/use_at_your_risk
# 3. run "chmod +x generateMakefile.sh"
# 4. run "./generateMakefile.sh"
# 5. run "make -f MakefileGenerated.testFortranInit"
# 6. run "./testFortranInit", there should be no output

#Modify these

path=../../../buildPolyMPO-CPU/test

#Leave these alone

libs=$(cat $path/CMakeFiles/testFortranInit.dir/link.txt)
flags=$(cat $path/CMakeFiles/testFortranInit.dir/flags.make)
file=MakefileGenerated.testFortranInit

#Remove extra
flags=$(echo $flags | sed "s@#.*with@@g")
flags=$(echo $flags | sed "s@Fortran_DEFINES = @@g")
flags=$(echo $flags | sed "s@Fortran_INCLUDES = @@g")
flags=$(echo $flags | sed "s@Fortran_FLAGS = @@g")

# break into multiple lines
libs=$(echo $libs | sed "s@\ @ \\\\\n\t@g")
flags=$(echo $flags | sed "s@\ @ \\\\\n\t@g")

printf "testFortranInit: testFortranInit.o \n\t$libs \n\ntestFortranInit.o: testFortranInit.f90 \n\t$flags -c testFortranInit.f90 -o testFortranInit.o" > $file

# replace ".." with path
sed -i "s@\.\.@$path/\.\.@g" $file

# set to build in current directory
sed -i "s@CMakeFiles/testFortranInit.dir/@@g" $file

# replace ".f90.o" with ".o" file
sed -i "s@\.f90\.o@\.o@g" $file

# add path to ".a" files without a path
sed -i "s@\t\([a-zA-Z]*\)\.a@\t$path/\1.a@g" $file

# remove netcdf and hdf5
sed -i "/netcdf/d" $file
sed -i "/NETCDF/d" $file
sed -i "/hdf5/d" $file

# run from use_at_your_own_risk
sed -i 's@testFortranInit.f90@../testFortranInit.f90@g' $file