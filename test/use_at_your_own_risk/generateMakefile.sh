
#Instructions
# 1.  Build polyMPO
# 2.  Modify polyMPO/test/testFortranInit.f90 so that it compiles again
# 3.  cd to buildPolyMPO and run "make VERBOSE=1"
# 4.  Copy the libs starting at compiler after "Linking Fortran executable testFortranInit"
# 5.  Copy the flags starting at compiler after "Building Fortran object test/CMakeFiles/testFortranInit.dir/testFortranInit.f90.o"
# 6.  Replace varibles below
# 7.  cd to polyMPO/test and run "chmod +x use_at_your_own_risk/generateMakefile.sh"
# 8.  cd to polyMPO/test and run "./use_at_your_own_risk/generateMakefile.sh"
# 9.  cd to polyMPO/test and run "make -f use_at_your_own_risk/MakefileGenerated.testFortranInit"
# 10. cd to polyMPO/test and run "./testFortranInit", there should be no output

#Modify these

path=/lore/castia5/cranium/cpu-pumi-pic/buildPolyMPO-CPU/test

#Leave these alone

libs=$(cat $path/CMakeFiles/testFortranInit.dir/link.txt)
flags=$(cat $path/CMakeFiles/testFortranInit.dir/flags.make)

#Remove extra
flags=$(echo $flags | sed "s@#.*with@@g")
flags=$(echo $flags | sed "s@Fortran_DEFINES = @@g")
flags=$(echo $flags | sed "s@Fortran_INCLUDES = @@g")
flags=$(echo $flags | sed "s@Fortran_FLAGS = @@g")
echo $flags

# break into multiple lines
libs=$(echo $libs | sed "s@\ @ \\\\\n\t@g")
flags=$(echo $flags | sed "s@\ @ \\\\\n\t@g")

file=use_at_your_own_risk/MakefileGenerated.testFortranInit
printf "testFortranInit: testFortranInit.o \n\t$libs \n\ntestFortranInit.o: testFortranInit.f90 \n\t$flags -c testFortranInit.f90 -o testFortranInit.o" > $file

# replace ".." with path
sed -i "s@\.\.@$path/\.\.@g" $file

# set to build in current directory
sed -i "s@CMakeFiles/testFortranInit.dir/@@g" $file

# replace ".f90.o" with ".o" file
sed -i "s@\.f90\.o@\.o@g" $file

# add path to ".a" files without a path
sed -i "s@\t\([a-zA-Z]*\)\.a@\t$path/\1.a@g" $file