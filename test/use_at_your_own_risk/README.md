The Makefiles in this directoery were generated manually by extracting the
compile and link line produced by CMake (`make VERBOSE=1`) for the GPU build.
It was tested by compiling one of our fortran tests that calls polyMPO APIs.

A big draw back of this approach, which has already bitten us once, is that if
we add new libraries or headers the Makefiles have to be recreated/updated
manually.

It is a crude approach, but since MPAS/MPM uses a Makefile based system and
isn't using CMake, autotools, pkgconfig etc. we don't have a significantly
better choice. 
