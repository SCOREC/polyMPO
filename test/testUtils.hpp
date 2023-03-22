#include "pmt_MPM.hpp"
#include "pmt_wachspressBasis.hpp"

using namespace polyMpmTest;

/*  \breif Init a simple test square mesh
 *  \details initialize a hard-coded unit square mash.
 *  There will be triangualr to Octahedral elements 
 *  the factor will make factors of duplicates to test as a larger mesh
 */
Mesh initTestMesh(int factor);

/*
 * \breif Initialize test particles to a given mesh
 * \details  Init particles at central position of the element
 * with a rand 4-6 duplicate of the particles and link them in an MPM object
 */
MPM initTestMPM(Mesh& meshObj);

/*
 * \breif Wachspress driver function
 * \details Wachspress driver function to calculate the Wachspress and gradF at each MP
 * (e.g. the functions in src/pmt_wachspressBasis.hpp)
 */
void interpolateWachspress(MPM& mpm);
