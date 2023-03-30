#ifndef POLYMPMTEST_TEST_TESTUTILS_H
#define POLYMPMTEST_TEST_TESTUTILS_H

#include "pmt_MPM.hpp"
#include "pmt_wachspressBasis.hpp"

using namespace polyMpmTest;

/*  \breif Init a simple test square mesh
 * 
 *  \details initialize a hard-coded unit square mash.
 *  There will be triangualr to Octahedral elements 
 *  element = 10*factor vertices = 19*factor
 *  
 *  \param factor: the function makes factors of duplicates to test as a larger mesh
 *
 *  \return return the hard-coded mesh. 
 */
Mesh initTestMesh(int factor);

/*
 * \breif Initialize test particles to a given mesh
 * 
 * \details  Init particles at central position of the element
 * with a rand 4-6 duplicate of the particles and link them in an MPM object
 * 
 * \param any resonable meshObj
 *
 * \return the MPM object with the meshObj and test MPs
 */
MPM initTestMPM(Mesh& meshObj);

/*
 * \breif Wachspress driver function
 *
 * \details Wachspress driver function to calculate the Wachspress and gradF at each MP
 * (e.g. the functions in src/pmt_wachspressBasis.hpp)
 *
 * \param any resonable mpm to test
 */
void interpolateWachspress(MPM& mpm);


/*
 * TODO: finish the comments
 */
Vector2View InitT2LDelta(int size);

/*
 * according to mp_init_poly_v5.py by Oncar
 */
MPM initMPMWithRandomMPs(Mesh& meshObj, int factor);

#endif
