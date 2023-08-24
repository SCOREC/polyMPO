#ifndef POLYMPO_TEST_UTIL 
#define POLYMPO_TEST_UTIL 

#include "pmpo_MPMesh.hpp"
#include "pmpo_wachspressBasis.hpp"

#define TEST_EPSILON 1e-6
#define TEST_PI 3.14159265358979323846
/*
 * \breif Wachspress driver function
 *
 * \details Wachspress driver function to calculate the Wachspress and gradF at each MP
 * (e.g. the functions in src/pmt_wachspressBasis.hpp)
 *
 * \param any resonable mpm to test
 */
void interpolateWachspress2DTest(polyMPO::MPMesh& mpMesh);

void interpolateWachspress3DTest(polyMPO::MPMesh& mpMesh);

extern "C" void setWithMPASMeshByFortran(void** p_mpMesh,
                                         const char* fileNameame,
                                         const int length);
void printVTP(polyMPO::MPMesh& mpMesh);
#endif
