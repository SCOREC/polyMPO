#include "pmpo_MPMesh.hpp"
#include "pmpo_wachspressBasis.hpp"

/*
 * \breif Wachspress driver function
 *
 * \details Wachspress driver function to calculate the Wachspress and gradF at each MP
 * (e.g. the functions in src/pmt_wachspressBasis.hpp)
 *
 * \param any resonable mpm to test
 */
void interpolateWachspress(polyMPO::MPMesh& mpMesh);
