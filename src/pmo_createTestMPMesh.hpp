#include "pmt_MPMesh.hpp"
#include "pmt_wachspressBasis.hpp"

namespace polyMpmTest {

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
polyMpmTest::Mesh initTestMesh(int factor);

/*
 * \brief Initialize test particles to a given mesh and return the MPMesh object
 *        that relates the mesh to the material points
 * 
 * \details  Init particles at central position of the element
 * with a rand 4-6 duplicate of the particles and link them in an MPMesh object
 * 
 * \param meshObj (in) a mesh
 * \param mpPerElement (in) number of material points per element
 *
 * \return the MPMesh object with the meshObj and test MPs
 */
polyMpmTest::MPMesh initTestMPMesh(polyMpmTest::Mesh& meshObj,
                                   std::vector<int>& mpPerElement);
polyMpmTest::MPMesh initTestMPMesh(polyMpmTest::Mesh& meshObj);
polyMpmTest::MaterialPoints* initTestMPs(polyMpmTest::Mesh& mesh,
                                         std::vector<int>& mpPerElement);

}
