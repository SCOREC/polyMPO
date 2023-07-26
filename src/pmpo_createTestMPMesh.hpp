#include "pmpo_MPMesh.hpp"
#include "pmpo_wachspressBasis.hpp"

namespace polyMPO {

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
polyMPO::Mesh* initTestMesh(int testMeshOption, int scaleFactor);

/*
 * \brief Initialize test particles to a given mesh and return the MPMesh object
 *        that relates the mesh to the material points
 * 
 * \details  Init particles at central position of the element
 * with a rand 4-6 duplicate of the particles and link them in an MPMesh object
 * 
 * \param meshObj (in) a mesh
 *
 * \param setMPOption (in) assign to every elements with 4,5,6,4,5,6...4,5,6 random MPs
 *
 * \return the MPMesh object with the meshObj and test MPs
 */
polyMPO::MPMesh initTestMPMesh(polyMPO::Mesh* meshObj,
                                   int setMPOption);
polyMPO::MaterialPoints* initTestMPs(polyMPO::Mesh* mesh,
                                         int testMPsOption);

}