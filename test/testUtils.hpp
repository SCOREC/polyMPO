#include "MPM.hpp"
#include "mesh.hpp"
#include "materialPoints.hpp"

using namespace polyMpmTest;

/* TODO: add some comments
 * Init a simple test square mesh
 */
Mesh initTestMesh(int factor);

/*
 * Init particles at central position with a rand 4-6 duplicate
 * put in an MPM file;
 */
MPM initTestMPM(Mesh& meshObj);
