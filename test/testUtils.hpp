#ifndef POLYMPMTEST_TEST_TESTUTILS_H
#define POLYMPMTEST_TEST_TESTUTILS_H

#include "pmt_MPM.hpp"
#include "pmt_wachspressBasis.hpp"

#define MPMTEST_PI 3.14159265358979323846
#define MPMTEST_EPSILON 1.19209e-07 //using the float epsilon
#define randSeed 12345
#define arbitraryInt 10

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
Vector2View initT2LDelta(const int size, const double range = -1, const int randomSeed = randSeed);

/*
 * according to mp_init_poly_v5.py by Oncar
 */
MPM initMPMWithRandomMPs(Mesh& meshObj, int factor, const int randomSeed = randSeed);

/*
 * TODO: comments
 */
MPM initMPMWithCenterMPs(Mesh& meshObj, int factor);

/*
 * TODO: comments
 */
Vector2View initT2LDeltaRankineVortex(MPM mpm, Vector2 center, const int numEdge, const double avgLength, const double Gamma);

/*
 * TODO: comments
 */
Vector2View initT2LTest1(MPM mpm, double percent1, double percent2, double percent3, double percent4, const int randomSeed = randSeed);

/*
 * TODO: comments
 */
Vector2View initT2LTest2(MPM mpm, const int MPAcross, const int randomSeed = randSeed);

/*
 * TODO: comments
 */
Vector2View initT2LTest3(MPM mpm, const int MPAcross = 0,const double percent = 1.0, const int randomSeed = randSeed);

/*
 * TODO: comments
 */
double calcAvgLengthOfEdge(Mesh mesh, int printOut = 0);

/*
 * TODO: comments
 */
void runT2LSimple(MPM mpm, int size = -1, double range = -1, const int loopTimes = 1, const int printVTP = -1, const int randomSeed = randSeed);

/*
 * TODO: comments
 */
void runT2LRankineVortex(MPM mpm, Vector2 center, const int numEdge = 1, double avgLength = -1, const double Gamma = 1, const int loopTimes = 1, const int printVTP = -1);

/*
 * TODO: comments
 */
void runT2LRandomWithProportion(MPM mpm, const double p0 = 0.25, const double p1 = 0.25, const double p2 = 0.25, const double p3 = 0.25, const int loopTimes = 1, const int printVTP = -1, const int randomSeed = randSeed);

/*
 * TODO: comments
 */
void runT2LWithOneMPAcross(MPM mpm, const int MPAcross = 0, const int loopTimes = 1, const int printVTP = -1, const int randomSeed = randSeed);

/*
 * TODO: comments
 */
void runCVTRandomWithProportion(MPM mpm, const double p0 = 0.25, const double p1 = 0.25, const double p2 = 0.25, const double p3 = 0.25, const int loopTimes = 1, const int printVTP = -1, const int randomSeed = randSeed);

/*
 * TODO: comments
 */
void runCVTElmCenterBasedRandomWithProportion(MPM mpm, const double p0 = 0.25, const double p1 = 0.25, const double p2 = 0.25, const double p3 = 0.25, const int loopTimes = 1, const int printVTP = -1, const int randomSeed = randSeed);

/*
 * TODO: comments
 */
void runCVTAllAcrossTest(MPM mpm, const int MPAcross = 0, const double percent = 1.0, const int loopTimes = 1, const int printVTP = -1, const int randomSeed = randSeed);

/*
 * TODO: comments
 */
void runCVTElmAllAcrossTest(MPM mpm, const int MPAcross = 0, const double percent = 1.0, const int loopTimes = 1, const int printVTP = -1, const int randomSeed = randSeed);

/*
 * TODO: comments
 */
void printMPs(MPM mpm);

#endif

