#include "pmt_MPM.hpp"
#include "pmt_assembly.hpp"
#include "testUtils.hpp"


int main() {
    Kokkos::initialize();

//  test Vector2
    auto v = polyMpmTest::Vector2();
    v[0] = 1;
    v[1] = 2;
    PMT_ALWAYS_ASSERT(v[0] == 1);
    PMT_ALWAYS_ASSERT(v[1] == 2);
//  test Vector operator
    auto v1 = polyMpmTest::Vector2(1,2);
    auto v2 = polyMpmTest::Vector2(3,4);
    auto v3 = v1 + v2;
    PMT_ALWAYS_ASSERT(v3[0] == 4);
    PMT_ALWAYS_ASSERT(v3[1] == 6);
    auto v4 = v1 - v2;
    PMT_ALWAYS_ASSERT(v4[0] == -2);
    PMT_ALWAYS_ASSERT(v4[1] == -2);
    auto v5 = v1 * 2;
    PMT_ALWAYS_ASSERT(v5[0] == 2);
    PMT_ALWAYS_ASSERT(v5[1] == 4);
    auto v6 = -v1;
    PMT_ALWAYS_ASSERT(v6[0] == -1);
    PMT_ALWAYS_ASSERT(v6[1] == -2);
    auto v7 = v1.dot(v2);
    PMT_ALWAYS_ASSERT(v7 == 11);
    auto v8 = v1.cross(v2);
    PMT_ALWAYS_ASSERT(v8 == -2);
    auto v9 = v1.magnitude();
    PMT_ALWAYS_ASSERT(v9 - sqrt(5) < 1e-6);

    //run assembly and test Wachspress
    {
        auto mesh = initTestMesh(1);
        auto mpm = initTestMPM(mesh);
        
        auto meshFromMPM = mpm.getMesh();
        PMT_ALWAYS_ASSERT(meshFromMPM.getNumVertices() == 19);
        PMT_ALWAYS_ASSERT(meshFromMPM.getNumElements() == 10);
        
        polyMpmTest::assembly(mpm);
        interpolateWachspress(mpm);              
    }
    Kokkos::finalize();

    return 0;
}
