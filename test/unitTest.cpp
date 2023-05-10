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
        auto mpPerElement = std::vector<int>({5,4,5,6,6,5,4,6,5,5});
        auto mpm = initTestMPM(mesh, mpPerElement); //creates test MPs
        
        auto meshFromMPM = mpm.getMesh();
        PMT_ALWAYS_ASSERT(meshFromMPM.getNumVertices() == 19);
        PMT_ALWAYS_ASSERT(meshFromMPM.getNumElements() == 10);
        
        //run non-physical assembly (mp -to- mesh vertex) kernel
        auto vtxField = polyMpmTest::assembly(mpm);
        //check the result
        auto vtxField_h = Kokkos::create_mirror_view(vtxField);
        Kokkos::deep_copy(vtxField_h, vtxField);
        const std::vector<double> vtxFieldExpected = {
          1.768750, 4.528750, 17.660000, 8.228750,
          8.406250, 2.818750, 4.978750, 5.708750,
          6.486250, 10.551786, 13.014286, 15.114286,
          9.714286, 8.631786, 4.990000, 1.050000,
          2.830000, 5.694286, 3.914286
        };
        PMT_ALWAYS_ASSERT(vtxField_h.size() == vtxFieldExpected.size());
        for(size_t i=0; i<vtxField_h.size(); i++) {
          auto res = polyMpmTest::isEqual(vtxField_h(i),vtxFieldExpected[i], 1e-6);
          if(!res) {
            fprintf(stderr, "computed value for vtx %ld, %.6f, does not match expected value %.6f\n",
                    i, vtxField_h(i), vtxFieldExpected[i]);
          }
          PMT_ALWAYS_ASSERT(res);
        }


        interpolateWachspress(mpm);              
    }
    Kokkos::finalize();

    return 0;
}
