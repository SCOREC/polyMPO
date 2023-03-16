

#include "MPM.hpp"
#include "mesh.hpp"
#include "materialPoints.hpp"
#include "utils.hpp"
#include "assembly.hpp"
#include "testUtils.hpp"


int main() {

//  test Vector2
    auto v = polyMpmTest::Vector2();
    v[0] = 1;
    v[1] = 2;
    assert(v[0] == 1);
    assert(v[1] == 2);
//  test Vector operator
    auto v1 = polyMpmTest::Vector2(1,2);
    auto v2 = polyMpmTest::Vector2(3,4);
    auto v3 = v1 + v2;
    assert(v3[0] == 4);
    assert(v3[1] == 6);
    auto v4 = v1 - v2;
    assert(v4[0] == -2);
    assert(v4[1] == -2);
    auto v5 = v1 * 2;
    assert(v5[0] == 2);
    assert(v5[1] == 4);
    auto v6 = -v1;
    assert(v6[0] == -1);
    assert(v6[1] == -2);
    auto v7 = v1.dot(v2);
    assert(v7 == 11);
    auto v8 = v1.cross(v2);
    assert(v8 == -2);
    auto v9 = v1.magnitude();
    assert(v9 - sqrt(5) < 1e-6);  


    return 0;
}
