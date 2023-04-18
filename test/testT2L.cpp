#include "pmt_MPM.hpp"

#include "testUtils.hpp"

int main(int argc, char* argv[]) {
    PMT_ALWAYS_ASSERT(argc == 2);
    Kokkos::initialize();
    int factor = atoi(argv[1]);
    //test T2L
    printf("TestvT2L with a factor of: %d\n",factor);
    {
        //auto mesh = Mesh::readMPASMesh("./mesh.QU.1920km.151026.nc");
        auto mesh = Mesh::readMPASMesh("./GIS.nc");
        
        //6122
        auto mpm = initMPMWithRandomMPs(mesh,factor); 
        //10,20,40,80,160,320,640
        //TODO: 10,40,160,640
        
        //Vector2View dx = initT2LDelta(mpm.getMPs().getCount(),4000,12345);
        //test T2L in materialPoints
        //mpm.T2LTracking(dx);

        //calcAvgLengthOfEdge(mesh,1); 
        //for(int i=0; i<100; i++){ 
        //    Vector2View dx = initT2LDeltaRankineVortex(mpm, Vector2(150000, -2000000), 15, 11500, 1);
        //    mpm.T2LTracking(dx,-1);
        //}
        
        //Test init with fractions:
        double p0 = 0.99;
        double p1 = 0.01;
        double p2 = 0.0;
        double p3 = 0.0;
        printf("\tfraction: %.2f %.2f %.2f %.2f\n",p0,p1,p2,p3);
//TODO: 1.00 0.00 0.00 0.00
//      1.00 1MP  0.00 0.00
//      1.00 0.00 1MP  0.00
//      1.00 0.00 0.00 1MP
//      0.90 0.10 0.00 0.00
//      0.90 0.00 0.10 0.00
//      0.90 0.00 0.00 0.10
//      0.90 0.05 0.05 0.00 //look the same??
//
        //Test1
        for(int i=0; i<5; i++){
            Vector2View dx1 = initT2LTest1(mpm,p0,p1,p2,p3);
            mpm.T2LTracking(dx1,-1);
        }
        //Test2
        //Vector2View dx2 = initT2LTest2(mpm.getMPs().getCount(),p0,p1,p2,p3,11500.0);
        //mpm.T2LTracking(dx2,2);
    }    
    Kokkos::finalize();
    return 0;
}
