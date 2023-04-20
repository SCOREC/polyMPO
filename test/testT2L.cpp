#include "pmt_MPM.hpp"

#include "testUtils.hpp"

int main(int argc, char* argv[]) {
    PMT_ALWAYS_ASSERT(argc == 7);
    Kokkos::initialize();
    int factor = atoi(argv[1]);
    double p0 = atof(argv[2]);
    double p1 = atof(argv[3]);
    double p2 = atof(argv[4]);
    double p3 = atof(argv[5]);
    int MPAcross = atoi(argv[6]);
    //test T2L
    printf("TestT2L with a factor of: %d\n",factor);
    {
        //auto mesh = Mesh::readMPASMesh("./mesh.QU.1920km.151026.nc");
        auto mesh = Mesh::readMPASMesh("./GIS.nc");//numElm == 6122
        
        auto mpm = initMPMWithRandomMPs(mesh,factor); 
        //10,20,40,80,160,320,640
        //TODO: 10,40,160,640
        
        //calcAvgLengthOfEdge(mesh,1); 
        //runT2LSimple(mpm);

        //Test Rankine
        //runT2LRankineVortex(mpm, Vector2(150000,-2000000), 15, 11500, 1, 100, -1);

        //Test1
        //Test init with fractions:
        //TODO: 1.00 0.00 0.00 0.00
        //      0.90 0.10 0.00 0.00
        //      0.90 0.00 0.10 0.00
        //      0.90 0.00 0.00 0.10
        //      0.90 0.05 0.05 0.00 //look the same??
        if (MPAcross == 0){
            for(int i=0; i<5; i++)
                runT2LRandomWithProportion(mpm, p0, p1, p2, p3, 1, -1);
        }else{
        //Test2
        //TODO: 1.00 1MP  0.00 0.00
        //      1.00 0.00 1MP  0.00
        //      1.00 0.00 0.00 1MP
            for(int i=0; i<5; i++)
                runT2LWithOneMPAcross(mpm, MPAcross, 1, -1);
        }
    }    
    Kokkos::finalize();
    return 0;
}
