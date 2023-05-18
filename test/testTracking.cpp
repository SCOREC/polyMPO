#include "pmt_MPM.hpp"
#include "testUtils.hpp"

//TODO: add usage comments
int main(int argc, char* argv[]) {
    PMT_ALWAYS_ASSERT(argc == 7);
    Kokkos::initialize();
    int factor = atoi(argv[1]);
    double p0 = atof(argv[2]);
    double p1 = atof(argv[3]);
    double p2 = atof(argv[4]);
    double p3 = atof(argv[5]);
    //MPAcross control which way we test
    //0:   T2L tracking with p0,p1,p2,p3 across centralized MPs TODO change naming
    //1-5: T2L tracking with one MP across randomMPs with one across
    //6:   CVT tracking with p0,p1,p2,p3 across centralized MPs
    int MPAcross = atoi(argv[6]);
    //testTracking
    printf("Test MPs tracking with a factor of: %d\n",factor);
    {
        //auto mesh = Mesh::readMPASMesh("./mesh.QU.1920km.151026.nc");
        //auto mesh = Mesh::readMPASMesh("./GIS.nc");//numElm == 6122 
        auto mesh = Mesh::readMPASMesh("./grid_mesh_out.nc");//numElm == 512 
        //auto mesh = Mesh::readMPASMesh("./grid_full.nc");//numElm == 
        MPM mpm = initMPMWithCenterMPs(mesh,factor); 
        //10,20,40,80,160,320,640
        //TODO: test 1M and 10M
        //      factor of 163, 1630
        
        //calcAvgLengthOfEdge(mesh,1); 
        //runT2LSimple(mpm);

        //Test Rankine
        //runT2LRankineVortex(mpm, Vector2(150000,-2000000), 15, 11500, 1, 100, -1);

        //test for timing
       // printMPs(mpm);
        if (MPAcross == 0){
            for(int i=0; i<5; i++){
                mpm = initMPMWithCenterMPs(mesh,factor); 
                runT2LRandomWithProportion(mpm, p0, p1, p2, p3, 1, -1);
            }
        }else if(MPAcross < 6){
            for(int i=0; i<5; i++){
                mpm = initMPMWithCenterMPs(mesh,factor); 
                runT2LWithOneMPAcross(mpm, MPAcross, 1, -1);
            }
        }else if(MPAcross == 6){
            for(int i=0; i<5; i++){
                mpm = initMPMWithCenterMPs(mesh,factor); 
                runCVTRandomWithProportion(mpm, p0, p1, p2, p3, 1, -1);
            }
        }else if(MPAcross == 7){
            for(int i=0; i<5; i++){
                mpm = initMPMWithCenterMPs(mesh,factor); 
                runCVTElmCenterBasedRandomWithProportion(mpm, p0, p1, p2, p3, 1, -1);
            }
        }else if(MPAcross == 8){
            for(int i=0; i<5; i++){
                mpm = initMPMWithCenterMPs(mesh,factor); 
                runCVTAllAcrossTest(mpm, 1, 1, -1);
            }
        }else if(MPAcross == 9){
            for(int i=0; i<5; i++){
                mpm = initMPMWithCenterMPs(mesh,factor); 
                runCVTAllAcrossTest(mpm, 2, 1, -1);
            }
        }else if(MPAcross == 10){
            for(int i=0; i<5; i++){
                mpm = initMPMWithCenterMPs(mesh,factor); 
                runCVTAllAcrossTest(mpm, 10, 0.4, 1, -1);
            }
        }else if(MPAcross == 11){
            for(int i=0; i<5; i++){
                mpm = initMPMWithCenterMPs(mesh,factor); 
                runCVTElmAllAcrossTest(mpm, 1, 1.0, 1, -1);
            }
        }else if(MPAcross == 12){
            for(int i=0; i<5; i++){
                mpm = initMPMWithCenterMPs(mesh,factor); 
                runCVTElmAllAcrossTest(mpm, 2, 1.0, 1, -1);
            }
        }else if(MPAcross == 13){
            for(int i=0; i<5; i++){
                mpm = initMPMWithCenterMPs(mesh,factor); 
                runCVTElmAllAcrossTest(mpm, 10, 0.4, 1, -1);
            }
        }else{
            mpm = initMPMWithRandomMPs(mesh,factor);
            auto dx = initDeltaRandomDir(mpm.getMesh().getNumElements()*factor,100000);
            MPM mpmOneMP = initMPMOneMP(mesh);
            auto dxOneMP = initOneDelta();
            if(MPAcross == 14)
                mpm.T2LTracking(dx,-1);
            else if(MPAcross == 15)
                mpm.CVTTrackingEdgeCenterBased(dx,-1);
            else if(MPAcross == 16)
                mpmOneMP.CVTTrackingElmCenterBased(dxOneMP,-1);
            //PMT_ALWAYS_ASSERT(false);
        }
        //printMPs(mpm);
    }    
    Kokkos::finalize();
    return 0;
}
