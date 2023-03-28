#include <Kokkos_Core.hpp>
#include "pmt_utils.hpp"
#include "pmt_MPM.hpp"

namespace polyMpmTest{

void MPM::T2LTracking(Vector2View dx){
    Kokkos::parallel_for("check",dx.size(),KOKKOS_LAMBDA(const int i){
        //printf("%d:(%.3f,%.3f)\n",i,dx(i)[0],dx(i)[1]);
    });   
} 

} 
