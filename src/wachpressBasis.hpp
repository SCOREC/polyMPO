#ifndef POLYMPMTEST_WACHPRESS_BASIS_H
#define POLYMPMTEST_WACHPRESS_BASIS_H

#include "utils.hpp"
#include "MPM.hpp"

KOKKOS_INLINE_FUNCTION
void getBasisAndGradByAreaGblForm(Vector2 MP,
                                  int numVtxs,
                                  Vector2* vtxCoords,
                                  double* basis,
                                  double* gradBasis){
    //TODO: cp from ByArea
    Vector2 e[maxVerti + 1];
    Vector2 p[maxVerti];
    double w[maxVerti];
    for (int i = 0; i < numVerti; i++){
        e[i + 1] = v[i + 1] - v[i];
        p[i] = v[i] - xp;
    }
    e[0] = e[numVerti];

    double c[maxVerti];
    double a[maxVerti];
    for (int i = 0; i < numVerti; i++){
        c[i] = e[i].cross(e[i + 1]);
        a[i] = p[i].cross(e[i + 1]);
    }
    double wSum = 0.0;

    double wdx[maxVerti];
    double wdy[maxVerti];
    initArrayWith(wdx, maxVerti, 0.0);
    initArrayWith(wdy, maxVerti, 0.0);
    double wdxSum = 0.0;
    double wdySum = 0.0;
    for (int i = 0; i < numVerti; i++){
        double aProduct = 1.0;
        for (int j = 0; j < numVerti - 2; j++){
            int index1 = (j + i + 1) % numVerti;
            aProduct *= a[index1];

            double productX = 1.0;
            double productY = 1.0;
            for (int k = 0; k < j; k++){
                int index2 = (i + k + 1) % numVerti;
                productX *= a[index2];
                productY *= a[index2];
            }
            productX *= -(v[index1 + 1][1] - v[index1][1]);
            productY *= v[index1 + 1][0] - v[index1][0];
            for (int k = j + 1; k < numVerti - 2; k++){
                int index2 = (i + k + 1) % numVerti;
                productX *= a[index2];
                productY *= a[index2];
            }
            wdx[i] += productX;
            wdy[i] += productY;
        }
        wdx[i] *= c[i];
        wdy[i] *= c[i];
        wdxSum += wdx[i];
        wdySum += wdy[i];
        w[i] = c[i] * aProduct;
        wSum += w[i];
    }

    double wSumInv = 1.0 / wSum;
    for (int i = 0; i < numVerti; i++){
        phi[i] = w[i] * wSumInv;
        gradientPhi[i] = Vector2(wdx[i] * wSumInv - w[i] * wSumInv * wSumInv * wdxSum, wdy[i] * wSumInv - w[i] * wSumInv * wSumInv * wdySum);
    }
}


#endif
