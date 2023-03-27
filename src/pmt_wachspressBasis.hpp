#ifndef POLYMPMTEST_WACHSPRESS_BASIS_H
#define POLYMPMTEST_WACHSPRESS_BASIS_H

#include "pmt_utils.hpp"
#include "pmt_MPM.hpp"

namespace polyMpmTest{

/** \brief calculate the basis and gradient of Basis for a give MP with its element Vtxs
 *
 *  \details based on the 4.1 section from:
 *  https://www.mn.uio.no/math/english/people/aca/michaelf/papers/gbc.pdf
 *
 *  \param MP: single Material Point
 *
 *  \param numVtxs: number of vertices of vtxCoords
 *
 *  \param vtxCoords: vertices of its corresponding element
 *
 *  \param basis && gadBasis: hold the return values
 *
 *  \return basis and gradient of basis
 * */
KOKKOS_INLINE_FUNCTION
void getBasisAndGradByAreaGblForm(Vector2 MP,
                                  int numVtxs,
                                  Vector2* vtxCoords,
                                  double* basis,
                                  Vector2* gradBasis){
    Vector2 e[maxVtxsPerElm + 1];
    Vector2 p[maxVtxsPerElm];
    double w[maxVtxsPerElm];
    for (int i = 0; i < numVtxs; i++){
        e[i + 1] = vtxCoords[i + 1] - vtxCoords[i];
        p[i] = vtxCoords[i] - MP;
    }
    e[0] = e[numVtxs];

    double c[maxVtxsPerElm];
    double a[maxVtxsPerElm];
    for (int i = 0; i < numVtxs; i++){
        c[i] = e[i].cross(e[i + 1]);
        a[i] = p[i].cross(e[i + 1]);
    }
    double wSum = 0.0;

    double wdx[maxVtxsPerElm];
    double wdy[maxVtxsPerElm];
    initArray(wdx, maxVtxsPerElm, 0.0);
    initArray(wdy, maxVtxsPerElm, 0.0);
    double wdxSum = 0.0;
    double wdySum = 0.0;
    for (int i = 0; i < numVtxs; i++){
        double aProduct = 1.0;
        for (int j = 0; j < numVtxs - 2; j++){
            int index1 = (j + i + 1) % numVtxs;
            aProduct *= a[index1];

            double product = 1.0;
            for (int k = 0; k < j; k++){
                int index2 = (i + k + 1) % numVtxs;
                product *= a[index2];
            }
            for (int k = j + 1; k < numVtxs - 2; k++){
                int index2 = (i + k + 1) % numVtxs;
                product *= a[index2];
            }
            wdx[i] += product * -(vtxCoords[index1 + 1][1] - vtxCoords[index1][1]);
            wdy[i] += product * (vtxCoords[index1 + 1][0] - vtxCoords[index1][0]);
        }
        wdx[i] *= c[i];
        wdy[i] *= c[i];
        wdxSum += wdx[i];
        wdySum += wdy[i];
        w[i] = c[i] * aProduct;
        wSum += w[i];
    }

    double wSumInv = 1.0 / wSum;
    double wSumInvSq = wSumInv * wSumInv;
    for (int i = 0; i < numVtxs; i++){
        basis[i] = w[i] * wSumInv;
        gradBasis[i] = Vector2(wdx[i] * wSumInv - w[i] * wSumInvSq * wdxSum, wdy[i] * wSumInv - w[i] * wSumInvSq * wdySum);
    }
}

} //namespace polyMpmTest end
#endif
