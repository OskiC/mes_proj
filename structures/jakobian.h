//
// Created by xoska on 23.10.2024.
//

#ifndef MES_PROJ_JAKOBIAN_H
#define MES_PROJ_JAKOBIAN_H

#include "elemuniv.h"

namespace oc {

    class Jakobian{
    private:
        double J[2][2];
        double J1[2][2]; // J^-1
        double detJ;


    public:
        std::vector<std::vector<double>> dN_dX;
        std::vector<std::vector<double>> dN_dY;

        void calcJakob(ElemUniv& elemUniv, double x[4], double y[4], int pointIndex);
        void calcJakobInver(ElemUniv& elemUniv, int pointIndex);
        void calcDetJ();
        double getDet();
        void calc_dN_dX_dN_dY(ElemUniv& elemUniv, int pointIndex);
        void printJakob();

        void computeDerivatives(){

        }
    };

} // oc

#endif //MES_PROJ_JAKOBIAN_H
