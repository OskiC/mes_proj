//
// Created by xoska on 23.10.2024.
//

#include "jakobian.h"

namespace oc {
    double Jakobian::getDet(){
        return detJ;
    }

    void Jakobian::calcDetJ() {
        detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    }

    void Jakobian::calcJakob(ElemUniv& elemUniv, double x[4], double y[4], int pointIndex) {
        J[0][0] = 0.0;
        J[0][1] = 0.0;
        J[1][0] = 0.0;
        J[1][1] = 0.0;

        for (int i = 0; i < 4; i++) {
            J[0][0] += elemUniv.dN_dXi[pointIndex][i] * x[i];  // dx/dξ
            J[0][1] += elemUniv.dN_dEta[pointIndex][i] * x[i]; // dx/dη
            J[1][0] += elemUniv.dN_dXi[pointIndex][i] * y[i];  // dy/dξ
            J[1][1] += elemUniv.dN_dEta[pointIndex][i] * y[i]; // dy/dη
        }
    }

    void Jakobian::calc_dN_dX_dN_dY(oc::ElemUniv &elemUniv, int pointIndex) {
        dN_dX.resize(4, std::vector<double>(4, 0.0));
        dN_dY.resize(4, std::vector<double>(4, 0.0));
        for (int i = 0; i < 4; i++) {
            dN_dX[pointIndex][i] = J1[0][0] * elemUniv.dN_dXi[pointIndex][i] + J1[0][1] * elemUniv.dN_dEta[pointIndex][i];
            dN_dY[pointIndex][i] = J1[1][0] * elemUniv.dN_dXi[pointIndex][i] + J1[1][1] * elemUniv.dN_dEta[pointIndex][i];
        }
    }

    void Jakobian::calcJakobInver(ElemUniv& elemUniv, int pointIndex) {
        if (detJ != 0) {
            double invDetJ = 1.0 / detJ;

            J1[0][0] = invDetJ * J[1][1];
            J1[0][1] = -invDetJ * J[0][1];
            J1[1][0] = -invDetJ * J[1][0];
            J1[1][1] = invDetJ * J[0][0];

            calc_dN_dX_dN_dY(elemUniv, pointIndex);
        }
    }



    void Jakobian::printJakob() {
        for(auto & i : J){
            for(double j : i) {
                std::cout << j << " ";
            }
            std::cout << "\n";
        }

        std::cout << "\nMacierz J^-1\n";
        for(auto & i : J1){
            for(double j : i){
                std::cout << j << " ";
            }
            std::cout << "\n";
        }
    }
} // oc