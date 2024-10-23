//
// Created by xoska on 12.10.2024.
//

#include "structs.h"

namespace oc {
    double Jakobian::getDet(){
        return detJ;
    }

    void ElemUniv::initialize(int num) {
        if (num == 4) {
            ksi = {-1.0 / sqrt(3), 1.0 / sqrt(3), 1.0 / sqrt(3), -1.0 / sqrt(3)};
            eta = {-1.0 / sqrt(3), -1.0 / sqrt(3), 1.0 / sqrt(3), 1.0 / sqrt(3)};
            dN_dXi.resize(4, std::vector<double>(4));
            dN_dEta.resize(4, std::vector<double>(4));
        }

        for (int i = 0; i < num; i++) {
            dN_dXi[i][0] = -0.25 * (1 - eta[i]);
            dN_dXi[i][1] = 0.25 * (1 - eta[i]);
            dN_dXi[i][2] = 0.25 * (1 + eta[i]);
            dN_dXi[i][3] = -0.25 * (1 + eta[i]);

            dN_dEta[i][0] = -0.25 * (1 - ksi[i]);
            dN_dEta[i][1] = -0.25 * (1 + ksi[i]);
            dN_dEta[i][2] = 0.25 * (1 + ksi[i]);
            dN_dEta[i][3] = 0.25 * (1 - ksi[i]);
        }
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

    void Jakobian::calcJakobInver() {
        if (detJ != 0) {
            double invDetJ = 1.0 / detJ;

            J1[0][0] = invDetJ * J[1][1];
            J1[0][1] = -invDetJ * J[0][1];
            J1[1][0] = -invDetJ * J[1][0];
            J1[1][1] = invDetJ * J[0][0];
        }
    }

    void Jakobian::printJakob() {
        for(auto & i : J){
            for(double j : i) {
                std::cout << j << " ";
            }
            std::cout << "\n";
        }
    }
} // oc