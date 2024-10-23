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
        }

        for (int i = 0; i < num; i++) {
            std::vector<double> dXi_row(4, 0.0);
            std::vector<double> dEta_row(4, 0.0);

            dXi_row[0] = -0.25 * (1 - eta[i]);
            dXi_row[1] =  0.25 * (1 - eta[i]);
            dXi_row[2] =  0.25 * (1 + eta[i]);
            dXi_row[3] = -0.25 * (1 + eta[i]);

            dEta_row[0] = -0.25 * (1 - ksi[i]);
            dEta_row[1] = -0.25 * (1 + ksi[i]);
            dEta_row[2] =  0.25 * (1 + ksi[i]);
            dEta_row[3] =  0.25 * (1 - ksi[i]);

            dN_dXi.push_back(dXi_row);
            dN_dEta.push_back(dEta_row);
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

        std::cout << "\nMacierz J^-1\n";
        for(auto & i : J1){
            for(double j : i){
                std::cout << j << " ";
            }
            std::cout << "\n";
        }
    }
} // oc