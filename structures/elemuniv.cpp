//
// Created by xoska on 23.10.2024.
//

#include "elemuniv.h"

namespace oc {
    void ElemUniv::initialize(int num) {
        if (num == 4) {
            ksi = {-1.0 / sqrt(3), 1.0 / sqrt(3), -1.0 / sqrt(3), 1.0 / sqrt(3)};
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
} // oc