#include "elemuniv.h"

namespace oc {
    void ElemUniv::initialize(int num) {
        // Clear any existing values to allow re-initialization
        ksi.clear();
        eta.clear();
        dN_dXi.clear();
        dN_dEta.clear();

        if (num == 4) {
            // 2x2 Gauss Quadrature Points
            ksi = {-1.0 / sqrt(3), 1.0 / sqrt(3), -1.0 / sqrt(3), 1.0 / sqrt(3)};
            eta = {-1.0 / sqrt(3), -1.0 / sqrt(3), 1.0 / sqrt(3), 1.0 / sqrt(3)};
        } else if (num == 9) {
            // 3x3 Gauss Quadrature Points
            double a = sqrt(3.0 / 5.0);
            ksi = {-a, 0, a, -a, 0, a, -a, 0, a};
            eta = {-a, -a, -a, 0, 0, 0, a, a, a};
        } else if (num == 16) {
            // 4x4 Gauss Quadrature Points
            double a = sqrt((3 + 2 * sqrt(6.0 / 5.0)) / 7.0);
            double b = sqrt((3 - 2 * sqrt(6.0 / 5.0)) / 7.0);
            ksi = {-a, -b, b, a, -a, -b, b, a, -a, -b, b, a, -a, -b, b, a};
            eta = {-a, -a, -a, -a, -b, -b, -b, -b, b, b, b, b, a, a, a, a};
        } else {
            std::cerr << "Unsupported number of points for Gauss quadrature." << std::endl;
            return;
        }

        // Resize and calculate dN_dXi and dN_dEta for each point
        dN_dXi.resize(num);
        dN_dEta.resize(num);

        for (int i = 0; i < num; i++) {
            std::vector<double> dXi_row(num, 0.0);
            std::vector<double> dEta_row(num, 0.0);

            if (num == 4 || num == 9 || num == 16) {
                // Compute derivatives based on the 4-node element for each integration point
                dXi_row[0] = -0.25 * (1 - eta[i]);
                dXi_row[1] = 0.25 * (1 - eta[i]);
                dXi_row[2] = 0.25 * (1 + eta[i]);
                dXi_row[3] = -0.25 * (1 + eta[i]);

                dEta_row[0] = -0.25 * (1 - ksi[i]);
                dEta_row[1] = -0.25 * (1 + ksi[i]);
                dEta_row[2] = 0.25 * (1 + ksi[i]);
                dEta_row[3] = 0.25 * (1 - ksi[i]);
            }

            dN_dXi[i] = dXi_row;
            dN_dEta[i] = dEta_row;
        }
    }
    void ElemUniv::print_dnx(){
        std::cout << "dN_dXi and dN_dEta values for each integration point:\n";
        for (int pointIndex = 0; pointIndex < dN_dXi.size(); ++pointIndex) {
            std::cout << "Integration Point " << pointIndex + 1 << ":\n";
            std::cout << "dN_dXi: ";
            for (double val : dN_dXi[pointIndex]) {
                std::cout << val << " ";
            }
            std::cout << "\ndN_dEta: ";
            for (double val : dN_dEta[pointIndex]) {
                std::cout << val << " ";
            }
            std::cout << "\n\n";
        }
    }
}
