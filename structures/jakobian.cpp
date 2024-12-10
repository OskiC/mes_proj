//
// Created by xoska on 23.10.2024.
//

#include "jakobian.h"

namespace oc {
    Jakobian::Jakobian(int num) : numPoints(num) {
        dN_dX.resize(numPoints, std::vector<double>(numPoints, 0.0));
        dN_dY.resize(numPoints, std::vector<double>(numPoints, 0.0));
    }

    double Jakobian::getDet(){
        return detJ;
    }

    void Jakobian::calcDetJ() {
        detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    }

    void Jakobian::calcJakob(ElemUniv& elemUniv, double x[], double y[], int pointIndex) {
        J[0][0] = 0.0;
        J[0][1] = 0.0;
        J[1][0] = 0.0;
        J[1][1] = 0.0;

        for (int i = 0; i < numPoints; i++) {
            J[0][0] += elemUniv.dN_dXi[pointIndex][i] * x[i];  // dx/dξ
            J[1][0] += elemUniv.dN_dEta[pointIndex][i] * x[i]; // dx/dη
            J[0][1] += elemUniv.dN_dXi[pointIndex][i] * y[i];  // dy/dξ
            J[1][1] += elemUniv.dN_dEta[pointIndex][i] * y[i]; // dy/dη
        }
    }

    void Jakobian::calc_dN_dX_dN_dY(ElemUniv& elemUniv, int pointIndex) {
        dN_dX.resize(numPoints, std::vector<double>(numPoints, 0.0));
        dN_dY.resize(numPoints, std::vector<double>(numPoints, 0.0));

        for (int i = 0; i < numPoints; i++) {
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

    void Jakobian::computeHpc(double k, int pointIndex) {
        std::vector<std::vector<double>> Hpc(numPoints, std::vector<double>(numPoints, 0.0));

        for (int i = 0; i < numPoints; i++) {
            for (int j = 0; j < numPoints; j++) {
                Hpc[i][j] = (dN_dX[pointIndex][i] * dN_dX[pointIndex][j] + dN_dY[pointIndex][i] * dN_dY[pointIndex][j]) * k * detJ;
            }
        }

        if (Hpc_list.size() <= pointIndex) {
            Hpc_list.resize(pointIndex + 1);
        }
        Hpc_list[pointIndex] = Hpc;
    }

    std::vector<std::vector<double>> Jakobian::computeTotalH(double k, double dV) {
        // Initialize the weights vector based on the number of points
        std::vector<double> weights;
        if (numPoints == 4) {
            weights = {1.0, 1.0, 1.0, 1.0}; // For 2x2 quadrature
        } else if (numPoints == 9) {
            weights = {5.0/9.0 * 5.0/9.0, 8.0/9.0 * 5.0/9.0, 5.0/9.0 * 5.0/9.0,
                       8.0/9.0 * 5.0/9.0, 8.0/9.0 * 8.0/9.0, 8.0/9.0 * 5.0/9.0,
                       5.0/9.0 * 5.0/9.0, 8.0/9.0 * 5.0/9.0, 5.0/9.0 * 5.0/9.0}; // For 3x3 quadrature
        } else if (numPoints == 16) {
            double w1 = (18 - sqrt(30)) / 36.0;
            double w2 = (18 + sqrt(30)) / 36.0;

            weights = {
                    w1 * w1, w1 * w2, w1 * w2, w1 * w1,
                    w2 * w1, w2 * w2, w2 * w2, w2 * w1,
                    w2 * w1, w2 * w2, w2 * w2, w2 * w1,
                    w1 * w1, w1 * w2, w1 * w2, w1 * w1
            }; // For 4x4 quadrature
        }

        // Initialize the total H matrix
        std::vector<std::vector<double>> H(numPoints, std::vector<double>(numPoints, 0.0));

        // Sum contributions from each integration point
        for (int pointIndex = 0; pointIndex < numPoints; pointIndex++) {
            const auto& Hpc = Hpc_list[pointIndex];
            double weightFactor = weights[pointIndex];

            for (int i = 0; i < numPoints; i++) {
                for (int j = 0; j < numPoints; j++) {
                    H[i][j] += Hpc[i][j] * weightFactor * dV;
                }
            }

            /*std::cout << "Macierz H dla punktu calkowania " << pointIndex + 1 << ":\n";
            for (const auto& row : H) {
                for (double val : row) {
                    if (val != 0.0) {
                        std::cout << val << " ";
                    } else {
                        std::cout << "    "; // Adds spaces for better alignment
                    }
                }
                std::cout << "\n";
            }*/
        }

        return H;
    }

    void Jakobian::printJakob() {
        std::cout << "Matrix J:\n";
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                std::cout << J[i][j] << " ";
            }
            std::cout << "\n";
        }

        std::cout << "\nMatrix J^-1:\n";
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                std::cout << J1[i][j] << " ";
            }
            std::cout << "\n";
        }
    }

}

// oc

    /*void Jakobian::computeH(double k, double dV) {
        H.resize(numPoints, std::vector<double>(numPoints, 0.0)); // Initialize H matrix

        std::vector<double> weights;
        if (numPoints == 4) {
            weights = {1.0, 1.0, 1.0, 1.0}; // Uniform for 2x2 quadrature
        } else if (numPoints == 9) {
            weights = {5.0/9.0 * 5.0/9.0, 8.0/9.0 * 5.0/9.0, 5.0/9.0 * 5.0/9.0,
                       8.0/9.0 * 5.0/9.0, 8.0/9.0 * 8.0/9.0, 8.0/9.0 * 5.0/9.0,
                       5.0/9.0 * 5.0/9.0, 8.0/9.0 * 5.0/9.0, 5.0/9.0 * 5.0/9.0}; // 3x3 quadrature
        } else if (numPoints == 16) {
            weights = {
                    0.3478548 * 0.3478548, 0.3478548 * 0.6521452, 0.3478548 * 0.6521452, 0.3478548 * 0.3478548,
                    0.6521452 * 0.3478548, 0.6521452 * 0.6521452, 0.6521452 * 0.6521452, 0.6521452 * 0.3478548,
                    0.6521452 * 0.3478548, 0.6521452 * 0.6521452, 0.6521452 * 0.6521452, 0.6521452 * 0.3478548,
                    0.3478548 * 0.3478548, 0.3478548 * 0.6521452, 0.3478548 * 0.6521452, 0.3478548 * 0.3478548
            }; // 4x4 quadrature
        }

        for (int pointIndex = 0; pointIndex < numPoints; pointIndex++) {
            std::vector<std::vector<double>> Hpc(numPoints, std::vector<double>(numPoints, 0.0));

            for (int i = 0; i < numPoints; i++) {
                for (int j = 0; j < numPoints; j++) {
                    Hpc[i][j] = (dN_dX[pointIndex][i] * dN_dX[pointIndex][j] +
                                 dN_dY[pointIndex][i] * dN_dY[pointIndex][j]) * getDet();
                }
            }

            double weightFactor = weights[pointIndex];
            for (int i = 0; i < numPoints; i++) {
                for (int j = 0; j < numPoints; j++) {
                    H[i][j] += Hpc[i][j] * weightFactor * dV;
                }
            }

            std::cout << "H matrix contribution at integration point " << pointIndex << ":\n";
            for (const auto &row : Hpc) {
                for (double val : row) {
                    std::cout << val << " ";
                }
                std::cout << "\n";
            }
            std::cout << "\n";
        }
    }
     */
