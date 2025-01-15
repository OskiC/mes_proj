//
// Created by xoska on 20.11.2024.
//

#include "solvSystem.h"

namespace oc {
    SolvSystem::SolvSystem(int numNodes) {
        globalMatrix.resize(numNodes, std::vector<double>(numNodes, 0.0));
        globalVector.resize(numNodes, 0.0);
        globalC.resize(numNodes, std::vector<double>(numNodes, 0.0));
    }

    void SolvSystem::addToGlobalMatrix(int i, int j, double value) {
        globalMatrix[i][j] += value;
    }

    void SolvSystem::addToGlobalC(const double localC[4][4], const int IDs[]) {
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                int global_i = IDs[i] - 1; // Map local row index to global row index
                int global_j = IDs[j] - 1; // Map local column index to global column index
                globalC[global_i][global_j] += localC[i][j]; // Add value to the correct position
            }
        }
    }

    void SolvSystem::printMatrix() {
        for (const auto& row : globalMatrix) {
            for (double val : row) {
                std::cout << val << " ";
            }
            std::cout << "\n";
        }
    }

    void SolvSystem::printC() {
        for (const auto& row : globalC) {
            for (double val : row) {
                std::cout << val << " ";
            }
            std::cout << "\n";
        }
    }

    void SolvSystem::addToGlobalVector(const std::vector<double>& localP, const int IDs[]) {
        for (int i = 0; i < 4; ++i) { // Loop through the 4 node IDs
            int globalIndex = IDs[i] - 1; // Convert 1-based ID to 0-based index
            globalVector[globalIndex] += localP[i];
        }
    }

    void SolvSystem::printVector(){
        std::cout << "Global Vector: [";
        for(auto val : globalVector){
            std::cout << val << ", ";
        }
        std::cout << "]\n";
    }

    void SolvSystem::resetGlobalSystem() {
        std::fill(globalMatrix.begin(), globalMatrix.end(), std::vector<double>(globalMatrix.size(), 0.0));
        std::fill(globalVector.begin(), globalVector.end(), 0.0);
        std::fill(globalC.begin(), globalC.end(), std::vector<double>(globalC.size(), 0.0));
    }

    void SolvSystem::solveSystem(const std::vector<std::vector<double>>& H_plus_C_over_Dt, const std::vector<double>& rightHandSide, std::vector<double>& solution) {
        int size = H_plus_C_over_Dt.size();

        // Check if dimensions match
        if (size != rightHandSide.size() || size != solution.size()) {
            throw std::runtime_error("Matrix and vector dimensions do not match");
        }

        // Solve the system using Gaussian elimination with partial pivoting

        // Copy the matrix to not modify the original
        std::vector<std::vector<double>> A = H_plus_C_over_Dt;
        std::vector<double> b = rightHandSide;

        // Forward elimination
        for (int i = 0; i < size; i++) {
            // Find pivot
            int maxRow = i;
            for (int k = i + 1; k < size; k++) {
                if (std::abs(A[k][i]) > std::abs(A[maxRow][i])) {
                    maxRow = k;
                }
            }

            // Swap maximum row with current row
            if (maxRow != i) {
                std::swap(A[i], A[maxRow]);
                std::swap(b[i], b[maxRow]);
            }

            // Make all rows below this one 0 in current column
            for (int k = i + 1; k < size; k++) {
                double factor = A[k][i] / A[i][i];
                for (int j = i; j < size; j++) {
                    A[k][j] -= factor * A[i][j];
                }
                b[k] -= factor * b[i];
            }
        }

        // Back substitution
        for (int i = size - 1; i >= 0; i--) {
            solution[i] = b[i];
            for (int j = i + 1; j < size; j++) {
                solution[i] -= A[i][j] * solution[j];
            }
            solution[i] /= A[i][i];
        }
    }

    void SolvSystem::solve(const std::vector<double>& t0, double timeStep, std::vector<double>& t1) {
        std::vector<std::vector<double>> H_plus_C_over_Dt = this->getGlobalMatrix();
        std::vector<std::vector<double>> C_over_Dt = this->getGlobalC();
        std::vector<double> P_global = this->getGlobalVector();

        int size = H_plus_C_over_Dt.size();

        // Form the left-hand side matrix [H] + [C]/Δτ
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                H_plus_C_over_Dt[i][j] += C_over_Dt[i][j] / timeStep;
            }
        }

        // Form the right-hand side vector ([C]/Δτ) * {t_0} + {P}
        std::vector<double> rightHandSide(size, 0.0);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                rightHandSide[i] += (C_over_Dt[i][j] / timeStep) * t0[j];
            }
            rightHandSide[i] += P_global[i];
        }

        // Solve the system
        this->solveSystem(H_plus_C_over_Dt, rightHandSide, t1);
    }

    const std::vector<std::vector<double>> &SolvSystem::getGlobalMatrix() const {
        return globalMatrix;
    }

    const std::vector<double> &SolvSystem::getGlobalVector() const {
        return globalVector;
    }

    const std::vector<std::vector<double>> &SolvSystem::getGlobalC() const {
        return globalC;
    }
} // oc