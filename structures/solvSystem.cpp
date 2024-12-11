//
// Created by xoska on 20.11.2024.
//

#include "solvSystem.h"

namespace oc {
    SolvSystem::SolvSystem(int numNodes) {
        globalMatrix.resize(numNodes, std::vector<double>(numNodes, 0.0));
    }

    void SolvSystem::addToGlobalMatrix(int i, int j, double value) {
        globalMatrix[i][j] += value;
    }

    void SolvSystem::printMatrix() {
        for (const auto& row : globalMatrix) {
            for (double val : row) {
                std::cout << val << " ";
            }
            std::cout << "\n";
        }
    }


} // oc