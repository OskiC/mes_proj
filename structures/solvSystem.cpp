//
// Created by xoska on 20.11.2024.
//

#include "solvSystem.h"

namespace oc {
    SolvSystem::SolvSystem(int numNodes) {
        globalMatrix.resize(numNodes, std::vector<double>(numNodes, 0.0));
    }

    void SolvSystem::assembleLocalToGlobal(const std::vector<std::vector<double>>& localH, const std::vector<int>& elementNodes){
        for (size_t i = 0; i < elementNodes.size(); ++i) {
            for (size_t j = 0; j < elementNodes.size(); ++j) {
                int globalRow = elementNodes[i] - 1;
                int globalCol = elementNodes[j] - 1;

                globalMatrix[globalRow][globalCol] += localH[i][j];
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


} // oc