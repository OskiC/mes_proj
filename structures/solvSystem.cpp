//
// Created by xoska on 20.11.2024.
//

#include "solvSystem.h"

namespace oc {
    SolvSystem::SolvSystem(int numNodes) {
        globalMatrix.resize(numNodes, std::vector<double>(numNodes, 0.0));
        globalVector.resize(numNodes, 0.0);
    }

    void SolvSystem::addToGlobalMatrix(int i, int j, double value) {
        globalMatrix[i][j] += value;
    }

    void SolvSystem::addHbcToGlobalMatrix(const std::vector<std::vector<double>>& Hbc_global) {
        for (int i = 0; i < globalMatrix.size(); ++i) {
            for (int j = 0; j < globalMatrix[i].size(); ++j) {
                globalMatrix[i][j] += Hbc_global[i][j];
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




} // oc