//
// Created by xoska on 20.11.2024.
//

#ifndef MES_PROJ_SOLVSYSTEM_H
#define MES_PROJ_SOLVSYSTEM_H

#include <vector>
#include <iostream>


namespace oc {

    class SolvSystem {
    private:
        std::vector<std::vector<double>> globalMatrix;
        std::vector<double> globalVector;

    public:
        SolvSystem(int numNodes);
        void assembleLocalToGlobal(const std::vector<std::vector<double>>& localH, const std::vector<int>& elementNodes);
        void printMatrix();
        void solve();
    };

} // oc

#endif //MES_PROJ_SOLVSYSTEM_H
