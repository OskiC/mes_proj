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
        std::vector<std::vector<double>> globalC;

    public:
        SolvSystem(int numNodes);
        void addToGlobalMatrix(int i, int j, double value);
        void addHbcToGlobalMatrix(const std::vector<std::vector<double>>& Hbc_global);
        void addToGlobalC(const std::vector<std::vector<double>>& localC, const int IDs[]);
        void addToGlobalVector(const std::vector<double>& localP, const int IDs[]);

        void printMatrix();
        void printC();
        void printVector();

        void resetGlobalSystem();
        void solveSystem(const std::vector<std::vector<double>>& H_plus_C_over_Dt, const std::vector<double>& rightHandSide, std::vector<double>& solution);
        void solve(const std::vector<double>& t0, double timeStep, std::vector<double>& t1);

        const std::vector<std::vector<double>> &getGlobalMatrix() const;
        const std::vector<double> &getGlobalVector() const;
        const std::vector<std::vector<double>> &getGlobalC() const;
    };

} // oc

#endif //MES_PROJ_SOLVSYSTEM_H
