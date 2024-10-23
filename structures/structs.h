//
// Created by xoska on 12.10.2024.
//

#ifndef MES_PROJ_STRUCTS_H
#define MES_PROJ_STRUCTS_H

#include <vector>
#include <cmath>
#include <iostream>

namespace oc {
    struct ElemUniv{
        std::vector<double> ksi;
        std::vector<double> eta;

        std::vector<std::vector<double>> dN_dXi;
        std::vector<std::vector<double>> dN_dEta;

        void initialize(int num);
    };

    class Jakobian{
    private:
        double J[2][2];
        double J1[2][2]; // J^-1
        double detJ;
        
    public:
        void calcJakob(ElemUniv& elemUniv, double x[4], double y[4], int pointIndex);
        void calcJakobInver();
        void calcDetJ();
        double getDet();
        void printJakob();
    };

    struct Node {
        double x, y;
    };

    struct Element {
        int ID[4];
        Jakobian jakobian;
    };


} // oc

#endif //MES_PROJ_STRUCTS_H
