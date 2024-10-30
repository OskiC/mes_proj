//
// Created by xoska on 23.10.2024.
//

#ifndef MES_PROJ_ELEMUNIV_H
#define MES_PROJ_ELEMUNIV_H

#include <cmath>
#include <iostream>
#include <vector>

namespace oc {

    struct ElemUniv{
        std::vector<double> ksi;
        std::vector<double> eta;

        std::vector<std::vector<double>> dN_dXi;
        std::vector<std::vector<double>> dN_dEta;

        void initialize(int num);
        void print_dnx();
    };

} // oc

#endif //MES_PROJ_ELEMUNIV_H
