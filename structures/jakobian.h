#ifndef MES_PROJ_JAKOBIAN_H
#define MES_PROJ_JAKOBIAN_H

#include "elemuniv.h"
#include "..\structures\GlobalData.h"
#include <iomanip>

namespace oc {

    class Jakobian{
    private:
        double J[2][2];
        double J1[2][2]; // J^-1
        double detJ;
        int numPoints;
        std::vector<std::vector<std::vector<double>>> Hpc_list;


    public:
        std::vector<std::vector<double>> dN_dX;
        std::vector<std::vector<double>> dN_dY;
        Jakobian() = default;
        Jakobian(int num);

        void calcJakob(ElemUniv& elemUniv, double x[4], double y[4], int pointIndex);
        void calcJakobInver(ElemUniv& elemUniv, int pointIndex);
        void calcDetJ();
        double getDet();
        void calc_dN_dX_dN_dY(ElemUniv& elemUniv, int pointIndex);
        void printJakob();

        void computeHpc(double k, int pointIndex);
        std::vector<std::vector<double>> computeTotalH(double k, double dV);

        void printH();
    };

} // oc

#endif //MES_PROJ_JAKOBIAN_H
