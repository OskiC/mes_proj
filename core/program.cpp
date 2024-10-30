#include "program.h"

namespace oc {

    void Program::zadanie1(std::string &fName){
        oc::GlobalData globalData{};
        globalData.parseFromFile(fName);

        oc::Grid grid;
        grid.parseNodes(fName);
        grid.parseElements(fName);

        globalData.printData();

        std::cout << "\n\n\n\n\n";

        grid.printNodes();
        grid.printElemets();
    }

    void Program::zadanie2(){
        double correctVal1 = 46.0 / 3.0 ; // 15.333
        double correctVal2 = 236.0 / 9.0; // 26.222

        for(int i = 1; i <= 3; i++) {
            double tempRes = oc::gaussQuadrature1d(i, oc::f1d);
            std::cout << "Wynik dla " << i << " wezlow: " << tempRes
                    << " Roznica z analitycznym: " << std::abs(tempRes - correctVal1) << std::endl;

        }

        for(int i = 1; i <= 3; i++) {
            double tempRes = oc::gaussQuadrature2d(i, oc::f2d);
            std::cout << "Wynik dla " << i << " wezlow: " << tempRes
                      << " Roznica z analitycznym: " << std::abs(tempRes - correctVal2) << std::endl;

        }
    }

    void Program::zadanie3(){
        int numPoints = 4; // Number of Gauss points (for 2D, typically 4 points)

        ElemUniv elemUniv;
        elemUniv.initialize(numPoints); // 2d - 4punkty
        //elemUniv.print_dnx();

        double x[4] = {0.01, 0.025, 0.025, 0};
        double y[4] = {-0.01, 0, 0.025, 0.025};

        double k = 30;   // Example thermal conductivity
        double dV = 1.0; // Example volume element

        std::vector<Jakobian> jakobians(numPoints);
        for(int i = 0; i < numPoints; i++){
            Jakobian jakobian(numPoints);
            jakobians.push_back(jakobian);
        }

        // Loop through each Gauss point to calculate Jacobians and derivatives
        for(int i = 0; i < numPoints; i++) {
            jakobians[i].calcJakob(elemUniv, x, y, i);  // Calculate Jacobian at point i
            jakobians[i].calcDetJ();                    // Calculate determinant of Jacobian
            jakobians[i].calcJakobInver(elemUniv, i);   // Calculate inverse Jacobian

            std::cout << "Wyznacznik Jakobiana w " << i + 1 << ": " << jakobians[i].getDet() << std::endl;
            jakobians[i].printJakob();

                // Print derivatives with respect to x and y
            std::cout << "\nPochodne gaussa w pkt " << i + 1 << ":\n";
            for (int j = 0; j < 4; j++) {
                std::cout << "dN" << j + 1 << "/dx: " << jakobians[i].dN_dX[i][j]
                          << ", dN" << j + 1 << "/dy: " << jakobians[i].dN_dY[i][j] << std::endl;
            }
            std::cout << "\n";
            jakobians[i].computeH(k, dV);
            std::cout << "Macierz H:" << std::endl;
            jakobians[i].printH();
        }


    }

    void Program::zadanie4(){
    }


} // oc