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

} // oc