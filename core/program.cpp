#include "program.h"

namespace oc {

    void Program::zadanie1(){
        std::string fName = "../data/test14.txt";
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
        for(int i = 1; i <= 3; i++)
            std::cout << oc::gaussQuadrature1d(i, oc::f1d) << std::endl;

        for(int i = 1; i <= 3; i++)
            std::cout << oc::gaussQuadrature2d(i, oc::f2d) << std::endl;

    }

} // oc