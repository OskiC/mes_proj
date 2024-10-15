#include "structures\GlobalData.h"
#include "structures\grid.h"
#include "helpers\gauss.h"

int main() {
    std::string fName = "../data/test14.txt";

    /*oc::GlobalData globalData{};
    globalData.parseFromFile(fName);

    oc::Grid grid;
    grid.parseNodes(fName);
    grid.parseElements(fName);

    globalData.printData();

    std::cout << "\n\n\n\n\n";

    grid.printNodes();
    grid.printElemets();
    */
    for(int i = 1; i <= 3; i++)
        std::cout << oc::gaussQuadrature1d(i, oc::f1d) << std::endl;

    for(int i = 1; i <= 3; i++)
        std::cout << oc::gaussQuadrature2d(i, oc::f2d) << std::endl;

    return 0;
}
