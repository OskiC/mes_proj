#include "structures\GlobalData.h"
#include "structures\grid.h"
#include <filesystem>

int main() {
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

    std::cout << std::filesystem::current_path();
    return 0;
}
