#ifndef MES_PROJ_GLOBALDATA_H
#define MES_PROJ_GLOBALDATA_H

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

namespace oc {

    class GlobalData {
    private:
        int simulationTime;
        int simulationStepTime;
        int conductivity;
        int alfa;
        int tot;
        int initialTemp;
        int density;
        int specificHeat;
        int nodesNumber;
        int elementsNumber;

    public:
        GlobalData() = default;

        void parseFromFile(const std::string& fileName);
        void printData();
    };

} // oc

#endif //MES_PROJ_GLOBALDATA_H
