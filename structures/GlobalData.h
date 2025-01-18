#ifndef MES_PROJ_GLOBALDATA_H
#define MES_PROJ_GLOBALDATA_H

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

namespace oc {

    class GlobalData {
    private:
        int simulationTime;
        int simulationStepTime;
    public:
        int getSimulationStepTime() const;

    private:
        int conductivity;
        int alfa;
        int tot;
        int initialTemp;
        int density;
        int specificHeat;
        int nodesNumber;
        int elementsNumber;
        std::vector<int> bc;

    public:
        GlobalData() = default;

        void parseFromFile(const std::string& fileName);
        void printData();

        int getSimulationTime() const;

        int getConductivity() const;

        int getAlfa() const;

        int getTot() const;

        int getInitialTemp() const;

        int getDensity() const;

        int getSpecificHeat() const;

        int getNodesNumber() const;

        int getElementsNumber() const;

        std::vector<int>& getBoundaryConditions();


    };

} // oc

#endif //MES_PROJ_GLOBALDATA_H
