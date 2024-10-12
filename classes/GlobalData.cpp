#include "GlobalData.h"

namespace oc {
    GlobalData::GlobalData(int simTime, int simStep, int cond, int a, int t, int init, int dens, int specHeat, int nN, int eN) {
        this->simulationTime = simTime;
        simulationStepTime = simStep;
        conductivity = cond;
        alfa = a;
        tot = t;
        initialTemp = init;
        density = dens;
        specificHeat = specHeat;
        nodesNumber = nN;
        elementsNumber = eN;
    }

    void GlobalData::parseFromFile(const std::string& fileName) {
        std::ifstream inputFile(fileName);
        if (!inputFile.is_open()) {
            std::cerr << "Error opening file!" << std::endl;
            return;
        }

        std::string line;
        while (std::getline(inputFile, line)) {
            std::istringstream iss(line);
            std::string key;
            if (line.find("SimulationTime") != std::string::npos) {
                iss >> key >> simulationTime;
            }
            else if (line.find("SimulationStepTime") != std::string::npos) {
                iss >> key >> simulationStepTime;
            }
            else if (line.find("Conductivity") != std::string::npos) {
                iss >> key >> conductivity;
            }
            else if (line.find("Alfa") != std::string::npos) {
                iss >> key >> alfa;
            }
            else if (line.find("Tot") != std::string::npos) {
                iss >> key >> tot;
            }
            else if (line.find("InitialTemp") != std::string::npos) {
                iss >> key >> initialTemp;
            }
            else if (line.find("Density") != std::string::npos) {
                iss >> key >> density;
            }
            else if (line.find("SpecificHeat") != std::string::npos) {
                iss >> key >> specificHeat;
            }
            else if (line.find("Nodesnumber") != std::string::npos) {
                iss >> key >> nodesNumber;
            }
            else if (line.find("Elementsnumber") != std::string::npos) {
                iss >> key >> elementsNumber;
            }
        }

        inputFile.close();
    }

    void GlobalData::printData()
    {
        std::cout << "SimulationTime: " << simulationTime << "\n";
        std::cout << "SimulationStepTime: " << simulationStepTime << "\n";
        std::cout << "Conductivity: " << conductivity << "\n";
        std::cout << "Alfa: " << alfa << "\n";
        std::cout << "Tot: " << tot << "\n";
        std::cout << "InitialTemp: " << initialTemp << "\n";
        std::cout << "Density: " << density << "\n";
        std::cout << "SpecificHeat: " << specificHeat << "\n";
        std::cout << "Nodes number: " << nodesNumber << "\n";
        std::cout << "Elements number: " << elementsNumber << "\n";
    }

} // oc