#include "GlobalData.h"

namespace oc {

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
            else if (line.find("*BC") != std::string::npos) {
                // Wczytanie węzłów warunków brzegowych
                std::string nodesLine;
                std::getline(inputFile, nodesLine); // wczytujemy linię z węzłami
                std::istringstream nodesStream(nodesLine);
                std::string node;
                while (std::getline(nodesStream, node, ',')) { // rozdzielamy po przecinkach
                    try {
                        bc.push_back(std::stoi(node));
                    } catch (const std::invalid_argument& e) {
                        std::cerr << "Invalid node in BC: " << node << std::endl;
                    }
                }
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
        std::cout << "Boundary Conditions (BC): ";
        for (const auto& bcx : bc) {
            std::cout << bcx << " ";
        }
        std::cout << "\n";
    }

    int GlobalData::getSimulationTime() const {
        return simulationTime;
    }

    int GlobalData::getConductivity() const {
        return conductivity;
    }

    int GlobalData::getAlfa() const {
        return alfa;
    }

    int GlobalData::getTot() const {
        return tot;
    }

    int GlobalData::getInitialTemp() const {
        return initialTemp;
    }

    int GlobalData::getDensity() const {
        return density;
    }

    int GlobalData::getSpecificHeat() const {
        return specificHeat;
    }

    int GlobalData::getNodesNumber() const {
        return nodesNumber;
    }

    int GlobalData::getElementsNumber() const {
        return elementsNumber;
    }

    std::vector<int>& GlobalData::getBoundaryConditions(){
        return bc;
    }

    int GlobalData::getSimulationStepTime() const {
        return simulationStepTime;
    }

} // oc