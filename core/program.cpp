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

    void Program::zadanie3() {
        int numPoints = 16; // Number of Gauss points (for 2D, typically 4 points)

        ElemUniv elemUniv;
        elemUniv.initialize(numPoints); // 2d - 4punkty
        //elemUniv.print_dnx();

        double x[4] = {0.01, 0.025, 0.025, 0};
        double y[4] = {-0.01, 0, 0.025, 0.025};

        double k = 30;   // Example thermal conductivity
        double dV = 1.0; // Example volume element

        Jakobian jakobian(numPoints);

        for (int i = 0; i < numPoints; i++) {
            jakobian.calcJakob(elemUniv, x, y, i);
            jakobian.calcDetJ();
            jakobian.calcJakobInver(elemUniv, i);
            jakobian.computeHpc(k, i);
        }

        auto H = jakobian.computeTotalH(k, dV);

        std::cout << "Koncowa macierz H:\n";
        for (const auto& row : H) {
            for (double val : row) {
                if (val != 0.0) {
                    std::cout << std::setw(10) << std::fixed << std::setprecision(4) << val << " ";
                } else {
                    std::cout << std::setw(10) << " ";
                }
            }
            std::cout << "\n";
        }
    }

    void Program::zadanie4(std::string &fName){
        oc::GlobalData globalData{};
        globalData.parseFromFile(fName);

        oc::Grid grid;
        grid.parseNodes(fName);
        grid.parseElements(fName);

        double k = globalData.getConductivity();
        double dV = 1.0;

        int numPoints = 16;

        const auto& elements = grid.getElements();
        const auto& nodes = grid.getNodes();

        for (const auto& element : elements) {
            double x[4], y[4];

            for (int i = 0; i < 4; ++i) {
                int nodeID = element.ID[i] - 1;
                x[i] = nodes[nodeID].x;
                y[i] = nodes[nodeID].y;
            }

            oc::Jakobian jakobian(numPoints);
            oc::ElemUniv elemUniv;
            elemUniv.initialize(numPoints);
            std::vector<std::vector<double>> H_total;

            for (int p = 0; p < numPoints; ++p) {
                // Calculate the Jacobian for the integration point
                jakobian.calcJakob(elemUniv, x, y, p);
                jakobian.calcDetJ();
                jakobian.calcJakobInver(elemUniv, p);
                jakobian.calc_dN_dX_dN_dY(elemUniv, p);
                jakobian.computeHpc(k, p);  // Calculate H matrix for the integration point
            }

            H_total = jakobian.computeTotalH(k, dV);

            std::cout << "Final H matrix for element " << &element - &elements[0] + 1 << ":\n";
            for (const auto& row : H_total) {
                bool firstValue = true;
                for (double val : row) {
                    if (val != 0) { // Skip zeros
                        if (!firstValue) {
                            std::cout << " ";
                        }
                        std::cout << std::setw(10) << val;
                        firstValue = false;
                    }
                }
                std::cout << "\n";
            }
            std::cout << "-----------------------------\n";
        }

    }

    void Program::zadanie5(std::string& fName) {
        oc::GlobalData globalData{};
        globalData.parseFromFile(fName);

        oc::Grid grid;
        grid.parseNodes(fName);
        grid.parseElements(fName);

        double k = globalData.getConductivity();
        double dV = 1.0;
        int numPoints = 4;

        const auto &elements = grid.getElements();
        const auto &nodes = grid.getNodes();

        // Tworzymy obiekt klasy SolvSystem z liczbą węzłów
        oc::SolvSystem solvSystem(nodes.size());

        // Przechodzimy przez wszystkie elementy i obliczamy macierz H
        for (const auto &element: elements) {
            double x[4], y[4];

            // Wczytujemy współrzędne węzłów
            for (int i = 0; i < 4; ++i) {
                int nodeID = element.ID[i] - 1;
                x[i] = nodes[nodeID].x;
                y[i] = nodes[nodeID].y;
            }

            // Inicjalizujemy obiekt Jakobianu i innych pomocniczych
            oc::Jakobian jakobian(numPoints);
            oc::ElemUniv elemUniv;
            elemUniv.initialize(numPoints);

            // Tablica na macierz H lokalną
            std::vector<std::vector<double>> H_local(4, std::vector<double>(4, 0.0));

            for (int p = 0; p < numPoints; ++p) {
                jakobian.calcJakob(elemUniv, x, y, p);
                jakobian.calcDetJ();
                jakobian.calcJakobInver(elemUniv, p);
                jakobian.calc_dN_dX_dN_dY(elemUniv, p);
                jakobian.computeHpc(k, p); // Obliczamy H dla punktu całkowania
            }

            H_local = jakobian.computeTotalH(k, dV); // Sumujemy wkład punktów całkowania

            // Agregacja do macierzy globalnej
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    int global_i = element.ID[i] - 1; // ID globalne wiersza
                    int global_j = element.ID[j] - 1; // ID globalne kolumny

                    solvSystem.addToGlobalMatrix(global_i, global_j, H_local[i][j]);
                }
            }
        }
        globalData.printData();
        solvSystem.printMatrix();
    }

} // oc