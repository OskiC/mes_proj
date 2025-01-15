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

    void Program::zadanie6(std::string& fName) {
        oc::GlobalData globalData{};
        globalData.parseFromFile(fName);

        oc::Grid grid;
        grid.parseNodes(fName);
        grid.parseElements(fName);

        double k = globalData.getConductivity();
        double alpha = globalData.getAlfa(); // Współczynnik wymiany ciepła
        double t_ot = globalData.getTot();  // Temperatura otoczenia
        double dV = 1.0;
        int numPoints = 4;

        const auto &elements = grid.getElements();
        const auto &nodes = grid.getNodes();

        oc::SolvSystem solvSystem(nodes.size());

        int elemNr = 1;
        for (Element element : elements) {
            double x[4], y[4];

            oc::Jakobian jakobian(numPoints);
            oc::ElemUniv elemUniv;
            elemUniv.initialize(numPoints);

            std::cout << "-------------------------------------------\n";
            std::cout << "Element nr: " << elemNr << "\n";

            // Wczytujemy współrzędne węzłów
            for (int i = 0; i < 4; ++i) {
                int nodeID = element.ID[i] - 1;
                x[i] = nodes[nodeID].x;
                y[i] = nodes[nodeID].y;
            }

            // Punkt całkowania
            std::cout << "Punkty calkowania:\n";
            for (int p = 0; p < numPoints; ++p) {
                std::cout << "Punkt nr " << p + 1
                          << " ksi = " << elemUniv.ksi[p]
                          << ", eta = " << elemUniv.eta[p] << "\n";
            }

            // Jakobian i jego detJ


            for (int p = 0; p < numPoints; ++p) {
                jakobian.calcJakob(elemUniv, x, y, p);
                jakobian.calcDetJ();
                jakobian.calcJakobInver(elemUniv, p);
                jakobian.calc_dN_dX_dN_dY(elemUniv, p);
                jakobian.computeHpc(k, p);

                std::cout << "Macierz Jakobianu dla punktu nr " << p + 1 << ":\n";
                jakobian.printJakob();
                std::cout << "det[J] = " << jakobian.getDet() << "\n";

                std::cout << "Wartosci dN/dx oraz dN/dy:\n";
                std::cout << "Wartosci dN/dx dla punktu nr " << p + 1 << ":\n";
                const auto& dN_dX = jakobian.get_dN_dX();
                for (const auto& row : dN_dX) {
                    for (double val : row) {
                        std::cout << val << " ";
                    }
                    std::cout << "\n";
                }

                std::cout << "Wartosci dN/dy dla punktu nr " << p + 1 << ":\n";
                const auto& dN_dY = jakobian.get_dN_dY();
                for (const auto& row : dN_dY) {
                    for (double val : row) {
                        std::cout << val << " ";
                    }
                    std::cout << "\n";
                }
            }

            // Macierz H
            std::vector<std::vector<double>> H_local(4, std::vector<double>(4, 0.0));
            H_local = jakobian.computeTotalH(k, dV);

            std::cout << "Macierz H dla elementu nr " << elemNr << ":\n";
            for (auto &row : H_local) {
                for (double val : row) {
                    std::cout << val << " ";
                }
                std::cout << "\n";

            }

            // Macierz Hbc
            element.calculateHbc(nodes, alpha, numPoints);
            std::cout << "Macierz Hbc dla elementu nr " << elemNr << ":\n";
            for (auto &row : element.Hbc) {
                for (double val : row) {
                    std::cout << val << " ";
                }
                std::cout << "\n";

            }

            for(int i = 0; i < 4; ++i){
                for(int j = 0; j < 4; ++j){
                    int global_i = element.ID[i] - 1;
                    int global_j = element.ID[j] - 1;

                    solvSystem.addToGlobalMatrix(global_i, global_j, H_local[i][j]);
                    solvSystem.addToGlobalMatrix(global_i, global_j, element.Hbc[i][j]);
                }
            }

            // Wektor P
            element.calculateP(nodes, alpha, t_ot, numPoints);
            solvSystem.addToGlobalVector(element.P, element.ID);
            std::cout << "Wektor {P} dla elementu nr " << elemNr << ":\n";
            for (double val : element.P) {
                std::cout << val << " ";
            }
            std::cout << "\n";
            elemNr++;
        }
        solvSystem.printMatrix();
        solvSystem.printVector();
    }

    void Program::zadanie7(std::string &fName) {
        oc::GlobalData globalData{};
        globalData.parseFromFile(fName);

        oc::Grid grid;
        grid.parseNodes(fName);
        grid.parseElements(fName);

        double k = globalData.getConductivity();
        double alpha = globalData.getAlfa(); // Współczynnik wymiany ciepła
        double t_ot = globalData.getTot();  // Temperatura otoczenia
        double density = globalData.getDensity();
        double specificHeat = globalData.getSpecificHeat();
        double dV = 1.0;
        int numPoints = 16;

        const auto &elements = grid.getElements();
        const auto &nodes = grid.getNodes();

        oc::SolvSystem solvSystem(nodes.size());

        int elemNr = 1;
        for (Element element : elements) {
            double x[4], y[4];

            oc::Jakobian jakobian(numPoints);
            oc::ElemUniv elemUniv;
            elemUniv.initialize(numPoints);

            std::cout << "-------------------------------------------\n";
            std::cout << "Element nr: " << elemNr << "\n";

            // Wczytujemy współrzędne węzłów
            for (int i = 0; i < 4; ++i) {
                int nodeID = element.ID[i] - 1;
                x[i] = nodes[nodeID].x;
                y[i] = nodes[nodeID].y;
            }

            // Punkt całkowania
            //std::cout << "Punkty calkowania:\n";
            for (int p = 0; p < numPoints; ++p) {
              //  std::cout << "Punkt nr " << p + 1
                //          << " ksi = " << elemUniv.ksi[p]
                  //        << ", eta = " << elemUniv.eta[p] << "\n";
            }

            // Jakobian i jego detJ


            for (int p = 0; p < numPoints; ++p) {
                jakobian.calcJakob(elemUniv, x, y, p);
                jakobian.calcDetJ();
                jakobian.calcJakobInver(elemUniv, p);
                jakobian.calc_dN_dX_dN_dY(elemUniv, p);
                jakobian.computeHpc(k, p);

                //std::cout << "Macierz Jakobianu dla punktu nr " << p + 1 << ":\n";
                //jakobian.printJakob();
                //std::cout << "det[J] = " << jakobian.getDet() << "\n";

                //std::cout << "Wartosci dN/dx oraz dN/dy:\n";
                //std::cout << "Wartosci dN/dx dla punktu nr " << p + 1 << ":\n";
                const auto& dN_dX = jakobian.get_dN_dX();
                for (const auto& row : dN_dX) {
                    for (double val : row) {
                  //      std::cout << val << " ";
                    }
                    //std::cout << "\n";
                }

                //std::cout << "Wartosci dN/dy dla punktu nr " << p + 1 << ":\n";
                const auto& dN_dY = jakobian.get_dN_dY();
                for (const auto& row : dN_dY) {
                    for (double val : row) {
                  //      std::cout << val << " ";
                    }
                    //std::cout << "\n";
                }
            }

            // Macierz H
            std::vector<std::vector<double>> H_local(4, std::vector<double>(4, 0.0));
            H_local = jakobian.computeTotalH(k, dV);

            //std::cout << "Macierz H dla elementu nr " << elemNr << ":\n";
            for (auto &row : H_local) {
                for (double val : row) {
                //    std::cout << val << " ";
                }
                //std::cout << "\n";

            }

            // Macierz Hbc
            element.calculateHbc(nodes, alpha, numPoints);
            //std::cout << "Macierz Hbc dla elementu nr " << elemNr << ":\n";
            for (auto &row : element.Hbc) {
                for (double val : row) {
                //    std::cout << val << " ";
                }
                //std::cout << "\n";

            }

            for(int i = 0; i < 4; ++i){
                for(int j = 0; j < 4; ++j){
                    int global_i = element.ID[i] - 1;
                    int global_j = element.ID[j] - 1;

                    solvSystem.addToGlobalMatrix(global_i, global_j, H_local[i][j]);
                    solvSystem.addToGlobalMatrix(global_i, global_j, element.Hbc[i][j]);
                }
            }

            // Wektor P
            element.calculateP(nodes, alpha, t_ot, numPoints);
            solvSystem.addToGlobalVector(element.P, element.ID);
            std::cout << "Wektor {P} dla elementu nr " << elemNr << ":\n";
            for (double val : element.P) {
                std::cout << val << " ";
            }
            std::cout << "\n";

            std::vector<std::vector<double>> C = element.calculateC(density, specificHeat, jakobian.getDet(), numPoints);
            solvSystem.addToGlobalC(C, element.ID);
            //std::cout << "Matrix C dla elementu nr: " << elemNr << std::endl;
            for (const auto& row : C) {
                for (double val : row) {
              //      std::cout << val << " ";
                }
                //std::cout << "\n";
            }

            elemNr++;
        }
        std::cout << "\n\n\n";
        solvSystem.printMatrix();
        solvSystem.printVector();
        solvSystem.printC();
    }

    void Program::zadanie8(std::string &fName) {
        oc::GlobalData globalData{};
        globalData.parseFromFile(fName);

        oc::Grid grid;
        grid.parseNodes(fName);
        grid.parseElements(fName);

        double k = globalData.getConductivity();
        double alpha = globalData.getAlfa(); // Heat transfer coefficient
        double t_ot = globalData.getTot();  // Ambient temperature
        double density = globalData.getDensity();
        double specificHeat = globalData.getSpecificHeat();
        double dV = 1.0;
        int numPoints = 16;
        double timeStep = globalData.getSimulationStepTime(); // Time step Δτ
        int numTimeSteps = globalData.getSimulationTime() / timeStep; // Number of time steps

        const auto &elements = grid.getElements();
        const auto &nodes = grid.getNodes();

        oc::SolvSystem solvSystem(nodes.size());

        // Initialize the temperature vector with initial conditions
        std::vector<double> t0(nodes.size(), globalData.getInitialTemp());
        std::vector<double> t1(nodes.size(), 0.0); // For storing the solution at next time step

        // Output header
        std::cout << "Time[s] MinTemp MaxTemp\n";

        // Time-stepping loop
        for (int timeStepIndex = 0; timeStepIndex < numTimeSteps; ++timeStepIndex) {
            std::cout << "-------------------------------------------\n";
            std::cout << "Time Step: " << timeStepIndex + 1 << "\n";

            // Reset global system for each time step
            solvSystem.resetGlobalSystem();

            for (Element element : elements) {
                double x[4], y[4];

                oc::Jakobian jakobian(numPoints);
                oc::ElemUniv elemUniv;
                elemUniv.initialize(numPoints);

                // Load node coordinates
                for (int i = 0; i < 4; ++i) {
                    int nodeID = element.ID[i] - 1;
                    x[i] = nodes[nodeID].x;
                    y[i] = nodes[nodeID].y;
                }

                // Compute Jacobian and related matrices for each integration point
                for (int p = 0; p < numPoints; ++p) {
                    jakobian.calcJakob(elemUniv, x, y, p);
                    jakobian.calcDetJ();
                    jakobian.calcJakobInver(elemUniv, p);
                    jakobian.calc_dN_dX_dN_dY(elemUniv, p);
                    jakobian.computeHpc(k, p);
                }

                // Compute local H matrix
                std::vector<std::vector<double>> H_local = jakobian.computeTotalH(k, dV);

                // Compute Hbc matrix
                element.calculateHbc(nodes, alpha, numPoints);

                // Aggregate H and Hbc into global system
                for(int i = 0; i < 4; ++i){
                    for(int j = 0; j < 4; ++j){
                        int global_i = element.ID[i] - 1;
                        int global_j = element.ID[j] - 1;
                        solvSystem.addToGlobalMatrix(global_i, global_j, H_local[i][j] + element.Hbc[i][j]);
                    }
                }

                // Compute P vector
                element.calculateP(nodes, alpha, t_ot, numPoints);
                solvSystem.addToGlobalVector(element.P, element.ID);

                // Compute C matrix
                std::vector<std::vector<double>> C = element.calculateC(density, specificHeat, jakobian.getDet(), numPoints);
                solvSystem.addToGlobalC(C, element.ID);
            }

            // Solve the system for each time step
            solvSystem.solve(t0, timeStep, t1);

            // Update temperatures for next iteration
            t0 = t1;

            // Find min and max temperatures
            double minTemp = *std::min_element(t0.begin(), t0.end());
            double maxTemp = *std::max_element(t0.begin(), t0.end());

            // Output the results in the desired format
            std::cout << (timeStepIndex + 1) * timeStep << " " << minTemp << " " << maxTemp << "\n";
        }

        // Final output (optional, since we're outputting at each step)
        std::cout << "\n\nSimulation Completed\n";
    }

} // oc