//
// Created by xoska on 12.10.2024.
//

#include "grid.h"

namespace oc {
    void Grid::parseNodes(const std::string& fileName) {
        std::ifstream inputFile(fileName);
        if (!inputFile.is_open()) {
            std::cerr << "Error opening file!" << std::endl;
            return;
        }

        std::string line;
        bool inNodeSection = false;

        while (std::getline(inputFile, line)) {
            if (line.find("*Node") != std::string::npos) {
                inNodeSection = true;
                continue;
            }

            if (inNodeSection) {
                if (line.empty() || line[0] == '*') {
                    inNodeSection = false;
                    break;
                }

                std::istringstream iss(line);
                int nodeNumber;
                double x, y;
                char comma;  // To handle commas in the input

                if (iss >> nodeNumber >> comma >> x >> comma >> y) {
                    Node newNode = { x, y };
                    nodes.push_back(newNode);
                }
            }
        }

        inputFile.close();
        nN = nodes.size();
    }

    void Grid::parseElements(const std::string& fileName) {
        std::ifstream inputFile(fileName);
        if (!inputFile.is_open()) {
            std::cerr << "Error opening file!" << std::endl;
            return;
        }

        std::string line;
        bool inElementSection = false;

        while (std::getline(inputFile, line)) {
            if (line.find("type=DC2D4") != std::string::npos) {
                inElementSection = true;
                continue;
            }

            if (inElementSection) {
                if (line.empty() || line[0] == '*') {
                    inElementSection = false;
                    break;
                }

                // Parse element data: elementID, node1, node2, node3, node4
                std::istringstream iss(line);
                int elementID;
                int node1, node2, node3, node4;
                char comma;

                if (iss >> elementID >> comma >> node1 >> comma >> node2 >> comma >> node3 >> comma >> node4) {
                    // Create a new Element struct and push it to the vector
                    Element newElement = { {node1, node2, node3, node4} };
                    elements.push_back(newElement);
                }
            }
        }

        inputFile.close();
        nE = elements.size();
    }

    void Grid::printNodes() {
        for (int i = 0; i < this->nodes.size(); ++i) {
            std::cout << "Node " << i + 1 << ": x = " << this->nodes[i].x << ", y = " << this->nodes[i].y << std::endl;
        }
    }

    void Grid::printElemets() {
        std::cout << "Elements:" << std::endl;
        for (int i = 0; i < this->elements.size(); ++i) {
            std::cout << "Element " << i + 1 << ": Nodes = ";
            for (int j = 0; j < 4; ++j) {
                std::cout << this->elements[i].ID[j] << " ";
            }
            std::cout << std::endl;
        }
    }

    const std::vector<Node> &Grid::getNodes() const {
        return nodes;
    }

    const std::vector<Element> &Grid::getElements() const {
        return elements;
    }

} // oc