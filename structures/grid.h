//
// Created by xoska on 12.10.2024.
//

#ifndef MES_PROJ_GRID_H
#define MES_PROJ_GRID_H

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include "structs.h"

namespace oc {

    class Grid {
        int nN;  // liczba wezlow
        int nE;  // liczba elementow
        std::vector<Node> nodes;
        std::vector<Element> elements;

    public:
        void parseNodes(const std::string& fileName);

        void parseElements(const std::string& fileName);

        void printNodes();
        void printElemets();

        void calculateAllHbc();

        const std::vector<Node> &getNodes() const;

        const std::vector<Element> &getElements() const;


    };

} // oc

#endif //MES_PROJ_GRID_H
