//
// Created by xoska on 12.10.2024.
//

#ifndef MES_PROJ_STRUCTS_H
#define MES_PROJ_STRUCTS_H

#include "jakobian.h"

namespace oc {
    struct Node {
        double x, y;
        bool bc = false;
    };

    struct Element {
        int ID[4];
        Jakobian jakobian;
        double Hbc[4][4] = {0};

        Element(int id1, int id2, int id3, int id4) {
            ID[0] = id1;
            ID[1] = id2;
            ID[2] = id3;
            ID[3] = id4;
        }

        void calculateHbc(const std::vector<Node>& nodes);
    };


} // oc

#endif //MES_PROJ_STRUCTS_H
