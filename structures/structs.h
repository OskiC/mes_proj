//
// Created by xoska on 12.10.2024.
//

#ifndef MES_PROJ_STRUCTS_H
#define MES_PROJ_STRUCTS_H

#include "jakobian.h"

namespace oc {
    struct Node {
        double x, y;
    };

    struct Element {
        int ID[4];
        Jakobian jakobian;
    };


} // oc

#endif //MES_PROJ_STRUCTS_H
