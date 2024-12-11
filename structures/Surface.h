//
// Created by xoska on 11.12.2024.
//

#ifndef MES_PROJ_SURFACE_H
#define MES_PROJ_SURFACE_H
#include <vector>

namespace oc {

    class Surface {
    public:
        std::vector<std::vector<double>> N;
        std::vector<double> SurfacePoints;

        void initializeSurface(int npc);
    };

} // oc

#endif //MES_PROJ_SURFACE_H
