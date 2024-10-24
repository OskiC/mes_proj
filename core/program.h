#ifndef MES_PROJ_PROGRAM_H
#define MES_PROJ_PROGRAM_H

#include "..\structures\GlobalData.h"
#include "..\structures\grid.h"
#include "..\structures\structs.h"
#include "..\helpers\gauss.h"

namespace oc {

    class Program {
    public:
        Program() = default;

        void zadanie1(std::string &fName);
        void zadanie2();
        void zadanie3();

    };

} // oc

#endif //MES_PROJ_PROGRAM_H
