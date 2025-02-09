#ifndef MES_PROJ_PROGRAM_H
#define MES_PROJ_PROGRAM_H


#include <algorithm>
#include "..\structures\GlobalData.h"
#include "..\structures\grid.h"
#include "..\structures\structs.h"
#include "..\helpers\gauss.h"
#include "..\structures\solvSystem.h"
#include "..\helpers\functions.h"

namespace oc {

    class Program {
    public:
        Program() = default;

        void zadanie1(std::string &fName);
        void zadanie2();
        void zadanie3();
        void zadanie4(std::string &fName);
        void zadanie5(std::string &fName);
        void zadanie6(std::string &fName);
        void zadanie7(std::string &fName);
        void zadanie8(std::string &fName);

    };

} // oc

#endif //MES_PROJ_PROGRAM_H
