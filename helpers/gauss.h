#ifndef MES_PROJ_GAUSS_H
#define MES_PROJ_GAUSS_H

#include <vector>
#include <cmath>

namespace oc {



    class GaussPoint {
    private:
        std::vector<double> x;
        std::vector<double> weight;

    public:
        explicit GaussPoint(int x);
        const std::vector<double>& getX();
        const std::vector<double>& getWeight();
    };

    double f1d(double x);
    double f2d(double x, double y);

    double gaussQuadrature1d(int size, double (*f)(double));
    double gaussQuadrature2d(int size, double (*f)(double, double));

} // oc

#endif //MES_PROJ_GAUSS_H
