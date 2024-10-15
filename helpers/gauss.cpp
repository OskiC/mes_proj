#include "gauss.h"

namespace oc {

    GaussPoint::GaussPoint(int x) {
        if(x == 1){
            this->x.push_back(0);
            this->weight.push_back(2);
        }
        else if(x == 2){
            this->x.push_back(-sqrt(1.0/3.0));
            this->x.push_back(sqrt(1.0/3.0));
            this->weight.push_back(1);
            this->weight.push_back(1);
        }
        else if(x == 3){
            this->x.push_back(-sqrt(3.0/5.0));
            this->x.push_back(0);
            this->x.push_back(sqrt(3.0/5.0));

            this->weight.push_back(5.0/9.0);
            this->weight.push_back(8.0/9.0);
            this->weight.push_back(5.0/9.0);
        }
    }

    const std::vector<double>& GaussPoint::getX(){
        return x;
    }

    const std::vector<double>& GaussPoint::getWeight() {
        return weight;
    }

    double gaussQuadrature1d(int size, double (*f)(double)){
        double integral = 0;
        GaussPoint gaussPoint(size);

        for(int i = 0; i < gaussPoint.getX().size(); ++i){
            double pointX = gaussPoint.getX()[i];
            double weight = gaussPoint.getWeight()[i];
            integral += weight * f(pointX);
        }

        return integral;
    }

    double gaussQuadrature2d(int size, double (*f)(double, double)){
        double integral = 0;
        GaussPoint gaussPoint(size);

        for(int i = 0; i < gaussPoint.getX().size(); ++i){
            for(int j = 0; j < gaussPoint.getX().size(); ++j){
                integral += gaussPoint.getWeight()[i] * gaussPoint.getWeight()[j] *
                            f(gaussPoint.getX()[i], gaussPoint.getX()[j]);
            }
        }

        return integral;
    }

    double f1d(double x){
        return 5 * x * x + 3 * x + 6;
    }

    double f2d(double x, double y) {
        return 5 * x * x * y * y + 3 * x * y + 6;
    }


} // oc