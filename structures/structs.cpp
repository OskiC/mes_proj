//
// Created by xoska on 12.10.2024.
//

#include "structs.h"

namespace oc {
    void Element::calculateHbc(const std::vector<Node> &nodes, double alfa) {
        double gauss_points[2] = { -1.0 / sqrt(3), 1.0 / sqrt(3) };
        double weights[2] = { 1.0, 1.0 };

        for (int edge = 0; edge < 4; edge++) {
            int n1 = ID[edge] - 1;
            int n2 = ID[(edge + 1) % 4] - 1;

            if (nodes[n1].bc && nodes[n2].bc) {
                double length = sqrt(pow(nodes[n2].x - nodes[n1].x, 2) +
                                     pow(nodes[n2].y - nodes[n1].y, 2));

                for (int gp = 0; gp < 2; gp++) {
                    double ksi = gauss_points[gp];

                    double N_edge[2] = { 0.5 * (1 - ksi), 0.5 * (1 + ksi) };

                    for (int i = 0; i < 2; i++) {
                        for (int j = 0; j < 2; j++) {
                            int local_i = (edge + i) % 4;
                            int local_j = (edge + j) % 4;

                            //std::cout << alfa * N_edge[i] * N_edge[j] * weights[gp] * (length / 2) << std::endl;
                            Hbc[local_i][local_j] += alfa * N_edge[i] * N_edge[j] * weights[gp] * (length / 2);

                        }
                    }
                }
            }
        }
    }

    void Element::calculateP(const std::vector<Node> &nodes, double alfa, double t_ot) {
        P.resize(4, 0.0); // Wektor P (4 węzły)

        const double w[2] = {1.0, 1.0}; // Wagi dla 2-punktowej kwadratury Gaussa
        const double ksi[2] = {-1.0 / sqrt(3), 1.0 / sqrt(3)}; // Współrzędne punktów Gaussa

        for (int edge = 0; edge < 4; ++edge) {
            int id1 = ID[edge] - 1;
            int id2 = ID[(edge + 1) % 4] - 1;

            if (nodes[id1].bc && nodes[id2].bc) { // Jeśli krawędź ma warunki brzegowe
                double x1 = nodes[id1].x, y1 = nodes[id1].y;
                double x2 = nodes[id2].x, y2 = nodes[id2].y;

                double length = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
                double detJ = length / 2.0;

                for (int i = 0; i < 2; ++i) { // Iteracja po punktach Gaussa
                    double N1 = (1.0 - ksi[i]) / 2.0;
                    double N2 = (1.0 + ksi[i]) / 2.0;

                    P[edge] += w[i] * alfa * t_ot * N1 * detJ;
                    P[(edge + 1) % 4] += w[i] * alfa * t_ot * N2 * detJ;
                }
            }
        }
    }

} // oc