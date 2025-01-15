//
// Created by xoska on 12.10.2024.
//

#include "structs.h"

namespace oc {
    void Element::calculateHbc(const std::vector<Node> &nodes, double alfa, int numPoints) {
        // Clear Hbc matrix
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                Hbc[i][j] = 0.0;
            }
        }

        IntegrationPoints ip(numPoints);

        for (int edge = 0; edge < 4; edge++) {
            int n1 = ID[edge] - 1;
            int n2 = ID[(edge + 1) % 4] - 1;

            if (nodes[n1].bc && nodes[n2].bc) {
                double length = std::sqrt(std::pow(nodes[n2].x - nodes[n1].x, 2) +
                                          std::pow(nodes[n2].y - nodes[n1].y, 2));

                //std::cout << "Edge " << edge << " length: " << length << ", Nodes: " << n1 << ", " << n2 << std::endl;

                for (int gp = 0; gp < ip.points.size(); gp++) {
                    double ksi = ip.points[gp].first;

                    double N_edge[2] = { 0.5 * (1 - ksi), 0.5 * (1 + ksi) };

                    for (int i = 0; i < 2; i++) {
                        for (int j = 0; j < 2; j++) {
                            int local_i = (edge + i) % 4;
                            int local_j = (edge + j) % 4;

                            // Adjusting weight calculation
                            long double weight = (ip.weights[gp] * (length / 2)) / 2.0L;

                            double contribution = alfa * N_edge[i] * N_edge[j] * weight;
                            //std::cout << "Contribution for (" << local_i << "," << local_j << "): " << contribution << std::endl;
                            Hbc[local_i][local_j] += contribution;
                        }
                    }
                }
            }
        }
        /*
        // Print final Hbc matrix
        std::cout << "Final Hbc matrix:" << std::endl;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                std::cout << Hbc[i][j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "----------------------------------------" << std::endl;
        */
    }

    void Element::calculateP(const std::vector<Node> &nodes, double alfa, double t_ot, int numPoints) {
        P.resize(4, 0.0); // Initialize vector P with 4 elements to 0.0

        IntegrationPoints ip(numPoints);

        for (int edge = 0; edge < 4; ++edge) {
            int id1 = ID[edge] - 1;
            int id2 = ID[(edge + 1) % 4] - 1;

            if (nodes[id1].bc && nodes[id2].bc) { // If the edge has boundary conditions
                double x1 = nodes[id1].x, y1 = nodes[id1].y;
                double x2 = nodes[id2].x, y2 = nodes[id2].y;

                long double length = std::sqrt(std::pow(x2 - x1, 2) + std::pow(y2 - y1, 2));
                long double detJ = length / 2.0L;

                for (int i = 0; i < ip.points.size(); ++i) { // Iterate over Gauss points
                    double ksi = ip.points[i].first; // For edge integration, we use only ksi

                    double N1 = (1.0 - ksi) / 2.0;
                    double N2 = (1.0 + ksi) / 2.0;

                    // Adjusting weight calculation as per previous discussion
                    long double weight = (ip.weights[i] * detJ) / 2.0L;

                    P[edge] += weight * alfa * t_ot * N1;
                    P[(edge + 1) % 4] += weight * alfa * t_ot * N2;
                }
            }
        }
    }

    std::vector<std::vector<double>> Element::calculateC(const double rho, const double temp, const double detJ, double numPoints) {

        IntegrationPoints integrationPoints(numPoints);
        // Initialize the C matrix (4x4 for 2D elements)
        std::vector<std::vector<double>> C(4, std::vector<double>(4, 0.0));

        // Loop through integration points
        for (size_t i = 0; i < integrationPoints.points.size(); ++i) {
            double ksi = integrationPoints.points[i].first;
            double eta = integrationPoints.points[i].second;
            double weight = integrationPoints.weights[i];

            // Shape functions at the integration point
            double N1 = 0.25 * (1 - ksi) * (1 - eta);
            double N2 = 0.25 * (1 + ksi) * (1 - eta);
            double N3 = 0.25 * (1 + ksi) * (1 + eta);
            double N4 = 0.25 * (1 - ksi) * (1 + eta);

            // Shape function matrix (vector form)
            std::vector<double> N = {N1, N2, N3, N4};

            // Compute Cpc = rho * c_p * (N * N^T)
            std::vector<std::vector<double>> Cpc(4, std::vector<double>(4, 0.0));
            for (int r = 0; r < 4; ++r) {
                for (int c = 0; c < 4; ++c) {
                    Cpc[r][c] = rho * temp * N[r] * N[c] * detJ;
                }
            }
            // Accumulate into the global C matrix with weights
            for (int r = 0; r < 4; ++r) {
                for (int c = 0; c < 4; ++c) {
                    C[r][c] += Cpc[r][c] * weight;
                }
            }
        }

        return C;
    }



    IntegrationPoints::IntegrationPoints(int num) {
        points.clear();
        weights.clear();

        if (num == 4) {
            // 2x2 Gauss-Legendre quadrature
            double sqrtRes = 1.0 / std::sqrt(3.0);
            points = {
                    {-sqrtRes, -sqrtRes},
                    { sqrtRes, -sqrtRes},
                    { sqrtRes,  sqrtRes},
                    {-sqrtRes,  sqrtRes}
            };
            weights = {1.0, 1.0, 1.0, 1.0}; // All weights are 1 for 2x2 quadrature
        } else if (num == 9) {
            // 3x3 Gauss-Legendre quadrature
            double sqrtRes = std::sqrt(3.0 / 5.0);
            points = {
                    {-sqrtRes, -sqrtRes},
                    { 0.0,         -sqrtRes},
                    { sqrtRes, -sqrtRes},
                    {-sqrtRes,  0.0},
                    { 0.0,          0.0},
                    { sqrtRes,  0.0},
                    {-sqrtRes,  sqrtRes},
                    { 0.0,          sqrtRes},
                    { sqrtRes,  sqrtRes}
            };
            weights = {
                    5.0 / 9 * 5.0 / 9, 8.0 / 9 * 5.0 / 9, 5.0 / 9 * 5.0 / 9,
                    5.0 / 9 * 8.0 / 9, 8.0 / 9 * 8.0 / 9, 5.0 / 9 * 8.0 / 9,
                    5.0 / 9 * 5.0 / 9, 8.0 / 9 * 5.0 / 9, 5.0 / 9 * 5.0 / 9
            };
        } else if (num == 16) {
            double a = sqrt((3 + 2 * sqrt(6.0 / 5.0)) / 7.0);
            double b = sqrt((3 - 2 * sqrt(6.0 / 5.0)) / 7.0);
            points = {
                    {-a, -a}, {-b, -a}, {b, -a}, {a, -a},
                    {-a, -b}, {-b, -b}, {b, -b}, {a, -b},
                    {-a, b}, {-b, b}, {b, b}, {a, b},
                    {-a, a}, {-b, a}, {b, a}, {a, a}
            };

            double w1 = (18 - sqrt(30)) / 36.0;
            double w2 = (18 + sqrt(30)) / 36.0;

            weights = {
                    w1 * w1, w1 * w2, w1 * w2, w1 * w1,
                    w2 * w1, w2 * w2, w2 * w2, w2 * w1,
                    w2 * w1, w2 * w2, w2 * w2, w2 * w1,
                    w1 * w1, w1 * w2, w1 * w2, w1 * w1
            };
        } else{
            throw std::invalid_argument("Unsupported number of integration points.");
        }
    }

} // oc