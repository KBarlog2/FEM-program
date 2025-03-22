#include "Elem4.h"
#include <iostream>
#include <iomanip>


//glowny plik z obliczeniami jakobianu, macierzy H, warunkow brzegowych
void Elem4::Obliczanie(int gaussOrder) {

    vector<double> gaussPoints, gaussWeights;
    wybierzPKT(gaussOrder, gaussPoints, gaussWeights);


    for (int e = 0; e < elements.size(); e++) {
        //cout << "\n\n\n\nElement " << e + 1 << ":\n\n";
        //cout << "\nElement " << e + 1 << ":\n";
        Element& elem = elements[e];
        global_H.resize(num_nodes, vector<double>(num_nodes, 0.0));
        global_HBC.resize(num_nodes, vector<double>(num_nodes, 0.0));
        global_P.resize(num_nodes, 0.0);
        double H_elem[4][4] = { 0 };
        double HBC_elem[4][4] = { 0 };
        double P_elem[4] = { 0 };

        int id = 0;
        for (double eta : gaussPoints) {
            for (double ksi : gaussPoints) {
                double dX_dE = 0, dY_dE = 0, dX_dN = 0, dY_dN = 0;

                //skladowe J
                for (int n = 0; n < 4; n++) {
                    int nodeIndex = elem.node_ids[n] - 1;
                    dX_dE += dN_dE[id][n] * nodes[nodeIndex].x;
                    dY_dE += dN_dE[id][n] * nodes[nodeIndex].y;
                    dX_dN += dN_dN[id][n] * nodes[nodeIndex].x;
                    dY_dN += dN_dN[id][n] * nodes[nodeIndex].y;
                }

                // macierz J
                J[id][0][0] = dX_dE;
                J[id][0][1] = dY_dE;
                J[id][1][0] = dX_dN;
                J[id][1][1] = dY_dN;


                detJ[id] = dX_dE * dY_dN - dX_dN * dY_dE;
                JacobianValues.push_back(detJ[id]);

                odwr_J[id][0][0] = dY_dN / detJ[id];
                odwr_J[id][0][1] = -dY_dE / detJ[id];
                odwr_J[id][1][0] = -dX_dN / detJ[id];
                odwr_J[id][1][1] = dX_dE / detJ[id];

                
                cout << fixed << setprecision(6);
                cout << "\nMacierz Jacobiego w punkcie " << id + 1 << ":\n";
                cout << "[" << J[id][0][0] << ", " << J[id][0][1] << "]\n";
                cout << "[" << J[id][1][0] << ", " << J[id][1][1] << "]\n";
                cout << "DetJ: " << detJ[id] << "\n";
                

                //obliczanie pochodnych wzgledem J
                for (int n = 0; n < 4; n++) {
                    dN_dX[id][n] = odwr_J[id][0][0] * dN_dE[id][n] + odwr_J[id][0][1] * dN_dN[id][n];
                    dN_dY[id][n] = odwr_J[id][1][0] * dN_dE[id][n] + odwr_J[id][1][1] * dN_dN[id][n];
            //        cout << "dN_dX[" << n + 1 << "] = " << dN_dX[id][n]
            //            << ", dN_dY[" << n + 1 << "] = " << dN_dY[id][n] << "\n";
                }

                // macierz H
                for (int i = 0; i < 4; i++) {
                    for (int j = 0; j < 4; j++) {
                        double dN_dX_i = dN_dX[id][i];
                        double dN_dX_j = dN_dX[id][j];
                        double dN_dY_i = dN_dY[id][i];
                        double dN_dY_j = dN_dY[id][j];

                        H_elem[i][j] += wsp_K * (dN_dX[id][i] * dN_dX[id][j] + dN_dY[id][i] * dN_dY[id][j]) * detJ[id] * gaussWeights[id % gaussOrder] * gaussWeights[id / gaussOrder];
                    }
                }
                id++;
            }
        }






        // Warunki brzegowe
        for (int edge = 0; edge < 4; edge++) {
            //n1 i n2 - indeksy wezlow krawedzi
            int n1 = elem.node_ids[edge] - 1;
            int n2 = elem.node_ids[(edge + 1) % 4] - 1;

            //sprawdzenie czy naleza do warunkow brzegowych
            if (find(boundary_conditions.begin(), boundary_conditions.end(), n1 + 1) != boundary_conditions.end() &&
                find(boundary_conditions.begin(), boundary_conditions.end(), n2 + 1) != boundary_conditions.end()) {






                for (int gp = 0; gp < gaussPoints.size(); gp++) {
                    double ksi = gaussPoints[gp];


                    double N_edge[2] = { 0.5 * (1 - ksi), 0.5 * (1 + ksi) };

                    //obliczenie dlugosci krawedzi
                    double length = sqrt(pow(nodes[n2].x - nodes[n1].x, 2) +
                        pow(nodes[n2].y - nodes[n1].y, 2));


                    for (int i = 0; i < 2; i++) {
                        for (int j = 0; j < 2; j++) {

                            //zmiana na indeksy macierzy 4x4 elementu
                            int local_i = (edge + i) % 4;
                            int local_j = (edge + j) % 4;


                            HBC_elem[local_i][local_j] +=
                                alfa * N_edge[i] * N_edge[j] * gaussWeights[gp] * (length / 2);
                        }
                    }

                    //Obliczanie wektora P

                    for (int i = 0; i < 2; i++) {
                        int local_i = (edge + i) % 4;
                        P_elem[local_i] += alfa * tot * N_edge[i] * gaussWeights[gp] * (length / 2);
                        //P_elem[local_i] += static_cast<int>(round(alfa * tot * N_edge[i] * gaussWeights[gp] * (length / 2)));
                    }
                }
            }
        }




        
        //macierz HBC
        cout << "\nMacierz HBC dla elementu " << e + 1 << ":\n";
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                cout << fixed << setprecision(4) << HBC_elem[i][j];
                if (j < 3) cout << " ";
            }
            cout << "\n";
        }
        





        //macierz H
        cout << "\nMacierz H dla elementu " << e + 1 << ":\n";
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                cout << fixed << setprecision(4) << H_elem[i][j];
                if (j < 3) cout << " ";
            }
            cout << "\n";
        }






        //Obliczanie dodanych macierzy
        double H_total[4][4] = { 0 };
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                H_total[i][j] = H_elem[i][j] + HBC_elem[i][j];

            }
        }




        
        //macierz final
        cout << "\nMacierz H + HBC dla elementu " << e + 1 << ":\n";
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                cout << fixed << setprecision(4) << H_total[i][j];
                if (j < 3) cout << " ";
            }
            cout << "\n";
        }
        




        // Globalna macierz H
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                int global_i = elem.node_ids[i] - 1;
                int global_j = elem.node_ids[j] - 1;
                global_H[global_i][global_j] += H_elem[i][j];//  //+hbc
                global_HBC[global_i][global_j] += H_elem[i][j] + HBC_elem[i][j];
            }
        }




        
        //Wektor P
        cout << "\nWektor P dla elementu " << e + 1 << ":\n";
        for (int i = 0; i < 4; i++) {
            cout << P_elem[i] << " ";
        }
        cout << endl;
        



        // P globalne
        for (int i = 0; i < 4; i++) {
            double global_i = elem.node_ids[i] - 1;
            global_P[global_i] += P_elem[i];
        }
    }
}

