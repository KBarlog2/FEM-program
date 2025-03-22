#include "Elem4.h"
#include <iostream>
#include <iomanip>

void Elem4::MacierzC(int gaussOrder) {
    vector<vector<double>> lokalna_C(4, vector<double>(4, 0.0));  // lokalne C
    global_H.resize(num_nodes, vector<double>(num_nodes, 0.0));
    global_C.resize(num_nodes, vector<double>(num_nodes, 0.0));  // globalne C

    vector<double> gaussPoints, gaussWeights;
    wybierzPKT(gaussOrder, gaussPoints, gaussWeights);


    for (int e = 0; e < elements.size(); e++) {
        Element& elem = elements[e];
        fill(lokalna_C.begin(), lokalna_C.end(), vector<double>(4, 0.0));
        double det_J;
        int jacobianIndex = e * gaussPoints.size() * gaussPoints.size(); 

        for (int gp = 0; gp < gaussPoints.size() * gaussPoints.size(); gp++) {
            double ksi = gaussPoints[gp % gaussPoints.size()];
            double eta = gaussPoints[gp / gaussPoints.size()];
            det_J = JacobianValues[jacobianIndex + gp]; 


            //funkcje kszta³tu
            double N[4] = {
                0.25 * (1 - ksi) * (1 - eta),
                0.25 * (1 + ksi) * (1 - eta),
                0.25 * (1 + ksi) * (1 + eta),
                0.25 * (1 - ksi) * (1 + eta)
            };

            //obliczanie lokalnej C
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    lokalna_C[i][j] += density * specific_heat * N[i] * N[j] * det_J *
                        gaussWeights[gp / gaussPoints.size()] * gaussWeights[gp % gaussPoints.size()];
                }
            }
        }

        //print local C 
        cout << "\nMacierz lokalna C dla elementu - " << e + 1 << ":\n";
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                cout << fixed << setprecision(4) << lokalna_C[i][j] << " ";
            }
            cout << endl;
        }

        // Agregacja macierzy lokalnej do globalnej
        for (int i = 0; i < 4; i++) {
            int global_i = elem.node_ids[i] - 1; 
            for (int j = 0; j < 4; j++) {
                int global_j = elem.node_ids[j] - 1; 
                if (global_i < num_nodes && global_j < num_nodes) { 
                    global_C[global_i][global_j] += lokalna_C[i][j];
                }
            }
        }
    }

    //print global C
    cout << "\nGlobalna macierz C:\n";
    for (int i = 0; i < num_nodes; i++) {
        for (int j = 0; j < num_nodes; j++) {
            cout << fixed << setprecision(4) << global_C[i][j] << " ";
        }
        cout << endl;
    }
}
