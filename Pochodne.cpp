#include "Elem4.h"
#include <iostream>
#include <iomanip>
//Obliczanie funkcji kszta³tu wzglêdem wspó³rzêdnych
void Elem4::Pochodne(int gaussOrder) {
    vector<double> gaussPoints, gaussWeights;
    wybierzPKT(gaussOrder, gaussPoints, gaussWeights);


    int id = 0;
    for (double eta : gaussPoints) {
        for (double ksi : gaussPoints) {
            dN_dE[id][0] = -0.25 * (1 - eta);
            dN_dE[id][1] = 0.25 * (1 - eta);
            dN_dE[id][2] = 0.25 * (1 + eta);
            dN_dE[id][3] = -0.25 * (1 + eta);

            dN_dN[id][0] = -0.25 * (1 - ksi);
            dN_dN[id][1] = -0.25 * (1 + ksi);
            dN_dN[id][2] = 0.25 * (1 + ksi);
            dN_dN[id][3] = 0.25 * (1 - ksi);


            id++;
        }
    }
}
