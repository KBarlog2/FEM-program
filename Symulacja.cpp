#include "Elem4.h"
#include <iostream>
#include <iomanip>
#include <algorithm> 
#include <vector>   

void Elem4::Symulacja() {

    cout << "Symulacja:" << endl;

    //obliczanie ilosci krokow
    int totalSteps = static_cast<int>(simulation_time / simulation_step_time);
    double dt = simulation_step_time;

    // wektory
    vector<vector<double>> H_new = global_HBC;
    vector<double> T(num_nodes, initial_temp);
    vector<double> T_new(num_nodes, 0.0);
    vector<double> P_new = global_P;


    for (int t = 0; t < totalSteps; t++) {
        // H_new = H + C / dt
        for (int i = 0; i < num_nodes; i++) {
            for (int j = 0; j < num_nodes; j++) {
                H_new[i][j] = global_HBC[i][j] + global_C[i][j] / dt;
            }
        }

        // P_new = P + C / dt * T
        for (int i = 0; i < num_nodes; i++) {
            P_new[i] = global_P[i];
            for (int j = 0; j < num_nodes; j++) {
                P_new[i] += global_C[i][j] / dt * T[j];
            }
        }

        // H_new * T_new = P_new
        RozszerzenieUkladu(H_new, T_new, P_new);


        T = T_new;

        //min/max
        double minTemp = *min_element(T.begin(), T.end());
        double maxTemp = *max_element(T.begin(), T.end());


        // wyswietlanie temp
        cout << fixed << setprecision(0);
        cout << (t + 1) * dt << " s: ";
        cout << "Min: " << fixed << setprecision(10) << minTemp
            << " | Max: " << maxTemp  << " | Temperatura: ";
        for (double temp : T) {
            cout << fixed << setprecision(2) << temp << " ";
        }
        cout << endl;
    }
}

void Elem4::RozszerzenieUkladu(const vector<vector<double>>& H, vector<double>& x, const vector<double>& P) {
    int n = H.size();
    // Rozszerzenie - macierz H + wektor P
    vector<vector<double>> rozszerzenie(n, vector<double>(n + 1, 0.0));

    // Tworzenie rozszerzonej macierzy
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            rozszerzenie[i][j] = H[i][j];
        }
        rozszerzenie[i][n] = P[i];
    }

    for (int i = 0; i < n; i++) {
        // Pivotowanie
        for (int k = i + 1; k < n; k++) {
            if (fabs(rozszerzenie[k][i]) > fabs(rozszerzenie[i][i])) {
                swap(rozszerzenie[i], rozszerzenie[k]);
            }
        }

        // Ustawienie elementu diagonalnego na 1 i eliminacja kolumny
        for (int k = i + 1; k < n; k++) {
            double wspó³czynnik = rozszerzenie[k][i] / rozszerzenie[i][i];
            for (int j = i; j <= n; j++) {
                rozszerzenie[k][j] -= wspó³czynnik * rozszerzenie[i][j];
            }
        }
    }

    // Zamiana wartoœci
    x.assign(n, 0.0);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = rozszerzenie[i][n] / rozszerzenie[i][i];
        for (int j = i - 1; j >= 0; j--) {
            rozszerzenie[j][n] -= rozszerzenie[j][i] * x[i];
        }
    }
}