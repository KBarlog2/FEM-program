// Elem4.h
#ifndef ELEM4_H
#define ELEM4_H

#include <vector>
#include <string>

using namespace std;

class Elem4 {
public:

    
    double dN_dE[16][4];
    double dN_dN[16][4];
    double J[16][2][2];
    double odwr_J[16][2][2];
    double detJ[16];
    double dN_dX[16][4];
    double dN_dY[16][4];
    double H[16][4][4];
    double suma_H[4][4] = { 0 };

 

    double wsp_K, alfa, tot, initial_temp, density, specific_heat;
    int simulation_time, simulation_step_time, num_nodes, num_elements;

    struct Node {
        int id;
        double x, y;
    };

    struct Element {
        int id;
        vector<int> node_ids;
    };

    vector<Node> nodes;
    vector<Element> elements;
    vector<int> boundary_conditions;
    vector<vector<double>> global_H;
    vector<vector<double>> global_HBC;
    vector<double> global_P;
    vector<vector<double>> global_C;
    vector<vector<double>> macierzC;  
    vector<double> JacobianValues;

    Elem4();
    void loadData(const string& filename);
    void printData();
    void PrintGlobalH();
    void PrintWektorP();
    void wybierzPKT(int order, vector<double>& points, vector<double>& weights);
    void Pochodne(int gaussOrder);
    void Obliczanie(int gaussOrder);
    void MacierzC(int gaussOrder);
    void Symulacja();
    void RozszerzenieUkladu(const vector<vector<double>>& H, vector<double>& x, const vector<double>& P);


};

#endif // ELEM4_H
