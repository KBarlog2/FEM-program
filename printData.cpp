#include "Elem4.h"
#include <iostream>
#include <iomanip>

//printowanie podstawowych danych
void Elem4::printData() {
    cout << "  Simulation Time: " << simulation_time << endl;
    cout << "  Simulation Step Time: " << simulation_step_time << endl;
    cout << "  Conductivity (K): " << wsp_K << endl;
    cout << "  Alfa: " << alfa << endl;
    cout << "  Tot: " << tot << endl;
    cout << "  Initial Temperature: " << initial_temp << endl;
    cout << "  Density: " << density << endl;
    cout << "  Specific Heat: " << specific_heat << endl;
    cout << "  Number of Nodes: " << num_nodes << endl;
    cout << "  Number of Elements: " << num_elements << endl;


    cout << "\nNodes:\n";
    for (const auto& node : nodes) {
        cout << "  ID: " << node.id << ", X: " << node.x << ", Y: " << node.y << endl;
    }


    cout << "\nElements:\n";
    for (const auto& element : elements) {
        cout << "  ID: " << element.id << ", Nodes: ";
        for (int nid : element.node_ids) {
            cout << nid << " ";
        }
        cout << endl;
    }


    cout << "\nBoundary Conditions (BC):\n";
    for (int bc : boundary_conditions) {
        cout << "  Node ID: " << bc << endl;
    }
}


//print wektor P
void Elem4::PrintWektorP() {
    cout << "\nGlobalny wektor P:\n";
    for (int i = 0; i < num_nodes; i++) {
        cout << global_P[i] << " ";
    }
    cout << "\n\n";
}

//print H globalne
void Elem4::PrintGlobalH() {
    cout << "\n\nGlobalna macierz H:\n";
    for (int i = 0; i < num_nodes; i++) {
        for (int j = 0; j < num_nodes; j++) {
            cout << fixed << setprecision(4) << global_H[i][j] << " ";
        }
        cout << endl;
    }

    cout << "\n\nGlobalna macierz H z warunkami brzegowymi BC:\n";
    for (int i = 0; i < num_nodes; i++) {
        for (int j = 0; j < num_nodes; j++) {
            cout << fixed << setprecision(4) << global_HBC[i][j] << " ";
        }
        cout << endl;
    }
}

