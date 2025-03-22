#include "Elem4.h"
#include <fstream>
#include <sstream>
#include <iostream>

void Elem4::loadData(const string& filename) {
    ifstream file(filename);
    if (!file) {
        cerr << "Nie mo¿na znaleŸæ pliku: " << filename << endl;
        exit(1);
    }

    string line;
    while (getline(file, line)) {
        istringstream iss(line);


        if (line.find("SimulationTime") != string::npos) {
            iss.ignore(15);
            iss >> simulation_time;
        }
        else if (line.find("SimulationStepTime") != string::npos) {
            iss.ignore(19);
            iss >> simulation_step_time;
        }
        else if (line.find("Conductivity") != string::npos) {
            iss.ignore(13);
            iss >> wsp_K;
        }
        else if (line.find("Alfa") != string::npos) {
            iss.ignore(5);
            iss >> alfa;
        }
        else if (line.find("Tot") != string::npos) {
            iss.ignore(4);
            iss >> tot;
        }
        else if (line.find("InitialTemp") != string::npos) {
            iss.ignore(12);
            iss >> initial_temp;
        }
        else if (line.find("Density") != string::npos) {
            iss.ignore(8);
            iss >> density;
        }
        else if (line.find("SpecificHeat") != string::npos) {
            iss.ignore(13);
            iss >> specific_heat;
        }
        else if (line.find("Nodes number") != string::npos) {
            iss.ignore(13);
            iss >> num_nodes;
        }
        else if (line.find("Elements number") != string::npos) {
            iss.ignore(16);
            iss >> num_elements;
        }
        else if (line.find("*Node") != string::npos) {

            for (int i = 0; i < num_nodes; i++) {
                getline(file, line);
                istringstream node_stream(line);
                int id;
                double x, y;
                char comma;
                node_stream >> id >> comma >> x >> comma >> y;
                nodes.push_back({ id, x, y });
            }
        }
        else if (line.find("*Element") != string::npos) {

            for (int i = 0; i < num_elements; i++) {
                getline(file, line);
                istringstream elem_stream(line);
                int id, n1, n2, n3, n4;
                char comma;
                elem_stream >> id >> comma >> n1 >> comma >> n2 >> comma >> n3 >> comma >> n4;
                elements.push_back({ id, {n1, n2, n3, n4} });
            }
        }
        else if (line.find("*BC") != string::npos) {

            getline(file, line);
            istringstream bc_stream(line);
            int bc_id;
            while (bc_stream >> bc_id) {
                boundary_conditions.push_back(bc_id);
                if (bc_stream.peek() == ',') bc_stream.ignore();
            }
        }
    }
}

