#include "Elem4.h"
#include <iostream>
#include <iomanip>

void Elem4::wybierzPKT(int order, vector<double>& points, vector<double>& weights) {
    if (order == 2) {
        points = { -1.0 / sqrt(3), 1.0 / sqrt(3) };
        weights = { 1.0, 1.0 };
    }
    else if (order == 3) {
        points = { -sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0) };
        weights = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
    }
    else if (order == 4) {
        points = { -0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053 };
        weights = { 0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454 };
    }
    else {
        throw invalid_argument("Zly wybor.");
    }
}
