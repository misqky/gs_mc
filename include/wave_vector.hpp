#pragma once
#define _USE_MATH_DEFINES
#include <utility>
#include <cmath>
#include "const_para.hpp"

using namespace std;

extern pair<double, double> wave_vector(int i)
{
    int x = i % L;
    int y = i / L;
    double kx = 2.0 * M_PI / L * (x - (L - 2) / 2);
    double ky = 2.0 * M_PI / L * (y - (L - 2) / 2);
    return make_pair(kx, ky);
}