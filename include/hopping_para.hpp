#pragma once
#define _USE_MATH_DEFINES
#include <cmath>

using namespace std;

extern double hopping(int distances)
{
    double h1 = 1.0;
    double h3 = 0.7;
    double h5 = 0.5;
    if(distances == 1){ return h1; }
    else if(distances == 3){  return h3; }
    else { return h5 * pow((distances), 5.0 / 3); }
}