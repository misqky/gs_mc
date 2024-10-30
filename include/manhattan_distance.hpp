#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include"const_para.hpp"

using namespace std;

extern int manhattan_distance(int i)
{
    int x = i % L;
    int y = i / L;
    int d1 = abs(x) + abs(y);
    int d2 = abs(x - L) + abs(y);
    int d3 = abs(x) + abs(y - L);
    int d4 = abs(x - L) + abs(y - L);

    return min({d1,d2,d3,d4});
}