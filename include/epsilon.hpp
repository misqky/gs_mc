#pragma once
#include "const_para.hpp"
#include "hopping_para.hpp"
#include "manhattan_distance.hpp"

using namespace std;

extern double epsilon(double kx, double ky)
{
    double eps = 0.0;
    for (int i = 0; i < site_number; ++i)
    {
        int x = i % L;
        int y = i / L;
        if((x + y) % 2 == 1)
        {
            eps += -2.0 * hopping(manhattan_distance(i)) * cos(kx * x + ky * y);
        }

    }

    if(abs(eps) > zero_judge)
    {
        return eps;
    }

    return 0;
}