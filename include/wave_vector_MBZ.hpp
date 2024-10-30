#pragma once
#define _USE_MATH_DEFINES
#include <utility>
#include <cmath>
#include <vector>
#include "const_para.hpp"

using namespace std;

void wave_vector(vector<pair<double, double>>& k_vector)
{
    int i = 0;
    double kx, ky;
    while(i < site_number / 2)
    {
        for (int n = - (L / 2); n < (L / 2); ++n)
        {
            kx = 2.0 * M_PI / L * n;

            if(kx < 0)  // MBZを考慮
            {
                for(int m = - n - (L / 2); m <= n + (L / 2); ++m)   // (-n-L/2<=m<=n+L/2 for kx<0)
                {
                    ky = 2.0 * M_PI / L * m;
                    k_vector[i] = make_pair(kx, ky);
                    ++i;
                }
            }

            else
            {
                for(int m = n - (L / 2) + 1; m < - n + (L / 2); ++m) // (n-L/2<m<-n+L/2 for kx>=0)
                {
                    ky = 2 * M_PI / L * m;
                    k_vector[i] = make_pair(kx, ky);
                    ++i;
                }
            }
        }
    }
}