#pragma once
#include <complex>
#include <iostream>
#include "eigen-3.4.0/Eigen/Dense"
#include"const_para.hpp"
using namespace std;
using namespace Eigen;

complex<double> determinant(1, 1);

extern complex<double> cal_determinant(MatrixXcd A_inverse)
{
    for (int i = 0; i < site_number; ++i )
    {
        determinant *= A_inverse(i, i);
    }
    return determinant;
}