#pragma once
#include "epsilon.hpp"
#include "eigen-3.4.0/Eigen/Dense"
using namespace std;
using namespace Eigen;

void set_matrix(double epsilon, double H, Matrix2d& A_up, Matrix2d& A_down)
{

        A_up << epsilon, -H,
                -H,    -epsilon;

        A_down <<   epsilon,  H,
                H,  -epsilon;
}