#pragma once
#include <iostream>
#include <complex>
#include <vector>
#include "eigen-3.4.0/Eigen/Dense"
#include "const_para.hpp"
#include "set_g.hpp"
#include "make_pair.hpp"

using namespace std;
using namespace Eigen;

extern complex<double> cal_flipweight(vector<int> spin, MatrixXcd g_up_normal, MatrixXcd g_up_flip, MatrixXcd g_down_normal, MatrixXcd g_down_flip, pair<complex<double>, complex<double>> det, pair<MatrixXcd, MatrixXcd> A_inverse)
{
    auto nearest_pairs = make_nearest_pair(); // 一度呼び出して結果を保存
    pair<int, int> nearest_pair;
    complex<double> result(0, 0);
    complex<double> d_up_tilde, B_up_buffer, d_up_prime, B_up[site_number];
    complex<double> d_down_tilde, B_down_buffer, d_down_prime, B_down[site_number];
    complex<double> d_prime;

    for (const auto& nearest_pair : nearest_pairs)
    {
        int i = nearest_pair.first;
        int j = nearest_pair.second;
        int alpha_i = spin[i];
        int alpha_j = spin[j];

        d_up_tilde = complex<double>(0, 0);
        B_up_buffer = complex<double>(0, 0);
        d_down_tilde = complex<double>(0, 0);
        B_down_buffer = complex<double>(0, 0);

        if (alpha_i == alpha_j) continue;



        for (int ik = 0; ik < site_number; ++ik)
        {
            d_up_tilde += g_up_flip(ik, i) * A_inverse.first(i, ik);
            B_up_buffer += A_inverse.first(j, ik) * (g_up_flip(ik, i) - g_up_normal(ik, i));
            d_down_tilde += g_down_flip(ik, i) * A_inverse.second(i, ik);
            B_down_buffer += A_inverse.second(j, ik) * (g_down_flip(ik, i) - g_down_normal(ik, i));
        }

        d_up_tilde *= det.first;
        B_up_buffer *= det.first / d_up_tilde;
        d_down_tilde *= det.second;
        B_down_buffer *= det.second / d_down_tilde;

        for (int ik = 0; ik < site_number; ++ik)
        {
            B_up[ik] = -B_up_buffer * A_inverse.first(i, ik);
            B_down[ik] = -B_down_buffer * A_inverse.second(i, ik);
        }

        d_up_prime = complex<double>(0, 0);
        d_down_prime = complex<double>(0, 0);
        for (int ik = 0; ik < site_number; ++ik)
        {
            d_up_prime += g_up_flip(ik, j) * (A_inverse.first(j, ik) + B_up[ik]);
            d_down_prime += g_down_flip(ik, j) * (A_inverse.second(j, ik) + B_down[ik]);
        }

        d_up_prime *= d_up_tilde;
        d_down_prime *= d_down_tilde;
        d_prime = d_up_prime * d_down_prime;
        result += d_prime;
    }

    return result;
}
