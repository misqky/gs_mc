#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <map>
#include <vector>
#include <complex>
#include <tuple>
#include <cstdlib>
#include <utility>

#include "eigen-3.4.0/Eigen/Dense"
#include "set_spin.hpp"
#include "cal_determinant.hpp"
#include "set_g.hpp"
#include "cal_g.hpp"
#include "cal_determinant_mc.hpp"
#include "lu_decomposition.hpp"
#include "cal_flipweight.hpp"

using namespace std;
using namespace Eigen;

void montecarlo
(
    const vector<tuple<double, pair<double, double>, int>> eigen_updata,
    const vector<tuple<double, pair<double, double>, int>> eigen_downdata,
    const map<pair<double, double>, pair<Matrix2d, Matrix2d>> Unitary_map,
    complex<double>& H_xy,
    complex<double>& H_zz
)
{
    vector<int> spin;
    vector<int> spin_temp;
    Matrix2d Unitary_up;
    Matrix2d Unitary_down;
    MatrixXcd g_up(site_number, site_number);
    MatrixXcd g_up_flip(site_number, site_number);
    MatrixXcd g_down(site_number, site_number );
    MatrixXcd g_down_flip(site_number, site_number);
    pair<MatrixXcd, MatrixXcd> A_inverse;
    pair<MatrixXcd, MatrixXcd> A_inv_after;
    pair<complex<double>, complex<double>> determinant_1;
    pair<complex<double>, complex<double>> determinant_2;
    pair<complex<double>, complex<double>> det_updown_mc;
    pair<complex<double>, complex<double>> determinant_temp;
    complex<double> determinant_mc;
    complex<double> det_1;
    complex<double> det_2;
    double weight_1;
    double weight_2;
    int flip_site;
    double weight_ratio;
    double transition_probably;
    vector<double> S_z(site_number);
    double ising;
    complex<double> flip_weight;
    // complex<double> H_xy(0,0), H_zz(0,0);



    srand((unsigned)time(NULL));
    for (int i = 0; i < site_number; ++i)
    {
        // set_spin();
        spin.push_back(set_spin());
    }


    // gを計算して行列として設定
    cal_g(eigen_updata, spin, g_up, g_up_flip, Unitary_map, true);
    cal_g(eigen_downdata, spin, g_down, g_down_flip, Unitary_map, false);


    lu_decomposition(g_up, g_down, A_inverse, determinant_1);

    det_1 = determinant_1.first * determinant_1.second;
    weight_1 = pow(abs(det_1), 2);

    for (int step = 0; step < step_max; ++step)
    {
        spin_temp = spin;
        flip_site = static_cast<int>((site_number - 1) * (rand() / static_cast<double>(RAND_MAX)));



        spin_temp[flip_site] = 1 - spin[flip_site];   // 元のスピンのうちフリップサイト番だけ反転
        weight_2 = 0.0;

        if(step % 100 == 0)
        {
            cal_g(eigen_updata, spin, g_up, g_up_flip, Unitary_map, true);
            cal_g(eigen_downdata, spin, g_down, g_down_flip, Unitary_map, false);


            lu_decomposition(g_up, g_down, A_inv_after, determinant_2);
            det_2 = determinant_2.first * determinant_2.second;
            determinant_temp = determinant_2;
            //cout << "det2; " << det_2 << endl;
        }
        else
        {
            cal_determinant_mc(eigen_updata, eigen_downdata, Unitary_map, flip_site, spin, determinant_1, A_inverse, determinant_mc, det_updown_mc, A_inv_after);
            det_2 = determinant_mc;
            determinant_temp = det_updown_mc;
        }
        weight_2 = pow(abs(det_2), 2);

        weight_ratio = weight_1 / weight_2;
        transition_probably = 1.0 / (1.0 + weight_ratio);


        if(transition_probably >= rand() / RAND_MAX)
        {
            spin = spin_temp;
            det_1 = det_2;
            determinant_1 = determinant_temp;
            weight_1 = weight_2;
            A_inverse = A_inv_after;

            cal_g(eigen_updata, spin, g_up, g_up_flip, Unitary_map, true);
            cal_g(eigen_downdata, spin, g_down, g_down_flip, Unitary_map, false);

        }
        for (int j = 0; j < site_number; ++j)
        {
            S_z[j] = (double) spin[j] - 0.5;
        }

        ising = 0.0;

        for (int i = 0; i < site_number; ++i)
        {
            for (int j = i; j < site_number; ++j)
            {
                ising += S_z[i] * S_z[j];
            }
        }
        // cout << A_inverse << endl << endl;
        /*cout << "g_up:" << endl << g_up << endl << endl;
        cout << "g_down:" << endl << g_down << endl << endl;
        cout << "A_inverse first (up):" << endl << A_inverse.first << endl << endl;
        cout << "A_inverse second (down):" << endl << A_inverse.second << endl << endl;*/

        flip_weight = cal_flipweight(spin, g_up, g_up_flip, g_down, g_down_flip, determinant_1, A_inverse);

        H_xy += 1./2 * conj(flip_weight) * det_1 / weight_1;
        H_zz += ising * conj(det_1) * det_1 / weight_1;
        // 続きはここから

        // cout << weight_1 << endl;
    }
    H_xy /= step_max;  //ステップ数で平均
    H_zz /= step_max;  //ステップ数で平均
}
