#pragma once
#include <vector>
#include <utility>
#include <complex>
#include <map>
#include "eigen-3.4.0/Eigen/Dense"
#include "const_para.hpp"
#include "cal_determinant.hpp"
#include "set_g.hpp"
#include "cal_g.hpp"


using namespace std;
using namespace Eigen;

void cal_determinant_mc(const vector<tuple<double, pair<double, double>, int>> eigen_updata, const vector<tuple<double, pair<double, double>, int>> eigen_downdata, const map<pair<double, double>, pair<Matrix2d, Matrix2d>> Unitary_map, int flip_site, vector<int> spin, pair<complex<double>, complex<double>> det_updown_before, pair<MatrixXcd, MatrixXcd> A_inverse, dcomplex& det_after, pair<complex<double>, complex<double>>& det_updown_after, pair<MatrixXcd, MatrixXcd>& A_inv_after)
// (out) det, A_inv_after (in)それ以外
{
    MatrixXcd g_up_normal(site_number, site_number);
    MatrixXcd g_up_flip(site_number, site_number);
    MatrixXcd g_down_normal(site_number, site_number);
    MatrixXcd g_down_flip(site_number, site_number);
    complex<double> det_up_after(0, 0);
    complex<double> det_down_after(0, 0);
    complex<double> det_before(0, 0);
    complex<double> buf_up;
    complex<double> buf_down;


    det_before = det_updown_before.first * det_updown_before.second;
    for (int kappa = 0; kappa < eigen_updata.size(); ++kappa)
    {
        // upの各パラメータ
        double eigenvalue_up = get<0>(eigen_updata[kappa]);
        pair<double, double> k_value_up = get<1>(eigen_updata[kappa]);
        int xi_up = get<2>(eigen_updata[kappa]);

        // downの各パラメータ
        double eigenvalue_down = get<0>(eigen_downdata[kappa]);
        pair<double, double> k_value_down = get<1>(eigen_downdata[kappa]);
        int xi_down = get<2>(eigen_downdata[kappa]);

        // 保存されている U(k) を取得
        Matrix2d Unitary_up;
        Matrix2d Unitary_down;
        try
        {
            Unitary_up = Unitary_map.at(k_value_up).first;
            Unitary_down = Unitary_map.at(k_value_down).second;
        }
        catch(const out_of_range& e)
        {
            cerr << "k_valueがUnitary_mapで見つかりません(cal_determinant_mc.hpp)" << endl;
            exit(-1);
        }


        int flipsite_x = flip_site % L;
        int flipsite_y = flip_site / L;

            if ((spin[flip_site] != 0) && (spin[flip_site] != 1))
        {
            cerr << "スピンの向きに誤りがあります(set_g,hpp)" << endl;
            exit(-1);
        }

        set_g set_g_up(Unitary_up, k_value_up, xi_up, flipsite_x, flipsite_y, spin[flip_site]);
        set_g set_g_down(Unitary_down, k_value_down, xi_down, flipsite_x, flipsite_y, spin[flip_site]);

        g_up_normal(kappa, flip_site) = set_g_up.normal;
        g_up_flip(kappa, flip_site) = set_g_up.flip;
        g_down_normal(kappa, flip_site) = set_g_down.normal;
        g_down_flip(kappa, flip_site) = set_g_down.flip;

        det_up_after += g_up_flip(kappa, flip_site) * A_inverse.first(flip_site, kappa); // 福田修論p46(c.9)式/d
        det_down_after += g_down_flip(kappa, flip_site) * A_inverse.second(flip_site, kappa);
    }
    det_after = det_up_after * det_down_after;
    det_after *= det_before; // (c.9)式
    det_updown_after = make_pair(det_up_after, det_down_after);

    for (int i = 0; i < site_number; ++i)
    {
        if(i == flip_site)
        {
            for (int j = 0; j < site_number; ++j)
            {
                A_inv_after.first(i, j) = (det_updown_before.first / det_up_after) * A_inverse.first(i, j);  // A'_up^-1(i, j) = det_up / det'_up * A_up^-1(i, j)
                A_inv_after.second(i, j) = (det_updown_before.second / det_down_after) * A_inverse.second(i, j);  // A'_down^-1(i, j) = det_down / det'_down * A_down^-1(i, j)
            }
        }
        else
        {
            buf_up = (0, 0);
            buf_down = (0, 0);
            for (int kappa = 0; kappa < site_number; ++kappa)
            {
                buf_up += A_inverse.first(i, kappa) * (g_up_flip(kappa, flip_site) - g_up_normal(kappa, flip_site));
                buf_down += A_inverse.second(i, kappa) * (g_down_flip(kappa, flip_site) - g_down_normal(kappa, flip_site));
            }
            buf_up = (det_updown_before.first / det_up_after) * buf_up;
            buf_down = (det_updown_before.second / det_down_after) * buf_down;
            for (int j = 0; j < site_number; ++j)
            {
                A_inv_after.first(i, j) = A_inverse.first(i, j) - buf_up * A_inverse.first(flip_site, j);
                A_inv_after.second(i, j) = A_inverse.second(i, j) - buf_down * A_inverse.second(flip_site, j);
            }
        }
    }
}