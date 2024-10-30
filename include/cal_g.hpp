#pragma once
#include <iostream>
#include <utility>
#include <vector>
#include <tuple>

#include "const_para.hpp"
#include "set_g.hpp"
#include "eigen-3.4.0/Eigen/Dense"

using namespace std;
using namespace Eigen;

void cal_g(const vector<tuple<double, pair<double, double>, int>> eigen_data, vector<int> spin, MatrixXcd& g_normal, MatrixXcd& g_flip, const map<pair<double, double>, pair<Matrix2d, Matrix2d>> Unitary_map, bool up_or_down)
{
    for (int kappa = 0; kappa < eigen_data.size(); ++kappa)
    {
        double eigenvalue = get<0>(eigen_data[kappa]);
        pair<double, double> k_value = get<1>(eigen_data[kappa]);
        int xi = get<2>(eigen_data[kappa]);

        // Uの計算
        Matrix2d Unitary;
        // const mapに[k_value]ではなく.at(k_value)しか使えないのでそれに伴ってtry&catchを用いる
        try
        {
            if (up_or_down)  //eigen_updata is true
            {
                Unitary = Unitary_map.at(k_value).first;  // eigen_updata の場合
            }
            else             //eigen_downdata is false
            {
                Unitary = Unitary_map.at(k_value).second; // eigen_downdata の場合
            }
        }
        catch(const out_of_range& e)
        {
            cerr << "k_valueがUnitary_mapで見つかりません: (" << k_value.first << ", " << k_value.second << ")(cal_g.hpp)" << endl;
            exit(-1);
        }


        // g の計算
        for (int i = 0; i < site_number; ++i)
        {
            int site_x = i % L;
            int site_y = i / L;

            set_g set_g(Unitary, k_value, xi, site_x, site_y, spin[i]); //クラスを設定(.normal(通常版).flip(スピン反転版))
            g_normal(kappa, i) = set_g.normal;
            g_flip(kappa, i) = set_g.flip;
        }
    }
}