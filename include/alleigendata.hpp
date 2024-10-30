#pragma once
#include <utility>
#include <map>
#include <algorithm>
#include <vector>
#include <tuple>
#include "eigen-3.4.0/Eigen/Dense"

#include "wave_vector_MBZ.hpp"
#include "epsilon.hpp"
#include "set_matrix.hpp"
#include "const_para.hpp"
#include "diagonalize_matrix.hpp"
#include "eigen_data.hpp"

using namespace std;
using namespace Eigen;

extern  map<pair<double, double>, pair<Matrix2d, Matrix2d>> all_eigendata
(
    vector<tuple<double, pair<double, double>, int>>& eigen_updata,
    vector<tuple<double, pair<double, double>, int>>& eigen_downdata
)
{
    double H_field = 0.11;
    vector<pair<double, double>> k_vector(site_number / 2);
    pair<double, double> k;
    double epsln;
    Matrix2d A_up, A_down;
    pair<Vector2d, Vector2d> eigenvalues; // (eigenvalue_up, eigenvalue_down)
    map<pair<double, double>, pair<Matrix2d, Matrix2d>> Unitary_map; // (k, U_up(k), U_down(k))

    // 波数点の生成
    wave_vector(k_vector);

    for (int i = 0; i < site_number / 2; ++i)
    {
        k = k_vector[i];
        epsln = epsilon(k.first, k.second);
        set_matrix(epsln, H_field, A_up, A_down); // A_up, A_downを生成

        // 固有値とユニタリ行列を計算
        diagonalize_matrix(k, A_up, A_down, eigenvalues, Unitary_map);

        // 固有データの保存
        eigen_data(eigenvalues, k, eigen_updata, eigen_downdata);
    }

    // 固有値を小さい順にソート
    auto cmd = [](const tuple<double, pair<double, double>, int>& a, const tuple<double, pair<double, double>, int>& b)
    {
        return get<0>(a) < get<0>(b); // 固有値が小さい順にソート
    };

    // ソート処理
    sort(eigen_updata.begin(), eigen_updata.end(), cmd);
    sort(eigen_downdata.begin(), eigen_downdata.end(), cmd);

    return Unitary_map;
}
