#pragma once
#include <iostream>
#include "eigen-3.4.0/Eigen/Dense"
#include "const_para.hpp"
#include "wave_vector.hpp"
#include "epsilon.hpp"
#include "set_matrix.hpp"
#include "diagonalize_matrix.hpp"
#include "eigen_data.hpp"

using namespace std;
using namespace Eigen;

extern
for (int i = 0; i < site_number / 2; ++i)
    {
        k = wave_vector(i);
        epsln = epsilon(wave_vector(i).first, wave_vector(i).second);
        A = set_matrix(epsln, H_field);

        eigenvalue_vec = diagonalize_matrix(A).first; // 固有値を取得
        Unitary_map[k] = diagonalize_matrix(A).second;// ユニタリ行列を取得

        // 各計算結果を all_eigendata に蓄積する
        vector<tuple<double, pair<double, double>, int>> current_eigendata = eigen_data(eigenvalue_vec, k);
        all_eigendata.insert(all_eigendata.end(), current_eigendata.begin(), current_eigendata.end());
    }