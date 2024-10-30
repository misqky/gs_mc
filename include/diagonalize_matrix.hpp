#pragma once
#include <iostream>
#include <utility>
#include <map>
#include <vector>
#include "eigen-3.4.0/Eigen/Dense"

#include "epsilon.hpp"
#include "set_matrix.hpp"
#include "const_para.hpp"

using namespace std;
using namespace Eigen;

void diagonalize_matrix
(
    const pair<double, double>& k,
    const Matrix2d& A_up,
    const Matrix2d& A_down,
    pair<Vector2d, Vector2d>& eigenvalues,           // (eigenvalue_up, eigenvalue_down)
    map<pair<double, double>, pair<Matrix2d, Matrix2d>>& Unitary_map // (k, U_up(k), U_down(k))
)
{
    // 固有値と固有ベクトルを計算
    SelfAdjointEigenSolver<Matrix2d> eigen_solver_up(A_up);
    SelfAdjointEigenSolver<Matrix2d> eigen_solver_down(A_down);

    // エラーがあれば終了
    if ((eigen_solver_up.info() != Success) || (eigen_solver_down.info() != Success))
    {
        cerr << "固有値分解に失敗しました！ (diagonalize_matrix/hpp)" << endl;
        exit(-1);
    }

    // 固有値と固有ベクトルを保存
    eigenvalues.first = eigen_solver_up.eigenvalues(); 
    eigenvalues.second = eigen_solver_down.eigenvalues(); 

    // 固有ベクトル（ユニタリ行列）をマップに保存
    Unitary_map[k] = make_pair(eigen_solver_up.eigenvectors(), eigen_solver_down.eigenvectors());
}
