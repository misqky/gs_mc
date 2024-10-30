#pragma once
#include "eigen-3.4.0/Eigen/Dense"
#include <complex>

using namespace Eigen;
using namespace std;

// LU分解を行い、L, U行列および行列式を返す関数
void lu_decomposition(const MatrixXcd& g_up, const MatrixXcd& g_down, /*MatrixXcd& L, MatrixXcd& U,*/ pair<MatrixXcd, MatrixXcd>& A_inverse, pair<complex<double>, complex<double>>& determinant) // (in)g,(out)A_inverse,(out)determinant
{
    // LU分解を行う Eigen のさクラスを使う
    FullPivLU<MatrixXcd> lu_decomp_up(g_up);
    FullPivLU<MatrixXcd> lu_decomp_down(g_down);

    // LとUを取得
    // L = lu_decomp.matrixLU().triangularView<Lower>();  // 下三角行列L
    // U = lu_decomp.matrixLU().triangularView<Upper>();  // 上三角行列U

    // Lの対角成分を1にする
    // L.diagonal().setOnes();

    // 行列式を計算
    determinant.first = lu_decomp_up.determinant();
    determinant.second = lu_decomp_down.determinant();

    // 逆行列を計算
    A_inverse.first = lu_decomp_up.inverse();
    A_inverse.second = lu_decomp_down.inverse();
}
