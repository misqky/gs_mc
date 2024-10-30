#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <map>
#include <algorithm>
#include <tuple>
#include <cstdlib>
#include <time.h>
//#include <eigen-3.4.0/Eigen/Dense>
#include "include/eigen-3.4.0/Eigen/Dense"

#include "include/const_para.hpp"
#include "include/manhattan_distance.hpp"
#include "include/hopping_para.hpp"
#include "include/wave_vector_MBZ.hpp"
#include "include/epsilon.hpp"
#include "include/set_matrix.hpp"
#include "include/diagonalize_matrix.hpp"
#include "include/eigen_data.hpp"
#include "include/alleigendata.hpp"
#include "include/set_g.hpp"
#include "include/set_spin.hpp"
#include "include/montecarlo.hpp"
#include "include/lu_decomposition.hpp"

using namespace std;
using namespace Eigen;

// 各種関数定義
int manhattan_distance();                                           // マンハッタン距離の計算
double hopping();                                                   // ホッピングパラメータの設定
pair<double, double> wave_vector();                                 // 波数点の設定
double epsilon();                                                   // 一粒子エネルギーの計算
Matrix4d set_matrix();                                              // ハミルトニアン行列の作成
pair <Vector4d, Matrix4d> diagonalize_matrix();                     // 行列の対角化(固有値・ユニタリ行列を取得)
vector<tuple<double, pair<double, double>, int>> eigen_data();      // 波数点に対応した固有値・ξを記憶
map<pair<double, double>, pair<Matrix2d, Matrix2d>> all_eigendata();//
int set_spin();                                                     // スピンの向きの設定
complex<double> cal_determinant();                                  // 行列式を計算
double montecarlo();                                                // monte carlo計算
complex<double> cal_flipweight();

int main()
{
    vector<pair<double, double>> k_vector(site_number / 2);
    pair<double, double> k;
    double epsln;
    Matrix2d A_up, A_down;
    Vector4d eigenvalue_vec;
    Matrix4d Unitary; // ユニタリ行列 U
    complex<double> det;
    complex<double> H_xy, H_zz, H_total;
    double Energy_ave;
    // kに対応するユニタリ行列 U を保存するマップ
    map<pair<double, double>, pair<Matrix2d, Matrix2d>> Unitary_map;
    vector<tuple<double, pair<double, double>, int>> eigen_updata;
    vector<tuple<double, pair<double, double>, int>> eigen_downdata;

    Unitary_map = all_eigendata(eigen_updata, eigen_downdata);

for(int it = 0; it < iterations; ++it)
    {
        H_xy = complex<double>(0,0);
        H_zz = complex<double>(0,0);
        montecarlo(eigen_updata, eigen_downdata, Unitary_map, H_xy, H_zz);

        cout << "H_xy: " << H_xy << endl << "H_zz: " << H_zz << endl;;
        H_total = H_xy + H_zz;
        cout << "H_total: " << H_total << endl;
        cout << "基底エネルギー: " << real(H_total) / site_number << endl;
        Energy_ave += real(H_total) / site_number;
    }
    Energy_ave /= iterations;
    cout << "平均基底エネルギー: " << Energy_ave << endl;
    return 0;
}