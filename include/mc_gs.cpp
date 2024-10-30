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
#include "include/set_g.hpp"
#include "include/set_spin.hpp"
#include "include/montecarlo.hpp"
#include "include/lu_decomposition.hpp"

using namespace std;
using namespace Eigen;

// 各種関数定義
int manhattan_distance();                                      // マンハッタン距離の計算
double hopping();                                              // ホッピングパラメータの設定
pair<double, double> wave_vector();                            // 波数点の設定
double epsilon();                                              // 一粒子エネルギーの計算
Matrix4d set_matrix();                                         // ハミルトニアン行列の作成
pair <Vector4d, Matrix4d> diagonalize_matrix();                // 行列の対角化(固有値・ユニタリ行列を取得)
vector<tuple<double, pair<double, double>, int>> eigen_data(); // 波数点に対応した固有値・ξを記憶
int set_spin();                                                // スピンの向きの設定
complex<double> cal_determinant();                             // 行列式を計算
double montecarlo();                                           // monte carlo計算

int main() {
    double H_field = 0.11;
    vector<pair<double, double>> k_vector(site_number / 2);
    pair<double, double> k;
    double epsln;
    Matrix4d A;
    Vector4d eigenvalue_vec;
    Matrix4d Unitary; // ユニタリ行列 U
    vector<tuple<double, pair<double, double>, int>> all_eigendata;
    complex<double> det;
    complex<double> H_xy, H_zz;
    
    

    // kに対応するユニタリ行列 U を保存するマップ
    map<pair<double, double>, Matrix4d> Unitary_map;

    wave_vector(k_vector);

    // 波数ベクトルに対して固有値を計算し、eigendataに蓄積
    for (int i = 0; i < site_number / 2; ++i)
    {

        k = k_vector[i];
        epsln = epsilon(k.first, k.second);
        A = set_matrix(epsln, H_field);
        
        eigenvalue_vec = diagonalize_matrix(A).first; // 固有値を取得
        Unitary_map[k] = diagonalize_matrix(A).second;// ユニタリ行列を取得

        // 各計算結果を all_eigendata に蓄積する
        vector<tuple<double, pair<double, double>, int>> current_eigendata = eigen_data(eigenvalue_vec, k);
        all_eigendata.insert(all_eigendata.end(), current_eigendata.begin(), current_eigendata.end());
    }

    // 固有値でソート
    sort(all_eigendata.begin(), all_eigendata.end(),
        [](const tuple<double, pair<double, double>, int>& a, const tuple<double, pair<double, double>, int>& b) {
            return get<0>(a) < get<0>(b);//固有値が小さい順
        });

    montecarlo(all_eigendata, Unitary_map, H_xy, H_zz);

    cout << H_xy << endl << H_zz ;

    return 0;
}