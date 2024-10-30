#pragma once
#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <vector>
#include <utility>
#include <algorithm>
#include <complex>
#include "eigen-3.4.0/Eigen/Dense"
#include "const_para.hpp"

using namespace std;
using namespace Eigen;

// set_g クラスの定義
class set_g
{
public:
    complex<double> normal;  // スピンに応じた g.normal
    complex<double> flip;    // スピン反転に応じた g.flip

    // コンストラクタ：gの計算に必要な情報を引数として受け取る
    set_g(const Matrix2d& U, pair<double, double> k, int xi, double x, double y, int spin)
    {
        // スピンの値が0または1でなければエラー
        if ((spin != 0) && (spin != 1))
        {
            cerr << "スピンの向きに誤りがあります(set_g,hpp)" << endl;
            exit(-1);
        }

        // g.normal と g.flip の計算をそれぞれのメソッドで行う
        normal = calculate_g(U, k, xi, x, y, spin);       // 通常のスピンでの g
        flip = calculate_g(U, k, xi, x, y, 1 - spin);     // スピンが反転した場合の g
    }

private:
    // gを計算する関数
    complex<double> calculate_g(const Matrix2d& U, pair<double, double> k, int xi, double x, double y, int spin)
    {
        // Uの要素とdoubleって g を計算
        complex<double> g;
        g = (1.0 / sqrt(site_number)) * (U(xi, 0) + U(xi, 1) * exp(I * M_PI * (x + y)))
            * exp(I * (k.first * x + k.second * y));



        // 許容誤差を用いて十分小さい値を 0 とみなす
        double real_part = real(g);
        double imag_part = imag(g);

        if (abs(real_part) < zero_judge)
        {
            real_part = 0;
        }
        if (abs(imag_part) < zero_judge)
        {
            imag_part = 0;
        }

        g = complex<double> (real_part, imag_part);
        return g;
    }
};
