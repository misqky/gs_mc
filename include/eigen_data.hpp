#pragma once
#include <vector>
#include <utility>
#include <tuple>
#include <algorithm>
#include "eigen-3.4.0/Eigen/Dense"

using namespace std;
using namespace Eigen;

void eigen_data
(
    const pair<Vector2d, Vector2d>& eigenvalues,
    const pair<double, double>& k,
    vector<tuple<double, pair<double, double>, int>>& eigen_updata,
    vector<tuple<double, pair<double, double>, int>>& eigen_downdata
)
{
    // 固有値をデータに追加
    for (int xi = 0; xi < 2; ++xi) {
        // eigenvalues.firstは上スピンの固有値ベクトル、eigenvalues.secondは下スピンの固有値ベクトル
        eigen_updata.push_back(make_tuple(eigenvalues.first[xi], k, xi));  // 上スピン
        eigen_downdata.push_back(make_tuple(eigenvalues.second[xi], k, xi));  // 下スピン
    }
}
