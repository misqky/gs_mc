#pragma once
#include <vector>

#include "const_para.hpp"

using namespace std;

extern vector<pair<int, int>> make_nearest_pair()
{
    vector<pair<int, int>> nearest_pair;
    for (int i; i < site_number; ++i)
    {
        // 右隣のサイト (周期境界条件を考慮)
        if ((i + 1) % L == 0) {
            nearest_pair.push_back({i, i - (L - 1)});  // 行の最後の場合は最初のサイトとペア
        } else {
            nearest_pair.push_back({i, i + 1});  // それ以外は右隣とペア
        }

        // 上隣のサイト (周期境界条件を考慮)
        if (i >= site_number - L) {
            nearest_pair.push_back({i, i - (site_number - L)});  // 最上段の場合は最下段の対応するサイトとペア
        } else {
            nearest_pair.push_back({i, i + L});  // それ以外は上のサイトとペア
        }
    }
    return nearest_pair;
}