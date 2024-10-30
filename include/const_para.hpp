#pragma once
#include <complex>

using namespace std;

extern const int L = 4; // 系のサイズ (L x L)
extern const int site_number = L * L;
extern const double zero_judge = 1e-14;
extern const complex<double> I(0, 1); // 複素数i
extern const int step_max = 10000;
extern const int iterations = 100;