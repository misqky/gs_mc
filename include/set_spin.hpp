#pragma once
#include <cstdlib>
#include <time.h>

using namespace std;

extern int set_spin()
{
    return rand() % 2;
}