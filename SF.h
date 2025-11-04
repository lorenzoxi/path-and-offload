#ifndef INC_005_CCNC_DRONE_ENERGY_OPT_SF_H
#define INC_005_CCNC_DRONE_ENERGY_OPT_SF_H

#include <vector>
#include "SpreadingFactor.h"

static const std::vector<SpreadingFactor> SF125 = {
        {7,  125.0, 5.468, 0,  6},
        {8,  125.0, 3.125, 4,  8},
        {9,  125.0, 1.758, 6,  11},
        {10, 125.0, 0.977, 8,  14},
        {11, 125.0, 0.537, 10, 17},
        {12, 125.0, 0.293, 12, 20}
};

static const std::vector<SpreadingFactor> SF250 = {
        {7,  250.0, 10.937, 0, 4},
        {8,  250.0, 6.25,   3, 6},
        {9,  250.0, 3.516,  5, 8},
        {10, 250.0, 1.953,  6, 11},
        {11, 250.0, 1.074,  8, 13},
        {12, 250.0, 0.586,  9, 14}
};

static const std::vector<SpreadingFactor> SF500 = {
        {7,  500.0, 21.875, 1, 3},
        {8,  500.0, 12.5,   2, 4},
        {9,  500.0, 7.031,  3, 5},
        {10, 500.0, 3.906,  4, 7},
        {11, 500.0, 2.148,  5, 9},
        {12, 500.0, 1.171,  6, 10}
};


#endif
