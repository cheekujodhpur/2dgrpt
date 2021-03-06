/*
    This file is part of 2dgrpt.

    2dgrpt is my program to test some fast geodesic algorithms.

    Copyright (c) 2018 by Kumar Ayush
*/

#pragma once

#include <string>
#include "helpers.hpp"

namespace grpt {
    /*
     * Metric properties and calculations
     */
    class Metric {

        // TODO: Maybe make the whole class general dimensional
        Vector2d origin;

    public:
        Metric(){origin=Vector2d(0,0);};
        Metric(const Vector2d _origin):origin(_origin){};
        Matrix2d compute_metric(const Vector2d position,
                const std::string type); 

        void operator()( const std::vector<double> &x, 
                std::vector<double> &dxdt, 
                const double /* t */);


    };
}
