/*
    This file is part of 2dgrpt.

    2dgrpt is my program to test some fast geodesic algorithms.

    Copyright (c) 2018 by Kumar Ayush
*/

#include "metric.hpp"

using namespace grpt;

//TODO: Overload () operator to do this compute metric
Matrix2d Metric::compute_metric(const Vector2d position,
        const std::string type) {

    Matrix2d g = Matrix2d::Zero();
    Vector2d v = position-origin;

    if (type=="polar") {
        g(0,0) = 1;
        g(1,1) = pow(sin(v[0]*M_PI),2.0);
    }
    else if (type=="a1") {
        g(0,0) = 1;
        g(1,1) = v[0];   
    }
    else if (type=="a2") {
        g(0,0) = 1;
        g(1,1) = 1./v[0];
    }

    return g;
}
