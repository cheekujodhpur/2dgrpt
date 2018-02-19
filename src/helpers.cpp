/*
    This file is part of 2dgrpt.

    2dgrpt is my program to test some fast geodesic algorithms.

    Copyright (c) 2018 by Kumar Ayush
*/

#include "helpers.hpp"

using namespace grpt;

Vector2d random_point(const Vector2d p1, const Vector2d p2, const Vector2d p3){
    unsigned short int seed[] = {(short unsigned int)(p1.x()*1e5), (short unsigned int)10005, (short unsigned int)99956};
    double c1 = erand48(seed);
    double c2 = erand48(seed);
    double c3 = erand48(seed);
    return (p1*c1 + p2*c2 + p3*c3)/(c1+c2+c3);
}
