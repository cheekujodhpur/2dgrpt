/*
    This file is part of 2dgrpt.

    2dgrpt is my program to test some fast geodesic algorithms.

    Copyright (c) 2018 by Kumar Ayush
*/

#pragma once

#include <stdlib.h>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Geometry>

using namespace Eigen;

namespace grpt{

/// Returns a random point inside the triangle formed by three points
Vector2d random_point(const Vector2d p1, const Vector2d p2, const Vector2d p3);

}
