/*
    This file is part of 2dgrpt.

    2dgrpt is my program to test some fast geodesic algorithms.

    Copyright (c) 2018 by Kumar Ayush
*/

#pragma once

#include <stdlib.h>
#include <vector>
#include <unordered_map>
#include <algorithm>

#include <Eigen/Dense>
#include <Eigen/Geometry>

using namespace Eigen;

namespace grpt{

    /// Hash function for a std::vector<int>
    class vector_int_hasher {
        public:
            std::size_t operator()(std::vector<int> const& vec) const {
                std::size_t seed = vec.size();
                for(auto& i : vec) {
                    seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
                }
                return seed;
            }
    };

    /// Cross product of 2d vectors
    inline double cross(const Vector2d a, const Vector2d b) {
        return a.x()*b.y()-a.y()*b.x();
    }

    /// Returns a random point inside the triangle formed by three points
    Vector2d random_point(const Vector2d p1, const Vector2d p2, 
            const Vector2d p3);

    /// Check if a point lies on a given triangle or not
    bool check_for_incidence(const std::vector<Vector2d> vertices, 
            const std::tuple<int> triangle,
            const Vector2d point_of_intersection, const double threshold);

    /// Return the distance from the closest edge in a triangle
    double return_incidence(const std::vector<Vector2d> vertices, 
            const std::vector<int> triangle,
            const Vector2d point_of_intersection);

    /// Find the closest edge, use overlay mesh (ovmesh) to use a 
    /// subset of edges
    /// rs is the refined size of the overlay mesh 
    bool find_closest_edge(const std::unordered_map<std::vector<int>, 
            std::vector<std::vector<int>>> ovmesh, 
            const Vector2d pt, const std::vector<Vector2d> vertices, 
            const double rs, 
            std::vector<int> &min_edge, double &distance);

    /// Find the edge tentatively in the direction specified, 
    /// use overlay mesh (ovmesh) to use a subset of edges
    /// rs is the refined size of the overlay mesh 
    bool find_shooting_edge(const std::unordered_map<std::vector<int>, 
            std::vector<std::vector<int>>, vector_int_hasher> ovmesh, 
            const Vector2d pt, const Vector2d direction, 
            const std::vector<Vector2d> vertices, const double rs, 
            std::vector<int> &min_edge, double &mint);

    /// ovm is the overlay mesh, a and b are vertex ids
    /// xi and yi are the starting and endpoints obviously
    /// This functions stores which edges lie in which overlay mesh cells 
    bool bresenham_and_mesh(const std::unordered_map<std::vector<int>, 
            std::vector<std::vector<int>>, vector_int_hasher> ovm, 
            const double x1, const double y1, 
            const double x2, const double y2, 
            const int a, const int b);

}
