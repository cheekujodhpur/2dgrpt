/*
    This file is part of 2dgrpt.

    2dgrpt is my program to test some fast geodesic algorithms.

    Copyright (c) 2018 by Kumar Ayush
*/

#pragma once

#include <algorithm>
#include <array>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <stdlib.h>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Geometry>

using namespace Eigen;

namespace grpt{

    /// Forward Declarations
    class Metric;
    class Mesh;

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

	#undef M_PI

	#define M_PI         3.14159265358979323846
    /// Cast angles to [0, 2pi)
    inline double anglemod(const double ang) {
        double modded = fmod(2*M_PI+ang, 2*M_PI);
        return modded;
    }

    /// To execute a command
    inline std::string exec(const char* cmd) {
        std::array<char, 128> buffer;
        std::string result;
        std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
        if (!pipe) throw std::runtime_error("popen() failed!");
        while (!feof(pipe.get())) {
            if (fgets(buffer.data(), 128, pipe.get()) != nullptr)
                result += buffer.data();
        }
        return result;
    } 

    /// Calculate area of a triangle
    double calculate_area(const Vector2d p1, const Vector2d p2, 
            const Vector2d p3);


    /// Returns a random point inside the triangle formed by three points
    Vector2d random_point(const Vector2d p1, const Vector2d p2, 
            const Vector2d p3);

    /// Return the distance from the closest edge in a triangle
    double return_incidence(const std::vector<Vector2d> vertices, 
            const std::vector<int> triangle,
            const Vector2d point_of_intersection);

    /// Find the closest edge, use overlay mesh (ovmesh) to use a 
    /// subset of edges
    /// rs is the refined size of the overlay mesh 
    bool find_closest_edge(const std::unordered_map<std::vector<int>, 
            std::vector<std::vector<int>>, vector_int_hasher> ovmesh, 
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
    bool bresenham_and_mesh(std::unordered_map<std::vector<int>, 
            std::vector<std::vector<int>>, vector_int_hasher> &ovm, 
            const double x1, const double y1, 
            const double x2, const double y2, 
            const int a, const int b);

    /// nelder mead optimization for 2-simplex
    template <typename Lambda>
        double nelder_mead(double x1, double x2, const double tau, 
                Lambda&& f) {

            // The parameters for nelder_mead
            double alpha = 1.0;
            double gamma = 2;
            double rho = 0.5;
            double sigma = 0.5;

            // TODO: define this globally, or maybe as a parameter
            uint MAX_ITER = 50;
            uint iter_idx = 0;

            double x0, xr, xe, xc; 

            while(iter_idx++ < MAX_ITER) {

                // Step 1: Order
                if (f(x1) > f(x2)){
                    std::swap(x1, x2);
                }

                // TODO: Figure out this number 1e-1 to be something better
                if (fabs(f(x2)-f(x1)) < tau && fabs(x2-x1) < 1e-1)
                    return 0.5*(x1+x2);

                // Step 2: Centroid
                x0 = x1;

                // Step 3: Reflection
                xr = x0 + alpha*(x0-x2);

                if (f(xr) < f(x1)) {
                    xe = x0 + gamma*(xr - x0);
                    if (f(xe) < f(xr)){
                        x2 = xe;
                        continue;
                    }
                    else {
                        x2 = xr;
                        continue;
                    }
                }

                xc = x0 + rho*(x2-x0);
                if (f(xc) < f(x2)) {
                    x2 = xc;
                    continue;
                }

                x2 = x1 + sigma*(x2 - x1);
            }

            //TODO: Add exxception "[WARN]: max iter reached in nelder mead"
            return 0.5*(x1+x2);

        }

}
