/*
    This file is part of 2dgrpt.

    2dgrpt is my program to test some fast geodesic algorithms.

    Copyright (c) 2018 by Kumar Ayush
*/

#pragma once

#include "helpers.hpp"
#include "metric.hpp"

namespace grpt {

    class Mesh {

        /// The list of vertices
        std::vector<Vector2d> vertices;

        /// The list of vertices, each element specifies three vertices
        std::vector<std::vector<int>> triangles;

        /// The edge list, each element specifies two vertices
        std::vector<std::vector<int>> edges;

        /// For CMU triangle, need to store a region center and a property
        /// This property was a scaling of max area of triangle for us
        std::vector<std::tuple<Vector2d, double>> regions;

        /// Metric class
        Metric metric;

        /// Indices of neighbouring triangles in the list of triangles
        std::vector<int> neighbours;

        /// Overlay mesh (effectively kd tree)
        std::unordered_map<std::vector<int>, std::vector<std::vector<int>>, vector_int_hasher> overlay_mesh;

        /// Storing corners of the mesh
        std::vector<Vector2d> corners;

        /// Storing the origin of the mesh, the corners and other coordinates
        /// are translated with respect to this
        Vector2d origin;

        /// Stores a list of samples for an edge
        /// Each sample is a list of 2 vectors
        /// First vector tells the deviation,
        /// Second tells the initial direction
        std::unordered_map<std::vector<int>, std::vector<std::vector<Vector2d>>,
            vector_int_hasher> edge_data;

        /// Stores a list of samples for an edge
        /// This is used for interpolation
        /// Each sample is a list of 3 doubles
        /// First vector tells the lower angle of angle interval
        /// Second and third are slope and intercept resp.
        /// You can use these values to interpolate the output for an angle
        /// that falls in the resp. interval
        std::unordered_map<std::vector<int>, std::vector<std::vector<double>>,
            vector_int_hasher> edge_slope_data;

    public:
        
        Mesh(const std::vector<Vector2d> _corners, const Vector2d _origin);
    };

}