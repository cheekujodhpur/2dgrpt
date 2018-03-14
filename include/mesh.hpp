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
        std::vector<std::vector<int>> neighbours;

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

        /// Fetch number of triangles
        int get_num_triangles(void){return triangles.size();};

        /// Subdivides the mesh
        void submesh(const int n);

        /// Write to a poly file which CMU Triangle can process
        void write_to_poly(const double k);

        /// Give a call to Triangle binary and refine our mesh
        /*
         * We could use triangle directly into our program, but it is written
         * in C, and uses different data structures than us, and it would
         * take time to convert all of it. It is not a huge difference if
         * we write to file, call triangle and read from file versus directly
         * casting our data into their structures
         */
        void refine_using_Triangle(const double min_angle,
                const double overlay_size);

        /// Converting deviation data to interpolators wrt input angle
        void churn_edge_data();

        /// Plot the mesh
        void draw();

        /// Print all information
        void print();

        /// Find trial error
        double find_trial_error(const double mod, const double dt, 
                const int t_id, const Vector2d startpoint, 
                const Vector2d covdir, bool save=false, bool dbg=false);

        /// Throw a geodesic and mark the data on edge
        void throw_geodesic_mark(const int seed, const double tau, 
                double dt=0.01);

        /// Actually throw a discrete geodesic to check state
        void throw_geodesic_discrete(int seed);


        /// Some helper functions
        bool find_intersection(std::vector<int> this_triangle, Vector2d startpoint, 
                Vector2d direction, std::vector<int> &local_edge, double &t);
        
        /// Check if a point lies on a given triangle or not
        bool check_for_incidence(const std::vector<int> triangle,
                const Vector2d point_of_intersection, const double threshold,
                std::vector<int> &local_edge);

        /// Sample from optimiser function
        void sample_optimiser(const double dt);

    };

}
