/*
    This file is part of 2dgrpt.

    2dgrpt is my program to test some fast geodesic algorithms.

    Copyright (c) 2018 by Kumar Ayush
*/
#include "helpers.hpp"
#include "mesh.hpp"
#include <iostream>

/// Number of samples on each edge
#define NOS_EDGE 25

using namespace grpt;

int main(int argc, char **argv) {

    std::vector<Vector2d> corners = {Vector2d(0,0),
                                     Vector2d(0,5),
                                     Vector2d(5,0),
                                     Vector2d(5,5)};
    Mesh myMesh = Mesh(corners, Vector2d(0,0));
    myMesh.submesh(3);
    myMesh.write_to_poly(0.01);
    myMesh.refine_using_Triangle(18, 0.05); 

    // myMesh.print();
    // myMesh.draw();
    // myMesh.find_trial_error(0,0.01,0,Vector2d(0.55,0.13), Vector2d(0.7,0.7));
    
#pragma omp parallel for schedule(dynamic, 1)
    for(int i=0;i<NOS_EDGE*3*myMesh.get_num_triangles();++i){
        std::cerr <<  "Marking" << i << "..." << std::endl;
        myMesh.throw_geodesic_mark(i, 1e-4, 1e-2);
        std::cerr << "done." << std::endl;
    }

    myMesh.churn_edge_data();

    std::cout << "Now throwing discrete geodesics..." << std::endl;
    for(int i=0;i<50;i++) {
        myMesh.throw_geodesic_discrete(i);
    }

    return 0;
}
