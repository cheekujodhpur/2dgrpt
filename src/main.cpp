/*
    This file is part of 2dgrpt.

    2dgrpt is my program to test some fast geodesic algorithms.

    Copyright (c) 2018 by Kumar Ayush
*/
#include "helpers.hpp"
#include "mesh.hpp"
#include <iostream>

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

    return 0;
}
