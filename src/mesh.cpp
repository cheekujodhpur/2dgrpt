/*
    This file is part of 2dgrpt.

    2dgrpt is my program to test some fast geodesic algorithms.

    Copyright (c) 2018 by Kumar Ayush
*/

#include "mesh.hpp"
#include <algorithm>
#include <iostream>
#include <fstream>

using namespace grpt;

Mesh::Mesh(const std::vector<Vector2d> _corners, const Vector2d _origin) {

    corners = _corners;
    origin = _origin;
        
    //TODO: Figure out how to deal with default dict behaviour
    //self.edge_data = defaultdict(list)
    //self.edge_slope_data = defaultdict(list)

    // NOTE: assumes a vector of dimension 2 
    Vector2d center = std::accumulate(corners.begin(), corners.end(), 
            Vector2d(0,0))/corners.size();
    std::vector<float> angles;
    for (auto corner : corners) {
        Vector2d v_to_corner = corner-center;
        angles.push_back(atan2(v_to_corner.y(), v_to_corner.x()));
    }

    std::vector<std::pair<Vector2d, float>> sorting_pair;
    // assumes angles and corners have same size
    for(uint i = 0;i < corners.size(); i++){
        sorting_pair.push_back(
                std::pair<Vector2d, float>(corners[i],angles[i]));
    }

    auto sorting_pair_op = [](std::pair<Vector2d, float> i,
            std::pair<Vector2d, float> j){ return (i.second<j.second); };

    std::sort(sorting_pair.begin(), sorting_pair.end(), sorting_pair_op);

    std::transform(sorting_pair.begin(), sorting_pair.end(),
            corners.begin(),
            [] (std::pair<Vector2d, float> p) -> Vector2d {
                return p.first;});


    //TODO: maybe find a better way to do this
    std::vector<int> tmp = {0,1,2};
    std::vector<int> tmp2 = {0,2,3};
    triangles.push_back(tmp);
    triangles.push_back(tmp2);

    // Writing triangles is easy, the main purpose was always to sort the corners
    vertices.insert(vertices.end(), corners.begin(), corners.end());

    metric = Metric(origin);
}

void Mesh::submesh(const int n) {
        Metric &m = metric;
        int index = vertices.size();
        std::vector<std::vector<int>> new_triangles;

        for (int nit = 0;nit < n;++nit) {
            new_triangles.clear();
            for (auto triangle : triangles) {
                Vector2d p1 = (vertices[triangle[0]] + 
                        vertices[triangle[1]])/2.;

                Vector2d p2 = (vertices[triangle[1]] + 
                        vertices[triangle[2]])/2.;

                Vector2d p3 = (vertices[triangle[2]] + 
                        vertices[triangle[0]])/2.;

                int i1, i2, i3;
                std::vector<Vector2d>::iterator pit;

                pit = std::find(vertices.begin(), vertices.end(), p1);
                if ( pit != vertices.end() ) {
                    i1 = (uint)(pit-vertices.begin());
                }
                else {
                    vertices.push_back(p1);
                    i1 = index;
                    index++;
                }

                pit = std::find(vertices.begin(), vertices.end(), p2);
                if ( pit != vertices.end() ) {
                    i2 = (uint)(pit-vertices.begin());
                }
                else {
                    vertices.push_back(p2);
                    i2 = index;
                    index++;
                }

                pit = std::find(vertices.begin(), vertices.end(), p3);
                if ( pit != vertices.end() ) {
                    i3 = (uint)(pit-vertices.begin());
                }
                else {
                    vertices.push_back(p3);
                    i3 = index;
                    index++;
                }
                std::vector<int> tmp = {i1, i2, i3};
                new_triangles.push_back(tmp);
                tmp.clear();
                tmp.push_back(i1);tmp.push_back(triangle[1]);tmp.push_back(i2);
                new_triangles.push_back(tmp);
                tmp.clear();
                tmp.push_back(i2);tmp.push_back(triangle[2]);tmp.push_back(i3);
                new_triangles.push_back(tmp);
                tmp.clear();
                tmp.push_back(i3);tmp.push_back(triangle[0]);tmp.push_back(i1);
                new_triangles.push_back(tmp);
            }
            triangles = new_triangles;
        }

        for (auto triangle : triangles) {
            std::vector<int> striangle = triangle;
            std::sort(striangle.begin(), striangle.end());
            
            std::vector<int> tmp;
            tmp.push_back(striangle[0]); tmp.push_back(striangle[1]);
            edges.push_back(tmp);
            tmp.clear();
            tmp.push_back(striangle[1]); tmp.push_back(striangle[2]);
            edges.push_back(tmp);
            tmp.clear();
            tmp.push_back(striangle[2]); tmp.push_back(striangle[0]);
            edges.push_back(tmp);

            Vector2d center =(vertices[striangle[0]]+
                              vertices[striangle[1]]+
                              vertices[striangle[2]])/3.;

            Matrix2d g = m.compute_metric(center, "a2");

            regions.push_back(std::make_tuple(center, sqrt(g.determinant())));
        }
}

void Mesh::write_to_poly(const double k) {
    std::fstream polyfile;
    //TODO: maybe this should not be hardcoded
    // maybe it does not matter
    polyfile.open("input.poly", std::fstream::out);

    polyfile << vertices.size();
    polyfile << " 2 0 0" << std::endl;
           
    for(uint i=0;i<vertices.size();++i) {
        polyfile << (i+1) << " "
                 << vertices[i].x() << " "
                 << vertices[i].y() << std::endl;
    } 

    polyfile << edges.size() << " 0" << std::endl;
           
    for(uint i=0;i<edges.size();++i) {
        polyfile << (i+1) << " "
                 << (edges[i][0]+1) << " "
                 << (edges[i][1]+1) << std::endl;
    }

    // Number of holes
    polyfile << 0 << std::endl;

    polyfile << regions.size() << std::endl;

    for(uint i=0;i<regions.size();++i) {
        polyfile << (i+1) << " "
                 << std::get<0>(regions[i]).x() << " "
                 << std::get<0>(regions[i]).y() << " "
                 << k*std::get<1>(regions[i]) << std::endl;
    }
}

void Mesh::refine_using_Triangle(const double min_angle,
        const double overlay_size) {


    // TODO: Maybe not assume input.poly, write_to_poly writes here rn
    std::string cmd = "./triangle -an -q"+std::to_string(min_angle)+
        " input.poly";
    exec(cmd.c_str());
    /*
        """
        As protocol, always read from input.poly
        """
        os.system("./triangle -an -q"+str(min_angle)+" input.poly")

        nodes = open("input.1.node", "r").readlines()
        all_vertices = []
        num_nodes = [int(s) for s in nodes[0].split() if s.isdigit()][0]
        for i in range(num_nodes):
            node_data = [float(s) for s in nodes[i+1].split() if s.replace('.','',1).isdigit()]
            all_vertices.append((node_data[1], node_data[2]))

        self.vertices = all_vertices

        eles = open("input.1.ele", "r").readlines()
        all_triangles = []
        num_eles = [int(s) for s in eles[0].split() if s.isdigit()][0]
        for i in range(num_eles):
            ele_data = [int(s) for s in eles[i+1].split() if s.isdigit()]
            all_triangles.append((ele_data[1]-1, ele_data[2]-1, ele_data[3]-1))

        self.triangles = all_triangles

        neighbours = open("input.1.neigh", "r").readlines()
        neighbours_dict = {}
        num_neighs = [int(s) for s in eles[0].split() if s.isdigit()][0]
        for i in range(num_neighs):
            neigh_data = [int(s) for s in neighbours[i+1].split() if s.isdigit()]
            neighbours_dict[(neigh_data[0]-1)] = filter(lambda x:x>=0,map(lambda x:x-1, neigh_data[1:]))

        self.neighbours = neighbours_dict
        */
}

void Mesh::print() {

    std::cout << "Vertices:" << std::endl;
    for (auto vertex : vertices) {
        std::cout << vertex.x() << " " << vertex.y() << std::endl;
    }
    std::cout << std::endl;
    std::cout << "Triangles:" << std::endl;
    for (auto triangle : triangles) {
        std::cout << triangle[0] << " " << triangle[1] << 
            " " << triangle[2] << std::endl;
    }
    std::cout << std::endl;
    //TODO: Add more information
}
