/*
    This file is part of 2dgrpt.

    2dgrpt is my program to test some fast geodesic algorithms.

    Copyright (c) 2018 by Kumar Ayush
*/

#include "mesh.hpp"
#include "helpers.hpp"
#include "matplotlibcpp.hpp"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <boost/numeric/odeint.hpp>

#define N_CHECKS 5

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
    polyfile.close();
}

void Mesh::refine_using_Triangle(const double min_angle,
        const double overlay_size) {


    // TODO: Maybe not assume input.poly, write_to_poly writes here rn
    std::string cmd = "./triangle -an -q"+std::to_string(min_angle)+
        " input.poly";
    exec(cmd.c_str());

    // Read nodes
    std::fstream nodefile;
    nodefile.open("input.1.node", std::fstream::in);
    int num_nodes;
    int tmp;
    std::vector<Vector2d> all_vertices;
    nodefile >> num_nodes >> tmp >> tmp >> tmp;
    for(int i = 0;i<num_nodes;++i) {
        nodefile >> tmp;
        double x,y;
        nodefile >> x >> y;
        all_vertices.push_back(Vector2d(x,y));
        nodefile >> tmp;
    }
    vertices = all_vertices;
    nodefile.close();

    // Read triangles
    std::fstream elefile;
    elefile.open("input.1.ele", std::fstream::in);
    int num_eles;
    std::vector<std::vector<int>> all_triangles;
    elefile >> num_eles >> tmp >> tmp;
    for(int i = 0;i<num_eles;++i) {
        elefile >> tmp;
        int a,b,c;
        elefile >> a >> b >> c;
        // Triangle orders from 1
        std::vector<int> tmp_vec = {a-1,b-1,c-1};
        all_triangles.push_back(tmp_vec);
    }
    triangles = all_triangles;
    elefile.close();


    // Read neighbours
    std::fstream neifile;
    neifile.open("input.1.neigh", std::fstream::in);
    int num_neighs;
    std::vector<std::vector<int>> all_neighs;
    neifile >> num_neighs >> tmp;
    for(int i = 0;i<num_neighs;++i) {
        neifile >> tmp;
        
        std::vector<int> tmp_vec;
        // three neighbours at max
        for(int j = 0;j<3;j++) {
            int tmp_x;
            neifile >> tmp_x;
            if (tmp_x != -1)
                tmp_vec.push_back(tmp_x);
        }

        all_neighs.push_back(tmp_vec);
    }
    neighbours = all_neighs;
    neifile.close();

    // Fill the overlay mesh
    for (auto tl : triangles) {

        double x1,x2,y1,y2;
        int a,b;
        x1 = vertices[tl[0]].x()/overlay_size;
        y1 = vertices[tl[0]].y()/overlay_size;

        x2 = vertices[tl[1]].x()/overlay_size;
        y2 = vertices[tl[1]].y()/overlay_size;

        a = tl[0] < tl[1] ? tl[0] : tl[1];
        b = tl[0] < tl[1] ? tl[1] : tl[0];
        bresenham_and_mesh(overlay_mesh, x1, y1, x2, y2, a, b);

        x1 = vertices[tl[1]].x()/overlay_size;
        y1 = vertices[tl[1]].y()/overlay_size;

        x2 = vertices[tl[2]].x()/overlay_size;
        y2 = vertices[tl[2]].y()/overlay_size;

        a = tl[1] < tl[2] ? tl[1] : tl[2];
        b = tl[1] < tl[2] ? tl[2] : tl[1];
        bresenham_and_mesh(overlay_mesh, x1, y1, x2, y2, a, b);

        x1 = vertices[tl[2]].x()/overlay_size;
        y1 = vertices[tl[2]].y()/overlay_size;

        x2 = vertices[tl[0]].x()/overlay_size;
        y2 = vertices[tl[0]].y()/overlay_size;

        a = tl[2] < tl[0] ? tl[2] : tl[0];
        b = tl[2] < tl[0] ? tl[0] : tl[2];
        bresenham_and_mesh(overlay_mesh, x1, y1, x2, y2, a, b);

    }

}


void Mesh::churn_edge_data() {

    for (auto one_edge_data : edge_data) {
        std::vector<std::vector<Vector2d>> all_samples = one_edge_data.second;
        std::vector<std::vector<double>> anglesamples;
        for (auto sample : all_samples) {
            double input_angle = anglemod(atan2(sample[1].y(), sample[1].x()));
            double output_angle = anglemod(atan2(sample[0].y()+sample[1].y(),
                        sample[0].x()+sample[1].x()));

            std::vector<double> tmp = {input_angle, output_angle};
            anglesamples.push_back(tmp);
        }
        auto sorting_angles = [](std::vector<double> i,
                std::vector<double> j){ return (i[0]<j[0]); };

        std::sort(anglesamples.begin(), anglesamples.end(), sorting_angles);
        std::vector<std::vector<double>> slopesamples; 

        int N = anglesamples.size();

        for(int i = 0;i<N-1;i++) {
            double m = (anglesamples[i+1][1]-anglesamples[i][1])/
                (anglesamples[i+1][0]-anglesamples[i][0]);
            double b = anglesamples[i][1] - m*anglesamples[i][0];
            // already sorted by input_angle
            std::vector<double> tmp = {anglesamples[i+1][0], m, b};
            slopesamples.push_back(tmp);
        }

        if (N > 1) {
            double m = (anglesamples[N-1][1]-anglesamples[0][1])/
                (anglesamples[N-1][0]-anglesamples[0][0]);
            double b = anglesamples[0][1] - m*anglesamples[0][0];
            std::vector<double> tmp = {anglesamples[0][0], m, b};
            slopesamples.push_back(tmp);
        }
        edge_slope_data[one_edge_data.first] = slopesamples;
    }
}

void Mesh::draw(){

    std::vector<double> X,Y;
    for (auto vertex : vertices) {
        X.push_back(vertex.x());
        Y.push_back(vertex.y());
    }
    namespace plt = matplotlibcpp;
    plt::plot(X,Y);
    plt::show();
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

double Mesh::find_trial_error(const double mod, const double dt, 
        const int t_id, const Vector2d startpoint, 
        const Vector2d covdir, bool save, bool dbg) {
    
    auto sorting_corners_1 = [](Vector2d i,
            Vector2d j){ return (i.x()<j.x()); };
    std::sort(corners.begin(), corners.end(), 
            sorting_corners_1);

    double _lx, _rx, _ly, _ry;
    _lx = corners[0].x();
    _rx = corners[corners.size()-1].x();

    auto sorting_corners_2 = [](Vector2d i,
            Vector2d j){ return (i.y()<j.y()); };
    std::sort(corners.begin(), corners.end(), 
            sorting_corners_2);
    _ly = corners[0].y();
    _ry = corners[corners.size()-1].y();

    std::vector<int> this_triangle = triangles[t_id];

    std::vector<std::vector<double>> result;

    boost::numeric::odeint::runge_kutta4 < std::vector<double> > rk;    
    std::vector<double> state(4);
    state[0] = startpoint.x();
    state[1] = startpoint.y();
    state[2] = covdir.x();
    state[3] = covdir.y();
    // Initial time
    int rk_t=0;

    result.push_back(state);

    // Take a step
    rk.do_step(metric, state, rk_t, dt);
    rk_t += dt;
    result.push_back(state);

    Vector2d direction = (Vector2d(result[1][0], result[1][1]) - startpoint).normalized();

    Vector2d a = vertices[this_triangle[0]];
    Vector2d b = vertices[this_triangle[1]];
    Vector2d c = vertices[this_triangle[2]];

    std::vector<int> local_edge;

    double t = cross(a-startpoint, b-a)/cross(direction,b-a);
    double s = cross(a-startpoint, direction)/cross(direction, b-a);
    if (t <= 1 && s >= 0) {
        local_edge.push_back(this_triangle[0]);
        local_edge.push_back(this_triangle[1]);
    }
    else {
        t = cross(b-startpoint, c-b)/cross(direction,c-b);
        s = cross(b-startpoint, direction)/cross(direction, c-b);
        if (t <= 1 && s >= 0) {
            local_edge.push_back(this_triangle[1]);
            local_edge.push_back(this_triangle[2]);
        }
        else {
            t = cross(c-startpoint, a-c)/cross(direction,a-c);
            s = cross(c-startpoint, direction)/cross(direction, a-c);
            if (t <= 1 && s >= 0) {
                local_edge.push_back(this_triangle[2]);
                local_edge.push_back(this_triangle[0]);
            }
            else {
                std::cout << "[WARN] Segment empty" << std::endl;
            }
        }
    }

    Vector2d point_of_intersection = startpoint + t*direction*(1.+1e-5);
    
    // modify now
    Vector2d new_dir(cos(mod), sin(mod));

    // Save if flag set
    if (save) {
        std::sort(local_edge.begin(), local_edge.end());
        std::vector<Vector2d> tmp = {new_dir-direction, direction};
        edge_data[local_edge].push_back(tmp);
    }
    
    if (dbg && save) {
        std::cout << "[DBG]" << " "
                  << mod*180./M_PI << " "
                  << new_dir-direction << " "
                  << direction << std::endl;
    }

    if (dbg) {
        std::cout << "[DBG]" << " "
                  << mod*180./M_PI;
    }

    std::vector<std::vector<double>> dev_result;

    boost::numeric::odeint::runge_kutta4 < std::vector<double> > dev_rk;    
    std::vector<double> dev_state(4);

    new_dir.x() = new_dir.x()*point_of_intersection.x();
    new_dir.normalize();

    dev_state[0] = point_of_intersection.x();
    dev_state[1] = point_of_intersection.y();
    dev_state[2] = new_dir.x();
    dev_state[3] = new_dir.y();
    // Initial time
    int dev_rk_t=0;

    dev_result.push_back(dev_state);

    int Ni = 0;
    while (Ni < N_CHECKS && _lx < dev_state[0]
            && dev_state[0] < _rx 
            && _ly < dev_state[1]
            && dev_state[1]  < _ry) {

        dev_rk.do_step(metric, dev_state, dev_rk_t, dt);
        dev_rk_t += dt;
        dev_result.push_back(dev_state);
        Ni++;
    }

    Ni = 0;
    while (Ni < N_CHECKS && _lx < state[0]
            && state[0] < _rx 
            && _ly < state[1]
            && state[1]  < _ry) {

        rk.do_step(metric, state, rk_t, dt);
        rk_t += dt;
        result.push_back(state);
        Ni++;
    }


    // Compare target angle versus calculated angle
    Vector2d targetf(result[result.size()-1][0], result[result.size()-1][1]);
    Vector2d targeti(result[0][0], result[0][1]);

    Vector2d target = targetf-targeti;
    double target_angle = atan2(target.y(), target.x());

    Vector2d calcuf(dev_result[dev_result.size()-1][0], 
            dev_result[dev_result.size()-1][1]);
    Vector2d calcui(dev_result[0][0], dev_result[0][1]);

    Vector2d calcu = calcuf-calcui;
    double calcu_angle = atan2(calcu.y(), calcu.x());

    return fabs(target_angle-calcu_angle);
}
