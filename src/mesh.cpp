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
#include <ctime>
#include <boost/numeric/odeint.hpp>

#define N_CHECKS 20
#define EPS 1e-10

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
                tmp_vec.push_back(tmp_x-1);
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

    std::vector<int> local_edge;

    double t;
    if (!find_intersection(this_triangle, startpoint, 
                direction, local_edge, t)) {
        std::cerr << "[warn] segment empty" << std::endl;
        std::cerr << "trying incidence check..." << std::endl;
        if (!check_for_incidence(this_triangle, startpoint, 1e-4, local_edge)) {
            std::cerr << "[warn] incidence check failed too" << std::endl;
            std::cerr << "checking neighbours..." << std::endl;
            bool lflag = false;
            for (auto neighbour : neighbours[t_id]) {
                if (find_intersection(triangles[neighbour], startpoint, 
                            direction, local_edge, t)) {
                    std::cerr << "found something on neighbours..." << std::endl;
                    lflag = true;
                    break;
                }
            }
            if (!lflag) {
                for (auto neighbour : neighbours[t_id]) {
                    if (check_for_incidence(triangles[neighbour], startpoint, 1e-4, 
                                local_edge)) {
                        std::cerr << "found something on neighbours..." << std::endl;
                        lflag = true;
                        break;
                    }
                }
                if(!lflag)
                    std::cerr << "[warn] nothing on neighbours either" << std::endl;
            }
        }
    }

    Vector2d point_of_intersection = startpoint + t*direction;
    
    // modify now
    Vector2d new_dir(cos(mod), sin(mod));

    // Save if flag set
    if (save) {
#pragma omp critical
        {
            std::sort(local_edge.begin(), local_edge.end());
            std::vector<Vector2d> tmp = {new_dir-direction, direction};
            edge_data[local_edge].push_back(tmp);
        }
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

    new_dir.y() = new_dir.y()*point_of_intersection.x();
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

void Mesh::throw_geodesic_mark(const int seed, const double tau, 
        double dt) {

    int num_triangle = triangles.size();
    //triangle_id
    int t_id = seed % num_triangle; 
    std::vector<int> this_triangle = triangles[t_id];

    // Start point randomly chosen
    Vector2d startpoint = random_point(
            vertices[this_triangle[0]],
            vertices[this_triangle[1]],
            vertices[this_triangle[2]]
            );

    // Choosing a random direction
    /*
     * Create three directions, randomly choosing a point on three edges
     * Choose randomly among these three directions
     * This will be contradir by the way
     * we now will need to construct the covdir the other way round
     */
    typedef unsigned short int usint;
    usint rand_seed[3] = {usint(seed), usint(seed*seed),
        usint(seed*seed*seed)};

    double r1 = erand48(rand_seed), r2=erand48(rand_seed), 
           r3 = erand48(rand_seed);
    Vector2d tmp = r1*vertices[this_triangle[0]] +
               (1-r1)*vertices[this_triangle[1]];
    Vector2d dir1 = tmp-startpoint;

    tmp = r2*vertices[this_triangle[1]] +
        (1-r2)*vertices[this_triangle[2]];
    Vector2d dir2 = tmp-startpoint;

    tmp = r3*vertices[this_triangle[2]] +
        (1-r3)*vertices[this_triangle[0]];
    Vector2d dir3 = tmp-startpoint;

    Vector2d contradir;
    double rr_eta = erand48(rand_seed);
    if (rr_eta < 0.33) contradir = dir1;
    else if (rr_eta < 0.67) contradir = dir2;
    else contradir = dir3;

    contradir.normalize();

    Vector2d covdir = contradir;
    covdir.y() = covdir.y()*startpoint.x();
    // Renormalizing
    covdir.normalize();

    auto f_minimizer = [=](double x)->double{ 
        return find_trial_error(x, dt, t_id, startpoint, covdir); 
    };

    double one_angle = atan2(contradir.y(), contradir.x());
    std::cerr << "[DBG] Input Angle: " << one_angle*180./M_PI << std::endl;
    double margin = 0.5;

    double result = nelder_mead(one_angle-margin*M_PI, 
            one_angle+margin*M_PI, tau, f_minimizer);

    while (find_trial_error(result, dt, t_id, startpoint, covdir)>tau && 
            margin>=0.2) {

        margin = margin-0.1;
        std::cerr <<  "[DBG] Margin reduced to" << margin << std::endl;
        result = nelder_mead(one_angle-margin*M_PI, one_angle+margin*M_PI, 
                tau, f_minimizer);

    }
    std::cerr << "[DBG] Saving with error " << 
        find_trial_error(result, dt, t_id, startpoint, covdir, true) 
        << std::endl;

}

void Mesh::throw_geodesic_discrete(int seed) {

    // Arrangements for plotting
    namespace plt = matplotlibcpp;

    // discrete and continuous X,Y
    std::vector<double> X_d, Y_d, X_c, Y_c;

    int num_triangle = triangles.size();
    typedef unsigned short int usint;
    usint rand_seed[3] = {usint(seed), usint(12021997), usint(std::time(nullptr))};
    // triangle id
    int t_id = int(erand48(rand_seed)*num_triangle);
    std::vector<int> this_triangle = triangles[t_id];

    // Start point randomly chosen
    Vector2d startpoint = random_point(
            vertices[this_triangle[0]],
            vertices[this_triangle[1]],
            vertices[this_triangle[2]]
            );

    Vector2d direction = Vector2d(erand48(rand_seed)*2-1, 
            erand48(rand_seed)*2-1);
    // Renormalizing
    direction.normalize();

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


    // The integrator
    std::vector<std::vector<double>> result;

    boost::numeric::odeint::runge_kutta4 < std::vector<double> > rk;    
    std::vector<double> state(4); state[0] = startpoint.x();
    state[1] = startpoint.y();
    state[2] = direction.x();
    state[3] = direction.y();

    // Initial time
    int rk_t=0;
    double dt = 0.01;

    result.push_back(state);

    while (_lx < state[0]
            && state[0] < _rx 
            && _ly < state[1]
            && state[1]  < _ry) {

        rk.do_step(metric, state, rk_t, dt);
        rk_t += dt;
        result.push_back(state);
    }

    for(int i = 0;i<result.size();++i){
        X_c.push_back(result[i][0]);
        Y_c.push_back(result[i][1]);
    }

    direction = (Vector2d(result[1][0], result[1][1])
            - startpoint).normalized();


    int count_checker = 0;
    while(true) {

        std::vector<int> local_edge;

        double t;
        if (!find_intersection(this_triangle, startpoint, 
                    direction, local_edge, t)) {
            std::cerr << "[warn] segment empty" << std::endl;
            std::cerr << "trying incidence check..." << std::endl;
            if (!check_for_incidence(this_triangle, startpoint, 1e-4, local_edge)) {
                std::cerr << "[warn] incidence check failed too" << std::endl;
                std::cerr << "checking neighbours..." << std::endl;
                bool lflag = false;
                for (auto neighbour : neighbours[t_id]) {
                    if (find_intersection(triangles[neighbour], startpoint, 
                                direction, local_edge, t)) {
                        std::cerr << "found something on neighbours..." << std::endl;
                        lflag = true;
                        break;
                    }
                }
                if (!lflag) {
                    for (auto neighbour : neighbours[t_id]) {
                        if (check_for_incidence(triangles[neighbour], startpoint, 
                                    1e-4, local_edge)) {
                            std::cerr << "found something on neighbours..." 
                                      << std::endl;
                            lflag = true;
                            break;
                        }
                    }
                    if(!lflag)
                        std::cerr << "[warn] nothing on neighbours either" 
                                  << std::endl;
                }
            }
        }
        if(!local_edge.size()) {
            std::cout << "Bad boy exit...@" << count_checker << std::endl;
            break;
        }

        Vector2d tmp_plot = startpoint;
        for(double i = 0;i<t;i+=0.01*t) {
            tmp_plot += 0.01*t*direction;   
            X_d.push_back(tmp_plot.x());
            Y_d.push_back(tmp_plot.y());
        }
        Vector2d point_of_intersection = startpoint + t*direction;
        Vector2d ndir;

        std::sort(local_edge.begin(), local_edge.end());

        if (edge_data[local_edge].size()<2) {

            std::cout << "[WARN] Sad you didn't sample enough boi..."
                << std::endl;

            // TODO: this whole section is a scam
            Vector2d cov_direction = Vector2d(direction.x(), 
                    direction.y()*startpoint.x()).normalized();

            boost::numeric::odeint::runge_kutta4 < std::vector<double> > rk_ad;    
            std::vector<double> state_ad(4);

            Vector2d trialpt = startpoint - 0.025*direction;

            state_ad[0] = trialpt.x();
            state_ad[1] = trialpt.y();
            state_ad[2] = cov_direction.x();
            state_ad[3] = cov_direction.y();

            int rk_ad_t=0;
            rk_ad.do_step(metric, state_ad, rk_ad_t, 0.025);

            ndir = Vector2d(state_ad[0], state_ad[1])-trialpt;
        }

        else {
            // Lagrange interpolation
            std::vector<std::vector<double>> entry = 
                edge_slope_data[local_edge];

            int _n = entry.size();

            double in_ang = anglemod(atan2(direction.y(), direction.x()));
            ndir = Vector2d(-2,-2);

            if (entry[_n-1][0] <= in_ang && in_ang < entry[0][0]) {
                double n_ang = anglemod(entry[0][1]*in_ang + entry[0][2]);
                ndir = Vector2d(cos(n_ang), sin(n_ang));
            }
            else {

                for(int i = 1;i<_n-1;++i) {
                    if (entry[i-1][0] <= in_ang && in_ang < entry[i][0]) {
                        double n_ang = anglemod(entry[i][1]*in_ang + 
                                entry[i][2]);
                        ndir = Vector2d(cos(n_ang), sin(n_ang));
                        break;
                    }
                }
            }
            if (ndir.x()==-2) {
                double n_ang = anglemod(entry[_n-1][1]*in_ang + 
                        entry[_n-1][2]);
                ndir = Vector2d(cos(n_ang), sin(n_ang));
            }

        }


        auto check_edge = [=](int x){
            std::vector<int> el = triangles[x];
            return ((std::find(el.begin(), el.end(), local_edge[0])!=el.end())
                    &&
                    (std::find(el.begin(), el.end(), local_edge[1])!=el.end())
                   );
        };
        std::vector<int>::iterator it = std::find_if(
                neighbours[t_id].begin(), 
                neighbours[t_id].end(),
                check_edge);

        if (it!= neighbours[t_id].end())
            t_id = *it;
        else {
            std::cout << "Good boy exit...@" << count_checker << std::endl;
            break;
        }

        startpoint = point_of_intersection;
        this_triangle = triangles[t_id];

        direction = ndir.normalized();

        count_checker++;
        if (count_checker>=500) {
            std::cerr << "Exit as loop exceeded threshold..." << std::endl;
            break;
        }
    }

    plt::plot(X_c,Y_c, "b-");
    plt::plot(X_d,Y_d, "r-");
    plt::show();
}

bool Mesh::find_intersection(std::vector<int> this_triangle, Vector2d startpoint, 
        Vector2d direction, std::vector<int> &local_edge, double &t) {

    Vector2d a = vertices[this_triangle[0]];
    Vector2d b = vertices[this_triangle[1]];
    Vector2d c = vertices[this_triangle[2]];

    t = cross(a-startpoint, b-a)/cross(direction,b-a);
    double s = cross(a-startpoint, direction)/cross(direction, b-a);
    if (EPS <= t && t <= 1 && EPS <= s && s <= 1) {
        local_edge.push_back(this_triangle[0]);
        local_edge.push_back(this_triangle[1]);
    }
    else {
        t = cross(b-startpoint, c-b)/cross(direction,c-b);
        s = cross(b-startpoint, direction)/cross(direction, c-b);
        if (EPS <= t && t <= 1 && EPS <= s && s <= 1) {
            local_edge.push_back(this_triangle[1]);
            local_edge.push_back(this_triangle[2]);
        }
        else {
            t = cross(c-startpoint, a-c)/cross(direction,a-c);
            s = cross(c-startpoint, direction)/cross(direction, a-c);
            if (EPS <= t && t <= 1 && EPS <= s && s <= 1) {
                local_edge.push_back(this_triangle[2]);
                local_edge.push_back(this_triangle[0]);
            }
            else {
                return false;
            }
        }
    }
    return true;
}

bool Mesh::check_for_incidence(const std::vector<int> triangle, 
        const Vector2d point_of_intersection, const double threshold,
        std::vector<int> &local_edge) {

    Vector2d a = vertices[triangle[0]];
    Vector2d b = vertices[triangle[1]];
    Vector2d c = vertices[triangle[2]];

    Vector2d pos = point_of_intersection;

    if( fabs((pos-a).norm() + (pos-b).norm() - (a-b).norm()) < threshold ) {
        local_edge.push_back(triangle[0]);
        local_edge.push_back(triangle[1]);
    }
    else if( fabs((pos-b).norm() + (pos-c).norm() - (b-c).norm()) < threshold ) {
        local_edge.push_back(triangle[1]);
        local_edge.push_back(triangle[2]);
    }
    else if( fabs((pos-c).norm() + (pos-a).norm() - (c-a).norm()) < threshold ) {
        local_edge.push_back(triangle[2]);
        local_edge.push_back(triangle[0]);
    }
    else
        return false;

    return true;
}
