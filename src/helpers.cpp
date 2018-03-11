/*
    This file is part of 2dgrpt.

    2dgrpt is my program to test some fast geodesic algorithms.

    Copyright (c) 2018 by Kumar Ayush
*/

#include "helpers.hpp"

using namespace grpt;

double calculate_area(const Vector2d p1, const Vector2d p2, 
        const Vector2d p3) {

    Matrix3d area_mat;
    area_mat << 1, p1.x(), p1.y(),
                1, p2.x(), p2.y(),
                1, p3.x(), p3.y();

    return 0.5*fabs(area_mat.determinant());
}

Vector2d random_point(const Vector2d p1, const Vector2d p2, const Vector2d p3) {
    unsigned short int seed[] = {(short unsigned int)(p1.x()*1e5), 
        (short unsigned int)10005, (short unsigned int)99956};
    double c1 = erand48(seed);
    double c2 = erand48(seed);
    double c3 = erand48(seed);
    return (p1*c1 + p2*c2 + p3*c3)/(c1+c2+c3);
}


bool check_for_incidence(const std::vector<Vector2d> vertices, 
        const std::vector<int> triangle, 
        const Vector2d point_of_intersection, const double threshold) {

    Vector2d a = vertices[triangle[0]];
    Vector2d b = vertices[triangle[1]];
    Vector2d c = vertices[triangle[2]];

    Vector2d pos = point_of_intersection;

    if( fabs((pos-a).norm() + (pos-b).norm() - (a-b).norm()) < threshold ) 
        return true;
    else if( fabs((pos-b).norm() + (pos-c).norm() - (b-c).norm()) < threshold ) 
        return true;
    else if( fabs((pos-c).norm() + (pos-a).norm() - (c-a).norm()) < threshold ) 
        return true;

    return false;
}


double return_incidence(const std::vector<Vector2d> vertices, 
        const std::vector<int> triangle, 
        const Vector2d point_of_intersection) {

    Vector2d a = vertices[triangle[0]];
    Vector2d b = vertices[triangle[1]];
    Vector2d c = vertices[triangle[2]];

    Vector2d pos = point_of_intersection;

    double one = fabs((pos-a).norm() + (pos-b).norm() - (a-b).norm());
    double two = fabs((pos-b).norm() + (pos-c).norm() - (b-c).norm());
    double tre = fabs((pos-c).norm() + (pos-a).norm() - (c-a).norm());

    return std::min(std::min(one, two), tre);
}

bool find_closest_edge(const std::unordered_map<std::vector<int>, 
        std::vector<std::vector<int>>, vector_int_hasher> ovmesh, 
        const Vector2d pt, const std::vector<Vector2d> vertices, 
        const double rs, 
        std::vector<int> &min_edge, double &distance) {

    int x = int(pt[0]/rs);
    int y = int(pt[1]/rs);
    std::vector<int> key = {x,y};
    if (!ovmesh.count(key)) {
        return false;
    }

    // TODO: Change this hard coding to either +INF or some better estimate
    distance = 10000;
    min_edge = std::vector<int>(0,0);

    for(
    std::vector<std::vector<int>>::const_iterator it = ovmesh.at(key).begin(); 
    it != ovmesh.at(key).end(); ++it) {

        Vector2d a = vertices[(*it)[0]];
        Vector2d b = vertices[(*it)[1]];

        double dist = fabs((pt-a).norm() + (pt-b).norm() - (a-b).norm()); 
        if (dist < distance) {
            distance = dist;
            min_edge = *it;
        }
    }

    // TODO: check if this threshold should be something else
    if (distance < rs)
        return true;
    else
        return false;

}

bool find_shooting_edge(const std::unordered_map<std::vector<int>, 
        std::vector<std::vector<int>>, vector_int_hasher> ovmesh, 
        const Vector2d pt, const Vector2d direction, 
        const std::vector<Vector2d> vertices, const double rs, 
        std::vector<int> &min_edge, double &mint) {

    int x = int(pt[0]/rs);
    int y = int(pt[1]/rs);
    std::vector<int> key = {x,y};
    if (!ovmesh.count(key)) {
        return false;
    }

    // TODO: Change this hard coding to either +INF or some better estimate
    mint = 10000;
    min_edge = std::vector<int>(0,0);

    for(
    std::vector<std::vector<int>>::const_iterator it = ovmesh.at(key).begin();
    it != ovmesh.at(key).end(); ++it) {

        Vector2d a = vertices[(*it)[0]];
        Vector2d b = vertices[(*it)[1]];

        double t = cross((a-pt), (b-a))/cross(direction, (b-a));
        double s = cross((a-pt), direction)/cross(direction, (b-a));

        if (fabs(t) < mint && 0 <= s && s <= 1){
            mint = t;
            min_edge = *it;
        } 
    }

    // TODO: check if this threshold should be something else
    if (mint < 10000)
        return true;
    else
        return false;

}

bool bresenham_and_mesh(std::unordered_map<std::vector<int>, 
        std::vector<std::vector<int>>, vector_int_hasher> &ovm, 
        const double x1, const double y1, 
        const double x2, const double y2, 
        const int a, const int b) {

    //TODO: make it conservative
    int deltax = x2-x1;
    int deltay = y2-y1;

    if (deltax == 0) {
        int x = int(x1);
        int y;
        for(y=int(y1);y<=int(y2);++y) {


            std::vector<int> key = {x,y};
            if (!ovm.count(key))
                ovm[key] = std::vector<std::vector<int>>();
            else
                ovm[key].push_back(std::vector<int>(a,b));

        }
    }

    else {
        double deltaerr = fabs(double(deltay)/deltax);
        double error = 0.0;
        int y = int(y1);
        int x;
        for(x=int(x1);x<=int(x2);++x) {
            std::vector<int> key = {x,y};
            if (!ovm.count(key))
                ovm[key] = std::vector<std::vector<int>>();
            else
                ovm[key].push_back(std::vector<int>(a,b));

            error = error + deltaerr;
            while (error>=0.5){
                y += std::copysign(1, deltay);
                error -= 1.0;
            }

        }
    }

    return true;
}

double nelder_mead(double x1, double x2, const double tau, 
        const double *f(const double)) {

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
