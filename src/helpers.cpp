/*
    This file is part of 2dgrpt.

    2dgrpt is my program to test some fast geodesic algorithms.

    Copyright (c) 2018 by Kumar Ayush
*/

#include "helpers.hpp"

using namespace grpt;

Vector2d random_point(const Vector2d p1, const Vector2d p2, const Vector2d p3) {
    unsigned short int seed[] = {(short unsigned int)(p1.x()*1e5), (short unsigned int)10005, (short unsigned int)99956};
    double c1 = erand48(seed);
    double c2 = erand48(seed);
    double c3 = erand48(seed);
    return (p1*c1 + p2*c2 + p3*c3)/(c1+c2+c3);
}


bool check_for_incidence(const std::vector<Vector2d> vertices, const std::vector<int> triangle, 
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


double return_incidence(const std::vector<Vector2d> vertices, const std::vector<int> triangle, 
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

bool find_closest_edge(const std::unordered_map<std::vector<int>, std::vector<std::vector<int>>, vector_int_hasher> ovmesh, 
        const Vector2d pt, const std::vector<Vector2d> vertices, const double rs, 
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
    for(std::vector<std::vector<int>>::const_iterator it = ovmesh.at(key).begin(); it != ovmesh.at(key).end(); ++it) {
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

bool find_shooting_edge(const std::unordered_map<std::vector<int>, std::vector<std::vector<int>>, vector_int_hasher> ovmesh, 
        const Vector2d pt, const Vector2d direction, const std::vector<Vector2d> vertices, const double rs, 
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
    for(std::vector<std::vector<int>>::const_iterator it = ovmesh.at(key).begin(); it != ovmesh.at(key).end(); ++it) {
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
