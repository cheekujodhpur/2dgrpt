#!/usr/bin/env python

# this is the file where we draw geodescis
from helpers import *
from matplotlib.lines import Line2D
from mesh import Mesh
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import ode
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np
import os

NOS_EDGE = 10 #number of samples on each edge

# Functions for integration
# def f(t, y):
#     dy_0 = y[2]
#     dy_1 = y[3]
#     # dy_2 = -0.5*y[0]*(y[3]**2)
#     dy_2 = 0
#     dy_3 = 0
#     return np.array([dy_0, dy_1, dy_2, dy_3])

# def jac(t, y):
#     r1 = np.array([0, 0, 1, 0])
#     r2 = np.array([0, 0, 0, 1])
#     r3 = np.array([0, 0, 0, 0])
#     r4 = np.array([0, 0, 0, 0])
#     return np.array([r1,r2,r3,r4])

def f(t, y):
    dy_0 = y[2]
    dy_1 = y[3]/y[0]
    # dy_2 = -0.5*y[0]*(y[3]**2)
    dy_2 = 0.5*(y[3]**2)/(y[0]**2)
    dy_3 = 0
    return np.array([dy_0, dy_1, dy_2, dy_3])

def jac(t, y):
    r1 = np.array([0, 0, 1, 0])
    r2 = np.array([-y[3]/(y[0]*y[0]), 0, 0, 1/y[0]])
    r3 = np.array([-y[3]*y[3]/(y[0]**3), 0, 0, y[3]/(y[0]**2)])
    r4 = np.array([0,0,0,0])
    return np.array([r1,r2,r3,r4])


myMesh = Mesh([(0,0),(0,5),(5,0),(5,5)], (0,0))
myMesh.submesh(3)
myMesh.write_to_poly(k=0.05)
refined_size = 0.05
myMesh.refine_using_Triangle(18, refined_size)

def throw_geodesic_integrating(mesh, dt=0.01):
    num_triangle = len(mesh.triangles)
    t_id = int(np.random.random()*num_triangle) #triangle_id
    this_triangle = mesh.triangles[t_id]

    # Start point randomly chosen
    startpoint = random_point(
            mesh.vertices[this_triangle[0]],
            mesh.vertices[this_triangle[1]],
            mesh.vertices[this_triangle[2]]
            )

    # Choosing a random direction
    direction = np.array((np.random.random()*2-1, np.random.random()*2-1))
    # Renormalizing
    direction = direction / np.linalg.norm(direction)

    _lx = sorted(mesh.corners, key=lambda x:x[0])[0][0]
    _rx = sorted(mesh.corners, key=lambda x:x[0])[-1][0]
    _ly = sorted(mesh.corners, key=lambda x:x[1])[0][1]
    _ry = sorted(mesh.corners, key=lambda x:x[1])[-1][1]

    # Geodesic integrator
    g_int = ode(f, jac).set_integrator("dopri5")
    y0 = np.array([startpoint[0], startpoint[1], direction[0], direction[1]])
    g_int.set_initial_value(y0,0)

    result = []
    while g_int.successful() and (_lx < g_int.y[0] < _rx) and (_ly < g_int.y[1] < _ry):
        result.append(g_int.integrate(g_int.t+dt))

    return np.array(result)


def throw_geodesic_mark(mesh, ax, dt=0.01):

    _lx = sorted(mesh.corners, key=lambda x:x[0])[0][0]
    _rx = sorted(mesh.corners, key=lambda x:x[0])[-1][0]
    _ly = sorted(mesh.corners, key=lambda x:x[1])[0][1]
    _ry = sorted(mesh.corners, key=lambda x:x[1])[-1][1]

    g_int = ode(f, jac).set_integrator("dopri5")

    dx = 2*dt

    count = 0
    # iterate over all triangles
    for triangle in mesh.triangles:
        count = count + 1
        print "Triangle no.", count, "..."
        for xx in range(3):
            # pick an edge and go on
            edge = tuple(sorted([triangle[xx], triangle[(xx+1)%3]]))
            if len(mesh.edge_data[edge])<NOS_EDGE:
                # pick a point a bit off the edge, how do we decide how off the edge
                edge_v = np.array(mesh.vertices[edge[1]])-np.array(mesh.vertices[edge[0]])
                edge_mid = (np.array(mesh.vertices[edge[1]])+np.array(mesh.vertices[edge[0]]))/2.
                perp_v = np.array([edge_v[1], -edge_v[0]])
                perp_v = perp_v / np.linalg.norm(perp_v)
                # assume sigma is one, that is ratio of dt to dx

                for sgnf in [1,-1]:

                    perp_v = sgnf*perp_v

                    startpoint = edge_mid - dx*perp_v
                    perp_ang = anglemod(np.arctan2(perp_v[1], perp_v[0]))
                    
                    start_ang = perp_ang - ((NOS_EDGE-1)/float(NOS_EDGE))*np.pi
                    end_ang = perp_ang + ((NOS_EDGE-1)/float(NOS_EDGE))*np.pi
                    while start_ang < end_ang:
                        direction = np.array([np.cos(start_ang), np.sin(start_ang)])
                        y0 = np.array([startpoint[0], startpoint[1], direction[0], direction[1]])
                        g_int.set_initial_value(y0,0)
                        res = g_int.integrate(g_int.t+dt)
                        old_dir = res[:2] - startpoint
                        old_dir = old_dir / np.linalg.norm(old_dir)

                        oldfinpos = res
                        finpos = np.copy(res)
                        for dumiter in range(4):
                            oldfinpos = np.copy(finpos)
                            finpos = g_int.integrate(g_int.t+dt)

                        new_dir = finpos[:2]-oldfinpos[:2]
                        new_dir = new_dir / np.linalg.norm(new_dir)
                        mesh.edge_data[edge].append([new_dir-old_dir, old_dir])
                        start_ang += (1/float(NOS_EDGE))*np.pi


def throw_geodesic_discrete(mesh, ax):
    num_triangle = len(mesh.triangles)
    t_id = int(np.random.random()*num_triangle) #triangle_id
    this_triangle = mesh.triangles[t_id]

    # Start point randomly chosen
    startpoint = random_point(
            mesh.vertices[this_triangle[0]],
            mesh.vertices[this_triangle[1]],
            mesh.vertices[this_triangle[2]]
            )

    # Choosing a random direction
    direction = np.array((np.random.random()*2-1, np.random.random()*2-1))
    # Renormalizing
    direction = direction / np.linalg.norm(direction)

    _lx = sorted(mesh.corners, key=lambda x:x[0])[0][0]
    _rx = sorted(mesh.corners, key=lambda x:x[0])[-1][0]
    _ly = sorted(mesh.corners, key=lambda x:x[1])[0][1]
    _ry = sorted(mesh.corners, key=lambda x:x[1])[-1][1]

    # Geodesic integrator
    g_int = ode(f, jac).set_integrator("dopri5")
    y0 = np.array([startpoint[0], startpoint[1], direction[0], direction[1]])
    g_int.set_initial_value(y0,0)

    result = []
    while g_int.successful() and (_lx < g_int.y[0] < _rx) and (_ly < g_int.y[1] < _ry):
        result.append(g_int.integrate(g_int.t+0.01))

    result = np.array(result)
    ax.plot(result[:,0], result[:,1], color="black", linewidth=1)

    direction = result[0,:2]-startpoint
    direction = direction / np.linalg.norm(direction)

    count_checker = 0

    while True:

        # now figure out the triangle intersection
        a = np.array(mesh.vertices[this_triangle[0]])
        b = np.array(mesh.vertices[this_triangle[1]])
        c = np.array(mesh.vertices[this_triangle[2]])
        intersections = {}
        # ab
        t = np.cross((a-startpoint), (b-a))/np.cross(direction, (b-a))
        s = np.cross((a-startpoint), direction)/np.cross(direction, (b-a))
        intersections['ab'] = (t,s)
        # bc
        t = np.cross((b-startpoint), (c-b))/np.cross(direction, (c-b))
        s = np.cross((b-startpoint), direction)/np.cross(direction, (c-b))
        intersections['bc'] = (t,s)
        # ca
        t = np.cross((c-startpoint), (a-c))/np.cross(direction, (a-c))
        s = np.cross((c-startpoint), direction)/np.cross(direction, (a-c))
        intersections['ca'] = (t,s)

        # segment = filter(lambda x:intersections[x][1] <= 1 and intersections[x][1] >=0 and np.abs(intersections[x][0])>=1e-10, intersections.keys())
        segment = filter(lambda x:intersections[x][1] <= 1 and intersections[x][1] >=0, intersections.keys())
        try:
            segment = max(segment, key=lambda x: intersections[x][0])
        except:
            print "Segment empty", intersections
            break

        local_edge = (0,0)
        if segment == 'ab':
            local_edge = (this_triangle[0], this_triangle[1])
        elif segment == 'bc':
            local_edge = (this_triangle[1], this_triangle[2])
        elif segment == 'ca':
            local_edge = (this_triangle[2], this_triangle[0])

        point_of_intersection = startpoint + intersections[segment][0]*direction*(1+1e-5)

        # use sorted tuple for uniqueness
        if len(mesh.edge_data[tuple(sorted(local_edge))])<2:
            print "Sad you didn't sample enough boi..."
            # TODO: this whole section is a scam
            cov_direction = np.array([direction[0], direction[1]*startpoint[0]])
            ndir = cov_direction + np.array([0.5 * cov_direction[1] * cov_direction[1]/(startpoint[0]*startpoint[0]), 0])
            ndir[1] = ndir[1]/startpoint[0]
            ax.add_line(Line2D([startpoint[0],point_of_intersection[0]] \
                        ,[startpoint[1], point_of_intersection[1]],color="red", lw=1))
            # ndir = np.copy(direction)

        else:
            # Lagrange interpolation
            entry = mesh.edge_slope_data[tuple(sorted(local_edge))]
            _n = len(entry)

            in_ang = anglemod(np.arctan2(direction[1], direction[0]))
            ndir = []
            if entry[-1][0] <= in_ang < entry[0][0]:
                n_ang = anglemod(entry[0][1]*in_ang + entry[0][2])
                ndir = np.array([np.cos(n_ang), np.sin(n_ang)])
            else:
                for i in range(1, _n-1):
                    if entry[i-1][0] <= in_ang < entry[i][0]:
                        n_ang = anglemod(entry[i][1]*in_ang + entry[i][2])
                        ndir = np.array([np.cos(n_ang), np.sin(n_ang)])
                        break
            if len(ndir)==0:
                n_ang = anglemod(entry[_n-1][1]*in_ang + entry[_n-1][2])
                ndir = np.array([np.cos(n_ang), np.sin(n_ang)])

            # m = (entry[1][0]-entry[0][0])/(entry[1][1]*entry[1][1]-entry[0][1]*entry[0][1])
            # b = entry[0][0] - m*entry[0][1]*entry[0][1]
            # ndir = direction + m*direction*direction + b

            ax.add_line(Line2D([startpoint[0],point_of_intersection[0]] \
                        ,[startpoint[1], point_of_intersection[1]],color="green", lw=1))
            # print "using data", ndir


        pos = point_of_intersection
        #TODO: Maybe parameterize this 1e-4
        t_id = filter(lambda x:check_for_incidence(mesh.vertices, mesh.triangles[x], point_of_intersection, 1e-4),
                mesh.neighbours[t_id])

        if len(t_id)==1:
            t_id = t_id[0]
        elif len(t_id)>1:
            print "Error. Too many neighbours..."
            t_id = sorted(t_id, key=lambda x:return_incidence(mesh.vertices, mesh.triangles[x], point_of_intersection))[0]
            # break
        else:
            print "Good Boy Exit"
            break

        startpoint = point_of_intersection
        this_triangle = mesh.triangles[t_id]

        direction = np.copy(ndir)
        direction = direction / np.linalg.norm(direction)

        count_checker = count_checker + 1
        if count_checker>=5000:
            print "Exit as loop exceeded threshold..."
            break

    return point_of_intersection, direction


fig = plt.figure()
ax = fig.add_subplot(111)
# ax.set_xlabel("x")
# ax.set_ylabel("y")
# _lx = sorted(myMesh.corners, key=lambda x:x[0])[0][0]
# _rx = sorted(myMesh.corners, key=lambda x:x[0])[-1][0]
# ax.set_xlim([_lx, _rx])
# _ly = sorted(myMesh.corners, key=lambda x:x[1])[0][1]
# _ry = sorted(myMesh.corners, key=lambda x:x[1])[-1][1]
# ax.set_ylim([_ly, _ry])

# res = throw_geodesic_integrating(myMesh)
# ax.plot(res[:,0], res[:,1], color="black", linewidth=1)

# Plot the mesh
# myMesh.draw(ax)

# Discrete geodesic
throw_geodesic_mark(myMesh, ax, dt=5e-3)

myMesh.churn_edge_data()
# [ax.add_line(Line2D([i*refined_size,i*refined_size],[0,5],color="red",lw=.2)) for i in range(1,int(5/refined_size))]
# [ax.add_line(Line2D([0,5],[i*refined_size,i*refined_size],color="red",lw=.2)) for i in range(1,int(5/refined_size))]

# Drawing the edge data
def draw_edge_data(myMesh):
    for edge in myMesh.edge_data.keys():
        for el in myMesh.edge_data[edge]:
            left = myMesh.vertices[edge[0]]
            right = myMesh.vertices[edge[1]]
            x = (left[0]+right[0])/2.
            y = (left[1]+right[1])/2.
            try:
                ax.add_line(Line2D([x,x+0.01*el[1][0]] \
                            ,[y, y+0.01*el[1][1]],color="green", lw=1))

                ax.add_line(Line2D([x+0.01*el[1][0],x+2*0.01*el[1][0]+0.01*el[0][0]] \
                            ,[y+0.01*el[1][1],y+2*0.01*el[1][1]+0.01*el[0][1]],color="blue", lw=1))

            except:
                print "Error attempting to plot empty edge data..."
            # print entry


import pickle
pickle.dump(myMesh, open("mesh_ns.pkl", "wb"))
# print "firing goedesics"
# N1 = 5
# for i in range(N1):
#     if not i%1:
#         print i, "out of", N, "..."
#     throw_geodesic_discrete(myMesh, ax)

print "Done!"

# draw_edge_data(myMesh)
# print myMesh.edge_slope_data

# plt.show()

