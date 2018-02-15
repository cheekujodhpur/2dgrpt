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

N_checks = 10 # number of steps to error control
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
myMesh.write_to_poly(k=0.5)
refined_size = 0.05
myMesh.refine_using_Triangle(18, refined_size)

def throw_geodesic_integrating(mesh, dt=0.01):
    num_triangle = len(mesh.triangles)
    t_id = int(np.random.random()*num_triangle) #triangle_id
    this_triangle = mesh.triangles[t_id]

    # Start point randomly chosen
    startpoint = random_point( mesh.vertices[this_triangle[0]], mesh.vertices[this_triangle[1]], mesh.vertices[this_triangle[2]]
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


def find_trial_error(mesh, ax, mod, dt, t_id, startpoint, covdir, save=False, dbg=False):

    _lx = sorted(mesh.corners, key=lambda x:x[0])[0][0]
    _rx = sorted(mesh.corners, key=lambda x:x[0])[-1][0]
    _ly = sorted(mesh.corners, key=lambda x:x[1])[0][1]
    _ry = sorted(mesh.corners, key=lambda x:x[1])[-1][1]

    num_triangle = len(mesh.triangles)
    this_triangle = mesh.triangles[t_id]

    g_int = ode(f, jac).set_integrator("dopri5")
    y0 = np.array([startpoint[0], startpoint[1], covdir[0], covdir[1]])
    g_int.set_initial_value(y0,0)
    result = [y0]
    result.append(g_int.integrate(g_int.t+dt))
    direction = result[1][:2] - startpoint
    direction = direction / np.linalg.norm(direction)

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
        print "[WARN] Segment empty", intersections

    local_edge = (0,0)
    if segment == 'ab':
        local_edge = (this_triangle[0], this_triangle[1])
    elif segment == 'bc':
        local_edge = (this_triangle[1], this_triangle[2])
    elif segment == 'ca':
        local_edge = (this_triangle[2], this_triangle[0])

    point_of_intersection = startpoint + intersections[segment][0]*direction*(1+1e-5)

    if ax is not None:
        ax.add_line(Line2D([startpoint[0],point_of_intersection[0]] \
                    ,[startpoint[1], point_of_intersection[1]],color="red", lw=1))

    new_dir = np.array([np.cos(mod), np.sin(mod)])

    if save:
        mesh.edge_data[tuple(sorted(local_edge))].append([new_dir-direction, direction])

    if dbg is not False:
        print dbg, mod*180./np.pi, new_dir-direction, direction

    dev_g_int = ode(f, jac).set_integrator("dopri5")

    new_dir[0] = new_dir[0]*point_of_intersection[0]
    new_dir = new_dir / np.linalg.norm(new_dir)
    # making sure to take covariant direction here
    dev_y0 = np.array([point_of_intersection[0], point_of_intersection[1], new_dir[0], new_dir[1]])
    # print "checking", dev_y0[2:], "..." 
    dev_g_int.set_initial_value(dev_y0,0)
    dev_result = [dev_y0]

    global N_checks
    N = N_checks 
    Ni = 0
    while Ni < N and g_int.successful() and (_lx < g_int.y[0] < _rx) and (_ly < g_int.y[1] < _ry):
        result.append(g_int.integrate(g_int.t+dt))
        Ni = Ni + 1

    Ni = 0
    while Ni < N and dev_g_int.successful() and (_lx < dev_g_int.y[0] < _rx) and (_ly < dev_g_int.y[1] < _ry):
        dev_result.append(dev_g_int.integrate(dev_g_int.t+dt))
        Ni = Ni + 1

    result = np.array(result)
    dev_result = np.array(dev_result)
    if ax is not None:
        ax.plot(result[:,0], result[:,1], color = "green", lw = 1)
        ax.plot(dev_result[:,0], dev_result[:,1], color = "blue", lw = 1)

    target = result[-1][:2]-result[0][:2]
    target = np.arctan2(target[1], target[0])
    calcu = dev_result[-1][:2]-dev_result[0][:2]
    calcu = np.arctan2(calcu[1], calcu[0])

    # if save:
    #     print "saving with an error of", np.abs(target-calcu)
    #     return
    return np.abs(target-calcu)


def throw_geodesic_mark(mesh, ax, seed, tau, dt=0.01):

    num_triangle = len(mesh.triangles)
    t_id = seed % num_triangle #triangle_id
    this_triangle = mesh.triangles[t_id]

    # Start point randomly chosen
    startpoint = random_point(
            mesh.vertices[this_triangle[0]],
            mesh.vertices[this_triangle[1]],
            mesh.vertices[this_triangle[2]]
            )

    # Choosing a random direction
    covdir = np.array((np.random.random()*2-1, np.random.random()*2-1))
    # Renormalizing
    covdir = covdir / np.linalg.norm(covdir)

    def f_minimizer(x):
    #     return find_trial_error(mesh, None, x, dt, t_id, startpoint, covdir)
        return find_trial_error(mesh, ax, x, dt, t_id, startpoint, covdir)

    contradir = np.copy(covdir)
    contradir[1] = contradir[1]*startpoint[0]
    contradir = contradir / np.linalg.norm(contradir)

    one_angle = np.arctan2(contradir[1], contradir[0])
    result = nelder_mead(one_angle-0.4*np.pi, one_angle+0.4*np.pi, tau, f_minimizer)

    # if find_trial_error(mesh, None, result, dt, t_id, startpoint, covdir) > tau:
    #     # print "all error details:"
    #     # print "startpoint", startpoint
    #     # print "direction", covdir
    #     # find_trial_error(mesh, ax, result, dt, t_id, startpoint, covdir)
    #     def f_minimizer_dbg(x):
    #         return find_trial_error(mesh, ax, x, dt, t_id, startpoint, covdir, dbg="Checking...")
    #     nelder_mead(one_angle-0.45*np.pi, one_angle+0.45*np.pi, tau, f_minimizer_dbg)

    find_trial_error(mesh, None, result, dt, t_id, startpoint, covdir, save=True)


fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel("x")
ax.set_ylabel("y")
_lx = sorted(myMesh.corners, key=lambda x:x[0])[0][0]
_rx = sorted(myMesh.corners, key=lambda x:x[0])[-1][0]
ax.set_xlim([_lx, _rx])
_ly = sorted(myMesh.corners, key=lambda x:x[1])[0][1]
_ry = sorted(myMesh.corners, key=lambda x:x[1])[-1][1]
ax.set_ylim([_ly, _ry])

# res = throw_geodesic_integrating(myMesh)
# ax.plot(res[:,0], res[:,1], color="black", linewidth=1)

# Plot the mesh
myMesh.draw(ax)

# Discrete geodesic
for i in range(NOS_EDGE*len(myMesh.triangles)):
    print "Marking", i, "...",
    throw_geodesic_mark(myMesh, ax, i, 1e-4, dt=1e-2)
    print "done."

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


# draw_edge_data(myMesh)

import pickle
pickle.dump(myMesh, open("feb15.pkl", "wb"))

print "Done!"

plt.show()
# plt.savefig('fig_feb9.png')

