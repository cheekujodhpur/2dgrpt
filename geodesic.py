# this is the file where we draw geodescis

import os
from mesh import Mesh
import numpy as np
from scipy.integrate import odeint
from helpers import *

# Functions for integration
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
myMesh.refine_using_Triangle(18, 0.125)

num_triangle = len(myMesh.triangles)

# find start point in a random triangle
t_id = int(np.random.random()*num_triangle) #triangle_id
this_triangle = myMesh.triangles[t_id]

direction = np.array((np.random.random(), np.random.random()))
# direction = np.array((0.0,1.0))
direction = direction / np.linalg.norm(direction)

startpoint = random_point(
        myMesh.vertices[this_triangle[0]],
        myMesh.vertices[this_triangle[1]],
        myMesh.vertices[this_triangle[2]]
        )


# This is the covariant direction, I have to correct it to contravariant direction
direction[0] = direction[0]
direction[1] = direction[1]/startpoint[0]
direction = direction / np.linalg.norm(direction)

saved_start = startpoint[:]
saved_direction = np.copy(direction)
stid = t_id

print "triangle id", t_id
print "start", saved_start
print "direction", saved_direction

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
ax = fig.add_subplot(111)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_xlim([0,5])
ax.set_ylim([0,5])
# ax.set_zlim([0,1])

import time
time_ode = 0
time_mesh = 0

st = time.time()
t = np.linspace(0,5,1000)

##########
# ODEINT #
##########
from scipy.integrate import ode
myOde = ode(f, jac).set_integrator("dopri5")
y0 = np.array([saved_start[0], saved_start[1], saved_direction[0], saved_direction[1]])
myOde.set_initial_value(y0,0)
t1 = 5
dt = 5/1000.

ode_pc = 0

result = []
while myOde.successful() and myOde.t < t1:
    if result:
        xx, dist = find_closest_edge(myMesh.overlay_mesh, np.array(result[-1][:2]), myMesh.vertices)
        if xx != -1:
            l = Line2D([myMesh.vertices[xx[0]][0],myMesh.vertices[xx[1]][0]], \
                    [myMesh.vertices[xx[0]][1],myMesh.vertices[xx[1]][1]], lw=2)
            ax.add_line(l)
    result.append(myOde.integrate(myOde.t+dt))
    ode_pc = ode_pc + 1

result = np.array(result)
# result = odeint(dy_dt, y0, t)
# _X = 1.05*np.sin(result[:,0]*np.pi)*np.cos(result[:,1]*2*np.pi)
# _Y = 1.05*np.sin(result[:,0]*np.pi)*np.sin(result[:,1]*2*np.pi)
# _Z = 1.05*np.cos(result[:,0]*np.pi) 
# ax.plot(_X, _Y,_Z, color="blue")
time_ode = time_ode + time.time() - st
st = time.time()
ax.plot(result[:,0], result[:,1], color="black", linewidth=1)
# result = np.concatenate((result,odeint(dy_dt, y0, t)))

myOde2 = ode(f, jac).set_integrator("dopri5")
y0 = np.array([saved_start[0], saved_start[1], -saved_direction[0], -saved_direction[1]])
myOde2.set_initial_value(y0,0)
t1 = 5
dt = 5/1000.

result = []
while myOde2.successful() and myOde2.t < t1:
    if result:
        xx, dist = find_closest_edge(myMesh.overlay_mesh, np.array(result[-1][:2]), myMesh.vertices)
        if xx != -1:
            ax.add_line(Line2D([myMesh.vertices[xx[0]][0],myMesh.vertices[xx[1]][0]], \
                    [myMesh.vertices[xx[0]][1],myMesh.vertices[xx[1]][1]], lw=1))
    result.append(myOde2.integrate(myOde2.t+dt))
    ode_pc = ode_pc + 1
# result2 = odeint(dy_dt, y0, t)
result = np.array(result)
# _X = 1.05*np.sin(result[:,0]*np.pi)*np.cos(result[:,1]*2*np.pi)
# _Y = 1.05*np.sin(result[:,0]*np.pi)*np.sin(result[:,1]*2*np.pi)
# _Z = 1.05*np.cos(result[:,0]*np.pi) 
# ax.plot(_X, _Y, _Z, color="green")
time_ode = time_ode + time.time() - st
ax.plot(result[:,0], result[:,1], color="black", linewidth=1)
##########
# ODEINT #
##########


first_time = True
def do_iteration(rate=1e-4,pcol="black"):

    st = time.time()

    startpoint = saved_start[:]
    direction = -np.copy(saved_direction)

    # coloring faces
    fcolors = np.zeros(num_triangle)

    covdir = np.copy(direction)
    covdirn = -np.copy(covdir)
    global mesh_pc
    mesh_pc = 0
    flipped = False
    start_points = []
    directions = []
    global t_id
    this_triangle = myMesh.triangles[t_id]

    while True:

        start_points.append(startpoint)
        directions.append(direction)
        mesh_pc = mesh_pc + 1

        fcolors[myMesh.triangles.index(this_triangle)] = 1

        # now figure out the triangle intersection
        a = np.array(myMesh.vertices[this_triangle[0]])
        b = np.array(myMesh.vertices[this_triangle[1]])
        c = np.array(myMesh.vertices[this_triangle[2]])
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

        segment = filter(lambda x:intersections[x][1] <= 1 and intersections[x][1] >=0, intersections.keys())
        segment = max(segment, key=lambda x: intersections[x][0])
        point_of_intersection = startpoint + intersections[segment][0]*direction
        ax.add_line(Line2D([startpoint[0],point_of_intersection[0]] \
                    ,[startpoint[1], point_of_intersection[1]],color="red", lw=1))

        pos = point_of_intersection
        t_id = filter(lambda x:check_for_incidence(myMesh.vertices, myMesh.triangles[x], point_of_intersection, 1e-4),
                myMesh.neighbours[t_id])

        if len(t_id)==1:
            t_id = t_id[0]
        elif len(t_id)>1:
            print "Error. Too many neighbours"
            t_id = sorted(t_id, key=lambda x:return_incidence(myMesh.vertices, myMesh.triangles[x], point_of_intersection))[0]
            # break
        else:
            if not flipped:
                startpoint = saved_start
                this_triangle = myMesh.triangles[stid]
                t_id = stid
                direction = -np.copy(saved_direction)
                covdir = np.copy(covdirn)
                flipped = True
                continue
            else:
                break

        startpoint = point_of_intersection
        this_triangle = myMesh.triangles[t_id]
        covdir[0] = covdir[0] + rate * 0.5 * covdir[1] * covdir[1]/(startpoint[0]*startpoint[0])
        covdir = covdir / np.linalg.norm(covdir)

        direction[0] = covdir[0]
        direction[1] = covdir[1]/startpoint[0]
        direction = direction / np.linalg.norm(direction)

    global time_mesh
    time_mesh = time_mesh + time.time() - st
    global first_time
    if first_time:
        X = np.array([vertex[0] for vertex in myMesh.vertices])
        Y = np.array([vertex[1] for vertex in myMesh.vertices])
        triangles = np.array([list(triangle) for triangle in myMesh.triangles])
        ax.tripcolor(X,Y,triangles=triangles, facecolors=np.zeros(num_triangle), edgecolors='k', cmap='Blues')
        first_time = False
    # ax.scatter(_X, _Y, color=pcol, s=2)

do_iteration(pcol="green", rate=1e-3)

plt.show()

print "mesh size", len(myMesh.triangles)
print "computation save", time_mesh/float(time_ode)

