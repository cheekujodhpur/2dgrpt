# this is the file where we draw geodescis

from mesh import Mesh
import numpy as np
from scipy.integrate import odeint

# assuming p1, p2, p3 are tupes
def random_point(p1, p2, p3):
    p1 = np.array(p1)
    p2 = np.array(p2)
    p3 = np.array(p3)
    c1 = np.random.random()
    c2 = np.random.random()
    c3 = np.random.random()
    return (p1*c1 + p2*c2 + p3*c3)/(c1+c2+c3)


def check_for_incidence(triangle, point_of_intersection):
    a = np.array(myMesh.vertices[triangle[0]])
    b = np.array(myMesh.vertices[triangle[1]])
    c = np.array(myMesh.vertices[triangle[2]])
    pos = np.array(point_of_intersection)
    if abs(np.linalg.norm(pos-a) + np.linalg.norm(pos-b) - np.linalg.norm(a-b)) < 1e-4: 
        return True
    elif abs(np.linalg.norm(pos-b) + np.linalg.norm(pos-c) - np.linalg.norm(b-c)) < 1e-4: 
        return True
    elif abs(np.linalg.norm(pos-c) + np.linalg.norm(pos-a) - np.linalg.norm(c-a)) < 1e-4: 
        return True
    return False

myMesh = Mesh([(0,0),(0,1),(1,0),(1,1)], (0,0))
myMesh.submesh(6)
num_triangle = len(myMesh.triangles)

# find start point in a random triangle
t_id = int(np.random.random()*num_triangle) #triangle_id
this_triangle = myMesh.triangles[t_id]

direction = np.array((np.random.random(), np.random.random()))
direction = direction / np.linalg.norm(direction)

_theta = np.array([])
_phi = np.array([])

startpoint = random_point(
        myMesh.vertices[this_triangle[0]],
        myMesh.vertices[this_triangle[1]],
        myMesh.vertices[this_triangle[2]]
        )
saved_start = startpoint
saved_direction = direction

flipped = False

while True:

    # print map(lambda x:myMesh.vertices[x], np.array(this_triangle))
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
    _theta = np.concatenate((_theta, np.linspace(startpoint[0], point_of_intersection[0], 10)))
    _phi = np.concatenate((_phi, np.linspace(startpoint[1], point_of_intersection[1], 10)))

    pos = point_of_intersection
    # print pos
    # print segment
    if segment=='ab':
        triangles = filter(lambda x:(this_triangle[0] in x) or (this_triangle[1] in x), myMesh.triangles)
    elif segment=='bc':
        triangles = filter(lambda x:(this_triangle[1] in x) or (this_triangle[2] in x), myMesh.triangles)
    elif segment=='ca':
        triangles = filter(lambda x:(this_triangle[2] in x) or (this_triangle[0] in x), myMesh.triangles)

    # print triangles
    this_triangle = filter(lambda x:x!=this_triangle, triangles)
    # print this_triangle
    this_triangle = filter(lambda x:check_for_incidence(x, point_of_intersection),
            this_triangle)

    if len(this_triangle)==1:
        this_triangle = this_triangle[0]
    elif len(this_triangle)>1:
        print "Error. Too many neighbours"
        break
    else:
        if not flipped:
            startpoint = saved_start
            this_triangle = myMesh.triangles[t_id]
            direction = -direction
            flipped = True
            continue
        else:
            break

    startpoint = point_of_intersection
    direction[0] = direction[0] - 1 * 0.5* np.sin(2*startpoint[0]*np.pi)  * direction[1] * direction[1]
    # direction = direction / np.linalg.norm(direction)


# _theta = _theta*np.pi
# _phi = _phi*2*np.pi
# _X = 1.05*np.sin(_theta)*np.cos(_phi)
# _Y = 1.05*np.sin(_theta)*np.sin(_phi)
# _Z = 1.05*np.cos(_theta) 
_X = _theta
_Y = _phi
_Z = np.ones(len(_theta))*0.05

# print _X
# print _Y

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
ax = fig.add_subplot(111)
plt.hold(True)

# theta = np.array([vertex[0] for vertex in myMesh.vertices])*np.pi
# phi = np.array([vertex[1] for vertex in myMesh.vertices])*2*np.pi
# X = np.sin(theta)*np.cos(phi)
# Y = np.sin(theta)*np.sin(phi)
# Z = np.cos(theta) 
# triangles = np.array([list(triangle) for triangle in myMesh.triangles])
# ax.plot_trisurf(X,Y,Z,triangles=triangles,color="gray",linewidth=0.2)
# ax.scatter(_X,_Y,_Z, color="red", s=1)

# X = np.array([vertex[0] for vertex in myMesh.vertices])
# Y = np.array([vertex[1] for vertex in myMesh.vertices])
# Z = np.zeros(len(X))
# triangles = np.array([list(triangle) for triangle in myMesh.triangles])
# ax.plot_trisurf(X,Y,Z,triangles=triangles,shade=True,color="gray",linewidth=1)
ax.scatter(_X,_Y, color="red",s=0.2)

t = np.linspace(0,1,1000)
def dy_dt(y,t):
    dy_0 = y[2]
    dy_1 = y[3]*np.sin(np.pi*y[0])**2
    dy_2 = -0.5*np.sin(2*np.pi*y[0])*(y[3]**2)
    dy_3 = 0
    return np.array([dy_0, dy_1, dy_2, dy_3])

y0 = np.array([saved_start[0], saved_start[1], saved_direction[0], saved_direction[1]])
result = odeint(dy_dt, y0, t)
y0 = np.array([saved_start[0], saved_start[1], -saved_direction[0], -saved_direction[1]])
result = np.concatenate((result,odeint(dy_dt, y0, t)))
ax.plot(result[:,0], result[:,1], color="blue")

plt.show()
