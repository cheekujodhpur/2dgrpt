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

def return_incidence(triangle, point_of_intersection):
    a = np.array(myMesh.vertices[triangle[0]])
    b = np.array(myMesh.vertices[triangle[1]])
    c = np.array(myMesh.vertices[triangle[2]])
    pos = np.array(point_of_intersection)
    one = abs(np.linalg.norm(pos-a) + np.linalg.norm(pos-b) - np.linalg.norm(a-b)) 
    two = abs(np.linalg.norm(pos-b) + np.linalg.norm(pos-c) - np.linalg.norm(b-c)) 
    tre = abs(np.linalg.norm(pos-c) + np.linalg.norm(pos-a) - np.linalg.norm(c-a)) 
    return min(one, two, tre)

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
myMesh.submesh(6)
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


saved_start = startpoint[:]
saved_direction = np.copy(direction)

print "triangle id", t_id
print "start", saved_start
print "direction", saved_direction

import numpy as np
import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_xlim([0,5])
ax.set_ylim([0,5])
plt.hold(True)

t = np.linspace(0,5,1000)

from scipy.integrate import ode
myOde = ode(f, jac).set_integrator("dopri5")
y0 = np.array([saved_start[0], saved_start[1], saved_direction[0], saved_direction[1]])
myOde.set_initial_value(y0,0)
t1 = 5
dt = 5/1000.

result = []
while myOde.successful() and myOde.t < t1:
    result.append(myOde.integrate(myOde.t+dt))

result = np.array(result)
ax.plot(result[:,0], result[:,1], color="black", linewidth=3)

myOde2 = ode(f, jac).set_integrator("dopri5")
y0 = np.array([saved_start[0], saved_start[1], -saved_direction[0], -saved_direction[1]])
myOde2.set_initial_value(y0,0)
t1 = 5
dt = 5/1000.

result = []
while myOde2.successful() and myOde2.t < t1:
    result.append(myOde2.integrate(myOde2.t+dt))

result = np.array(result)
ax.plot(result[:,0], result[:,1], color="black", linewidth=3)

plt.show()
