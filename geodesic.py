# this is the file where we draw geodescis

from mesh import Mesh
import numpy as np

# assuming p1, p2, p3 are tupes
def random_point(p1, p2, p3):
    p1 = np.array(p1)
    p2 = np.array(p2)
    p3 = np.array(p3)
    return reduce(lambda x,y:x+y*np.random.random(),[p1, p2, p3])/3.

myMesh = Mesh([(0,0),(0,1),(1,0),(1,1)], (0,0))
myMesh.submesh(5)
num_triangle = len(myMesh.triangles)

# find start point in a random triangle
t_id = int(np.random.random()*num_triangle) #triangle_id

startpoint = random_point(
        myMesh.vertices[myMesh.triangles[t_id][0]],
        myMesh.vertices[myMesh.triangles[t_id][1]],
        myMesh.vertices[myMesh.triangles[t_id][2]]
        )

direction = np.array((np.random.random(), np.random.random()))
direction = direction / np.linalg.norm(direction)

#TODO: Find the direction we are in

print startpoint, direction

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
theta = np.array([vertex[0] for vertex in myMesh.vertices])*np.pi
phi = np.array([vertex[1] for vertex in myMesh.vertices])*2*np.pi
X = np.sin(theta)*np.cos(phi)
Y = np.sin(theta)*np.sin(phi)
Z = np.cos(theta)

triangles = np.array([list(triangle) for triangle in myMesh.triangles])
ax.plot_trisurf(X,Y,Z,triangles=triangles,color="gray",linewidth=0.2)

plt.show()
