# this is the file where we draw geodescis

from mesh import Mesh

myMesh = Mesh([(0,0),(0,1),(1,0),(1,1)], (0,0))
myMesh.submesh(5)

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
ax.plot_trisurf(X,Y,Z,triangles=triangles,shade=True,color="gray",linewidth=2)
plt.show()
