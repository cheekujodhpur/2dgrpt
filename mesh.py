import math
import numpy as np

from metric import Metric

# Defining the Mesh
class Mesh:
    """
    The class which defines a mesh, given the corner
    coordinates and the metric. The triangles are denser
    in regions where curvature is higher
    """
    vertices = [] # the list of vertices, each a tuple, 2-tuple
    triangles = [] # 3-tuples of ids of vertices making a triangle
    metric = None

    def __init__(self,corners,origin):
        """
        Needs:
        corners : list of four 2-tuples which describe the corners of the mesh
        """

        # NOTE: assumes a 2-tuple
        center = tuple(map(lambda x:float(x)/len(corners),reduce(lambda x,y:(x[0]+y[0],x[1]+y[1]),corners)))
        angles = []
        for corner in corners:
            vector = tuple(map(lambda x,y:x-y,corner,center))
            angles.append(math.atan2(vector[0],vector[1]))

        sorted_corners = [x for _,x in sorted(zip(angles,corners))][::-1]

        triangles = [(0,1,2),(0,2,3)]

        # Writing triangles is easy, the main purpose was always to sort the corners
        self.vertices = self.vertices + sorted_corners
        self.triangles = self.triangles + triangles

        self.metric = Metric(origin)


    def submesh(self,n=1):
        """
        Create one level of submesh
        """
        m = self.metric
        index = len(self.vertices)

        for nn in range(n):
            newTriangles = []
            for triangle in self.triangles:
                p1 = ((self.vertices[triangle[0]][0] + self.vertices[triangle[1]][0])/2.,(self.vertices[triangle[0]][1] + self.vertices[triangle[1]][1])/2.)
                p2 = ((self.vertices[triangle[1]][0] + self.vertices[triangle[2]][0])/2.,(self.vertices[triangle[1]][1] + self.vertices[triangle[2]][1])/2.)
                p3 = ((self.vertices[triangle[2]][0] + self.vertices[triangle[0]][0])/2.,(self.vertices[triangle[2]][1] + self.vertices[triangle[0]][1])/2.)

                center = tuple(np.mean(np.array(map(lambda x:np.array(self.vertices[x]),triangle)),axis=0))
                print p1, p2, p3, center

                mean_metric = np.mean([m.compute_metric(p1), m.compute_metric(p2), m.compute_metric(p3)], axis=0)
                center_metric = m.compute_metric(center)
                diff = np.linalg.norm(mean_metric-center_metric) 

                if diff < 1e-1:
                    continue

                self.vertices.append(p1)
                self.vertices.append(p2)
                self.vertices.append(p3)

                newTriangles.append((index,index+1,index+2))
                newTriangles.append((index,triangle[1],index+1))
                newTriangles.append((index+1,triangle[2],index+2))
                newTriangles.append((index+2,triangle[0],index))

                index = index+3
            self.triangles = newTriangles



myMesh = Mesh([(0,0),(0,1),(1,0),(1,1)], (0,0))
myMesh.submesh(1)

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X = np.array([vertex[0] for vertex in myMesh.vertices])
Y = np.array([vertex[1] for vertex in myMesh.vertices])
Z = np.zeros(len(X))
triangles = np.array([list(triangle) for triangle in myMesh.triangles])
ax.plot_trisurf(X,Y,Z,triangles=triangles,shade=True,color="gray",linewidth=2)
plt.show()
