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
    neighbours = {}

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


    def calculateArea(self, p1, p2, p3):
        p1 = np.array(p1)
        p2 = np.array(p2)
        p3 = np.array(p3)

        if len(p1)==2:
            p1 = np.insert(p1,-1,1)
        if len(p2)==2:
            p2 = np.insert(p2,-1,1)
        if len(p3)==2:
            p3 = np.insert(p3,-1,1)

        return 0.5*np.abs(np.linalg.det(np.array([p1,p2,p3])))


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

                three_metrics = np.array([m.compute_metric(p1), m.compute_metric(p2), m.compute_metric(p3)])
                diff = np.linalg.norm(three_metrics.std(0))
                # mean_metric = np.mean(three_metrics, axis=0)
                # center_metric = m.compute_metric(center)
                # diff = np.linalg.norm(mean_metric-center_metric) 

                area = self.calculateArea(p1,p2,p3)
                # TODO: replace this by determinantof metric
                # sint = np.sin(center[0]*np.pi)
                measure = area/(np.sqrt(center[0]))
                # call this 1e-2 `scale*scale`, scale=1e-1
                # relevance will be clear later
                if measure < 1e-4:
                    newTriangles.append(triangle)
                    continue

                if self.vertices.count(p1) > 0:
                    i1 = self.vertices.index(p1)
                else:
                    self.vertices.append(p1)
                    i1 = index
                    index = index + 1

                if self.vertices.count(p2) > 0:
                    i2 = self.vertices.index(p2)
                else:
                    self.vertices.append(p2)
                    i2 = index
                    index = index + 1

                if self.vertices.count(p3) > 0:
                    i3 = self.vertices.index(p3)
                else:
                    self.vertices.append(p3)
                    i3 = index
                    index = index + 1

                newTriangles.append((i1,i2,i3))
                newTriangles.append((i1,triangle[1],i2))
                newTriangles.append((i2,triangle[2],i3))
                newTriangles.append((i3,triangle[0],i1))

            self.triangles = newTriangles


# myMesh = Mesh([(0,0),(0,1),(1,0),(1,1)], (0,0))
# myMesh.submesh(5)
# 
# import numpy as np
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# # X = np.array([vertex[0] for vertex in myMesh.vertices])
# # Y = np.array([vertex[1] for vertex in myMesh.vertices])
# # Z = np.zeros(len(X))
# # triangles = np.array([list(triangle) for triangle in myMesh.triangles])
# # ax.plot_trisurf(X,Y,Z,triangles=triangles,shade=True,color="gray",linewidth=2)
# 
# theta = np.array([vertex[0] for vertex in myMesh.vertices])*np.pi
# phi = np.array([vertex[1] for vertex in myMesh.vertices])*2*np.pi
# X = np.sin(theta)*np.cos(phi)
# Y = np.sin(theta)*np.sin(phi)
# Z = np.cos(theta)
# 
# triangles = np.array([list(triangle) for triangle in myMesh.triangles])
# ax.plot_trisurf(X,Y,Z,triangles=triangles,shade=True,color="gray",linewidth=2)
# 
# plt.show()
