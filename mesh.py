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
    edges = set([]) # 2-tuples holding ids of vertices which are edges
    regions = []
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

        for triangle in self.triangles:
            striangle = sorted(triangle)
            self.edges.add((striangle[0], striangle[1]))
            self.edges.add((striangle[1], striangle[2]))
            self.edges.add((striangle[0], striangle[2]))
            center = (np.array(self.vertices[striangle[0]]) + np.array(self.vertices[striangle[1]]) + np.array(self.vertices[striangle[2]]))/3.
            g = m.compute_metric(center, t="a1")
            self.regions.append((center[0], center[1], np.sqrt(np.linalg.det(g))))
        self.edges = list(self.edges)


    def write_to_poly(self, k=1.0):
        with open("input.poly", "w") as outfile:
            outfile.write(str(len(self.vertices)))
            outfile.write(" 2 0 0\n")
            
            for i in range(len(self.vertices)):
                outfile.write(str(i+1) + " ")
                outfile.write(str(self.vertices[i][0]) + " " + \
                        str(self.vertices[i][1]) + "\n")
           
            outfile.write(str(len(self.edges)) + " 0\n")
            for i in range(len(self.edges)):
                outfile.write(str(i+1) + " ")
                outfile.write(str(self.edges[i][0]+1) + " " + \
                        str(self.edges[i][1]+1) + "\n")

            outfile.write("0\n") # no holes
            outfile.write(str(len(self.regions)) + "\n")
            for i in range(len(self.regions)):
                outfile.write(str(i+1) + " ")
                outfile.write(str(self.regions[i][0]) + " " + \
                        str(self.regions[i][1]) + " " + \
                        str(k*self.regions[i][2]) + "\n")


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
