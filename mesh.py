import math
from collections import defaultdict

# Defining the Mesh
class Mesh:
    """
    The class which defines a mesh, given the corner
    coordinates and the metric. The triangles are denser
    in regions where curvature is higher
    """
    vertices = [] # the list of vertices, each a tuple, 2-tuple
    triangles = [] # 3-tuples of ids of vertices making a triangle

    def __init__(self,corners):
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
        triangles = [(0,1,2),(0,2,3)] # Writing triangles is easy, the main purpose was always to sort the corners
        
        self.vertices = self.vertices + sorted_corners
        self.triangles = self.triangles + triangles

    def sortedOrder(self, order):
        # sort points counter-clockwise
        print order

    def submesh(self,n=1):
        """
        Create one level of submesh
        Catmull-Clark Subdivision Surface
        """
        #TODO: Make it Catmull-Clark; naive right now
        
        index = len(self.vertices)
        for nn in range(n):

            vertex_ids = [i for i in range(len(self.vertices))]
            edges = [(x,y) for y in vertex_ids for x in vertex_ids[:y]]
            faces = {}

            for edge in edges:
                anyfaces = filter(lambda x:x.count(edge[0]) and x.count(edge[1]), self.triangles)
                if len(anyfaces)==2:
                    faces[frozenset(edge)] = anyfaces

            facepoints = {}

            for triangle in self.triangles:
                center = tuple(map(lambda x:float(x)/len(self.vertices),\
                        reduce(lambda x,y:(x[0]+y[0],x[1]+y[1]),\
                        map(lambda x:self.vertices[x],triangle))))
                facepoints[frozenset(triangle)] = center
            
            print "facepoints", facepoints

            edgepoints = {}

            for edge in faces.keys():
                facefuredge = faces[frozenset(edge)]
                p1 = facepoints[frozenset(facefuredge[0])]
                p3 = self.vertices[list(edge)[0]]
                p4 = self.vertices[list(edge)[1]]
                p2 = facepoints[frozenset(facefuredge[1])]
                p = ((p1[0]+p2[0]+p3[0]+p4[0])/4., (p1[1]+p2[1]+p3[1]+p4[1])/4.)

                edgepoints[frozenset(edge)] = p

            print "edgepoints", edgepoints

            new_vertices = self.vertices[:]
            for i in vertex_ids:
                pFaces = filter(lambda x:x.count(i),self.triangles)
                if pFaces:
                    F = tuple(map(lambda x:float(x)/len(pFaces),\
                            reduce(lambda x,y:(x[0]+y[0],x[1]+y[1]),\
                        map(lambda x:facepoints[frozenset(x)],pFaces))))

                pEdges = filter(lambda x:x.count(i),[list(kk) for kk in faces.keys()])
                if pEdges:
                    R = tuple(map(lambda x:float(x)/len(pEdges),\
                            reduce(lambda x,y:(x[0]+y[0],x[1]+y[1]),\
                        map(lambda x:((self.vertices[x[0]][0]+self.vertices[x[1]][0])/2.,\
                        (self.vertices[x[0]][1]+self.vertices[x[1]][1])/2.),\
                        pEdges))))

                if pEdges and pFaces:
                    print F,R
                    little_x = (F[0] + 2*R[0] + (len(pEdges)-3)*self.vertices[i][0])/float(len(pEdges))
                    little_y = (F[1] + 2*R[1] + (len(pEdges)-3)*self.vertices[i][1])/float(len(pEdges))

                    new_vertices[i] = (little_x,little_y)

            self.vertices = new_vertices
            print "vertices", self.vertices

            newTriangles = []

            for frozen_face in facepoints.keys():
                face = list(frozen_face)

                p1 = facepoints[frozen_face]
                self.vertices.append(p1) #index update
                i1 = index
                index = index+1

                p2 = edgepoints[frozenset((face[0],face[1]))]
                self.vertices.append(p2) #index update
                i2 = index
                index = index+1

                p3 = self.vertices[face[0]]
                i3 = face[0]
                
                newTriangles.append((i1,i2,i3))
                ##

                p3 = self.vertices[face[1]]
                i3 = face[1]
                
                newTriangles.append((i1,i2,i3))
                ####

                p2 = edgepoints[frozenset((face[1],face[2]))]
                self.vertices.append(p2) #index update
                i2 = index
                index = index+1

                p3 = self.vertices[face[1]]
                i3 = face[1]
                
                newTriangles.append((i1,i2,i3))
                ##

                p3 = self.vertices[face[2]]
                i3 = face[2]
                
                newTriangles.append((i1,i2,i3))
                ####
                
                p2 = edgepoints[frozenset((face[2],face[0]))]
                self.vertices.append(p2) #index update
                i2 = index
                index = index+1

                p3 = self.vertices[face[2]]
                i3 = face[2]
                
                newTriangles.append((i1,i2,i3))
                ##

                p3 = self.vertices[face[0]]
                i3 = face[0]
                
                newTriangles.append((i1,i2,i3))
                ##

            self.triangles = newTriangles



myMesh = Mesh([(0,0),(0,1),(1,0),(1,1)])
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
