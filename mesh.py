from collections import defaultdict
import math
import numpy as np
import os
from helpers import bresenham_and_mesh
from helpers import anglemod
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
    overlay_mesh = {}
    corners = []
    origin = ()
    edge_data = {}
    edge_slope_data = {}

    def __init__(self,corners,origin):
        """
        Needs:
        corners : list of four 2-tuples which describe the corners of the mesh
        """

        self.corners = corners
        self.origin = origin
        self.edge_data = defaultdict(list)
        self.edge_slope_data = defaultdict(list)

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
            g = m.compute_metric(center, t="a2")
            self.regions.append((center[0], center[1], np.sqrt(np.linalg.det(g))))
        self.edges = list(self.edges)


    def write_to_poly(self, k=1.0):
        """
        As protocol, always write to input.poly
        """
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


    def refine_using_Triangle(self, min_angle, overlay_size):
        """
        As protocol, always read from input.poly
        """
        os.system("./triangle -an -q"+str(min_angle)+" input.poly")

        nodes = open("input.1.node", "r").readlines()
        all_vertices = []
        num_nodes = [int(s) for s in nodes[0].split() if s.isdigit()][0]
        for i in range(num_nodes):
            node_data = [float(s) for s in nodes[i+1].split() if s.replace('.','',1).isdigit()]
            all_vertices.append((node_data[1], node_data[2]))

        self.vertices = all_vertices

        eles = open("input.1.ele", "r").readlines()
        all_triangles = []
        num_eles = [int(s) for s in eles[0].split() if s.isdigit()][0]
        for i in range(num_eles):
            ele_data = [int(s) for s in eles[i+1].split() if s.isdigit()]
            all_triangles.append((ele_data[1]-1, ele_data[2]-1, ele_data[3]-1))

        self.triangles = all_triangles

        neighbours = open("input.1.neigh", "r").readlines()
        neighbours_dict = {}
        num_neighs = [int(s) for s in eles[0].split() if s.isdigit()][0]
        for i in range(num_neighs):
            neigh_data = [int(s) for s in neighbours[i+1].split() if s.isdigit()]
            neighbours_dict[(neigh_data[0]-1)] = filter(lambda x:x>=0,map(lambda x:x-1, neigh_data[1:]))

        self.neighbours = neighbours_dict


        # At this point, we hash everything
        # (1x1) -> 64 parts
        overlay_mesh = {}
        for tl in self.triangles:
            x1 = int(self.vertices[tl[0]][0]/overlay_size)
            y1 = int(self.vertices[tl[0]][1]/overlay_size)

            x2 = int(self.vertices[tl[1]][0]/overlay_size)
            y2 = int(self.vertices[tl[1]][1]/overlay_size)

            #TODO: Maybe I am doing double sorting
            # There is sorting later as well while geodesic
            a,b = tuple(sorted((tl[0],tl[1])))
            bresenham_and_mesh(overlay_mesh, x1, y1, x2, y2, a, b)

            x1 = int(self.vertices[tl[2]][0]/overlay_size)
            y1 = int(self.vertices[tl[2]][1]/overlay_size)

            x2 = int(self.vertices[tl[1]][0]/overlay_size)
            y2 = int(self.vertices[tl[1]][1]/overlay_size)

            a,b = tuple(sorted((tl[2],tl[1])))
            bresenham_and_mesh(overlay_mesh, x1, y1, x2, y2, a, b)

            x1 = int(self.vertices[tl[0]][0]/overlay_size)
            y1 = int(self.vertices[tl[0]][1]/overlay_size)

            x2 = int(self.vertices[tl[2]][0]/overlay_size)
            y2 = int(self.vertices[tl[2]][1]/overlay_size)

            a,b = tuple(sorted((tl[0],tl[2])))
            bresenham_and_mesh(overlay_mesh, x1, y1, x2, y2, a, b)

        self.overlay_mesh = overlay_mesh

    def churn_edge_data(self):
        print len(self.edge_data), "edges available..."
        print "Begin churning..."

        for edge in self.edge_data:
            all_samples = self.edge_data[edge]
            anglesamples = []
            for sample in all_samples:
                input_angle = anglemod(np.arctan2(sample[1][1], sample[1][0]))
                output_angle = anglemod(np.arctan2(sample[0][1] + sample[1][1], sample[1][0] + sample[0][0]))

                anglesamples.append((input_angle, output_angle))
                # print input_angle, output_angle

            anglesamples = sorted(anglesamples, key=lambda x:x[0])
            slopesamples = []
            N = len(anglesamples)
            for i in range(N-1):
                m = (anglesamples[i+1][1]-anglesamples[i][1])/(anglesamples[i+1][0]-anglesamples[i][0])
                b = anglesamples[i][1] - m*anglesamples[i][0]
                # already sorted by input_angle
                slopesamples.append((anglesamples[i+1][0], m, b))

            if N>1:
                m = (anglesamples[N-1][1]-anglesamples[0][1])/(anglesamples[N-1][0]-anglesamples[0][0])
                b = anglesamples[0][1] - m*anglesamples[0][0]
                slopesamples.append((anglesamples[0][0], m, b))

            self.edge_slope_data[edge] = slopesamples

        print "End churning..."

    def draw(self, ax):
        X = np.array([vertex[0] for vertex in self.vertices])
        Y = np.array([vertex[1] for vertex in self.vertices])
        triangles = np.array([list(triangle) for triangle in self.triangles])
        ax.tripcolor(X,Y,triangles=triangles, facecolors=np.zeros(len(self.triangles)), edgecolors='k', cmap='Blues')

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
