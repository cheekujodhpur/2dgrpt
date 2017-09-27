# this is the file where we draw geodescis

from mesh import Mesh
import numpy as np

# assuming p1, p2, p3 are tupes
def random_point(p1, p2, p3):
    p1 = np.array(p1)
    p2 = np.array(p2)
    p3 = np.array(p3)
    c1 = np.random.random()
    c2 = np.random.random()
    c3 = np.random.random()
    return (p1*c1 + p2*c2 + p3*c3)/(c1+c2+c3)

myMesh = Mesh([(0,0),(0,1),(1,0),(1,1)], (0,0))
myMesh.submesh(6)
num_triangle = len(myMesh.triangles)

# find start point in a random triangle
t_id = int(np.random.random()*num_triangle) #triangle_id
this_triangle = myMesh.triangles[t_id]

direction = np.array((np.random.random(), np.random.random()))
direction = direction / np.linalg.norm(direction)

while True:

    startpoint = random_point(
            myMesh.vertices[this_triangle[0]],
            myMesh.vertices[this_triangle[1]],
            myMesh.vertices[this_triangle[2]]
            )

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

    segment = filter(lambda x:intersections[x][1] <= 1 and intersections[x][1] >=0 and intersections[x][0]>=0, intersections.keys())[0]
    point_of_intersection = startpoint + intersections[segment][0]*direction

    print segment
    print this_triangle
    print point_of_intersection

    if segment=='ab':
        triangles = filter(lambda x:(this_triangle[0] in x) and (this_triangle[1] in x), myMesh.triangles)
    elif segment=='bc':
        triangles = filter(lambda x:(this_triangle[1] in x) and (this_triangle[2] in x), myMesh.triangles)
    elif segment=='ca':
        triangles = filter(lambda x:(this_triangle[2] in x) and (this_triangle[0] in x), myMesh.triangles)

    this_triangle = filter(lambda x:x!=this_triangle, triangles)
    if(this_triangle):
        this_triangle = this_triangle[0]
    else:
        break


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# theta = np.array([vertex[0] for vertex in myMesh.vertices])*np.pi
# phi = np.array([vertex[1] for vertex in myMesh.vertices])*2*np.pi
# X = np.sin(theta)*np.cos(phi)
# Y = np.sin(theta)*np.sin(phi)
# Z = np.cos(theta) 
# triangles = np.array([list(triangle) for triangle in myMesh.triangles])
# ax.plot_trisurf(X,Y,Z,triangles=triangles,color="gray",linewidth=0.2)
# ax.scatter(_X,_Y,_Z, color="red", s=1)

X = np.array([vertex[0] for vertex in myMesh.vertices])
Y = np.array([vertex[1] for vertex in myMesh.vertices])
Z = np.zeros(len(X))
triangles = np.array([list(triangle) for triangle in myMesh.triangles])
ax.plot_trisurf(X,Y,Z,triangles=triangles,shade=True,color="gray",linewidth=2)

plt.show()
