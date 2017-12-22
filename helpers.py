import numpy as np

# Cast to ndarray, would work
# if pi is tuple OR list
def random_point(p1, p2, p3):
    p1 = np.array(p1)
    p2 = np.array(p2)
    p3 = np.array(p3)
    c1 = np.random.random()
    c2 = np.random.random()
    c3 = np.random.random()
    return (p1*c1 + p2*c2 + p3*c3)/(c1+c2+c3)


# again, we simply cast everything
# for safety
def check_for_incidence(vertices, triangle, point_of_intersection, threshold):
    a = np.array(vertices[triangle[0]])
    b = np.array(vertices[triangle[1]])
    c = np.array(vertices[triangle[2]])
    pos = np.array(point_of_intersection)
    if abs(np.linalg.norm(pos-a) + np.linalg.norm(pos-b) - np.linalg.norm(a-b)) < threshold: 
        return True
    elif abs(np.linalg.norm(pos-b) + np.linalg.norm(pos-c) - np.linalg.norm(b-c)) < threshold: 
        return True
    elif abs(np.linalg.norm(pos-c) + np.linalg.norm(pos-a) - np.linalg.norm(c-a)) < threshold: 
        return True
    return False


def return_incidence(vertices, triangle, point_of_intersection):
    a = np.array(vertices[triangle[0]])
    b = np.array(vertices[triangle[1]])
    c = np.array(vertices[triangle[2]])
    pos = np.array(point_of_intersection)
    one = abs(np.linalg.norm(pos-a) + np.linalg.norm(pos-b) - np.linalg.norm(a-b)) 
    two = abs(np.linalg.norm(pos-b) + np.linalg.norm(pos-c) - np.linalg.norm(b-c)) 
    tre = abs(np.linalg.norm(pos-c) + np.linalg.norm(pos-a) - np.linalg.norm(c-a)) 
    return min(one, two, tre)


def find_closest_edge(ovmesh, pt, vertices):
    x = int(pt[0]/0.125)
    y = int(pt[1]/0.125) 
    if (x,y) not in ovmesh:
        return -1, 0
    distance = 10000
    minEdge = (0,0)
    for edge in ovmesh[(x,y)]:
        a = np.array(vertices[edge[0]])
        b = np.array(vertices[edge[1]])
        dist = abs(np.linalg.norm(pt-a) + np.linalg.norm(pt-b) - np.linalg.norm(a-b)) 
        if dist < distance:
            distance = dist
            minEdge = edge

    if distance<1e-2:
        return minEdge, distance
    else:
        return -1, 0
        
