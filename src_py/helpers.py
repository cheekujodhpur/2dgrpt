import numpy as np

phi = (1. + np.sqrt(5))/2.
resphi = 2. - phi

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


def find_closest_edge(ovmesh, pt, vertices, rs):
    x = int(pt[0]/rs)
    y = int(pt[1]/rs) 
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

    if distance<rs:
        return minEdge, distance
    else:
        return -1, 0


def find_shooting_edge(ovmesh, pt, direction, vertices, rs):
    x = int(pt[0]/rs)
    y = int(pt[1]/rs) 
    if (x,y) not in ovmesh:
        return -1, 0
    mint = 10000
    minEdge = (0,0)
    for edge in ovmesh[(x,y)]:
        a = np.array(vertices[edge[0]])
        b = np.array(vertices[edge[1]])

        t = np.cross((a-pt), (b-a))/np.cross(direction, (b-a))
        s = np.cross((a-pt), direction)/np.cross(direction, (b-a))

        if abs(t) < mint and 0 <= s <= 1:
            mint = t
            minEdge = edge

    if mint<10000:
        return minEdge, mint
    else:
        return -1, 0


def bresenham_and_mesh(ovm, x1, y1, x2, y2, a, b):
    """
    ovm is the overlay mesh, a and b are vertex ids
    xi and yi are the starting and endpoints obviously
    """
    #TODO: make it conservative
    deltax = x2-x1
    deltay = y2-y1

    if deltax == 0:
        x = int(x1)
        y = int(y1)
        for y in range(y1, y2+1):
            if (x,y) not in ovm:
                ovm[(x,y)] = set([])
            else:
                ovm[(x,y)].add((a, b))

    else:
        deltaerr = np.abs(float(deltay)/deltax)
        error = 0.0
        y = int(y1)
        for x in range(x1, x2+1):
            if (x,y) not in ovm:
                ovm[(x,y)] = set([])
            else:
                ovm[(x,y)].add((a, b))

            error = error + deltaerr
            while error>=0.5:
                y = y + np.sign(deltay) * 1
                error = error - 1.0


def anglemod(ang):
    modded = (2*np.pi+ang)%(2*np.pi)
    return modded

# Drawing the edge data
def draw_edge_data(myMesh):
    for edge in myMesh.edge_data.keys():
        for el in myMesh.edge_data[edge]:
            left = myMesh.vertices[edge[0]]
            right = myMesh.vertices[edge[1]]
            x = (left[0]+right[0])/2.
            y = (left[1]+right[1])/2.
            try:
                ax.add_line(Line2D([x,x+0.01*el[1][0]] \
                            ,[y, y+0.01*el[1][1]],color="green", lw=1))

                ax.add_line(Line2D([x+0.01*el[1][0],x+2*0.01*el[1][0]+0.01*el[0][0]] \
                            ,[y+0.01*el[1][1],y+2*0.01*el[1][1]+0.01*el[0][1]],color="blue", lw=1))

            except:
                print "Error attempting to plot empty edge data..."
            # print entry

# a, b, c are three points, tau is the error margin
def golden_section_search(a, b, c, tau, f):

    if b < c:
        x = b + resphi * (c-b)
    else:
        x = b - resphi * (b-a)

    if np.abs(c-a) < tau * (np.abs(b)+np.abs(x)):
        return (c+a)/2.

    # assert f(x) != f(b)

    if f(x) < f(b):
        if b < c:
            return golden_section_search(b, x, c, tau, f)
        else:
            return golden_section_search(a, x, b, tau, f)

    else:
        if b < c:
            return golden_section_search(a, b, x, tau, f)
        else:
            return golden_section_search(x, b, c, tau, f)


def nelder_mead(x1, x2, tau, f):

    alpha = 1.0
    gamma = 2
    rho = 0.5
    sigma = 0.5

    MAX_ITER = 50
    iter_idx = 0

    while True and iter_idx < MAX_ITER:
        iter_idx = iter_idx + 1

        # Step 1: Order
        if f(x1) > f(x2):
            x1, x2 = x2, x1

        # print x1, x2, f(x1), f(x2)

        if np.abs(f(x2)-f(x1)) < tau and np.abs(x2-x1) < 1e-1:
            return 0.5*(x1+x2)

        # Step 2: Centroid
        x0 = x1

        # Step 3: Reflection
        xr = x0 + alpha*(x0-x2)

        if f(xr) < f(x1):
            xe = x0 + gamma*(xr - x0)
            if f(xe) < f(xr):
                x2 = xe
                continue
            else:
                x2 = xr
                continue

        xc = x0 + rho*(x2-x0)
        if f(xc) < f(x2):
            x2 = xc
            continue

        x2 = x1 + sigma*(x2 - x1)

    print "[WARN]: max iter reached in nelder mead"
    return 0.5*(x1+x2)
