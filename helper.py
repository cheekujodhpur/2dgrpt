
import math

def triangulateQuad(corners):
    # NOTE: assumes a 2-tuple
    center = tuple(map(lambda x:float(x)/len(corners),reduce(lambda x,y:(x[0]+y[0],x[1]+y[1]),corners)))
    angles = []
    for corner in corners:
        vector = tuple(map(lambda x,y:x-y,corner,center))
        angles.append(math.atan2(vector[0],vector[1]))
        
    sorted_corners = [x for _,x in sorted(zip(angles,corners))][::-1]

    triangles = [(0,1,2),(0,2,3)]

    # Writing triangles is easy, the main purpose was always to sort the corners
    return triangles,sorted_corners

