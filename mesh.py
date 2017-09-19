
from helper import triangulateQuad

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
        triangles, newCorners = triangulateQuad(corners)
        self.vertices = self.vertices + newCorners
        self.triangles = self.triangles + triangles
