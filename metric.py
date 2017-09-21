import numpy as np

class Metric:
    """
    Metric properties and calculation
    polar metric special
    """
    
    origin = (0,0)

    def __init__(self, origin):
       self.origin = origin

    def compute_metric(self, position):
        o = np.array(self.origin)
        p = np.array(position)

        g = np.zeros((3,3))
        g[0,0] = 1.
        g[1,1] = np.sum((p-o)**2)
        v = p-o
        g[2,2] = g[1,1]*float(v[1]*v[1])/(v[0]*v[0]+v[1]*v[1])
        
        return g
        
