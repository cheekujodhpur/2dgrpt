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

        g = np.zeros((2,2))
        g[0,0] = 1
        v = p-o
        g[1,1] = (np.sin(v[0]*np.pi)**2)
        
        return g
        
