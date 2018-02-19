#!/usr/bin/env python

from mesh import *
x = Mesh([(0,0),(1,0),(1,1),(0,1)],(0,0))
x.submesh(3)
x.write_to_poly(k=0.001)

