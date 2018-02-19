from helpers import *
from matplotlib.lines import Line2D
from mesh import Mesh
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import ode
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np
import os
import pickle

myMesh = pickle.load(open("meshlin.pkl","rb"))
myMesh.churn_edge_data()
