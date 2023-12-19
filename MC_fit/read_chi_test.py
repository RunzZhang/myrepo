import matplotlib.pyplot as plt
import sys
import numpy as np
import math
import scipy.optimize
import random

try:
    fileprefix=sys.argv[1]
except:
    fileprefix="./Test_Dump/test2"