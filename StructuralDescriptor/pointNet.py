import math
import numpy as np
import types
import env.GlobalSetting

elements = [e for e in env.GlobalSetting._ELEMENTS_]
def pointNet(positions, types):
    points = []
    positions = positions
    for i in range(len(types)):
        type = list(np.zeros(len(elements)+1))
        type[elements.index(types[i])] = 1
        points.append(list(positions[i])+type)
    points = np.array(points)
    return points