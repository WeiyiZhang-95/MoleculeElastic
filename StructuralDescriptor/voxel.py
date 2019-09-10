import env.GlobalSetting
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import numpy as np
from sklearn.decomposition import PCA
import scipy.sparse

#VOXESCALE = 0.2
VOXESCALE = env.GlobalSetting.VOXESCALE
#THETA = 0.6
THETA = env.GlobalSetting.THETA
#CUTOFF = 3
CUTOFF = env.GlobalSetting.CUTOFF

__plotting__ = env.GlobalSetting.voxPlot
__channel__ = env.GlobalSetting.voxChannel

#elements = {'id': class element = (id,index,mass,radii), ...}
elements = env.GlobalSetting._ELEMENTS_
elements['unidentified'] = env.GlobalSetting._element_(id='unidentified',index=0,mass=0,radius=1)
for e in elements:
    elements[e].theta = THETA * elements[e].radius
    elements[e].omega = 1.0 / elements[e].theta
    elements[e].cutoff = math.ceil(CUTOFF * elements[e].theta / VOXESCALE)
systemCutoff = max([elements[e].cutoff for e in elements])

color = [[0, 0, 1], [1, 0.5, 0], [0, 1, 0], [1, 0, 0]]
color_scatter = ["blue", "orange", "green", "red"]
# H: blue, C: orange N: green O: red


cmaps1 = ["PiYG", "PRGn", "PuOr", "RdGy"]
cmaps2 = ['Purples','Blues','Greens','Reds']

def calcVoxels(deltaX, deltaY, deltaZ, ele):
    if deltaX == 0 and deltaY == 0 and deltaZ == 0 and elements[ele].theta == 0:
        return 1
    if elements[ele].theta == 0:
        return 0
    if (deltaX ** 2 + deltaY ** 2 + deltaZ ** 2) > elements[ele].cutoff ** 2:
        return 0

    return math.exp(-(deltaX ** 2 + deltaY ** 2 + deltaZ ** 2) * (VOXESCALE ** 2) / (2 * (elements[ele].theta ** 2)))

def waveTransform(deltaX, deltaY, deltaZ, ele):
    return math.cos(2 * math.pi * elements[ele].omega * VOXESCALE * ((deltaX ** 2 + deltaY ** 2 + deltaZ ** 2) ** 0.5))\
           * calcVoxels(deltaX, deltaY, deltaZ, ele)

def channel(position,type,waveTrans = True):
    species = []
    for e in type:
        if e not in species:
            species.append(e)
    size = [0, 0, 0]
    limit = np.array([[0, 0], [0, 0], [0, 0]])
    for i in range(len(size)):
        limit[i][-1] = int(max(position[:, i]) / VOXESCALE) + systemCutoff + 1
        limit[i][0] = int(min(position[:, i]) / VOXESCALE) - systemCutoff - 1
        size[i] = limit[i][-1]-limit[i][0]
    dataMatrix = np.zeros(shape=(len(species), size[0], size[1], size[2]), dtype=np.float32)
    env.GlobalSetting.__statusOutput__("voxel: Meshgrid Size: "  +  str(size))
    env.GlobalSetting.__statusOutput__("voxel: Computing Grid...")

    for p,t in zip(position,type):
        x = int((p[0]) / VOXESCALE) - limit[0][0]
        y = int((p[1]) / VOXESCALE) - limit[1][0]
        z = int((p[2]) / VOXESCALE) - limit[2][0]
        for m in range(x - elements[t].cutoff, x + elements[t].cutoff + 1):
            for n in range(y - elements[t].cutoff, y + elements[t].cutoff + 1):
                for l in range(z - elements[t].cutoff, z + elements[t].cutoff + 1):
                    if waveTrans:
                        dataMatrix[species.index(t)][m][n][l] += waveTransform(x - m, y - n, z - l, t)
                        #print(species.index(t),m,n,l,dataMatrix[species.index(t)][m][n][l])
                    else:
                        dataMatrix[species.index(t)][m][n][l] += calcVoxels(x - m, y - n, z - l, t)
    return dataMatrix, limit

def whole(position, waveTrans = True):
    l = ['unidentified' for i in range(len(position))]

    dataMatrix, limit = channel(position,l,waveTrans = waveTrans)
    return dataMatrix, limit

def plot(dataMatrix, limit,waveTrans = True):
    size = [0, 0, 0]
    for i in range(len(size)):
        size[i] = limit[i][-1]-limit[i][0]
    stat = dataMatrix.reshape(-1)
    # Just for Plotting.
    env.GlobalSetting.__statusOutput__("Plotting")
    fig = plt.figure(figsize=[10,10])
    #fig.suptitle(lab+'-'+ args[0]+'\nCUTOFF='+str(CUTOFF)+'\nGaussian Parameter='+str(THETA)+'\nVoxel Size(Ang)='+str(VOXESCALE))

    '''
    ax = fig.add_subplot(2, 2, 1, projection='3d')
    ax.title.set_text('Voxels')


    colors = np.zeros(dataMatrix[0].shape + (4,))
    #should work for at most 4 element types
    for i in range(len(species)):
        colors[..., 0] = color[i][0]
        colors[..., 1] = color[i][1]
        colors[..., 2] = color[i][2]
        colors[..., 3] = dataMatrix[i]
        ax.voxels(dataMatrix[i], facecolors=colors)
    '''
    ax = fig.add_subplot(2, 2, 3)
    ax.title.set_text('Histogram')
    ax.hist(stat, bins=10)

    ax = fig.add_subplot(2, 2, 4)
    ax.title.set_text('Histogram (0.2~max)')
    ax.hist(stat, bins=100, range=(0.2, stat.max()))

    ax = fig.add_subplot(2, 2, 2)
    ax.title.set_text('xoy Section')
    for i in range(0,len(dataMatrix)):
        if waveTrans:
            cmLimit = max([abs(min(stat)),max(stat)])
            ax.imshow(dataMatrix[i,:,:,(-limit[2][0])], cmap=cmaps1[i],vmin=-cmLimit, vmax=cmLimit, alpha=0.5)
        else:
            ax.imshow(dataMatrix[i,:,:,(-limit[2][0])], cmap=cmaps1[i], alpha=0.5)

    ax = fig.add_subplot(2, 2, 1)
    overAll = np.zeros((size[0], size[1], size[2]), dtype=np.float32)
    for i in range(size[0]):
        for j in range(size[1]):
            for k in range(size[2]):
                overAll[i,j,k] = sum(dataMatrix[:,i,j,k])
    cmLimit = max([abs(min(stat)),max(stat)])
    ax.imshow(sum(dataMatrix[:,:,:,(-limit[2][0])]), cmap=cmaps1[3],vmin=-cmLimit, vmax=cmLimit, alpha=0.5)
    #plt.savefig('CUTOFF=' +str(CUTOFF) + '-' +str(THETA)+ '-' + str(VOXESCALE) + '.jpg')
    return plt

