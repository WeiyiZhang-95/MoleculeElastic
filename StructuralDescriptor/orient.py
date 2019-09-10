import math
import numpy as np
import sys

def rotation(position,matix):
    if type(position[0]) == type(np.array([])):
        return np.array(np.mat(matix)*position.T).T
    else:
        return np.array(np.mat(matix)*position.T).T.flatten()

def trans(position,transVector):
    return position+transVector

def massCenter(positions,weight = []):
    if weight == []:
        weight = np.ones(len(positions))
    massCenter = [sum(x) for x in np.matrix.transpose(np.array([i * j for i, j in zip(np.array(weight), np.array(positions))]))] / sum(np.array(weight))
    return massCenter

def carbox(node1, node2):
    nodes2 = np.array(node2)
    nodes1 = np.array(node1)
    centers = np.array([sum(node1)/2,sum(node2)/2])
    x = centers[0]-centers[1]
    x = x / np.dot(x,x)**0.5
    directions = np.array([nodes1[0]-nodes1[1],nodes2[0]-nodes2[1]])
    if np.dot(directions[0],directions[1]) >= 0:
        y = directions[0] + directions[1]
    else:
        y = directions[0] - directions[1]
    z = np.cross(x,y)
    z = z / np.dot(z,z)**0.5
    y = np.cross(z,x)
    transMat = np.mat(np.matrix([x,y,z]))
    return transMat

def centering(positions,weight = []):
    if weight == []:
        weight = np.ones(len(positions))
    center = massCenter(positions,weight = weight)
    positions = trans(positions,-center)
    return positions

def orienting(positions):
    node1 = np.array([positions[0],positions[2]])
    node2 = np.array(positions[-2:])
    transM = carbox(node1, node2)
    positions = rotation(positions, transM)
    return positions