import numpy as np
import os

def addArray(arr, data, toEnd = True):
    '''
    add new data points to an array
    :param arr: old array
    :param data: new data
    :param toEnd: with order or to the end
    :return: new array
    '''
    if len(arr) == 0:
        return np.array([data])
    l = np.zeros(len(data))
    if len(data) != len(arr[0,:]):
        ValueError("Importing element format not match")
        return False
    if toEnd:
        p = len(data.flatten())
    else:
        p = max(np.where(arr[:,0]<data[0])[0]) * len(data)
    arr = np.insert(arr, p+l , data)
    arr = arr.reshape((int(len(arr)/len(data)),len(data)))
    return arr

def addNpy(file, data, toEnd = True, dir = os.getcwd()):
    '''
    add new data tp .npy file
    :param file: file name
    :param data: new data
    :param toEnd: with order or to the end
    :param dir: directory
    :return:
    '''
    if file in os.listdir(dir):
        f = np.load(file)
    else:
        f = np.array([])
    f = addArray(f,data, toEnd=toEnd)
    np.save(file,f)
    return

def writeFile(file, data, toEnd = True, dir = os.getcwd()):
    '''
    write new coming data to data files
    :param file: file name
    :param data: new data
    :param toEnd: with order or to the end
    :param dir: directory
    :return:
    '''
    line = '\t'.join([str(e) for e in data])
    line = line+'\n'
    f = open(file,'a')
    f = open(file,'r+')
    l = f.readlines()
    if toEnd:
        p = len(l)-1
    else:
        ids = np.array([float(e.split()[0]) for e in l])
        if not len(np.where(ids < float(data[0]))[0]):
            p = -1
        else:
            p = max(np.where(ids < float(data[0]))[0])
    l.insert(p+1,line)
    f.seek(0)
    f.writelines(l)
