# Pull tht factors representing the force-strain curve showing molecules behavior while being deformed

import types
import math
import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt
import scipy.interpolate as inter
import sys
import os
import numpy as np
import env.GlobalSetting


#import parameters in the calculation
__delta__ = env.GlobalSetting._fdDelta_
__points__ = env.GlobalSetting._fdPoints_
__plotting__ = env.GlobalSetting._fdPlot_

jumpDepthThreshold = env.GlobalSetting._jumpDepthThreshold_
jumpWidthThreshold = env.GlobalSetting._jumpWidthThreshold_
jumpSlopeThreshold = env.GlobalSetting._jumpSlopeThreshold_
aveFracLimit = env.GlobalSetting._aveFracLimit_
aveThreshold = env.GlobalSetting._aveThreshold_


def derivation(func, values, delta):
    '''
    calculate derivation of one function with given list of x values
    :param func: the function to be derivated
    :param values: the x points
    :param delta: delta used to calculate the derivation
    :return: results showing the slope of func at x = values
    '''
    return (func(values+delta)-func(values))/delta

def findMin(_from, _to, _function_,points = __points__):
    '''
    find points with minimum value in a given range in given function.
    :param _from: the start of considering range
    :param _to: end of the range
    :param _function_: aiming function
    :param points: Difference between points
    :return: the x and y of the minimum point
    '''
    step = (_to - _from)/points
    knots = [x*step + _from for x in range(points+1)]
    minNo=_function_(knots[0])
    minXNo=knots[0]
    yvalue = []
    for e in knots:
        yvalue.append(_function_(e))
        if yvalue[-1] < minNo:
            minNo = yvalue[-1]
            minXNo = e
    return minXNo,minNo

def curveMatch(x, y):
    """
    This function is designed to find the best s value for spline
    :param x:
    :param y:
    :return: the min and max of best s range
            if the spline cannot be smooth enough even with s = 5000(maxS), minS will be returned as False
    """
    threshold = 50
    def findMatchRange(x, y, step, start, end):
        trailList = np.linspace(start, end, (end-start)/step + 1)
        _from = min(x)
        _to = max(x)
        testPointsX = np.array([e * (_to - _from) / 500 + _from for e in range(500 + 1)])
        for trail in trailList:
            tempFunc = inter.UnivariateSpline (x, y, s = trail)
            if max(abs(derivation(tempFunc, testPointsX, 0.0001))) < threshold:
                if trail == 0:
                    return trail, trail+step
                else:
                    return trail-step, trail
        if threshold == 50:
            return 'False', trail
        return trail-step, trail
    minS = 200
    maxS = 5000
    while maxS-minS > 5:
        delta = (maxS-minS)/32.0
        minS, maxS = findMatchRange(x, y, delta, minS, maxS)
        if minS == 'False':
            return minS, maxS

    while minS < 1000 and threshold > 20:
        threshold = threshold / 1.2
        minS = 200
        maxS = 2000
        while maxS - minS > 5:
            delta = (maxS - minS) / 32.0
            minS, maxS = findMatchRange(x, y, delta, minS, maxS)
    while minS < 500 and threshold > 10:
        threshold = threshold / 1.1
        minS = 200
        maxS = 1000
        while maxS - minS > 5:
            delta = (maxS - minS) / 32.0
            minS, maxS = findMatchRange(x, y, delta, minS, maxS)
    return minS, maxS

def initCorrection(x, y):
    '''
    calculate the needed shift of the curve to make sure lowest energy occurs with no strain.
    (calculate the needed shift having minimum at x = 0)
    :param x: list of x values of the curve
    :param y: list of y values of the curve
    :return: x and y list after shifted and the shifts being applied
    '''
    coeffs = poly.polyfit(x,y,5)
    pfit = poly.Polynomial(coeffs)
    eo, PEo = findMin(min(x),max(x), pfit)
    x -= eo
    y -= PEo
    return x, y, eo, PEo

def plot(mcn, mtn, stiffC, stiffT, funcC,funcT,lab,plotting = __plotting__,dir = os.getcwd()):
    '''
    Plotting stiffness-strain curve
    :param mcn: Force-Strain data of compression process
    :param mtn: Force-Strain data of tension process
    :param stiffC: Stiffness-Strain data of compression process
    :param stiffT: Stiffness-Strain data of tension process
    :param funcC: function of smoothed Force-Strain curve of compression process
    :param funcT: function of smoothed Force-Strain curve of tension process
    :param lab: label of the deformed molecule
    :param plotting: whether the plots should be plotted
    :param dir: the directory to save the plots
    :return: N/A
    '''
    if plotting == 'yes' or plotting == 'Yes':
        fig = plt.figure(figsize=(8.5, 12))
        ax1 = fig.add_subplot(311)
        ax1.plot(mcn[:, 0], mcn[:, 2], label='Force - Compression', linewidth=0.8)
        ax1.plot(mtn[:, 0], mtn[:, 2], label='Force - Tension', linewidth=0.8)
        ax1.plot(mcn[:, 0], funcC(mcn[:, 0]), label='Force - Compression - smooth', linewidth=0.4)
        ax1.plot(mtn[:, 0], funcT(mtn[:, 0]), label='Force - Tension - smooth', linewidth=0.4)
        ax1.legend()
        ax1.set_xlabel('Strain (%)', family='Helvetica', fontsize=10, color='black')
        ax1.set_ylabel('Force (Kcal/mol/A)', family='Helvetica', fontsize=10, color='black')
        ax21 = fig.add_subplot(312)
        ax21.plot(mcn[:, 0], mcn[:, -3], label='Potential Energy - Compression', linewidth=0.5)
        ax21.plot(mtn[:, 0], mtn[:, -3], label='Potential Energy - Tension', linewidth=0.5)
        ax21.legend(loc=1)
        ax21.set_ylabel('Potential Energy (Kcal/mol)', family='Helvetica', fontsize=10, color='black')
        ax22 = ax21.twinx()
        ax22.plot(mcn[:, 0], mcn[:, -2], 'r', label='Total Energy - Compression', linewidth=0.5)
        ax22.plot(mtn[:, 0], mtn[:, -2], 'g', label='Total Energy - Tension', linewidth=0.5)
        ax22.legend(loc=2)
        ax22.set_ylabel('Total Energy (Kcal/mol/A)', family='Helvetica', fontsize=10, color='black')
        ax22.set_title('')
        ax31 = fig.add_subplot(313)
        ax31.plot(stiffC[:, 0], stiffC[:, 1], label='Stiffness - Compression', linewidth=0.5)
        ax31.plot(stiffT[:, 0], stiffT[:, 1], label='Stiffness - Tension', linewidth=0.5)
        ax31.legend()
        ax31.set_ylabel('Flexibility (Kcal/mol/A/%)', family='Helvetica', fontsize=10, color='black')
        ax31.set_title('')
        if dir[-1]!='/':
            dir.append('/')
        plt.savefig(dir + lab + "-FD-TE-curve.jpeg", transparent=True)

def stiffness(_from, _to, forceFunc, points = __points__):
    '''
    Calculate stiffness changing in the process
    :param _from: start strain of the deformation
    :param _to: end strain of the deformation
    :param forceFunc: function of smoothed force-strain data
    :param points: difference between points showing stiffness changing tendency
    :return: a list shows stiffness-strain data ([[strain1, stiffness @ strain1], ...])
    '''
    stiff = []
    step = (_to - _from) / points
    xPoints = [x * step + _from for x in range(points + 1)]
    deltaX = __delta__
    for e in xPoints:
        stiff.append(np.array([e, (forceFunc(e + deltaX) - forceFunc(e)) / deltaX * ((100 + e) / 100.0)]))
    stiff = np.array(stiff)
    return stiff

def jump(strain, stiffness, forceFunc, _depth_ = jumpDepthThreshold, _slope_ = jumpSlopeThreshold):
    '''
    Identify extraordinary force changes.
    :param strain: an array of strain values
    :param stiffness: stiffness value with strain in the strain array
    :param forceFunc: function of smoothed force-strain data
    :param _depth_: the least depth value as threshold to identify if jump phenomena occurs
    :param _slope_: the least slope value as threshold to identify if jump phenomena occurs
    :return: a list of jumps with the format [[starting strain, width of the jump, and depth of jump], ...]
    '''
    _start = [float('-inf'),0,0]
    _end = [float('-inf'),0,0]
    jumps = []
    for x,y in zip(strain,stiffness):
        if y < -1:   #The Threshold was set randomly, in need of more trials.
            if _start[0]==float('-inf'):
                _start = [x,y,forceFunc(x)]
            else:
                _end = [x,y,forceFunc(x)]
        else:
            if _start[0]!=float('-inf') and _end[0]!=float('-inf'):
                if abs(_start[2]-_end[2])>_depth_ and abs(_start[2]-_end[2])/abs(_start[0]-_end[0])>_slope_:
                    if len(jumps)>1 and _start[0] - (jumps[-1][0] + jumps[-1][1]) < 2:
                        if forceFunc(jumps[-1][0] + jumps[-1][0]) - _end[2] > 0:
                            jumps.append(np.array([jumps[-1][0], abs(jumps[-1][0] - _end[0]), abs(forceFunc(jumps[-1][0]) - _end[2])]))
                    else:
                        jumps.append(np.array([_start[0],abs(_start[0]-_end[0]),abs(_start[2]-_end[2])]))
            _start = [float('-inf'),0,0]
            _end = [float('-inf'),0,0]
    jumps = np.array(jumps)
    return jumps

def ave(stiff, noise, threshold = aveThreshold):
    '''
    Calculate average stiffness with a given curve from a starting points before y-value changes significantly.
    :param stiff: a list of stiffness
    :param noise: shows how noisy the curve was before smoothed
    :param threshold: threshold to identify where y-value changes significantly
    :return: average stiffness
    '''

    count = 1
    aveStiff = stiff[threshold]
    r = [threshold]
    for i in range(threshold,len(stiff) - 1):
        if i<200 and (count < threshold or abs(stiff[i + 1] - (aveStiff/count)) < aveFracLimit/noise):    #The Threshold was set randomly, in need of more trials.
            aveStiff += stiff[i + 1]
            r.append(i + 1)
            count += 1
    aveStiff = aveStiff / count

    i=0
    while abs(aveStiff- stiff[threshold+i]) > aveFracLimit/noise and i < count - 20:
        aveStiff = (aveStiff * (count - i) - stiff[threshold+i])/(count-1-i)
        i = i + 1
        r.pop(0)
    x = 0
    y = 0
    for e in r:
        x += stiff[e]
    x = x/len(r)

    return aveStiff, x

def calcCurveFeatures(mc, mt, plotting = __plotting__, lab='',dir = os.getcwd()):
    '''
    calculate the factors being used to represent force-strain curve
    :param mc: Force-Strain data of compression process
    :param mt: Force-Strain data of tension process
    :param plotting: whether plot the stiffness curve figures
    :param lab: molecule label
    :param dir: directory of files
    :return: the representing factors
    '''

    x ,y, eo, PEo = initCorrection(mc[:,0],mc[:,-2])

    mc[:,0] -= eo
    mt[:,0] -= eo
    mc[:,-2] -= PEo
    mt[:,-2] -= PEo

    minSC, maxS = curveMatch(mc[:, 0], mc[:, 2])
    fdFunctionC = inter.UnivariateSpline(mc[:, 0], mc[:, 2], s = maxS)
    fdNoiseC = maxS
    minST, maxS = curveMatch(mt[:, 0], mt[:, 2])
    fdFunctionT = inter.UnivariateSpline(mt[:, 0], mt[:, 2], s = maxS)
    fdNoiseT = maxS


    # Calculate the flexibility of the molecule with 'flexibility = PresentLength*Delta(Force)/Delta(Length)'
    # plot the flexibility information and save the data
    stiffC = stiffness(min(mc[:,0]), max(mc[:,0]), fdFunctionC)
    stiffT = stiffness(min(mt[:,0]), max(mt[:,0]), fdFunctionT)
    jumpC = jump(stiffC[:,0],stiffC[:,1], fdFunctionC)
    jumpT = jump(stiffT[:,0],stiffT[:,1], fdFunctionT)
    # Get the differences of tension and compression process and evaluate the jump abnormal situation
    _diff_ = max([ abs(fdFunctionC(e) - fdFunctionT(e)) for e in mc[:,0]])

    # Calculate the average of flexibility before and after bulk

    _aveCF_,refxCF = ave(stiffC[:,1],fdNoiseC)
    _aveCB_,refxCB = ave(stiffC[:,1][::-1],fdNoiseT)
    _aveTF_,refxTF = ave(stiffT[:,1],fdNoiseC)
    _aveTB_,refxTB = ave(stiffT[:,1][::-1],fdNoiseT)
    intersectionC = (fdFunctionC(refxCF)-fdFunctionC(refxCB)+_aveCB_*refxCB-_aveCF_*refxCF)/(_aveCB_-_aveCF_)
    intersectionT = (fdFunctionT(refxTF)-fdFunctionT(refxTB)+_aveTB_*refxTB-_aveTF_*refxTF)/(_aveTB_-_aveTF_)

    if plotting:
        plot(mc, mt, stiffC, stiffT, fdFunctionC, fdFunctionT, lab, plotting=__plotting__, dir = dir)

    return {"jumpC": jumpC, "jumpT": jumpT,
            "diff": _diff_, "aveCB": _aveCB_,
            "aveCF": _aveCF_, "aveTB": _aveTB_,
            "aveTF": _aveTF_,"intersectionC": intersectionC,"intersectionT": intersectionT}

def preAnalysis(oriData):
    '''
    preperation of data from '-ave-force.d' file
    :param oriData: data from '-ave-force.d' file
    :return: prepared data
    '''
    l = len(oriData)
    nav = 5
    giveUpNo =l-l//nav*nav
    if giveUpNo != 0:
        oriData = oriData[:-(giveUpNo)]
    len0 = oriData[0][2]
    dataReorg = []
    for e in oriData:
        f = [0.5 * (e[3] - e[6]), 0.5 * (e[4] - e[7]), 0.5 * (e[5] - e[8])]
        v = e[15:18]
        mod = (v[0]**2+v[1]**2+v[2]**2)**0.5
        d = [x / mod for x in v]
        fax = f[0]*d[0]+f[1]*d[1]+f[2]*d[2]
        fn = [f[0]-[x * fax for x in d][0],f[1]-[x * fax for x in d][1],f[2]-[x * fax for x in d][2]]
        mod = (fn[0]**2+fn[1]**2+fn[2]**2)**0.5
        dataReorg.append([100*(e[2]-len0)/len0,e[2]-len0,fax,mod,e[18],e[19],e[20]])
    data = []
    for i in range(int(len(dataReorg))):
        if i/nav == int(i/nav):
            data.append([])
            for j in range(len(dataReorg[0])):
                data[i//nav].append(0)
                if i//nav != 0:
                    data[i // nav - 1][j] = data[i // nav - 1][j]/nav
        for j in range(len(dataReorg[0])):
            data[i//nav][j]+=dataReorg[i][j]

    for j in range(len(dataReorg[0])):
        data[-1][j] = data[-1][j]/nav
    return data

def splitCT(data):
    '''
    seperate data to compression process and tnesion process
    :param data: raw data
    :return: splitted data
    '''
    mc = []
    mt = []
    for i in range(len(data)):
        if i == 0:
            mc.append(data[i])
            continue
        if (data[i-1][1] - data[i][1]) < 0:
            mt.append(data[i])
        else:
            mc.append(data[i])
    mc = sorted(mc, key=lambda x : x[0])
    mt = sorted(mt, key=lambda x : x[0])
    return mc, mt

def trim(data):
    '''
    Combine duplicated data (with same strain)
    :param data: raw data
    :return: trimmed data
    '''
    i = 0
    while i in range(len(data) - 1):
        count = 1.0
        newData = data[i]
        while i < (len(data) - 1) and data[i][0] == data[i + 1][0]:
            count = count + 1
            for j in range(len(data[i])):
                newData[j] = newData[j] + data[i + 1][j]
            data.pop(i)
        for j in range(len(data[i])):
            newData[j] = newData[j] / count
        data[i] = newData
        i = i + 1
    return data

def calcFileFeatures(res, plotting = __plotting__, lab='',dir = os.getcwd()):
    '''
    Calculate factors describing the carve data from '-ave-force.d' file.
    :param res: The opened file
    :param plotting: whether plot the results
    :param lab: label of given molecule
    :param dir: directory of used files
    :return: data describing factors
    '''
    # Read file in a list
    oriData = []
    for line in res:
        words = line.split()
        if len(words)==0:
            continue
        if words[0]!='#':
            tempt = []
            for word in words:
                tempt.append(float(word))
            oriData.append(tempt)

    data = preAnalysis(oriData)
    mc, mt = splitCT(data)

    #Get rid of repeated xvalue:
    mc = trim(mc)
    mt = trim(mt)

    return calcCurveFeatures(np.array(mc),np.array(mt), plotting = plotting, lab=lab, dir = dir)

class fd:
    def __init__(self,filename,dir = os.getcwd(), lab = '', plotting = __plotting__):
        if lab == '':
            lab = filename
        if dir[-1]!='/':
            dir.append('/')
        file = open(dir+filename,'r')
        self.factors = calcFileFeatures(file, plotting = plotting, lab=lab, dir = dir)

def test():
    test_1 = fd('-ave-force.d')
    print(test_1.factors)


if __name__ == '__main__':
    test()