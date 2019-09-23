# calculate moments of given molecule
import math
import numpy as np
import types
import env.GlobalSetting

#_coorID_ = np.array([['x','1',0],['y','2',1],['z','3',2]],dtype='object')
_coorID_ = env.GlobalSetting._coorID_

def toStr(n,base):
    '''
    convert a number to a string
    :param n: the number
    :param base: signals used to represent the number
    :return: string
    '''
    convertString = "0123456789ABCDEF"
    if n < base:
       return convertString[n]
    else:
       return toStr(n//base,base) + convertString[n%base]

def momentCalc(position,index,weight = []):
    '''
    calculate specific moment of given index
    :param position: molecule positions
    :param index: index of calculated moment
    :param weight: weight of each atoms
    :return:
    '''
    data = [0]
    if type(index) == type(1.0) or type(index) == type(1):
        index = str(int(index))
    if type(index) != type(''):
        raise ValueError('Moment index parameter wrong')
    if weight == []:
        weight = np.ones(len(position))
    if len(position)!=len(weight):
        raise ValueError('Wrong Weight')
    if len(index) == 0 or index =='0':
        return sum(weight)
    temp = 'data[0] = ( sum([(w '
    for e in index:
        if e not in _coorID_.reshape((1,_coorID_.size)):
            raise ValueError('Moment index parameter wrong')
        i,j= np.where(_coorID_ == e)
        temp+= '* p['+str(_coorID_[i[0],2])+'] '
    temp+=') for w, p in zip(weight, position)]))'
    exec(temp)
    return data[0]/(sum(weight))

def momOrder(positions,mass,order):
    '''
    calculate moments of one given order
    :param position: molecule positions
    :param mass: masses of atoms
    :param order: order of request moment
    :return: results
    '''
    num = int(3**order)
    base = int(int(toStr((num-1),3))/2)
    moment = {}
    for i in range(num):
        temp = list(str(int(toStr(i,3).zfill(order))+base))
        temp.sort(reverse = True)
        temp = sum([float(temp[r])*10**r for r in range(len(temp))])
        if temp not in moment:
            moment[int(temp)] = momentCalc(positions,temp,weight=mass)
    return moment

def gmomOrder(positions,order):
    '''
    geometry moments
    :param positions: atom positions
    :param order: inquiry order of moments
    :return: moments results
    '''
    return momOrder(positions,np.ones(len(positions)),order)

def mom(positions,masses,orders):
    '''
    calculate moments of demanded orders
    :param position: molecule positions
    :param mass: masses of atoms
    :param order: orders of request moments
    :return: results
    '''
    moment = {}
    for e in orders:
        temp = momOrder(positions,masses,e)
        for ele in temp:
            moment[ele]=temp[ele]
    return moment

def gmom(positions,orders):
    '''
    calculate geometry moments of demanded orders
    :param position: molecule positions
    :param order: orders of request moments
    :return: results
    '''
    return mom(positions,np.ones(len(positions)),orders)

def momAll(positions,masses,upto):
    '''
    calculate moments of demanded orders
    :param position: molecule positions
    :param masses: masses of atoms
    :param upto: demanded orders of moments will be [0, upto]
    :return: results
    '''
    return mom(positions,masses,range(upto+1))

def gmomAll(positions,upto):
    '''
    calculate geometry moments of demanded orders
    :param position: molecule positions
    :param upto: demanded orders of moments will be [0, upto]
    :return: results
    '''
    return gmom(positions,range(upto+1))

def returnMom(moments,ifIndex):
    '''
    return moments with requested format
    :param moments: moments
    :param ifIndex: if return moments with their index
    :return:
    '''
    if ifIndex:
        return [moments[e] for e in moments]
    else:
        return moments

def moment(positions, weight = [], momType = 'gmom',upto=-1,order=-1,orders = [],
                 index = '-1', ifIndex = False):
    '''

    :param positions:
    :param weight:
    :param momType: geometry or physical (using mass as weight)
    :param upto:
    :param order:
    :param orders:
    :param index:
    :param ifIndex:
    :return:
    '''

    moments = {}
    if index != '-1':
        if type(index) == type([]):
            for e in index:
                moments[e] = momentCalc(positions, e, weight=weight)
        else:
            moments[index] = momentCalc(positions, index, weight=weight)
        result = returnMom(moments, ifIndex)
        return result
    if order != -1:
        orders = [order]
    if upto != -1:
        orders = range(upto + 1)
    if momType == 'gmom':
        moments = gmom(positions, orders)
    else:
        moments = mom(positions, weight, orders)

    result = returnMom(moments,ifIndex)
    return result
