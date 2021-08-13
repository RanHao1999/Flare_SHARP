# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 13:12:50 2021

@author: 11730
"""
import numpy as np
from datetime import  datetime
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import warnings
from sklearn.datasets import make_moons
from sklearn import preprocessing
#functions that will be repeatly used in the whole program.
def string2float(a):
    lst = a.split(' ')
    while '' in lst:
        lst.remove('')
    for i in range(len(lst)):
        lst[i] = float(lst[i])
    return lst


def dataframe2matrix(a):  # turn dataframe, which is directly read from .txt, to matrix
    a = a.values
    res = []
    for i in range(len(a)):
        res.append(string2float(a[i,0]))
    return np.matrix(res)

def lstabs(lst): #get all the values' abs in a list
    for i in range(len(lst)):
        if lst[i] < 0:
            lst[i] = abs(lst[i])
    return lst

def linear_regression(train_x,train_y,title):
    #train_x and train_y are lists
    warnings.filterwarnings('ignore')
    train_x = np.array(train_x)
    train_y = np.array(train_y)
    
    times = 1000      #times for iteration
    lrate = 0.01      #learning rate, small value 
    w0,w1 = [1],[1]   #original model parameters, record every gradient drops
    losses = []       #record the loss of the value every iteration
    epoches = []      #record index for every iteration
    
    for i in range(times):
        # output allw0,w1 and loss in every drop
        epoches.append(i-1)
        loss = ((w0[-1]+w1[-1]*train_x-train_y)**2).sum()/2.0
        losses.append(loss)
#        print('{:4}> w0={:.6f},w1={:.6f},loss={:.6f}'.format(epoches[-1], w0[-1], w1[-1], losses[-1]))
        
        #correlation of w0 and w1, d0 and d1 are partial derivatives for loss in w0 and w1 direction
        d0 = (w0[-1] + w1[-1]*train_x -train_y).sum()
        d1 = ((w0[-1]+w1[-1]*train_x-train_y)*train_x).sum()
        w0.append(w0[-1]-lrate*d0)
        w1.append(w1[-1]-lrate*d1)
    pred_y = w0[-1] + w1[-1]*train_x   #regression
    plt.figure('Linear Regression',facecolor = 'lightgray')
    plt.title(title)
    plt.grid(linestyle = ':')
    plt.scatter(train_x,train_y,s=30,color='blue',label='Points')
    plt.plot(train_x,pred_y,color='orange',label = 'Regression Line')
    plt.legend()
    plt.show()
        
    return pred_y
        
def normalization(data):
    # for regression, we normalize the data with max_min normalization method
    # data should be numpy.ndarray
    res = np.ones(len(data))
    d_max = np.max(data)
    d_min = np.min(data)
    for i in range(len(data)):
        res[i] = float(data[i] - d_min) / float(d_max - d_min)
    return res


def average(list):
    return sum(list)/len(list)

def dict2list(dict):
    """
    Turn the dicts from .csv files to list for convenience. Normally, 'dict' is the 'row' in reader.
    :param dict: dict data directly read from .csv files
    :return:a list containing all the dict.values() in the form of float
    """
    values = dict.values()
    try:
        list = [float(x) for x in values]
    except:
        list = [x for x in values]
    return list

def list2date0(tlist):
    # Turn a list of strings which mean date into a datetime object
    for i in range(len(tlist)):
        tlist[i] = datetime.strptime(tlist[i], '%Y-%m-%dT%H:%M:%S.%f')

    return tlist

def list2date1(tlist):
    # Turn a list of strings which mean date into a datetime object
    for i in range(len(tlist)):
        tlist[i] = datetime.strptime(tlist[i], '%Y/%m/%d %H:%M')

    return tlist

def SHARP2float(SHARP_parameters):
    for i in range(len(SHARP_parameters)):
        for j in range(len(SHARP_parameters[0])):
            if SHARP_parameters[i][j] == '':
                SHARP_parameters[i][j] = (float(SHARP_parameters[i][j-1]) + float(SHARP_parameters[i][j+1]))/2.0
            else:
                SHARP_parameters[i][j] = float(SHARP_parameters[i][j])
    return SHARP_parameters


