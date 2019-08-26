# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 10:45:56 2019

@author: Noel
"""
import numpy as np
import random

def f(r, sigma):
    coef = 1
    #coef = 1/np.sqrt(2*np.pi*sigma**2)
    exp = -r**2/(2*sigma**2)
    return coef * np.exp(exp)



Xs = []
Ys = []
    
def makeImage(sigma):
    image = []
    for i in range(240):
        image.append([0]*240)
    for i in range(240):
        for j in range(240):
            x = i - 120.
            y = j - 120.
            r = np.sqrt(x**2 + y **2)
            image[i][j] = [f(r, sigma)]*4
#get shape right, filters            
            
    return image

def main():   
    global Xs
    global Ys         
    for i in range(100):        
        if random.random() < 0.5:
            Xs.append(makeImage(5.))
            Ys.append(0)
        else:
            Xs.append(makeImage(50.))
            Ys.append(1)
    Xs = np.array(Xs)
    print(Xs.shape)
    Ys = np.array(Ys)
    print(Ys.shape)
    np.save('test_gaussian_xs', Xs)
    np.save('test_gaussian_ys', Ys)
    
if __name__ == '__main__':
    main()
    