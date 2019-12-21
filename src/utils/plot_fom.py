#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 10:40:52 2019

@author: nabalone
"""

import numpy as np
import matplotlib.pyplot as plt
pred = np.load('/mnt/c/Users/Noel/Desktop/summer2019/dev/onOdyssey/ia_hazel/y_pred_from__cnn_aardvark_aug_concat.py_b_200_mask_c_mp_5_n_50_a_kfold_0_l_0.0025.npy')
test = np.load('/mnt/c/Users/Noel/Desktop/summer2019/dev/onOdyssey/ia_hazel/testy_pred_from__cnn_aardvark_aug_concat.py_b_200_mask_c_mp_5_n_50_a_kfold_0_l_0.0025.npy')
actual_ias = []
actual_ccs = []

for i in range(len(test)):
    if test[i] == 0:
        actual_ias.append(pred[i][0])
    else:
        actual_ccs.append(pred[i][0])
ntotia = len(actual_ias)
fomx = []
fomy = []
effs = []
pps = []
for i in range(200):
    cutoff = i/100
    
    nsubia = 0
    for n in actual_ias:
        if n > cutoff:
            nsubia += 1
            
    nsubnonia = 0
    for n in actual_ccs:
        if n > cutoff:
            nsubnonia += 1
    if nsubia==0 and nsubnonia==0:
        continue
    eff = nsubia/ntotia
    pp = nsubia/(nsubia + 5*nsubnonia)
    effs.append(eff)
    pps.append(pp)
    fomy.append(eff*pp)
    fomx.append(cutoff)
    
f = open('../dev/foleymandel_fom_plus.txt')
ff = f.readlines()
xs = []
ys = []
exs = []
eys = []
pxs = []
pys = []

for i in ff:
    try:
        px, py, ex, ey, x, y = i.split()
    except:
        break
    xs.append(float(x))
    ys.append(float(y))
    exs.append(float(ex))
    eys.append(float(ey))
    pxs.append(float(px))
    pys.append(float(py))
f.close()
plt.plot(fomx, fomy, color='k', lw=3, label='Chou+')
plt.plot(fomx, effs, color='r', lw=3, label='Chou+ efficiency')
plt.plot(fomx, pps, color='m', lw=3, label='Chou+ psuedopurity')
plt.plot(exs, eys, color='b', lw=3, label='Foley efficiency')
plt.plot(pxs, pys, color='c', lw=3, label='Foley psuedopurity')
plt.plot(xs, ys, color='g', lw=3, label='Foley-Mandel')
plt.plot(np.arange(0,1.1,0.1),[np.max(ys)]*11, linestyle = "--", color='g')
plt.plot(np.arange(0,1.1,0.1),[np.max(fomy)]*11, linestyle = "--", color='k')
plt.axis([0, 1, 0, 1])
plt.ylabel('Figure of Merit')
plt.xlabel('Classification Weighting')
plt.legend(loc=2)
plt.savefig('fom_good_with_effandpp.png')
            
    

a = '''
[[9.75498497e-01 2.45015509e-02]
 [1.87610015e-01 8.12389970e-01]
 [1.00000000e+00 9.20312093e-09]
 [9.68505085e-01 3.14949527e-02]
 [9.80335236e-01 1.96648184e-02]
 [2.94532895e-01 7.05467105e-01]
 [6.11041614e-05 9.99938846e-01]
 [1.00000000e+00 3.86906598e-13]
 [9.98797297e-01 1.20270671e-03]
 [1.00000000e+00 6.30463650e-12]
 [9.98951316e-01 1.04863709e-03]
 [8.87457550e-01 1.12542436e-01]
 [2.47553110e-01 7.52446949e-01]
 [3.90692294e-01 6.09307706e-01]
 [3.61451358e-01 6.38548613e-01]
 [9.99397993e-01 6.02020998e-04]
 [5.73698133e-02 9.42630172e-01]
 [4.17809635e-01 5.82190394e-01]
 [9.92366135e-01 7.63390446e-03]
 [9.01885724e-05 9.99909759e-01]
 [6.11041614e-05 9.99938846e-01]
 [5.73698133e-02 9.42630172e-01]
 [3.75418772e-08 1.00000000e+00]
 [9.80335236e-01 1.96648184e-02]
 [3.90692294e-01 6.09307706e-01]
 [4.94037539e-01 5.05962431e-01]
 [2.83390105e-01 7.16609895e-01]
 [9.99117196e-01 8.82767374e-04]
 [5.44878721e-01 4.55121309e-01]
 [3.61451358e-01 6.38548613e-01]
 [8.43645275e-01 1.56354696e-01]
 [2.69169658e-01 7.30830312e-01]
 [8.77990723e-02 9.12200928e-01]
 [3.14717799e-01 6.85282230e-01]
 [1.00000000e+00 8.16231851e-17]
 [6.77544813e-05 9.99932289e-01]
 [1.00000000e+00 1.96200121e-16]
 [8.87457550e-01 1.12542436e-01]
 [1.97676882e-01 8.02323103e-01]
 [8.43645275e-01 1.56354696e-01]
 [9.60534811e-01 3.94652374e-02]
 [2.25613620e-02 9.77438629e-01]
 [1.56118146e-10 1.00000000e+00]
 [1.00000000e+00 9.20312093e-09]
 [1.19422756e-01 8.80577266e-01]
 [9.97588873e-01 2.41110264e-03]
 [9.71022069e-01 2.89779603e-02]
 [6.79777741e-01 3.20222318e-01]
 [9.18254629e-03 9.90817428e-01]
 [1.00000000e+00 1.96200121e-16]
 [8.56707265e-06 9.99991417e-01]
 [1.00000000e+00 1.96200121e-16]
 [9.99397993e-01 6.02020998e-04]
 [9.98484433e-01 1.51553657e-03]
 [3.26889098e-01 6.73110962e-01]
 [9.73966420e-01 2.60336511e-02]
 [3.25903669e-02 9.67409611e-01]
 [8.43645275e-01 1.56354696e-01]
 [2.79215574e-01 7.20784485e-01]
 [9.80335236e-01 1.96648184e-02]
 [1.00000000e+00 3.86906598e-13]
 [4.94037539e-01 5.05962431e-01]
 [8.43645275e-01 1.56354696e-01]
 [7.58696616e-01 2.41303340e-01]
 [6.11041614e-05 9.99938846e-01]
 [3.92647879e-03 9.96073484e-01]
 [7.58696616e-01 2.41303340e-01]
 [9.18254629e-03 9.90817428e-01]
 [4.60357249e-01 5.39642811e-01]
 [4.10838068e-01 5.89161932e-01]
 [3.78712284e-10 1.00000000e+00]
 [2.33722001e-01 7.66278028e-01]
 [9.18254629e-03 9.90817428e-01]
 [5.47179130e-26 1.00000000e+00]
 [2.25613620e-02 9.77438629e-01]
 [3.58529687e-01 6.41470253e-01]
 [1.35919571e-01 8.64080429e-01]
 [1.00000000e+00 3.86906598e-13]
 [8.87457550e-01 1.12542436e-01]
 [8.87457550e-01 1.12542436e-01]
 [9.99397993e-01 6.02020998e-04]
 [7.18386114e-01 2.81613916e-01]
 [6.44599438e-01 3.55400622e-01]
 [6.79777741e-01 3.20222318e-01]
 [5.45264930e-02 9.45473552e-01]
 [1.00000000e+00 6.30463650e-12]
 [1.00000000e+00 2.57019124e-13]
 [2.69169658e-01 7.30830312e-01]
 [1.00000000e+00 2.58657723e-17]
 [1.54152215e-01 8.45847845e-01]
 [8.77990723e-02 9.12200928e-01]
 [1.94440067e-01 8.05559933e-01]
 [1.00000000e+00 6.89141277e-09]
 [1.00000000e+00 4.27379174e-11]
 [9.99801576e-01 1.98485286e-04]
 [9.68505085e-01 3.14949527e-02]
 [1.00000000e+00 3.86906598e-13]
 [2.98221022e-01 7.01778948e-01]
 [8.77990723e-02 9.12200928e-01]
 [4.11270261e-01 5.88729799e-01]
 [7.28696525e-01 2.71303475e-01]
 [1.56118146e-10 1.00000000e+00]
 [1.00000000e+00 9.20312093e-09]
 [5.09423852e-01 4.90576118e-01]
 [3.26907337e-01 6.73092663e-01]
 [1.00000000e+00 6.30463650e-12]
 [3.64735693e-04 9.99635220e-01]
 [4.45694804e-01 5.54305196e-01]
 [1.00000000e+00 1.68278656e-29]
 [3.26889098e-01 6.73110962e-01]
 [8.43645275e-01 1.56354696e-01]
 [8.96089852e-01 1.03910178e-01]
 [8.43645275e-01 1.56354696e-01]
 [1.00000000e+00 1.68278656e-29]
 [5.28947175e-01 4.71052885e-01]
 [3.08457404e-01 6.91542566e-01]]
'''
#import numpy as np
#import matplotlib.pyplot as plt
#nums = a.split('[')
#print(len(nums))
#nums = nums[2:]
#aa = []
#for num in nums:
#    num = num.replace(']', '')
#    num = num.strip()
#    a, b = num.split()
#    aa.append(float(a))
#    
#aa.sort()
#tot = len(aa)
#y = []
#for i in range(tot):
#    y.append(i/tot)
#plt.plot(aa, y)
