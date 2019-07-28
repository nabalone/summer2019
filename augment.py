import numpy as np
from scipy import ndimage
import math
import matplotlib.pyplot as plt

import random

X = np.load('../x_all.npy')[:3]

# for naming plots files
namecount = 0
def namegen():
    global namecount
    namecount += 1
    return 'deletable\Ximage' + str(namecount) + ".png"
    
def augment(images, num):
    aug_images = list(images)
    rots = np.array(range(1, 40))*360/40
    random.shuffle(rots)
    m = 39/len(images)
    for i in range(39):
        for j in range(len(images)):
            if len(aug_images) >= num:
                return aug_images
            else:
                aug_images.append(ndimage.rotate(images[j], rots[int((i*m +j) % 39)], reshape=False))
                
                
X_aug = augment(X, 11)
for data in X_aug:
    plt.imshow(data, interpolation='nearest', cmap='gray')
    plt.savefig(namegen())
    plt.close()
    