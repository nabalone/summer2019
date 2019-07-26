    
def augment(images, num):
    aug_images = images[:]
    rots = numpy.array(range(1, 40))*360/40
    random.shuffle(rots)
    m = len(images)/39
    for i in range(39):
        for j in range(len(images)):
            if len(aug_images) >= num:
                return aug_images
            else:
                aug_images.append(scipy.ndimage.rotate(image, rots[(i*m +j) mod 39], reshape=False))