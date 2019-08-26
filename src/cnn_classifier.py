
from __future__ import print_function
import keras
from keras.datasets import cifar10
from keras.preprocessing.image import ImageDataGenerator
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation, Flatten
from keras.layers import Conv2D, MaxPooling2D
from keras.callbacks import EarlyStopping
import os
import numpy as np
from sklearn.utils.class_weight import compute_class_weight
import math
import random
from sklearn.metrics import confusion_matrix
from scipy import ndimage

def augment(images, num):
    print("augmenting")
    print(len(images))
    aug_images = list(images)
    rots = np.array(range(1, 40))*360/40
    random.shuffle(rots)
    m = 39/len(images)
    for i in range(100):
        for j in range(len(images)):
            if len(aug_images) >= num:
                print("returning")
                return aug_images
            else:
                aug_images.append(ndimage.rotate(images[j], rots[int((i*m +j) % 39)], reshape=False))

def shuffle(X, y):
    X = np.array(X)
    y = np.array(y)
    if len(X) != len(y):
        raise Exception("shuffle received unequal length arguments")
    new_ind = np.random.permutation(len(y))
    return (X[new_ind, :], y[new_ind])
#TODO check axes

def load():
    NUM = 200
    X_test = []
    X_train = []
    X_val = []
    y_test = []
    y_train = []
    y_val = []
    for i in range(5):
        print(i)
        raw = np.load("x_all_%s.npy" % i).astype('float32')/1000000.
        print(len(raw))
        random.seed(i)
        random.shuffle(raw)
#TODO handpick which goes in which
        #div = math.floor(len(raw)/3)
        dif = math.floor(0.2*len(raw))
	dif = min(dif, 2)
	print(div)
        testi = augment(raw[:div], NUM)
        print(type(testi))
        print(len(testi))
        vali = augment(raw[div:2*div], NUM)
        traini = augment(raw[2*div:], NUM)
        X_test.extend(testi)
        y_test.extend([i]*len(testi))
        X_train.extend(traini)
        y_train.extend([i]*len(traini))
        X_val.extend(vali)
        y_val.extend([i]*(len(vali)))
    X_test, y_test = shuffle(X_test, y_test)
    X_train, y_train = shuffle(X_train, y_train)    
    X_val, y_val = shuffle(X_val, y_val)       
    np.save("y_test_aardvark_aug", y_test)
    return( X_test, y_test, X_train, y_train, X_val, y_val)   
USE_SAVED = False

num_classes = 5
epochs = 25
#num_predictions = 20
save_dir = os.path.join(os.getcwd(), 'saved_models')
model_name = 'aardvark_aug.h5'
def main():
    # The data, split between train and test sets:
    #(x_train, y_train), (x_test, y_test) = cifar10.load_data()
    #print('x_train shape:', x_train.shape)
    #print(x_train.shape[0], 'train samples')
    #print(x_test.shape[0], 'test samples')

    #x_train = np.load("x_all.npy")
    #y_train = np.load("y_all.npy")

    X_test, y_test, X_train, y_train, X_val, y_val = load()
    y_test_orig = y_test
    y_train = keras.utils.to_categorical(y_train, num_classes)
    y_test = keras.utils.to_categorical(y_test, num_classes)
    y_val = keras.utils.to_categorical(y_val, num_classes)
    batch_size = len(X_train)
    # Convert class vectors to binary class matrices.
    #y_train = keras.utils.to_categorical(y_train, num_classes)
    #y_test = keras.utils.to_categorical(y_test, num_classes)

    model_path = os.path.join(save_dir, model_name)
    if USE_SAVED:
        model = keras.models.load_model(model_path)
    else:
        model = Sequential()
        model.add(Conv2D(32, (3, 3), padding='same',
                 input_shape=X_train.shape[1:]))
        model.add(Activation('relu'))
        model.add(MaxPooling2D(pool_size=(2, 2)))
        #model.add(Dropout(0.25))

        model.add(Flatten())
        model.add(Dense(512))
        model.add(Activation('relu'))
        #model.add(Dropout(0.5))
        model.add(Dense(num_classes))
        model.add(Activation('softmax'))

        # initiate RMSprop optimizer
        opt = keras.optimizers.Adam(lr=0.0005)
        #TODO check above

        # Let's train the model using RMSprop
        model.compile(loss='categorical_crossentropy',
              optimizer=opt,
              metrics=['accuracy'])

        
         #x_train = x_train.astype('float32')
        #x_test = x_test.astype('float32')
        #x_train /= 1000000
        #x_test /= 255

        #make balanced class weights
        #a = compute_class_weight('balanced', np.unique(y_orig), y_orig)
    #    sampleweights = []
    #    for sample in y_train:
    #        sampleweights.append(a[sample])
    #    sampleweights = np.array(sampleweights)
        es = EarlyStopping(monitor='val_acc', patience=5, verbose=1, baseline=0.35, restore_best_weights=True)
        model.fit(X_train, y_train,
              batch_size=batch_size,
              epochs=epochs,
              validation_data=(X_val, y_val),
              callbacks = [es],
              #sample_weight=sampleweights,
              shuffle=True)
        # Save model and weights
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)
            model_path = os.path.join(save_dir, model_name)
            model.save(model_path)
            print('Saved trained model at %s ' % model_path)
    y_pred = model.predict(X_test)
    y_pred2 = np.argmax(y_pred, 1)
    np.save("y_pred_aardvark25weighted_aug", y_pred2)
    print(y_pred)
    print(y_pred2)
    # Score trained model.
    scores = model.evaluate(X_test, y_test, verbose=1)
    print('Test loss:', scores[0])
    print('Test accuracy:', scores[1])

    cm = confusion_matrix(y_test_orig, y_pred2)
    print(cm)

if __name__ == "__main__":
    main()
