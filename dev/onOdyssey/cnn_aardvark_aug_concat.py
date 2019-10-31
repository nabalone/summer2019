
from __future__ import print_function
#import keras
#from keras.datasets import cifar10
#from keras.preprocessing.image import ImageDataGenerator
from keras.models import Model
from keras.layers import Input, Dense, Dropout, Activation, Flatten
from keras.layers import Conv2D, MaxPooling2D, concatenate
from keras.callbacks import EarlyStopping
from keras.optimizers import Adam
from keras.utils import to_categorical
from keras.models import load_model
import os
import numpy as np
#from sklearn.utils.class_weight import compute_class_weight
import math
import random
from sklearn.metrics import confusion_matrix
from scipy import ndimage
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-e', action='store_true', dest='use_extracted')
args = parser.parse_args() 


model_num = '10'

def augment_old(images, num):
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


def augment(images, corresponding_properties, num):
    aug_images = list(images)
    aug_props = list(corresponding_properties)
    l = len(images)
    if l > num:
        aug_images = aug_images[:num]
        aug_props = aug_props[:num]
    for i in range(num - l):
        image = images[i % l]
        image = ndimage.rotate(image, 360*random.random(), reshape=False)
        if random.random() > 0.5:
            image = np.flipud(image)
        if random.random() > 0.5:
            image = np.fliplr(image)
        aug_images.append(image)
        aug_props.append(corresponding_properties[i%l])
    return (aug_images, aug_props)
            

def shuffle(X, y, X_sep=None):
    X = np.array(X)
    y = np.array(y)
    if len(X) != len(y):
        raise Exception("shuffle received unequal length arguments")
    new_ind = np.random.permutation(len(y))
    
    if X_sep:
        X_sep = np.array(X_sep)
        if len(X) != len(X_sep):
            raise Exception("shuffle x_sep uneqal lengths")
        return (X[new_ind, :], y[new_ind], X_sep[new_ind, :])
    else:
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
    Xsep_test = []
    Xsep_train = []
    Xsep_val = []
    for i in range(5):
        print(i)
        raw = np.load("x_all_%s.npy" % i).astype('float32')/1000000.
        raw_sep = np.load("x_ans_%s.npy" % i).astype('float32')/1000000.
        print(len(raw))
        random.seed(i)
        random.shuffle(raw)
        random.seed(i)
        random.shuffle(raw_sep)
#TODO handpick which goes in which
        #div = math.floor(len(raw)/3)
        div = math.floor(0.2*len(raw))
        div = max(div, 3)
        print(div)
        testi, testp = augment(raw[:div], raw_sep[:div], NUM)
        vali, valp = augment(raw[div:2*div], raw_sep[div:2*div], NUM)
        traini, trainp = augment(raw[2*div:], raw_sep[2*div:], NUM)
        X_test.extend(testi)
        y_test.extend([i]*len(testi))
        Xsep_test.extend(testp)
        X_train.extend(traini)
        y_train.extend([i]*len(traini))
        Xsep_train.extend(valp)
        X_val.extend(vali)
        y_val.extend([i]*(len(vali)))
        Xsep_val.extend(valp)
    X_test, y_test, Xsep_test = shuffle(X_test, y_test, Xsep_test)
    X_train, y_train, Xsep_train = shuffle(X_train, y_train, Xsep_train)    
    X_val, y_val, Xsep_val = shuffle(X_val, y_val, Xsep_val)       
    np.save("y_test_aardvark_aug"+model_num, y_test)
    np.savez("aug_load_output", X_test, y_test, Xsep_test, X_train, y_train, Xsep_train, X_val, y_val, Xsep_val)   

    return(X_test, y_test, Xsep_test, X_train, y_train, Xsep_train, X_val, y_val, Xsep_val)   
USE_SAVED = False

num_classes = 5
epochs = 25
#num_predictions = 20
save_dir = os.path.join(os.getcwd(), 'saved_models')
model_name = 'aardvark_aug' + model_num + '.h5'
def main():
    # The data, split between train and test sets:
    #(x_train, y_train), (x_test, y_test) = cifar10.load_data()
    #print('x_train shape:', x_train.shape)
    #print(x_train.shape[0], 'train samples')
    #print(x_test.shape[0], 'test samples')

    #x_train = np.load("x_all.npy")
    #y_train = np.load("y_all.npy")

    all_data =  np.load("aug_load_output.npz") #load()
    X_test = all_data['arr_0'] 
    y_test= all_data['arr_1']
    Xsep_test = all_data['arr_2']
    X_train= all_data['arr_3']
    y_train= all_data['arr_4']
    Xsep_train= all_data['arr_5']
    X_val= all_data['arr_6']
    y_val= all_data['arr_7']
    Xsep_val= all_data['arr_8']

    y_test_orig = y_test
    y_train = to_categorical(y_train, num_classes)
    y_test = to_categorical(y_test, num_classes)
    y_val = to_categorical(y_val, num_classes)
    batch_size = len(X_train)
    # Convert class vectors to binary class matrices.
    #y_train = keras.utils.to_categorical(y_train, num_classes)
    #y_test = keras.utils.to_categorical(y_test, num_classes)

    model_path = os.path.join(save_dir, model_name)
    if USE_SAVED:
        model = load_model(model_path)
    else:
        input1 = Input(shape=X_train.shape[1:])
        input2 = Input(shape=Xsep_train.shape[1:])
        conv = Conv2D(32, (3, 3), padding='same')(input1)
        act1 = Activation('relu')(conv)
        pool1 = MaxPooling2D(pool_size=(2, 2))(act1)
        #model.add(Dropout(0.25))

        flatten = Flatten()(pool1)
        dense1 = Dense(512, activation=('relu'))(flatten)
        #model.add(Dropout(0.5))
        concat = concatenate([dense1, input2], axis=-1)
        dense2 = Dense(num_classes, activation=('softmax'))(concat)

        # initiate RMSprop optimizer
        opt = Adam(lr=0.0005)
        #TODO check above

        # Let's train the model using RMSprop
        model = Model(input=[input1, input2], output=dense2)
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
        es = EarlyStopping(monitor='val_acc', patience=10, verbose=1, baseline=0.4, restore_best_weights=True)
        model.fit(x=[X_train, Xsep_train], y=y_train,
              batch_size=batch_size,
              epochs=epochs,
              validation_data=([X_val, Xsep_val], y_val),
              callbacks = [es],
              #sample_weight=sampleweights,
              shuffle=True)
        # Save model and weights
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)
        model_path = os.path.join(save_dir, model_name)
        model.save(model_path)
        print('Saved trained model at %s ' % model_path)
    y_pred = model.predict([X_test, Xsep_test])
    y_pred2 = np.argmax(y_pred, 1)
    np.save("y_pred_aardvark25weighted_aug"+model_num, y_pred2)
    np.save("y_test", y_test)
    print(y_pred)
    print(y_pred2)
    # Score trained model.
    scores = model.evaluate([X_test, Xsep_test], y_test, verbose=1)
    print('Test loss:', scores[0])
    print('Test accuracy:', scores[1])

    cm = confusion_matrix(y_test_orig, y_pred2)
    print(cm)
    np.save("cm_aa"+model_num, cm, allow_pickle=True, fix_imports=True)

if __name__ == "__main__":
    main()
