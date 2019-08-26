from __future__ import print_function
import keras
from keras.datasets import cifar10
from keras.preprocessing.image import ImageDataGenerator
from keras.models import Sequential, load_model
from keras.layers import Dense, Dropout, Activation, Flatten
from keras.layers import Conv2D, MaxPooling2D
import os
from pad_and_sort import load_data
import plot_cm as c

from keras.wrappers.scikit_learn import KerasClassifier
from sklearn.model_selection import StratifiedKFold, LeaveOneOut
from sklearn.model_selection import cross_val_score
import numpy as np

USE_SAVED = False
TRAIN = True

num_classes = 5
epochs = 50

save_dir = os.path.join(os.getcwd(), 'saved_models')
model_name = 'one.h5'
model_path = os.path.join(save_dir, model_name)

# Convert class vectors to binary class matrices.
X = np.load("x_all.npy")
Y = np.load("y_all.npy")
#Y = keras.utils.to_categorical(Y, num_classes)

if USE_SAVED:
        model = load_model(model_path)
else:
        model = Sequential()
        model.add(Conv2D(32, (3, 3), padding='same',
                         input_shape=x_train.shape[1:]))
        model.add(Activation('relu'))
        model.add(Conv2D(32, (3, 3)))
        model.add(Activation('relu'))
        model.add(MaxPooling2D(pool_size=(2, 2)))
        #model.add(Dropout(0.25))

        model.add(Conv2D(64, (3, 3), padding='same'))
        model.add(Activation('relu'))
        model.add(Conv2D(64, (3, 3)))
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
        opt = keras.optimizers.rmsprop(lr=0.0001, decay=1e-6)

        # Let's train the model using RMSprop
        model.compile(loss='categorical_crossentropy',
                      optimizer=opt,
                      metrics=['accuracy'])

        x_train = x_train.astype('float32')
        x_test = x_test.astype('float32')
        x_train /= 1000000#255
        x_test /= 1000000#255

    if TRAIN:
        model.fit(X, Y,
                  batch_size=len(X),
                  epochs=epochs,
                  shuffle=True,
                  class_weight='balanced'
                  )

# Save model and weights
if not os.path.isdir(save_dir):
    os.makedirs(save_dir)
model_path = os.path.join(save_dir, model_name)
model.save(model_path)
print('Saved trained model at %s ' % model_path)


y_pred = model.predict(X)
np.save("ypred_all", y_pred)
print("y_pred:")
print(y_pred)
print("y:")
print(Y)
c.plot_confusion_matrix(Y, y_pred, name_extension="_all")
