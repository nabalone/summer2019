
from __future__ import print_function
from keras.models import Model
from keras.layers import Input, Dense, Dropout, Activation, Flatten
from keras.layers import Conv2D, MaxPooling2D, AveragePooling2D, concatenate
from keras.callbacks import EarlyStopping
from keras.optimizers import Adam
from keras.utils import to_categorical
from keras.models import load_model
import os
import numpy as np
from sklearn.model_selection import StratifiedKFold
import random
from sklearn.metrics import confusion_matrix
from scipy import ndimage
import argparse
import json
import glob

PROJ_HOME = os.environ['DATA_SRCDIR']
OUTPUT_DIR = PROJ_HOME + '/src/outputs/'

parser = argparse.ArgumentParser()
parser.add_argument('-e', action='store_true', dest='use_extracted', 
    help='also input the mannually extracted properties from SEP')
parser.add_argument('-t', action='store_true', dest='test',
    help='for testing purposes, pass in the answers in place of the manually extracted properties')
parser.add_argument('-b', nargs=1, type=int, 
    help='batch size')
parser.add_argument('-p', nargs=1, type=int, dest='pooling',
    help='add extra pooling layer of this factor')
parser.add_argument('-d', nargs=1, type=float, dest='dropout',
    help='add extra dropout layer of this factor, default 0.25')
parser.add_argument('-c', action='store_true', dest='conv',
    help='add extra conv layer?')
parser.add_argument('-m', nargs=1, type=int,
    help='model number for naming.')
parser.add_argument('-n', nargs=1, type=int,
    help='number of epochs')
parser.add_argument('-s', action='store_true', dest='save',
    help='save model and results? DO NOT USE WITH SIMULTANEOUS RUNS')
parser.add_argument('-a', action='store_true', dest='ia_only',
    help='classify ia vs. other only')
parser.add_argument('-l', nargs=1, type=float,
    help='learning rate')
parser.add_argument('--no_sls', action='store_true', dest='no_sls',
    help='remove sls from sample')
parser.add_argument('--alt', action='store_true', dest='use_alt',
    help='MIGHT BE ALTERNATE VALIDATION SET use alternate shuffle (seed_offset 10)')
parser.add_argument('-3', action='store_true', dest='three_categories',
    help='classify between three categories: ia; ii and iin; ibc and sls')
parser.add_argument('--mask', action='store_true',
    help='include most likely host mask layer')
parser.add_argument('--mp', nargs=1, type=int, dest='mean_pooling',
    help='add an extra mean pooling layer of this factor')
parser.add_argument('--kfold', action='store', default='-1', dest='k_fold',
    help='use kfold crossval with 4 folds')
parser.add_argument('--all', action='store_true',
    help='Train on entire sample and save model, no validation')
args, _remaining = parser.parse_known_args() # parser.parse_args() 

#TODO move utils to separate file


def augment(images, corresponding_properties, num, rotate=True):
  
    aug_images = []
    aug_props = list(corresponding_properties)
    l = len(images)

    for i in range(num):
        image = images[i % l]
        if rotate:
            image = ndimage.rotate(image, 360*random.random(), reshape=False)
            if random.random() > 0.5:
                image = np.flipud(image)
            if random.random() > 0.5:
                image = np.fliplr(image)
        aug_images.append(image)
    if args.use_extracted:
        raise
#TODO Restore!!!!!
        #aug_props.append(corresponding_properties[i%l])
    return (aug_images, aug_props)
            
#UNUSED
#TODO: use instead of random.seed random.shuffle
def shuffle(X, y, X_sep=None):
    X = np.array(X)
    y = np.array(y)
    if len(X) != len(y):
        raise Exception("shuffle received unequal length arguments")
    new_ind = np.random.permutation(len(y))
    
    if X_sep:
        X_sep = np.array(X_sep)
        if len(X) != len(X_sep):
            print(len(X))
            print(len(X_sep))
            raise Exception("shuffle x_sep uneqal lengths")
        if args.use_extracted:
            raise
        #TODO restore!!!!!!!
        return (X[new_ind, :], y[new_ind], [])#X_sep[new_ind, :])
    else:
        return (X[new_ind, :], y[new_ind])
#TODO check axes
        
def load_fixed_kfold(ia_only=False, three=False, mask=False, num_splits=12, 
                     seed_offset=0):
    
    DATASET_DIR = OUTPUT_DIR + 'cnn_datasets'
    if ia_only:
        DATASET_DIR = DATASET_DIR + '_ia'
    elif three:
        DATASET_DIR = DATASET_DIR + '_three_cats'
    if mask:
        DATASET_DIR = DATASET_DIR + '_mask'
        
    # clear leftover datasets of this type from previous runs
    old_datasets = glob.glob(DATASET_DIR + '/*')
    for old_dataset in old_datasets:
        os.remove(old_dataset)
    #make dataset directory if it does not exist
    if not os.path.isdir(DATASET_DIR): 
        os.mkdir(DATASET_DIR)
    DATASET_DIR = DATASET_DIR + "/"
    
    n_ia = len(np.load(OUTPUT_DIR + "x_all2_0.npy"))
    n_ibc = len(np.load(OUTPUT_DIR + "x_all2_1.npy"))
    n_ii = len(np.load(OUTPUT_DIR + "x_all2_2.npy"))
    n_iin = len(np.load(OUTPUT_DIR + "x_all2_3.npy"))
    n_sls = len(np.load(OUTPUT_DIR + "x_all2_4.npy"))
    n_tot = n_ia + n_ibc + n_ii + n_iin + n_sls 
    
    if num_splits != 'loo' and num_splits != 'all':
        if type(num_splits)!= int or num_splits < 1:
            raise Exception("num_splits argument must be 'loo', 'all', or a  \
                            positive integer")
        elif num_splits > np.min([n_ia, n_ibc, n_ii, n_iin, n_sls]):
            raise Exception("num_splits must be smaller than the smallest \
                            class size, %s. use Leave One Out instead with \
                            num_splits='loo" \
                            % np.min([n_ia, n_ibc, n_ii, n_iin, n_sls]))
    
    extrastring = str(seed_offset) if seed_offset>0 else ''

    filname = "aug_all"

    if ia_only:
        filname = filname + "_ia"
    elif three:
        filname = filname + "_three_cats"
    if mask:
        filname = filname + "_mask"   
    
    # the numbers to which to augment each type
    # so we end up with equal total numbers in the groups we are classifying 
    # between (e.g. as many Ia as non-Ia if doing ia vs. other) and within 
    # a group (such as non-Ia), numbers are proportional to type prevalence
    # e.g. the ratio of ii to ibc in the augmented non-Ia sample is the same as 
    # the ratio of ia to non-ia in the original sample
    if ia_only:
        aug_to = [n_ia, 
                    np.round(n_ia*n_ibc/(n_tot - n_ia)), 
                    np.round(n_ia*n_ii/(n_tot - n_ia)),
                    np.round(n_ia*n_iin/(n_tot - n_ia)),
                    np.round(n_ia*n_sls/(n_tot - n_ia))]
    elif three:
        aug_to = [n_ia, 
                    np.round(n_ia*n_ibc/(n_ibc+n_sls)),
                    np.round(n_ia*n_ii/(n_ii+n_iin)),
                    np.round(n_ia*n_iin/(n_ii+n_iin)),
                    np.round(n_ia*n_sls/(n_ibc+n_sls))]
    else:
        largest_sample_size = np.max([n_ia, n_ibc, n_ii, n_iin, n_sls])
        aug_to = [largest_sample_size]*5
    
    #load raw data
    raw = [[],[],[],[],[]] 
    for i in range(5):#for each type
        #load all data of that type
        #dividing by arbitrary 1000000 to avoid overflows
        if mask:
            raw[i] = np.load(OUTPUT_DIR + "x_all2_%s.npy" % i).astype('float32')/1000000.
            #raw_sep = np.load("x_ans_%s.npy" % i).astype('float32')/1000000.
        else:
            raw[i] = np.load(OUTPUT_DIR + "x_all_%s.npy" % i).astype('float32')/1000000.
            #raw_sep = np.load("x_ans_%s.npy" % i).astype('float32')/1000000.
            
        #shuffle the data for this type
        random.seed(i + seed_offset)
        random.shuffle(raw[i])
        #random.seed(i+seed_offset)
        #random.shuffle(raw_sep)
    
    
    
    '''Normal k folding:'''
    if num_splits != 'all' and num_splits != 'loo':
        
        # for each of the 5 types, we split that type into num_splits splits 
        #folds[i] is a list of (train, test) splits of sn type i
        folds = [[],[],[],[],[]]
        for i in range(5):
            #TODO reove stratification, its useless since we're just doing 1 type  
            #TODO remove one of the shuffles. we already shuffled all raw[i] above
            folds[i] = list(StratifiedKFold(n_splits=num_splits, shuffle=True, 
                                         random_state=1+seed_offset).split(raw[i], [i]*len(raw[i])))
        for j in range(num_splits): #create the jth fold
            print('j: %s' %j)
            jth_X_train = []
            jth_y_train = []
            jth_X_test = []
            jth_y_test = []
            for i in range(5): #iterate through the types
                print('i: %s' %i)
                train_of_i, test_of_i = folds[i][j]
    #TODO this has been changed from the sep
                #augment data so we have aug_to[i], the right number of samples 
                train_aug, _train_aug_sep = augment(raw[i][train_of_i], 
                                                    [0]*len(train_of_i), 
                                                    aug_to[i])
                #crop all images from 240x240 to 160x160 to eliminate the blank
                #corners caused by rotation in augmentation
                # distribute among folds
                jth_X_train.extend(crop(train_aug))
                jth_y_train.extend([i]*len(train_aug))
                
                #only "augmented" up to original length, so no rotations added
                test_aug, _test_aug_sep = augment(raw[i][test_of_i], 
                                                  [0]*len(train_of_i), 
                                                  len(raw[i][test_of_i])) 
                jth_X_test.extend(crop(test_aug))
                jth_y_test.extend([i]*len(test_aug))   
          
            # after we have added the samples from all 5 types to each fold,
            # we can shuffle and save
        
            #TODO make sure shuffle works
            print("saving")
            random.seed(100+j+ seed_offset)
            random.shuffle(jth_X_train)
            random.seed(100+j+ seed_offset)
            random.shuffle(jth_y_train)
            np.savez(DATASET_DIR + filname+extrastring+'_fold_%s'%j, 
                 jth_X_train, jth_y_train, 
                 jth_X_test, jth_y_test)       
                     
#TODO I DON'T KNOW IF ANY OF THE FOLLOWING IS CORRECT
            
    # if not plain kfolding, make "all_of" arrays: augmented of entire set of each type
    else:
        for i in range(5):
            X_all_of = [[],[],[],[],[]]    
            train_aug_all, _train_aug_sep_all = augment(raw[i], 
                                                        [0]*len(raw[i]), 
                                                        aug_to[i])
            X_all_of[i] = list(crop(train_aug_all)) 
        
    
    '''LOO Folding'''
    if num_splits == 'loo':
        count = 0
        for i in range(5): 
            for ind in range(len(raw[i])): #for each SN of that type
                X_test_fold = [raw[i][ind]] #set that SN aside as test
                y_test_fold = [i]
                
                X_train_fold = [] #create the training set
                y_train_fold = [] #fill answers swet in parallel
                for j in range(5):
                    # for each type, add all of that type to training set 
                    # except for i, which is the type of the test SN
                    # in which case add all except for the test SN which is at ind
                    if j==i:
                        #omit the set-aside SN, augment from the rest
                        raw_without_ind = np.delete(raw[j], ind, 0)
                        aug_without_ind, _aug_without_ind_sep = augment(raw_without_ind, 
                                                                [0]*len(raw_without_ind), 
                                                                aug_to[j])#X_all_of[j][:ind] + X_all_of[j][ind+1:] 
                        X_train_fold.extend(list(crop(aug_without_ind)))
                        y_train_fold.extend([j]*len(aug_without_ind))
                    else:
                        X_train_fold.extend(X_all_of[j])
                        y_train_fold.extend([j]*len(X_all_of[j]))
                 
                # shuffle training set and save
                random.seed(1000*i+ind+ seed_offset)
                random.shuffle(X_train_fold)
                random.seed(1000*i+ind+ seed_offset)
                random.shuffle(y_train_fold)
                
                print(count)
                print(np.array(X_train_fold).shape)
                print(np.array(y_train_fold).shape)
                print('\n')
                
                np.savez(DATASET_DIR + filname+extrastring+'_fold_%s'%count, 
                         X_train_fold, y_train_fold, 
                         X_test_fold, y_test_fold)  
                count+=1
                
    elif num_splits == 'all':
        X_train_fold = []
        y_train_fold = []
        for i in range(5):
            X_train_fold.extend(X_all_of[i])
            y_train_fold.extend([i]*len(X_all_of[i]))
        
        random.seed(23+ seed_offset)
        random.shuffle(X_train_fold)
        random.seed(23+ seed_offset)
        random.shuffle(y_train_fold)
        
        np.savez(DATASET_DIR + filname+extrastring+'_all', 
                X_train_fold, y_train_fold, 
                [], [])   
        print(np.array(X_train_fold).shape)
        print(np.array(y_train_fold).shape)
            
                

            
def crop(ls):
    a = np.array(ls)
    b = a[:, 40:-40, 40:-40]
    return(list(b))


def main():

    if args.test:
        print("TEST NOT YET IMPLEMENTED")

    no_sls = args.no_sls

    if int(args.k_fold) >= 0:
        k_folded=True
        fold_num = int(args.k_fold)
    else:
        k_folded=False
        
    if args.m:
        model_num = str(args.m[0])
    else:
        model_num= '11'

    if args.l:
        LR = float(args.l[0])
    else:
        LR = 0.0005

    ia_only = args.ia_only
    three_categories = args.three_categories


    USE_SAVED = False
    #TODO restore
    if args.ia_only:
        num_classes = 2
    elif three_categories:
        num_classes = 3
    else:
        num_classes = 5

    if args.n:
        epochs = args.n[0]
    else:
    #TODO restore
        epochs = 35

    save_dir = OUTPUT_DIR 
    model_name = 'aardvark_aug' + model_num + '.h5'

    if (args.use_extracted) and ia_only:
        raise(ValueError, "NOT YET IMPLEMENTED IA WITH EXTRACTED")
        exit(1)

#TODO restore


    filname = "aug_all"
    DATASET_DIR = OUTPUT_DIR + "cnn_datasets"
    if ia_only:
        filname = filname + "_ia"
        DATASET_DIR = DATASET_DIR + "_ia"
    elif three_categories:
        filname = filname + "_three_cats"
        DATASET_DIR = DATASET_DIR + "_three_cats"
    if args.mask:
        filname = filname + "_mask"
        DATASET_DIR = DATASET_DIR + "_mask"
    if k_folded:
        filname = filname + "_fold_%s"%fold_num
    
    if args.all:
        filname = filname + "_all"
        
    if (k_folded and args.all) or (not k_folded and not args.all):
        raise Exception("Must use either --k_fold or --all")
    DATASET_DIR = DATASET_DIR + "/"
    all_data = np.load(DATASET_DIR + filname + '.npz')

    X_full_train = all_data['arr_0']
    y_full_train = all_data['arr_1']
    X_full_test = all_data['arr_2']
    y_full_test = all_data['arr_3']

    if ia_only: 
        y_full_train = np.where(y_full_train==0, 0, 1)
        y_full_test = np.where(y_full_test==0, 0, 1)

    if three_categories:
        raise Exception("three categories not yet implemented")
        #make all sls into ibc
        #y_full = np.where(y_full==4, 1, y_full)
        
        #make all iin into ii
        #y_full = np.where(y_full==3, 2, y_full)

    if no_sls:
        raise Exception("no_sls not yet implemented")
        #y_full_orig = y_full
        #y_full = y_full[y_full_orig != 4]
        #X_full = X_full[y_full_orig != 4]

    y_full_test_orig = y_full_test
    y_full_train = to_categorical(y_full_train, num_classes)
    y_full_test = to_categorical(y_full_test, num_classes)

   
        
    if args.b:
        batch_size = args.b[0]
    else:
        batch_size = len(X_full_train)
        
    model_path = os.path.join(save_dir, model_name)
    def get_model(in1_shape, in2_shape):
        if USE_SAVED:
            return(load_model(model_path))
        else:
            input1 = Input(shape=in1_shape)
            if args.use_extracted:
                input2 = Input(shape=in2_shape)
            if args.conv:
                conv = Conv2D(32, (3, 3), padding='same')(input1)
                if args.dropout:
                    dropout_conv=Dropout(args.dropout[0])(conv)
                    prev_for_act1 = dropout_conv
                else:
                    prev_for_act1 = conv
                
                act1 = Activation('relu')(prev_for_act1)
                if args.dropout:
                    dropout_act1=Dropout(args.dropout[0])(act1)
                    prev_for_pool1 = dropout_act1
                else:
                    prev_for_pool1 = act1
                if args.mean_pooling:
                    pool1 = AveragePooling2D(pool_size=(2, 2))(prev_for_pool1)
                else:
                    pool1 = MaxPooling2D(pool_size=(2, 2))(prev_for_pool1)
                if args.dropout:
                    dropout_pool1=Dropout(args.dropout[0])(pool1)
                    prev_for_convb = dropout_pool1
                else:
                    prev_for_convb = pool1
    
    #TODO proofread this whole section
                convb = Conv2D(32, (3, 3), padding='same')(prev_for_convb)
                if args.dropout:
                    dropout_convb=Dropout(args.dropout[0])(convb)
                    prev_for_act1b = dropout_convb
                else:
                    prev_for_act1b = convb
                
                act1b = Activation('relu')(prev_for_act1b)
                if args.dropout:
                    dropout_act1b=Dropout(args.dropout[0])(act1b)
                    prev_for_pool1b = dropout_act1b
                else:
                    prev_for_pool1b = act1b
                if args.mean_pooling:
                    pool1b = AveragePooling2D(pool_size=(2, 2))(prev_for_pool1b)
                else:
                    pool1b = MaxPooling2D(pool_size=(2, 2))(prev_for_pool1b)
                if args.dropout:
                    dropout_pool1b=Dropout(args.dropout[0])(pool1b)
                    prev_for_conv2 = dropout_pool1b
                else:
                    prev_for_conv2 = pool1b
    
    
            else:
                prev_for_conv2=input1
            
            conv2 = Conv2D(32, (3, 3), padding='same')(prev_for_conv2)
            if args.dropout:
                dropout_conv2=Dropout(args.dropout[0])(conv2)
                prev_for_act2 = dropout_conv2
            else:
                prev_for_act2 = conv2
    
            act2 = Activation('relu')(prev_for_act2)
            if args.dropout:
                dropout_act2=Dropout(args.dropout[0])(act2)
                prev_for_pool2= dropout_act2
            else:
                prev_for_pool2 = act2
            if args.mean_pooling:
                pool2 = AveragePooling2D(pool_size=(2, 2))(prev_for_pool2)
            else:
                pool2 = MaxPooling2D(pool_size=(2, 2))(prev_for_pool2)
            
            if args.pooling:
                pool3 = MaxPooling2D(pool_size=(args.pooling[0], args.pooling[0]))(pool2)
                flatten=Flatten()(pool3)
            elif args.mean_pooling:
                pool3 = AveragePooling2D(pool_size=(args.mean_pooling[0], args.mean_pooling[0]))(pool2)
                flatten=Flatten()(pool3)
            else:
                flatten = Flatten()(pool2)
            if args.dropout:
                dropout = Dropout(args.dropout[0])(flatten)
                dense1 = Dense(512, activation=('relu'))(dropout)
            else:
                dense1 = Dense(512, activation=('relu'))(flatten)
            if args.use_extracted:
                concat = concatenate([dense1, input2], axis=-1)
                dense2 = Dense(num_classes, activation=('softmax'))(concat)
            else:
                dense2 = Dense(num_classes, activation=('softmax'))(dense1)
            # initiate RMSprop optimizer
            opt = Adam(lr=LR)
            #TODO check above
    
            # Let's train the model using RMSprop
            if args.use_extracted:
                ins = [input1, input2]
            else:
                ins = input1
            model = Model(inputs=ins, outputs=dense2)
            model.compile(loss='categorical_crossentropy',
                  optimizer=opt,
                  metrics=['accuracy'])
            return model

    # Early stopping not used because validation accuracy does not decrease 
    # with overfitting but just continues to fluctuate with noise after 
    # plateauing
    #es = EarlyStopping(monitor='val_acc', patience=50, verbose=1, baseline=0.4, restore_best_weights=True)


    if args.use_extracted:
        raise("not yet implemented")

#TODO fix in2_shape
    model = get_model(in1_shape=X_full_train.shape[1:], in2_shape=X_full_train.shape[1:])
        
    if k_folded:        
        model.fit(x=X_full_train, y=y_full_train,
          batch_size=batch_size,
          epochs=epochs,
          validation_data=(X_full_test, y_full_test),
          shuffle=True)
            
        y_pred = model.predict(X_full_test)
        y_pred_2 = np.argmax(y_pred, 1)
          
        if args.ia_only:
            if not os.path.isdir(OUTPUT_DIR + 'cnn_kfold_results_ia'):
                os.mkdir(OUTPUT_DIR + 'cnn_kfold_results_ia')
            np.save(OUTPUT_DIR + 'cnn_kfold_results_ia/y_pred_ia_fold%s' % args.k_fold, y_pred_2)
            np.save(OUTPUT_DIR + 'cnn_kfold_results_ia/y_true_ia_fold%s' % args.k_fold, y_full_test_orig)
        else:
            if not os.path.isdir(OUTPUT_DIR + 'cnn_kfold_results'):
                os.mkdir(OUTPUT_DIR + 'cnn_kfold_results')
            np.save(OUTPUT_DIR + 'cnn_kfold_results/y_pred_fold%s' % args.k_fold, y_pred_2)
            np.save(OUTPUT_DIR + 'cnn_kfold_results/y_true_fold%s' % args.k_fold, y_full_test_orig)
        cm = confusion_matrix(y_full_test_orig, y_pred_2)
        print(cm)
    
    elif args.all:
        model.fit(x=X_full_train, y=y_full_train,
          batch_size=batch_size,
          epochs=epochs,
          shuffle=True)
        
    model_path = OUTPUT_DIR + 'final_trained_cnn'
    if args.ia_only:
        model_path = model_path + "_ia_only"
    if not args.mask:
        model_path = model_path + "_no_mask"
    model_path = model_path + ".h5"
    model.save(model_path)
    print('Saved trained model at %s ' % model_path)
    json_string = model.to_json()
    with open(OUTPUT_DIR + 'cnn_architecture.json', 'w+') as f:
        json.dump(json_string, f)
        

if __name__ == "__main__":
    main()
