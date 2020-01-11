from itertools import chain

import keras
from keras.layers import Dropout, Flatten, Dense, Input, Concatenate, Conv2D, MaxPooling2D
from collections import defaultdict
import numpy as np
import json
import sys
import random
import pickle
import os


def datafetcher(path, imgpath, classification, imageclass, splitratio):
    imgfiles = os.listdir(imgpath)
    with open(f'{imgpath}{imgfiles[0]}', "rb") as pa:
        image = pickle.load(pa)
    imagelen = len(image)
    pixellen = len(image[0])

    if not classification:
        names = []
        labels = {}
        if os.path.exists(f'{path}subimage_filtered.json'):
            for line in open(f'{path}subimage_filtered.json'):
                data = json.loads(line)
                name = f'{data["image"]}.txt'
                names.append(name)
                labels[name] = data[imageclass]
        elif os.path.exists(f'{path}subimage.json'):
            for line in open(f'{path}subimage.json'):
                data = json.loads(line)
                name = f'{data["image"]}.txt'
                names.append(name)
                labels[name] = data[imageclass]
        else:
            print('No metadata for images exists')
        splits = round(len(names) * float(splitratio))
        random.shuffle(names)
        trainlist = names[0:splits]
        vallist = names[splits:]
        partition = {'train': trainlist, 'validation': vallist}

    else:
        labels = {}
        if os.path.exists(f'{path}subimage_filtered.json'):
            for line in open(f'{path}subimage_filtered.json'):
                data = json.loads(line)
                name = f'{data["image"]}.txt'
                labels[name] = data[imageclass]
        elif os.path.exists(f'{path}subimage.json'):
            for line in open(f'{path}subimage.json'):
                data = json.loads(line)
                name = f'{data["image"]}.txt'
                labels[name] = data[imageclass]
        else:
            print('No metadata for images exists')

        labels2 = defaultdict(list)
        for k, v in labels.items():
            labels2[v].append(k)

        partition = defaultdict(list)
        for f in labels2:
            random.shuffle(labels2[f])  # Shuffles to get random images into training and validation
            splits = round(len(labels2[f]) * float(splitratio))
            trainlist = (labels2[f][0:splits])
            vallist = (labels2[f][splits:])
            partition['train'].append(trainlist)
            partition['validation'].append(vallist)

        partition['train'] = list(chain.from_iterable(partition['train']))
        partition['validation'] = list(chain.from_iterable(partition['validation']))

    return partition, labels, imagelen, pixellen


# Developing the data generator
class DataGenerator(keras.utils.Sequence):
    'Generates data for Keras'

    def __init__(self, path, list_IDs, labels, batch_size, size, n_channels, n_classes, shuffle):
        'Initialization'
        self.size = size
        self.path = path
        self.batch_size = batch_size
        self.labels = labels
        self.list_IDs = list_IDs
        self.n_channels = n_channels
        self.n_classes = n_classes
        self.shuffle = shuffle
        self.on_epoch_end()

    def __len__(self):
        'Denotes the number of batches per epoch'
        return int(np.floor(len(self.list_IDs) / self.batch_size))

    def __getitem__(self, index):
        'Generate one batch of data'
        # Generate indexes of the batch
        indexes = self.indexes[index * self.batch_size:(index + 1) * self.batch_size]

        # Find list of IDs
        list_IDs_temp = [self.list_IDs[k] for k in indexes]
        # Generate data
        X, y = self.__data_generation(list_IDs_temp)
        return X, y

    def on_epoch_end(self):
        'Updates indexes after each epoch'
        self.indexes = np.arange(len(self.list_IDs))
        if self.shuffle:
            np.random.shuffle(self.indexes)

    def __data_generation(self, list_IDs_temp):
        'Generates data containing batch_size samples'
        # Initialization
        X = np.empty((self.batch_size, self.size[1], self.size[0], self.n_channels))
        y = np.empty((self.batch_size), dtype="S20")

        # Generate data
        for i, ID in enumerate(list_IDs_temp):
            # Store sample
            with open(f'{imagepath}{ID}', "rb") as pa:
                image = pickle.load(pa)
            image = np.array(image)
            image = image[:, :, 0:self.n_channels]
            X[i,] = image
            y[i] = self.labels[ID]

        if not classification:
            return X, y
        else:
            return X, keras.utils.to_categorical(y, num_classes=self.n_classes)


# Developing the neural network
def nnmodel(imglen, pixellen, classification, n_channels, n_classes, imageclass):
    input = Input(shape=(imglen, pixellen, n_channels,))
    x = Conv2D(16, kernel_size=(1, 1), activation='relu', padding='same')(input)
    x = MaxPooling2D(pool_size=(2, 2))(x)

    x1 = Conv2D(16, kernel_size=(3, 3), activation='relu', padding='same')(input)
    x1 = MaxPooling2D(pool_size=(2, 2))(x1)

    x2 = Conv2D(16, kernel_size=(5, 5), activation='relu', padding='same')(input)
    x2 = MaxPooling2D(pool_size=(2, 2))(x2)

    x = Concatenate()([x, x1, x2])
    x = Flatten()(x)
    x = Dense(64, activation='relu')(x)
    x = Dropout(rate=0.25)(x)
    if not classification:
        output = Dense(1, activation='linear')(x)
    else:
        output = Dense(n_classes, activation='softmax')(x)
    model = keras.Model(input, output)

    if classification:
        if n_classes == 2:
            model.compile(loss='binary_crossentropy', metrics=['accuracy'], optimizer='adam')
        else:
            model.compile(loss='categorical_crossentropy', metrics=['accuracy'], optimizer='adam')
    else:
        model.compile(loss='mse', metrics=['mse'], optimizer='rmsprop')

    # Create callbacks
    if classification:
        checkpoint = keras.callbacks.ModelCheckpoint(f'Best-{imageclass}.h5', monitor='val_accuracy',
                                                     save_best_only=True)
    else:
        checkpoint = keras.callbacks.ModelCheckpoint(f'Best-{imageclass}.h5', monitor='val_mse', save_best_only=True)
    early_stopping = keras.callbacks.EarlyStopping(monitor='val_loss', patience=4)
    callbacks_list = [checkpoint, early_stopping]

    return model, callbacks_list


if __name__ == '__main__':
    # Read datapath from config file
    with open('config.json') as json_file:
        data = json.load(json_file)

    datapath = data['path']
    metapath = f'{datapath}metadata/'
    imagepath = f'{datapath}images/'

    # Cmd inputs
    classification = sys.argv[1] == 'T'
    imageclass = sys.argv[2]
    splitratio = sys.argv[3]

    nameofclass = imageclass.replace('/', '')

    output = datafetcher(metapath, imagepath, classification, imageclass, splitratio)
    partition = output[0]
    labels = output[1]
    imglen = output[2]
    pixellen = output[3]

    ###Classes and Channels:
    if classification:
        n_classes = len(labels)
    else:
        n_classes = 1

    n_channels = 4
    #######################

    params = {'size': (pixellen, imglen),
              'batch_size': 128,
              'n_classes': n_classes,
              'n_channels': n_channels,
              'shuffle': True}

    training_generator = DataGenerator(imagepath, partition['train'], labels, **params)
    validation_generator = DataGenerator(imagepath, partition['validation'], labels, **params)

    output = nnmodel(imglen, pixellen, classification, n_channels, n_classes, nameofclass)
    model = output[0]
    callbacks_list = output[1]
    model.fit_generator(generator=training_generator, validation_data=validation_generator, epochs=50,
                        callbacks=callbacks_list)

# python3 network.py F m/z 0.8
# python3 network.py T Length 0.8
