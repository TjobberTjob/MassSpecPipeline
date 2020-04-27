import gzip
import json
import math
import os
import random
import re
import sys
from collections import defaultdict
from itertools import chain
import keras
import matplotlib.pyplot as plt
import numpy as np
from keras.engine.saving import load_model
from keras.layers import Dense, Input, Flatten, Conv2D, MaxPooling2D, Concatenate
# from keras.utils import plot_model
from simplejson import loads


def getclass(word, string):
    if word == '"size"':
        output = string[re.search('\[', string).span()[0]:re.search(']', string).span()[1]]
    else:
        ab = [f for f in [m.start() for m in re.finditer('"', string)] if f > re.search(word, string).span()[1]]
        output = string[ab[0] + 1: ab[1]]
    return output


def createnetworkfile(lenMS2):
    if whichMS == 'both' or whichMS == 'ms2':
        print('Creating MS2 data structure')

        outfile = open(f'{metapath}subimage_filtered_network.json', 'w')

        if lenMS2 == 'max':
            lenMS2 = min([json.loads(line)['ms2arraylength'] for line in open(f'{metapath}{filetosuse}')
                          if 'ms2arraylength' in json.loads(line)])

        if whichMS == 'both':
            for line in open(f'{metapath}{filetosuse}', 'r'):
                data = json.loads(line)

                if data['ms2arraylength'] >= lenMS2 and data['datacollected'] == 'both':
                    outfile.write(json.dumps(data) + '\n')
            outfile.close()
        else:
            for line in open(f'{metapath}{filetosuse}', 'r'):
                data = json.loads(line)

                if data['ms2arraylength'] >= lenMS2 and (
                        data['datacollected'] == 'ms2' or data['datacollected'] == 'both'):
                    outfile.write(json.dumps(data) + '\n')
            outfile.close()


def datafetcher(path, imgpath, imageclass, test_accessions, whichMS):
    print('Data fetching')
    if whichMS == 'both' and os.path.exists(f'{path}subimage_filtered_network.json'):
        filetouse = 'subimage_filtered_network.json'
    elif os.path.exists(f'{path}subimage_filtered.json'):
        filetouse = 'subimage_filtered.json'
    else:
        filetouse = 'subimage.json'

    for lines in open(f'{path}{filetouse}', 'r'):
        imgname = lines.split(', "')[0][11:-1]
        break
    with gzip.GzipFile(f'{imgpath}{imgname}', 'r') as fin:
        fullinfoimage = json.loads(fin.read().decode('utf-8'))
    image = fullinfoimage['ms1']
    imagelen = len(image)
    pixellen = len(image[0])

    accs = [acc.split(', "')[1][-10:-1] for acc in open(f'{path}{filetouse}', 'r')]
    accs = np.unique(accs)
    random.shuffle(accs)
    test_accs = accs[0:test_accessions]

    testfiles = []
    trainvalfiles = []
    for acc in open(f'{path}{filetosuse}'):
        accession = acc.split(', "')[1][-10:-1]
        name = acc.split(', "')[0][11:-1]
        if accession in test_accs:
            testfiles.append(name)
        else:
            trainvalfiles.append(name)

    with open('config.json') as json_file:
        config = json.load(json_file)['networkattributes']
    splitratio = config['training_percentage'] / 100
    random.shuffle(trainvalfiles)
    splits = round(len(trainvalfiles) * float(splitratio))
    trainfiles = trainvalfiles[:splits]
    validationfiles = trainvalfiles[splits:]

    partition = {'train': trainfiles, 'validation': validationfiles, 'test': testfiles}

    tests = defaultdict(list)
    for f in testfiles:
        tests[f.split('-')[-1][:-5]].append(f)

    labels = {}
    testlabels = {}
    for line in open(f'{path}{filetosuse}'):
        data = loads(line)
        name = data['image']
        label = data[imageclass]

        if name in tests[name.split('-')[-1][:-5]]:
            testlabels[name] = label
        else:
            labels[name] = label

    for f in partition:
        print(f'Datapoint in {f}: {len(partition[f])}')

    return partition, labels, imagelen, pixellen, testlabels


# Developing the data generator
class DataGenerator(keras.utils.Sequence):
    'Generates data for Keras'

    def __init__(self, path, list_IDs, labels, batch_size, size, n_channels, n_classes, shuffle, MS, mslen):
        'Initialization'
        self.size = size
        self.path = path
        self.batch_size = batch_size
        self.labels = labels
        self.list_IDs = list_IDs
        self.n_channels = n_channels
        self.n_classes = n_classes
        self.shuffle = shuffle
        self.MS = MS
        self.mslen = mslen
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
        if self.MS == 'ms1':
            X = np.empty((self.batch_size, self.size[1], self.size[0], self.n_channels))

        elif self.MS == 'ms2':
            X = np.empty((self.batch_size, self.mslen))

        else:
            X = np.empty((self.batch_size, self.size[1], self.size[0], self.n_channels))
            X2 = np.empty((self.batch_size, self.mslen))

        y = np.empty(self.batch_size, dtype=float)

        # Generate data
        for i, ID in enumerate(list_IDs_temp):
            # Store sample
            with gzip.GzipFile(f'{imagepath}{ID}', 'r') as fin:
                fullinfoimage = json.loads(fin.read().decode('utf-8'))
            ms1 = fullinfoimage['ms1']
            ms2 = fullinfoimage['ms2']
            if self.MS == 'ms1':
                image = np.array(ms1)
                image = image[:, :, 0:self.n_channels]
                X[i,] = image

            elif self.MS == 'ms2':
                mz_array = ms2[0]
                # rt_array = ms2[1]
                rt_array = [math.log(intval) for intval in ms2[1]]
                X[i,] = list(chain.from_iterable([mz_array, rt_array]))

            else:
                image = np.array(ms1)
                image = image[:, :, 0:self.n_channels]
                X[i,] = image

                mz_array = ms2[0]
                # rt_array = ms2[1]
                rt_array = [math.log(intval) for intval in ms2[1]]
                X2[i,] = list(chain.from_iterable([mz_array, rt_array]))

            y[i] = self.labels[ID]

        if classification:
            try:
                y = keras.utils.to_categorical(y, num_classes=self.n_classes)
            except:
                print(y)
                quit()

        if self.MS == 'both':
            return [X, X2], y

        else:
            return X, y


def r2(y_true, y_pred):
    from keras import backend as K
    SS_res = K.sum(K.square(y_true - y_pred))
    SS_tot = K.sum(K.square(y_true - K.mean(y_true)))
    return 1 - SS_res / (SS_tot + K.epsilon())


# Developing the neural network
def nnmodel(imglen, pixellen, classification, n_channels, n_classes, imageclass, metapath, patience, whichMS, lenMS2):
    # HIDDEN LAYERS
    if whichMS == 'ms1':  # MS1
        input = Input(shape=(imglen, pixellen, n_channels,))
        x = Conv2D(124, kernel_size=(5, 5), activation='relu', padding='same')(input)
        x = MaxPooling2D(pool_size=(2, 2))(x)
        x1 = Conv2D(124, kernel_size=(3, 3), activation='relu', padding='same')(input)
        x1 = MaxPooling2D(pool_size=(2, 2))(x1)
        x2 = Conv2D(124, kernel_size=(2, 2), activation='relu', padding='same')(input)
        x2 = MaxPooling2D(pool_size=(2, 2))(x2)
        x = Concatenate()([x, x1, x2])
        x = Flatten()(x)
        x = Dense(64, activation='relu')(x)
        x = Dense(32, activation='relu')(x)

    elif whichMS == 'ms2':  # MS2
        input = Input(shape=(lenMS2 * 2,))
        x = Dense(128, activation='relu')(input)
        x = Dense(64, activation='relu')(x)
        x = Dense(32, activation='relu')(x)

    else:  # BOTH
        input = Input(shape=(imglen, pixellen, n_channels,))  # MS1
        x = Conv2D(124, kernel_size=(5, 5), activation='relu', padding='same')(input)
        x = MaxPooling2D(pool_size=(2, 2))(x)
        x1 = Conv2D(124, kernel_size=(3, 3), activation='relu', padding='same')(input)
        x1 = MaxPooling2D(pool_size=(2, 2))(x1)
        x2 = Conv2D(124, kernel_size=(2, 2), activation='relu', padding='same')(input)
        x2 = MaxPooling2D(pool_size=(2, 2))(x2)
        x = Concatenate()([x, x1, x2])
        x = Flatten()(x)

        input2 = Input(shape=(lenMS2,))  # MS2
        x1 = Dense(128, activation='relu')(input2)
        x1 = Dense(64, activation='relu')(x1)
        x1 = Dense(32, activation='relu')(x1)

        x = Concatenate()([x, x1])  # Combine
        x = Dense(64, activation='relu')(x)
        x = Dense(32, activation='relu')(x)

    if whichMS == 'ms1' or whichMS == 'both':  # After combining
        x = Dense(128, activation='relu')(x)
        x = Dense(64, activation='relu')(x)
    # OUTPUT LAYER
    if classification:
        if n_classes == 2:
            output = Dense(1, activation='sigmoid')(x)
        else:
            output = Dense(n_classes, activation='linear')(x)
    else:
        output = Dense(n_classes, activation='softmax')(x)
    # COMBINE MODEL
    if whichMS == 'ms1' or whichMS == 'ms2':
        model = keras.Model(inputs=input, outputs=output)
    else:
        model = keras.Model(inputs=[input, input2], outputs=output)
    # print(model.summary())
    # plot_model(model, to_file="model.png")

    if classification:
        if n_classes == 2:
            model.compile(loss='binary_crossentropy', metrics=['accuracy'], optimizer='adam')
        else:
            model.compile(loss='categorical_crossentropy', metrics=['accuracy'], optimizer='adam')
    else:
        model.compile(loss='mse', metrics=[r2], optimizer='rmsprop')

    # Create callbacks
    if classification:
        checkpoint = keras.callbacks.ModelCheckpoint(f'{metapath}Best-{imageclass}.h5', monitor='val_accuracy',
                                                     save_best_only=True)
    else:
        checkpoint = keras.callbacks.ModelCheckpoint(f'{metapath}Best-{imageclass}.h5', monitor='val_r2',
                                                     save_best_only=True)
    early_stopping = keras.callbacks.EarlyStopping(monitor='val_loss', patience=patience)
    callbacks_list = [checkpoint, early_stopping]

    return model, callbacks_list


if __name__ == '__main__':
    # Read datapath from config file
    with open('config.json') as json_file:
        data = json.load(json_file)

    datapath = data['path']
    metapath = f'{datapath}metadata/'
    imagepath = f'{datapath}images/'

    with open('config.json') as json_file:
        config = json.load(json_file)['networkattributes']
    test_accessions = config['test_accessions']
    n_channels = config['n_channels']
    batch_size = config['batch_size']
    epochs = config['epochs']
    patience = config['early_stopping']
    setseed = config['setseed'] == 'True'
    whichMS = config['MS']
    lenMS2 = config['lenms2']

    if os.path.exists(f'{metapath}subimage_filtered.json'):
        filetosuse = 'subimage_filtered.json'
    elif os.path.exists(f'{metapath}subimage_filtered.json'):
        filetosuse = 'subimage.json'
    else:
        print('No datafile exists')
        quit()

    # Creating network files
    createnetworkfile(lenMS2)

    if setseed:
        random.seed(1)

    classification = sys.argv[1].lower()
    if not (classification == 'c' or classification == 'r'):
        print('classification or regression problem not input correctly.')
        quit()
    imageclass = sys.argv[2]
    classification = classification == 'c'

    nameofclass = imageclass.replace('/', '')
    output = datafetcher(metapath, imagepath, imageclass, test_accessions, whichMS)
    partition = output[0]
    labels = output[1]
    imglen = output[2]
    pixellen = output[3]
    testlabels = output[4]

    if classification:
        classes = [json.loads(line)[imageclass] for line in open(f'{metapath}subimage_filtered.json', 'r') if
                   str(imageclass) in json.loads(line)]
        n_classes = len(np.unique(classes))
    else:
        n_classes = 1
    print(f'# of classes: {n_classes}')

    params = {'size': (pixellen, imglen),
              'batch_size': batch_size,
              'n_classes': n_classes,
              'n_channels': n_channels,
              'shuffle': True,
              'MS': whichMS,
              'mslen': lenMS2}

    training_generator = DataGenerator(imagepath, partition['train'], labels, **params)
    validation_generator = DataGenerator(imagepath, partition['validation'], labels, **params)

    output = nnmodel(imglen, pixellen, classification, n_channels, n_classes, nameofclass, metapath, patience, whichMS,
                     lenMS2)
    model = output[0]
    callbacks_list = output[1]
    history = model.fit_generator(generator=training_generator, validation_data=validation_generator, epochs=epochs,
                                  callbacks=callbacks_list)

    if classification:
        plt.plot(history.history['accuracy'])
        plt.plot(history.history['val_accuracy'])
        plt.title('model accuracy')
        plt.ylabel('accuracy')
        plt.xlabel('epoch')
        plt.legend(['train', 'validation'], loc='upper left')
        plt.savefig(f'{metapath}{imageclass}.png')
    else:
        plt.plot(history.history['loss'])
        plt.plot(history.history['val_loss'])
        plt.title('MSE over time')
        plt.ylabel('MSE')
        plt.xlabel('epoch')
        plt.legend(['train', 'validation'], loc='upper left')
        plt.savefig(f'{metapath}{imageclass}.png')

    print('Creating and running model')
    model = load_model(f'{metapath}Best-{imageclass}.h5')
    test_generator = DataGenerator(imagepath, partition['test'], testlabels, **params)
    testaccuracy = model.evaluate_generator(test_generator)
    if classification:
        print(f'Accuracy on test data. Loss: {testaccuracy[0]}. Accuracy: {testaccuracy[1]}')
    else:
        print(f'Accuracy on test data. Loss: {testaccuracy[0]}. R^2: {testaccuracy[1]}')

# python3 network.py r m/z
# python3 network.py c Length_class
# python3 network.py c Seq_class
# python3 network.py c Modi_class
