import json
import os
import pickle
import random
import sys
from collections import defaultdict
from itertools import chain
import keras
import matplotlib.pyplot as plt
import numpy as np
from keras.engine.saving import load_model
from keras.layers import Dropout, Dense, Input, Flatten, Conv2D, MaxPooling2D
from keras.utils import plot_model


def datafetcher(path, imgpath, classification, imageclass, splitratio, test_accessions):
    print('Data fetching')
    imgfiles = os.listdir(imgpath)
    with open(f'{imgpath}{imgfiles[0]}', "rb") as pa:
        image = pickle.load(pa)
    imagelen = len(image)
    pixellen = len(image[0])

    if not classification:  # Single output networks
        names = []
        labels = {}
        if os.path.exists(f'{path}subimage_filtered.json'):
            for line in open(f'{path}subimage_filtered.json'):
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

    else:  # Multiple output classifier
        accessions = [json.loads(line)['accession'] for line in open(f'{path}subimage_filtered.json') if
                      'accession' in json.loads(line)]
        random.shuffle(accessions)
        testlist = [f'{json.loads(line)["image"]}.txt' for line in open(f'{path}subimage_filtered.json') if
                    json.loads(line)['accession'] in accessions[0:test_accessions]]

        labels = {}
        testlabels = {}
        if os.path.exists(f'{path}subimage_filtered.json'):
            for line in open(f'{path}subimage_filtered.json'):
                data = json.loads(line)

                if f'{data["image"]}.txt' not in testlist:
                    name = f'{data["image"]}.txt'
                    labels[name] = data[imageclass]
                else:
                    name = f'{data["image"]}.txt'
                    testlabels[name] = data[imageclass]
        else:
            print('No metadata for images exists')

        # Make train and validation data
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

        # Make test data
        testlabels2 = defaultdict(list)
        for k, v in testlabels.items():
            testlabels2[v].append(k)

        for f in testlabels2:
            testlist = testlabels2[f]
            partition['test'].append(testlist)

        partition['train'] = list(chain.from_iterable(partition['train']))
        partition['validation'] = list(chain.from_iterable(partition['validation']))
        partition['test'] = list(chain.from_iterable(partition['test']))

        for f in partition:
            print(f'Datapoint in {f}: {len(partition[f])}')
        # for f in testlabels2:
        #     print(f'Test data - {f}: {len(testlabels2[f])}')

    return partition, labels, imagelen, pixellen, testlabels


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
        y = np.empty((self.batch_size), dtype=float)

        # Generate data
        for i, ID in enumerate(list_IDs_temp):
            # Store sample
            with open(f'{imagepath}{ID}', "rb") as pa:
                image = pickle.load(pa)
            image = np.array(image)
            image = image[:, :, 0:self.n_channels]
            X[i,] = image
            y[i] = self.labels[ID]

        if classification:
            y = keras.utils.to_categorical(y, num_classes=self.n_classes)
            return X, y
        else:
            return X, y


# Developing the neural network
def nnmodel(imglen, pixellen, classification, n_channels, n_classes, imageclass, metapath, patience):
    input = Input(shape=(imglen, pixellen, n_channels,))
    x = Conv2D(16, kernel_size=(5, 5), activation='relu', padding='same')(input)
    x = MaxPooling2D(pool_size=(2, 2))(x)
    # x = Concatenate()([x, x1, x2])
    x = Flatten()(x)
    x = Dropout(rate=0.25)(x)
    x = Dense(64, activation='relu')(x)
    if not classification:
        output = Dense(1, activation='linear')(x)
    else:
        output = Dense(n_classes, activation='softmax')(x)
    model = keras.Model(input, output)
    plot_model(model, to_file="model.png")

    if classification:
        if n_classes == 2:
            model.compile(loss='binary_crossentropy', metrics=['accuracy'], optimizer='adam')
        else:
            model.compile(loss='categorical_crossentropy', metrics=['accuracy'], optimizer='adam')
    else:
        model.compile(loss='mse', metrics=['mse'], optimizer='rmsprop')

    # Create callbacks
    if classification:
        checkpoint = keras.callbacks.ModelCheckpoint(f'{metapath}Best-{imageclass}.h5', monitor='val_accuracy',
                                                     save_best_only=True)
    else:
        checkpoint = keras.callbacks.ModelCheckpoint(f'{metapath}Best-{imageclass}.h5', monitor='val_mse',
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
    TFversion = config['TF']


    # Cmd inputs
    classification = sys.argv[1] == 'C'
    imageclass = sys.argv[2]
    splitratio = sys.argv[3]

    if not classification == 'C' or classification == 'R':
        print('classification or regression problem not input correctly.')
        quit()

    nameofclass = imageclass.replace('/', '')

    output = datafetcher(metapath, imagepath, classification, imageclass, splitratio, test_accessions)
    partition = output[0]
    labels = output[1]
    imglen = output[2]
    pixellen = output[3]
    testlabels = output[4]

    if classification:
        classes = [json.loads(line)[imageclass] for line in open(f'{metapath}subimage_filtered.json', 'r') if str(imageclass) in json.loads(line)]
        n_classes = len(np.unique(classes))
    else:
        n_classes = 1

    params = {'size': (pixellen, imglen),
              'batch_size': batch_size,
              'n_classes': n_classes,
              'n_channels': n_channels,
              'shuffle': True}

    training_generator = DataGenerator(imagepath, partition['train'], labels, **params)
    validation_generator = DataGenerator(imagepath, partition['validation'], labels, **params)

    output = nnmodel(imglen, pixellen, classification, n_channels, n_classes, nameofclass, metapath, patience)
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
        plt.plot(history.history['mse'])
        plt.plot(history.history['val_mse'])
        plt.title('model accuracy')
        plt.ylabel('Mean squared error')
        plt.xlabel('epoch')
        plt.legend(['train', 'validation'], loc='upper left')
        plt.savefig(f'{metapath}{imageclass}.png')

    print('Creating and running model')
    model = load_model(f'{metapath}Best-{imageclass}.h5')
    test_generator = DataGenerator(imagepath, partition['test'], testlabels, **params)
    testaccuracy = model.evaluate_generator(test_generator)
    print(f'Accuracy on test data. Loss: {testaccuracy[0]}. Accuracy: {testaccuracy[1]}')

# python3 network.py R m/z 0.8
# python3 network.py C Length 0.8
# python3 network.py C Seq_class 0.8
# python3 network.py C Modi_class 0.8
