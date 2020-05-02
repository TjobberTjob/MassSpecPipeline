import numpy as np
from simplejson import loads
from keras.utils import Sequence, to_categorical
import gzip

class DataGenerator(Sequence):
    'Generates data for Keras'

    def __init__(self, path, list_IDs, labels, batch_size, ms1size, ms2size, n_channels, n_classes, shuffle, MS,
                 minMS2, classification, imagepath):
        'Initialization'
        self.imagepath = imagepath
        self.classification = classification
        self.ms1size = ms1size
        self.ms2size = ms2size
        self.path = path
        self.batch_size = batch_size
        self.labels = labels
        self.list_IDs = list_IDs
        self.n_channels = n_channels
        self.n_classes = n_classes
        self.shuffle = shuffle
        self.MS = MS
        self.minMS2 = minMS2
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
            X = np.empty((self.batch_size, self.ms1size[0], self.ms1size[1], self.n_channels))

        elif self.MS == 'ms2':
            X = np.empty((self.batch_size, self.minMS2, self.ms2size[1]))

        else:
            X = np.empty((self.batch_size, self.ms1size[0], self.ms1size[1], self.n_channels))
            X2 = np.empty((self.batch_size, self.minMS2, self.ms2size[1]))

        if self.classification:
            y = np.empty(self.batch_size, dtype=int)
        else:
            y = np.empty(self.batch_size, dtype=float)

        # Generate data
        for i, ID in enumerate(list_IDs_temp):
            # Store sample
            with gzip.GzipFile(f'{self.imagepath}{ID}', 'r') as fin:
                fullinfoimage = loads(fin.read().decode('utf-8'))
            ms1 = fullinfoimage['ms1']
            ms2 = fullinfoimage['ms2']
            if self.MS == 'ms1':
                X[i,] = np.array(ms1)[:, :, 0:self.n_channels]

            elif self.MS == 'ms2':
                X[i,] = np.array(sorted(ms2, key=lambda x: x[1], reverse=True)[:self.minMS2])

            else:
                X[i,] = np.array(ms1)[:, :, 0:self.n_channels]
                X2[i,] = np.array(sorted(ms2, key=lambda x: x[1], reverse=True)[:self.minMS2])

            y[i] = self.labels[ID]

        if self.classification:
            y = to_categorical(y, num_classes=self.n_classes)

        if self.MS == 'both':
            return [X, X2], y

        else:
            return X, y