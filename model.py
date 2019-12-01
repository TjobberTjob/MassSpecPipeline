if 'import' == 'import':
	import keras
	import tensorflow
	import tensorboard
	from keras.utils import plot_model
	from keras.preprocessing.image import ImageDataGenerator
	from keras.layers import Activation, Dropout, Flatten, Dense, Input, Concatenate, BatchNormalization, Conv2D, MaxPooling2D
	from keras import backend as K
	from keras.models import Model
	from keras.callbacks import ModelCheckpoint
	from keras import regularizers
	import glob
	import re
	import pandas as pd
	import numpy as np
	import json
	import sys
	import random
	import pickle
	import os
	import operator

def datafetcher(path, impath,classorreg, imageclass, splitratio):
	imgfiles = glob.glob(impath+'/*')
	with open(imgfiles[0], "rb") as pa:
		image = pickle.load(pa)
	imagelen = len(image)

	if classorreg == 'regression':
		names = []
		labels = {}
		for line in open(path+'subimage_filtered.json'):
			data   = json.loads(line)
			name  = data['image']+".txt"
			names.append(name)
			labels[name] = data[imageclass]

		random.shuffle(names)
		print(len(names))
		print(splitratio)
		splits = round(len(names) * float(splitratio))
		trainlist   = names[0:splits]
		vallist 	= names[splits:]
		partition = {'train': trainlist, 'validation': vallist}
	
	elif classorreg == 'classification':
		labels = {}
		for line in open(path+'subimage_filtered.txt'):
			data   = json.loads(line)
			name  = data['image']+".txt"
			labels[name] = data[imageclass]

		partition = {'train': '', 'validation': ''}
		for f in labels:
			random.shuffle(labels[f])
			splits = round(len(labels[f])*(int(splitratio)/100))
			trainlist   = labels[0:splits]
			vallist 	= labels[splits:]
			partition[train.append(trainlist)]
			partition[validation.append(vallist)]

	return partition, labels, imagelen

# Developing the data generator
class DataGenerator(keras.utils.Sequence):
	'Generates data for Keras'
	def __init__(self, list_IDs, labels, batch_size, dim, n_channels, n_classes, shuffle):
		'Initialization'
		self.dim = dim
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
		indexes = self.indexes[index*self.batch_size:(index+1)*self.batch_size]

		# Find list of IDs
		list_IDs_temp = [self.list_IDs[k] for k in indexes]

		# Generate data
		X, y = self.__data_generation(list_IDs_temp)

		return X, y

	def on_epoch_end(self):
		'Updates indexes after each epoch'
		self.indexes = np.arange(len(self.list_IDs))
		if self.shuffle == True:
			np.random.shuffle(self.indexes)

	def __data_generation(self, list_IDs_temp):
		'Generates data containing batch_size samples' # X : (n_samples, *dim, n_channels)
		# Initialization
		X = np.empty((self.batch_size, *self.dim, self.n_channels))
		y = np.empty((self.batch_size), dtype=int)

		# Generate data
		for i, ID in enumerate(list_IDs_temp):
			# Store sample
			X[i,] = np.load('data/' + ID + '.npy')

			# Store class
			y[i] = self.labels[ID]

		return X, keras.utils.to_categorical(y, num_classes=self.n_classes)


# def nnmodel():
# 	inputs = Input(shape = (200,200,3))
# 	x  = Conv2D(16, kernel_size=(1,1), activation = 'relu')(inputs)
# 	x1 = Conv2D(16, kernel_size=(1,1), activation = 'relu')(inputs)
# 	x1 = Conv2D(8, kernel_size=(3,3), activation = 'relu', padding = "same")(x1)
# 	x2 = Conv2D(8, kernel_size=(1,1), activation = 'relu')(inputs)
# 	x2 = Conv2D(4, kernel_size=(5,5), activation = 'relu', padding = "same")(x2)
# 	x  = Concatenate()([x, x1, x2])
# 	x  = Flatten()(x)
# 	x  = Dropout(rate = 0.25)(x)
# 	outputs = Dense(64, activation = 'relu')(x)
# 	if 	classorreg == 'classify':
# 		outputs = Dense(classes, activation = 'softmax')(x)
# 	elif classorreg == 'regression':
# 		outputs = Dense(1, activation = 'linear')(x)
# 	model   = keras.Model(inputs,outputs)

# 	if 	classorreg == 'classify':
# 		model.compile(loss = 'categorical_crossentropy', metrics=['accuracy'], optimizer = 'adam') 
# 	elif classorreg == 'regression':
# 		model.compile(loss = 'mean_absolute_percentage_error', metrics=['accuracy'], optimizer = 'adam') 


	#Create callbacks
	# checkpoint  = keras.callbacks.ModelCheckpoint('bestweight.h5', monitor='val_acc', save_best_only=True)
	# early_stopping = keras.callbacks.EarlyStopping(monitor='val_loss', patience=10)
	# callbacks_list = [checkpoint, early_stopping]
	# return model, callbacks_list


if 'run' == 'run':

	classorreg = sys.argv[1]
	imageclass = sys.argv[2]
	splitratio = sys.argv[3]

	#Pre-information and folderhandling
	# datapath = '/data/ProteomeToolsRaw/images/'
	datapath = 'Data/'
	metapath = datapath+'metadata/'
	imagepath = datapath+'images/'

	output = datafetcher(metapath, imagepath, classorreg, imageclass, splitratio)
	partition = output[0]
	labels = output[1] 
	imglen = output[2]

	params = {'dim': imglen,
			  'batch_size': 16,
			  'n_classes': 1,
			  'n_channels': 4,
			  'shuffle': True}

	training_generator = DataGenerator(partition['train'], labels, **params)
	validation_generator = DataGenerator(partition['validation'], labels, **params)

	output = nnmodel()
	model = output[0]
	callbacks_list = [output[1]]
	model.fit_generator(generator=training_generator,
						validation_data=validation_generator,
						use_multiprocessing=True,
						workers=6)

#python3 model.py regression m/z 0.8
