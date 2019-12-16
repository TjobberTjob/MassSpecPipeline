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
	import re
	import pandas as pd
	from collections import defaultdict
	import numpy as np
	import json
	import sys
	import random
	import pickle
	import os
	import operator

def datafetcher(path, imgpath,classification, imageclass, splitratio):
	imgfiles = os.listdir(imgpath)
	for f in imgfiles:
		print(os.path.exists(imgpath+imgfiles[0:10]))
	quit()
	with open(imgfiles[0], "rb") as pa:
		image = pickle.load(pa)
	imagelen = len(image)
	pixellen = len(image[0])

	if not classification:
		names = []
		labels = {}
		if os.path.exists(path+'subimage_filtered.json'):
			for line in open(path+'subimage_filtered.json'):
				data  = json.loads(line)
				name  = data['image']+".txt"
				names.append(name)
				labels[name] = data[imageclass]
		elif os.path.exists(path+'subimage.json'):
			for line in open(path+'subimage.json'):
				data  = json.loads(line)
				name  = data['image']+".txt"
				names.append(name)
				labels[name] = data[imageclass]
		else:
			print('No metadata for images exists')
		splits = round(len(names) * float(splitratio))
		trainlist   = names[0:splits]
		vallist 	= names[splits:]
		partition = {'train': trainlist, 'validation': vallist}
	
	else:
		labels = {}
		if os.path.exists(path+'subimage_filtered.json'):
			for line in open(path+'subimage_filtered.json'):
				data   = json.loads(line)
				name  = data['image']+".txt"
				labels[name] = data[imageclass]
		elif os.path.exists(path+'subimage.json'):
			for line in open(path+'subimage.json'):
				data   = json.loads(line)
				name  = data['image']+".txt"
				labels[name] = data[imageclass]
		else:
			print('No metadata for images exists')

		partition = {'train': [], 'validation': []}
		labels2 = defaultdict(list)
		for k, v in labels.items():
			labels2[v].append(k)
		for f in labels2:
			random.shuffle(labels2[f])
			splits = round(len(labels2[f]) * float(splitratio))
			trainlist   = (labels2[f][0:splits])
			vallist 	= (labels2[f][splits:])
			for x in trainlist:
				partition['train'].append(x)
			for x in vallist:
				partition['validation'].append(x)

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
		indexes = self.indexes[index*self.batch_size:(index+1)*self.batch_size]

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
		y = np.empty((self.batch_size), dtype = float)

		# Generate data
		for i, ID in enumerate(list_IDs_temp):
			# Store sample
			with open(imagepath+ID, "rb") as pa:
				image = pickle.load(pa)
			image = np.array(image)
			image = image[:,:,0:self.n_channels]
			X[i,] = image

			y[i] = self.labels[ID]
		if not classification:
			return X, y
		else:
			return X, keras.utils.to_categorical(y, num_classes=self.n_classes)


def nnmodel(imglen, pixellen, classification, n_channels, n_classes, imageclass):
	input = Input(shape=(imglen, pixellen, n_channels,))
	x  = Dense(128, activation = 'relu')(input)
	x  = Dense(128, activation = 'relu')(x)
	x  = Dense(128, activation = 'relu')(x)
	x  = Flatten()(x)
	x  = Dropout(rate = 0.25)(x)		
	if 	not classification:
		output = Dense(1, activation = 'linear')(x)
	else:
		output = Dense(n_classes, activation = 'softmax')(x)
	model   = keras.Model(input, output)

	if classification:
		model.compile(loss = 'categorical_crossentropy', metrics=['accuracy'], optimizer = 'adam') 
	else:
		model.compile(loss='mse', metrics=['mse'], optimizer='rmsprop') 


	#Create callbacks
	if classification:
		checkpoint  = keras.callbacks.ModelCheckpoint('Best-'+imageclass+'.h5', monitor='val_accuracy', save_best_only=True)
	else:
		checkpoint  = keras.callbacks.ModelCheckpoint('Best-'+imageclass+'.h5', monitor='val_mse', save_best_only=True)
	early_stopping = keras.callbacks.EarlyStopping(monitor='val_loss', patience=4)
	callbacks_list = [checkpoint, early_stopping]
	
	return model, callbacks_list


if __name__ == '__main__':
	#Read datapath from config file
	with open('config.json') as json_file:
   		data = json.load(json_file)
   	
	datapath = data['path']
	metapath 	= datapath+'metadata/'
	imagepath 	= datapath+'images/'

	#Cmd inputs
	classification 	= sys.argv[1] == 'T'
	imageclass 		= sys.argv[2]
	splitratio		= sys.argv[3]

	nameofclass = imageclass.replace('/','')

	output = datafetcher(metapath, imagepath, classification, imageclass, splitratio)
	partition = output[0]
	labels   = output[1]
	imglen   = output[2]
	pixellen = output[3]

	if classification:
		n_classes = len(labels)
	else:
		n_classes = 1

	n_channels = 4

	params = {'size': (pixellen,imglen),
			  'batch_size': 32,
			  'n_classes' : n_classes,
			  'n_channels': n_channels,
			  'shuffle': True}

	training_generator = DataGenerator(imagepath,partition['train'], labels, **params)
	validation_generator = DataGenerator(imagepath,partition['validation'], labels, **params)

	output 	= nnmodel(imglen, pixellen, classification, n_channels, n_classes, nameofclass)
	model 	= output[0]
	callbacks_list = output[1]
	model.fit_generator(generator=training_generator,validation_data=validation_generator, epochs = 50, callbacks = callbacks_list)

# python3 network.py F m/z 0.8
# python3 network.py T Charge 0.8
