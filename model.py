if 'import' == 'import':
	import keras
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
	import sys
	import os
	import operator

classorreg = sys.argv[1]

#Pre-information and folderhandling
# datapath = '/data/ProteomeToolsRaw/images/'
datapath = "Data/images/"
trainpath = datapath+'training/'
valpath = datapath+'validation/'

files = [f for f in glob.glob(trainpath + '*', recursive=True)]
classnum = {}
for f in files:
	folderclass = f[[m.start() for m in re.finditer('/', f)][-1]+1:]
	classnum.update({"class: "+str(folderclass) : len([f for f in glob.glob(trainpath+folderclass + '/*.png', recursive=True)])})
# print('classes are distributed thusly: \n'+str(classnum))
# print('Guessing only the most abundant would result in '+str((round(max(classnum.items(), key=operator.itemgetter(1))[1] / sum(classnum.values())*100)))+"% accuracy")

dirs = [os.path.dirname(p) for p in glob.glob(trainpath+'/*/*')]
classes = len(np.unique(dirs))

#Developing the imagegenerator
trainImageDataGen = ImageDataGenerator(rescale=1/255.)
validationImageDataGen = ImageDataGenerator(rescale=1/255.)
datapath = 'Data/Images/'
if 	classorreg == 'classify':
	trainGen = trainImageDataGen.flow_from_directory(trainpath,
	target_size=(200,200),
	batch_size=128,
	class_mode='categorical')

	valGen = validationImageDataGen.flow_from_directory(valpath,
	target_size=(200,200),
	batch_size=64,
	class_mode='categorical')

elif classorreg == 'regression':
	traindata = pd.read_pickle(trainpath+'data.txt')
	trainGen = trainImageDataGen.flow_from_dataframe(dataframe=traindata, directory = trainpath,
	x_col='image', y_col='class', class_mode="None",
	target_size=(200,200), batch_size=128)

	valdata = pd.read_pickle(valpath+'data.txt')
	valGen = validationImageDataGen.flow_from_dataframe(dataframe=valdata, directory = trainpath,
	x_col='image', y_col='class', class_mode="None",
	target_size=(200,200), batch_size=128)
	

#Create model
inputs = Input(shape = (200,200,3))
x  = Conv2D(16, kernel_size=(1,1), activation = 'relu')(inputs)
x1 = Conv2D(16, kernel_size=(1,1), activation = 'relu')(inputs)
x1 = Conv2D(8, kernel_size=(3,3), activation = 'relu', padding = "same")(x1)
x2 = Conv2D(8, kernel_size=(1,1), activation = 'relu')(inputs)
x2 = Conv2D(4, kernel_size=(5,5), activation = 'relu', padding = "same")(x2)
x  = Concatenate()([x, x1, x2])
x  = Flatten()(x)
x  = Dropout(rate = 0.25)(x)
outputs = Dense(64, activation = 'relu')(x)
if 	classorreg == 'classify':
	outputs = Dense(classes, activation = 'softmax')(x)
elif classorreg == 'regression':
	outputs = Dense(1, activation = 'linear')(x)
model   = keras.Model(inputs,outputs)

if 	classorreg == 'classify':
	model.compile(loss = 'categorical_crossentropy', metrics=['accuracy'], optimizer = 'adam') 
elif classorreg == 'regression':
	model.compile(loss = 'mean_absolute_percentage_error', metrics=['accuracy'], optimizer = 'adam') 


#Create callbacks
checkpoint  = keras.callbacks.ModelCheckpoint('bestweight.h5', monitor='val_acc', save_best_only=True)
early_stopping = keras.callbacks.EarlyStopping(monitor='val_loss', patience=10)
callbacks_list = [checkpoint, early_stopping]


trainfiles = [f for f in glob.glob(trainpath + '/*/*', recursive=True)]
valfiles = [f for f in glob.glob(valpath + '/*/*', recursive=True)]
#Run model
model.fit_generator(generator = trainGen, steps_per_epoch = len(trainfiles)//128, epochs = 50, callbacks = callbacks_list, validation_data = valGen, validation_steps = len(valfiles)//64)


