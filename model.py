if 'import' == 'import':
	import keras
	import tensorboard
	#import pillow
	from keras.utils import plot_model
	from keras.preprocessing.image import ImageDataGenerator
	from keras.layers import Activation, Dropout, Flatten, Dense, Input, Concatenate, BatchNormalization, Conv2D, MaxPooling2D
	from keras import backend as K
	from keras.models import Model
	from keras.callbacks import ModelCheckpoint
	from keras import regularizers
	import glob
	import re
	import numpy as np
	import os
	import operator

#Pre-information and folderhandling
datapath = "/data/ProteomeToolsRaw/Images/"
#datapath = "Data/Images"
trainpath = datapath+'training/'
valpath = datapath+'validation/'

files = [f for f in glob.glob(trainpath + "*", recursive=True)]
classnum = {}
for f in files:
	folderclass = f[[m.start() for m in re.finditer('/', f)][-1]+1:]
	classnum.update({"class: "+str(folderclass) : len([f for f in glob.glob(trainpath+folderclass + "/*.png", recursive=True)])})
print("classes are distributed thusly: \n"+str(classnum))
print("Guessing only the most abundant would result in "+str((round(max(classnum.items(), key=operator.itemgetter(1))[1] / sum(classnum.values())*100)))+"% accuracy")

dirs = [os.path.dirname(p) for p in glob.glob(trainpath+"/*/*")]
classes = len(np.unique(dirs))

#Developing the imagegenerator
trainImageDataGen = ImageDataGenerator(rescale=1/255.)
validationImageDataGen = ImageDataGenerator(rescale=1/255.)
datapath = "Data/Images/"
trainGen = trainImageDataGen.flow_from_directory(trainpath,
	target_size=(200,200),
	batch_size=128,
	class_mode="categorical")

valGen = validationImageDataGen.flow_from_directory(valpath,
	target_size=(200,200),
	batch_size=64,
	class_mode="categorical")	

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
outputs = Dense(classes, activation = 'softmax')(x)
model   = keras.Model(inputs,outputs)

model.compile(loss = 'categorical_crossentropy', metrics=['accuracy'], optimizer = 'adam') 
#Create callbacks
checkpoint  = keras.callbacks.ModelCheckpoint("bestweight.h5", monitor='val_acc', save_best_only=True)
early_stopping = keras.callbacks.EarlyStopping(monitor='val_loss', patience=10)
callbacks_list = [checkpoint, early_stopping]


trainfiles = [f for f in glob.glob(trainpath + "/*/*", recursive=True)]
valfiles = [f for f in glob.glob(valpath + "/*/*", recursive=True)]
print(len(trainfiles))
#Run model
model.fit_generator(generator = trainGen, steps_per_epoch = len(trainfiles)//128, epochs = 50, callbacks = callbacks_list, validation_data = valGen, validation_steps = len(valfiles)//64)
