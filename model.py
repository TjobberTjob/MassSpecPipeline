import keras
from keras.utils import plot_model
from keras.preprocessing.image import ImageDataGenerator
from keras.layers import Activation, Dropout, Flatten, Dense, Input, Concatenate, BatchNormalization, Conv2D, MaxPooling2D
from keras import backend as K
from keras.models import Model
from keras.callbacks import ModelCheckpoint
from keras import regularizers
import glob
import re
import operator

trainpath = "Data/Images/training/"
valpath = "Data/Images/validation/"
files = [f for f in glob.glob(trainpath + "*", recursive=True)]
classnum = {}
for f in files:
	folderclass = f[[m.start() for m in re.finditer('/', f)][-1]+1:]
	classnum.update({"class: "+str(folderclass) : len([f for f in glob.glob(trainpath+folderclass + "/*.png", recursive=True)])})
#print("classes are distributed thusly: \n"+str(classnum))
#print("Guessing only the most abundant would result in "+str((round(max(classnum.items(), key=operator.itemgetter(1))[1] / sum(classnum.values())*100,2)))+"% accuracy")

dirs = [os.path.dirname(p) for p in glob.glob(trainpath+"/*/*")]
classes = [] 
for x in dirs:  
	if x not in udirs: 
		classes.append(x)
classes = len(classes)
print(classes)
quit()

trainImageDataGen = ImageDataGenerator(rescale=1/255.)
validationImageDataGen = ImageDataGenerator(rescale=1/255.)
datapath = "Data/Images/"
trainGen = trainImageDataGen.flow_from_directory(trainpath,
	target_size=(300,300),
	batch_size=16,
	class_mode="categorical")

valGen = validationImageDataGen.flow_from_directory(valpath,
	target_size=(300,300),
	batch_size=16,
	class_mode="categorical")	

#Create model
inputs = Input(shape = (300,300,3))
x  = Conv2D(16, kernel_size=(1,1), activation = 'relu')(inputs)
x1 = Conv2D(16, kernel_size=(1,1), activation = 'relu')(inputs)
x1 = Conv2D(8, kernel_size=(3,3), activation = 'relu', padding = "same")(x1)
x2 = Conv2D(8, kernel_size=(1,1), activation = 'relu')(inputs)
x2 = Conv2D(4, kernel_size=(5,5), activation = 'relu', padding = "same")(x2)
x  = Concatenate()([x, x1, x2])
x  = Flatten()(x)
x  = Dropout(rate = 0.25)(x)
outputs = Dense(3, activation = 'softmax')(x)
model   = keras.Model(inputs,outputs)

model.compile(loss = 'categorical_crossentropy', metrics=['accuracy'], optimizer = 'adam') 
#Create callbacks
# checkpoint  = keras.callbacks.ModelCheckpoint("bestweight_C.h5", monitor='val_acc', save_best_only=True)
# tensorboard = keras.callbacks.TensorBoard('logs', update_freq='batch')
early_stopping = keras.callbacks.EarlyStopping(monitor='val_loss', patience=4)
callbacks_list = [checkpoint, tensorboard, early_stopping]

#Run model
model.fit_generator(generator = trainGen, steps_per_epoch = 100, epochs = 50, callbacks = callbacks_list, validation_data = valGen, validation_steps = 100)

# model.save('model_C_run.h5')