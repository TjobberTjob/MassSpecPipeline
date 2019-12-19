import os
import numpy as np
import pickle
with open('config.json') as json_file:
   	data = json.load(json_file)
   	
imgpath = data['path']+'images/'
imgfiles = os.listdir(imgpath)
for f in imgfiles:
	with open(imgpath+f, "rb") as pa:
		image = pickle.load(pa)

	image = np.array(image)
	print(image.shape)