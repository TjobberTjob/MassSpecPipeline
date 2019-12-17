import os
import numpy as np
import pickle
imgpath = '/data/ProteomeToolsRaw/images/'
imgfiles = os.listdir(imgpath)
imgshape = (165,66,4)
for f in imgfiles:
	with open(imgpath+f, "rb") as pa:
		image = pickle.load(pa)
	image = np.array(image)
	if image.shape != imgshape:
		print('hey') 