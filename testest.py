import os
import numpy as np
imgfiles = os.listdir('/data/ProteomeToolsRaw/images/')
imgshape = (165,66,4)
for f in imgfiles:
	ff = np.array(f) 
	print(ff.shape)
	quit()