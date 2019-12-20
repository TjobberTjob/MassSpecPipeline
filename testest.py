import os
import numpy as np
import json
import pickle
with open('config.json') as json_file:
   	data = json.load(json_file)
path = data['path']+'metadata/'
   	
sizedict = {}
for line in open(path+'subimage.json'):
	data  = json.loads(line)
	name = str(data['size'])
	if name in sizedict:
		sizedict[name] = sizedict[name]+1
	else:
		sizedict[name] = 1
print(sizedict)

