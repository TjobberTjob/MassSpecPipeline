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
	try:
		sizedict[str(data['size'])] = str(int(sizedict[data['size']])+1)
	except Exception:
		sizedict[str(data['size'])] = 1

print(sizedict)

