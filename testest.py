import json
import os
from collections import Counter
import numpy as np

with open('config.json') as json_file:
	data = json.load(json_file)
path = data['path'] + 'metadata/'

# names, length = [(json.loads(line)['image']+".txt") for line in open(path + 'subimage.json') if 'image' in json.loads(line)]
seen = []
leng = []
for line in open(path + 'subimage.json'):
	data = json.loads(line)
	if 'Sequence' in data:
		seen.append(str(data['Sequence']))
		leng.append(len(data['Sequence']))
Seen = np.unique(seen)
Leng = np.unique(leng)
a = {}
for f in Seen:
	a[str(f)] = seen.count(f)
print(Counter(a).most_common(4))
a = {}
for f in Leng:
	a[str(f)] = leng.count(f)
print(Counter(a).most_common(4))
print('files in '+str(len(Seen))+' different classes')
