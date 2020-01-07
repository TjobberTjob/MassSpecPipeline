import json
from collections import Counter
import time

with open('config.json') as json_file:
	data = json.load(json_file)
path = data['path'] + 'metadata/'

# names = [(json.loads(line)['image']+".txt") for line in open(path + 'subimage.json') if 'image' in json.loads(line)]
Seen = []
lendict = {}
for line in open(path + 'subimage.json'):
	try:
		data = json.loads(line)
		name = str(data['Sequence'])
		Seen.append(name)
		try:
			lendict[len(name)] = int(lendict[len(name)]) + 1
		except:
			lendict[len(name)] = 1
	except:
		pass
a = {}
for f in Seen:
	a[str(f)] = Seen.count(f)
print(Counter(a).most_common(4))
# print('files in '+str(i)+' different classes')
# print(lendict)
