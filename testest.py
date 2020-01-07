import json
from collections import Counter
import time

with open('config.json') as json_file:
	data = json.load(json_file)
path = data['path'] + 'metadata/'

start = time.time()
defdict = {}
names = [(json.loads(line)['image']+".txt") for line in open(path + 'subimage.json') if 'image' in json.loads(line)]
imgclass = [(json.loads(line)['m/z']) for line in open(path + 'subimage.json') if 'm/z' in json.loads(line)]
if len(names) == len(imgclass):
	for i in range(len(names)):
		defdict[names[i]] = imgclass[i]
end = time.time()
print(end - start)

start = time.time()
names = []
labels = {}
for line in open(path + 'subimage.json'):
	data = json.loads(line)
	name = data['image'] + ".txt"
	names.append(name)
	labels[name] = data['m/z']
end = time.time()
print(end - start)

quit()
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
