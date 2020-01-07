import json
from collections import Counter, defaultdict

with open('config.json') as json_file:
	data = json.load(json_file)
path = data['path'] + 'metadata/'

# names = [(json.loads(line)['image']+".txt") for line in open(path + 'subimage.json') if 'image' in json.loads(line)]
Seen = []
for line in open(path + 'subimage.json'):
	try:
		data = json.loads(line)
		name = str(data['Sequence'])
		Seen.append(name)
		lendict = defaultdict(len(name): (int(lendict[len(name)]) + 1))
a = {}
for f in Seen:
	a[str(f)] = Seen.count(f)
print(Counter(a).most_common(4))
# print('files in '+str(i)+' different classes')
print(lendict)
