import json
from collections import Counter

with open('config.json') as json_file:
	data = json.load(json_file)
path = data['path'] + 'metadata/'

# names, length = [(json.loads(line)['image']+".txt") for line in open(path + 'subimage.json') if 'image' in json.loads(line)]
seen = []
Seen = []
leng = []
Leng = []
for line in open(path + 'subimage.json'):
	data = json.loads(line)
	if 'Sequence' in data:
		Seen.append(str(data['Sequence']))
		leng.append(len(data['Sequence']))
		if str(data['Sequence']) not in Seen:
			Seen.append(str(data['Sequence']))
		if len(data['Sequence']) not in Leng:
			Leng.append(len(data['Sequence']))
a = {}
for f in Seen:
	a[str(f)] = seen.count(f)
print(Counter(a).most_common(4))
a = {}
for f in Leng:
	a[str(f)] = leng.count(f)
print(Counter(a).most_common(4))
# print('files in '+str(i)+' different classes')
# print(lendict)
