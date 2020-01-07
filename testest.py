import json
from collections import Counter, defaultdict

with open('config.json') as json_file:
	data = json.load(json_file)
path = data['path'] + 'metadata/'

Seen = []
lendict = defaultdict(list)
# Seen = [json.loads(line)['Sequence'] for line in open(path + 'subimage.json') if 'Sequence' in json.loads(line)]

for line in open(path + 'subimage.json'):
	try:
		data = json.loads(line)
		name = str(data['Sequence'])
		Seen.append(name)
		lendict[len(name)] = lendict[len(name)] + 1
	except:
		pass
a = {}
for f in Seen:
	a[str(f)] = Seen.count(f)
print(Counter(a).most_common(4))
# print('files in '+str(i)+' different classes')
print(lendict)
