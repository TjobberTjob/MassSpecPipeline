import json

with open('config.json') as json_file:
	data = json.load(json_file)
path = data['path'] + 'metadata/'

Seen = []
i = 0
for line in open(path + 'subimage.json'):
	data = json.loads(line)
	name = str(data['Sequence'])
	if name not in Seen:
		Seen.append(name)
		i += 1
print(Seen, i)
