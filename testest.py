import json
import os

with open('config.json') as json_file:
	data = json.load(json_file)
path = data['path'] + 'metadata/'

Seen = []
lendict = {}
i = 0
for line in open(path + 'subimage.json'):
	data = json.loads(line)
	name = str(data['Sequence'])
	if name not in Seen:
		Seen.append(name)
		i += 1
	try:
		lendict[len(name)] = int(lendict[len(name)]) + 1
	except:
		lendict[len(name)] = 1
print(Seen)

print(str(os.system('lsxl'))+' files in '+str(i)+' different classes')
print(lendict)
