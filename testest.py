import json
from collections import Counter

with open('config.json') as json_file:
	data = json.load(json_file)
path = data['path'] + 'metadata/'

for line in open(path + 'subimage.json'):
	data = json.loads(line)
	try:
		data['Sequence']
	except:
		print(data)
		quit()
Seen = []
lendict = {}
strdict = {}
i = 0
Seen = [json.loads(line) for line in open(path + 'subimage.json')]
print(Seen)
quit()
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
