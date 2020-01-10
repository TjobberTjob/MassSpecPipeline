import json

path = '/data/ProteomeToolsRaw/metadata/'
accessions = []
i = 0
for line in open(path + 'subimage_filtered.json'):
    data = json.loads(line)
    if 'Length' in data and data['Length'] == 11:
        i += 1
print(i)
# print(accessions)
