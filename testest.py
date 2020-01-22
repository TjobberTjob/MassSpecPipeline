import json
from collections import Counter
from statistics import mean

import numpy as np

with open('config.json') as json_file:
    data = json.load(json_file)
path = f'{data["path"]}metadata/'

scores = []
for line in open(f'{path}subimage.json'):
    try:
        if 'PIF' in json.loads(line):
            scores.append(float(json.loads(line)['PIF']))
    except:
        pass
meanpif = mean(scores)

seen = []
leng = []
for line in open(f'{path}subimage.json'):
    data = json.loads(line)
    if 'Sequence' in data and 'PIF' in data and data['PIF'] > meanpif:
        seen.append(str(data['Sequence']))
        leng.append(len(data['Sequence']))

Seen = np.unique(seen)
Leng = np.unique(leng)

a = {}
for f in Seen:
    a[str(f)] = seen.count(f)
print(f"{Counter(a).most_common(10)}  Total amount of classes: {len(Seen)}")

a = {}
for f in Leng:
    a[str(f)] = leng.count(f)
print(Counter(a).most_common(10))


# find /data/ProteomeToolsRaw/ -name file.mzML -exec rm -f {} \;
# find /data/ProteomeToolsRaw/ -name file-metadata.txt -exec rm -f {} \;
# find /data/ProteomeToolsRaw/ -name 1250x1000.txt -exec rm -f {} \;
