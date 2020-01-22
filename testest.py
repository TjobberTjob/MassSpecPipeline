import json
import os
import subprocess
from collections import Counter
from glob import glob
from statistics import mean

import numpy as np

# with open('config.json') as json_file:
#     data = json.load(json_file)
# path = f'{data["path"]}metadata/'
#
# seen = []
# leng = []
# for line in open(f'{path}subimage.json'):
#     data = json.loads(line)
#     if 'Sequence' in data:# and 'PIF' in data and data['PIF'] != 'Non Numérique' and float(data['PIF']) > 0.25:
#         seen.append(str(data['Sequence']))
#         leng.append(len(data['Sequence']))
#
# Seen = np.unique(seen)
# Leng = np.unique(leng)
#
# a = {}
# for f in Seen:
#     a[str(f)] = seen.count(f)
# print(f"{Counter(a).most_common(10)}  Total amount of classes: {len(Seen)}")
#
# a = {}
# for f in Leng:
#     a[str(f)] = leng.count(f)
# print(Counter(a).most_common(10))
listofshit = glob('/data/ProteomeToolsRaw/P*/*')
dictofshit = {}
for f in listofshit:
    dictofshit[f[33:]] = f[23:32]
outfile = open(f'/data/ProteomeToolsRaw/metadata/subimage2.json', 'a')
for line in open(f'/data/ProteomeToolsRaw/metadata/subimage.json', 'r'):
    data = json.loads(line)
    data['accession'] = dictofshit[data['Raw file']]
    outfile.write(json.dumps(data) + '\n')
outfile.close()


# find /data/ProteomeToolsRaw/ -name file.mzML -exec rm -f {} \;
# find /data/ProteomeToolsRaw/ -name file-metadata.txt -exec rm -f {} \;
# find /data/ProteomeToolsRaw/ -name 1250x1000.txt -exec rm -f {} \;
