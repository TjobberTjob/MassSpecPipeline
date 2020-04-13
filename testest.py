import json
import os
import subprocess
from collections import Counter
from glob import glob
from statistics import mean
import numpy as np



with open('config.json') as json_file:
    data = json.load(json_file)
path = f'{data["path"]}metadata/'
i = 0
lines_seen =  []
for line in open(f'{path}accessions.json'):
    data = json.loads(line)
    if data['maxquant'] and 'filetypes' in data and 'raw' in data['filetypes'] and line not in lines_seen:
        lines_seen.append(line)
        i += 1
print(i)
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
# print(f'{Counter(a).most_common(10)} Total images: {len(os.listdir(f"{path[:-9]}images/"))}')


# find /data/ProteomeToolsRaw/ -name file.mzML -exec rm -f {} \;
# find /data/ProteomeToolsRaw/ -name file-metadata.txt -exec rm -f {} \;
# find /data/ProteomeToolsRaw/ -name 1250x1000.txt -exec rm -f {} \;
