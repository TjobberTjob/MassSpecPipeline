import glob
import json
import os
import re
import time
import numpy as np
from itertools import chain
import random
import glob
import sys
from collections import defaultdict, Counter
import pandas as pd
from simplejson import loads


def getsizeandscore(path, scorecheck):
    getsizes = [lines[re.search('\[', lines).span()[0]:re.search(']', lines).span()[1]] for lines in
                open(f'{path}subimage.json')]
    uniquesizes = np.unique(getsizes)

    sizedict = defaultdict()
    for sizes in uniquesizes:
        sizedict[str(sizes)] = getsizes.count(sizes)

    for sizes in Counter(sizedict).most_common():
        size = str(sizes[0]).replace(" ", "")[1:-1].split(',')
        if len(size) == 3:
            ms1size = sizes[0]
            break

    if scorecheck[0]:
        getscores = [float(line[9:-1]) for lines in open(f'{path}subimage.json') for line in lines.split(', "') if
                     'score' in line.lower() and 'dp' not in line.lower() and str(ms1size) in lines.lower()]
        getabovehere = np.percentile(getscores, scorecheck[1])
    else:
        getabovehere = 'Fuckthis'
    return ms1size, getabovehere


def boomfilter(path, scorecheck, amountcheck, ms1size, getabovehere, filterclass, xmostfrequent):
    outfile = open(f'{path}subimage_filtered.json', 'w')
    seen = defaultdict(list)
    for line in open(f'{path}subimage.json'):
        jsonlist = line.split(', "')
        if scorecheck[0]:
            keys = ['image', 'size', filterclass, 'Score']
        else:
            keys = ['image', 'size', filterclass]
        values = [key.split('"')[-2] for key in jsonlist for part in keys if part in key and 'DP' not in key]

        if scorecheck[0] and len(values) == 4 and values[1] == ms1size and float(values[3]) > getabovehere:
            seen[values[2]].append(values[0])
        elif not scorecheck[0] and len(values) == 3 and values[1] == ms1size:
            seen[values[2]].append(values[0])

    amountdict = defaultdict()
    for classes in seen.keys():
        amountdict[classes] = len(seen[classes])

    if xmostfrequent == 'max':
        mostfrequent = len(seen.keys())
    else:
        mostfrequent = xmostfrequent

    if amountcheck[0]:
        mostcommon = [f[0] for f in Counter(amountdict).most_common(mostfrequent) if f[1] > int(amountcheck[1] / 100 * sum(amountdict.values()))]
        minamount = min([f[1] for f in Counter(amountdict).most_common(mostfrequent) if f[1] > int(amountcheck[1] / 100 * sum(amountdict.values()))])
    else:
        mostcommon = [f[0] for f in Counter(amountdict).most_common(mostfrequent)]
        minamount = min([f[1] for f in Counter(amountdict).most_common(mostfrequent) if f[0] == mostcommon[-1]])

    namelist = []
    for seq in seen:
        if seq in mostcommon:
            random.shuffle(seen[seq])
            for names in seen[seq][0:minamount]:
                namelist.append(names)

    quickcheckdict = defaultdict(list)
    for names in namelist:
        quickcheckdict[names.split('-')[-1][:-5]].append(names)

    i = 0
    for line in open(f'{path}subimage.json', 'r'):
        name = line.split(', "')[0][11:-1]

        if name in quickcheckdict[name.split('-')[-1][:-5]]:
            data = loads(line)
            data[f'{filterclass}_class'] = str([str(index) for index in mostcommon].index(data[filterclass]))
            outfile.write(json.dumps(data) + '\n')
            i += 1
    print(f'{i} lines written to filtered version')
    outfile.close()


def boomaccessions(path):
    outfile = open(f'{path}accessions_filtered.json', 'w')

    for line in open(f'{path}accessions.json', 'r'):
        data = json.loads(line)
        if 'allpeptides' in data and data['allpeptides'] and 'filetypes' in data and 'raw' in data['filetypes']:
            outfile.write(line)
            outfile.write(json.dumps(data) + '\n')
    outfile.close()


def combine(path):
    if os.path.exists(f'{path}subimage.json'):
        allimgs = [json.loads(line)['image'] for line in open(f'{path}subimage.json') if
                   'image' in json.loads(line)]
    else:
        allimgs = []

    with open(f'{path}subimage.json', 'a') as outfile:
        for imagejson in glob.glob(f'{path}subimage-*.json'):
            for line in open(imagejson, 'r'):
                data = json.loads(line)
                if data['image'] not in allimgs:
                    outfile.write(json.dumps(data) + '\n')
            os.remove(imagejson)
    quit()



if __name__ == '__main__':
    # Read datapath from config file
    with open('config.json') as json_file:
        data = json.load(json_file)

    path = f'{data["path"]}metadata/'
    scorecheck = data['filterscore']
    scorecheck[0] = scorecheck[0] == 'True'
    amountcheck = data['filteramount']
    amountcheck[0] = amountcheck[0] == 'True'

    if sys.argv[1] == 'combine':
        combine(path)

    elif sys.argv[1] == 'accessions':
        boomaccessions(path)

    else:
        start = time.time()
        print('Getting sizes and scores', end = '\r')
        output = getsizeandscore(path, scorecheck)
        size = output[0]
        scorepercentile = output[1]
        stop = time.time()
        print(f'Getting sizes and scores - {round(stop-start,5)} seconds elapsed')

        start = time.time()
        print('Creating filtered version', end='\r')
        boomfilter(path, scorecheck, amountcheck, size, scorepercentile, sys.argv[1], 10)
        stop = time.time()
        print(f'Creating filtered version complete - {round(stop-start,5)} seconds elapsed')

# python3 filehandler.py filter accessions
# python3 filehandler.py filter subimage PTM/Charge/Sequence/Length
