import glob
import gzip
import json
import os
import random
import re
import sys
import time
from collections import defaultdict, Counter
import numpy as np
from simplejson import loads


def getsizeandscore(path, scorecheck):
    getsizes = [lines[re.search('\[', lines).span()[0]:re.search(']', lines).span()[1]] for lines in
                open(f'{path}subimage.json')]
    ms1size = max(set(getsizes), key=getsizes.count)

    if scorecheck[0]:
        getscores = [float(line[9:-1]) for lines in open(f'{path}subimage.json') for line in lines.split(', "') if
                     'score' in line.lower() and 'dp' not in line.lower() and str(ms1size) in lines.lower()]
        getabovehere = np.percentile(getscores, scorecheck[1])
    else:
        getabovehere = 'Fuckthis'
    return ms1size, getabovehere


def boomfilter(path, scorecheck, amountcheck, ms1size, getabovehere, filterclass, xmostfrequent, minbinary):
    binary = False
    if len(sys.argv) > 2:
        binary = True
        binaryclass = str(sys.argv[2])

    outfile = open(f'{path}subimage_filtered.json', 'w')
    seen = defaultdict(list)
    for line in open(f'{path}subimage.json'):
        jsonlist = line.split(', "')
        if scorecheck[0]:
            keys = ['image', 'ms1size', filterclass, 'Score']
        else:
            keys = ['image', 'ms1size', filterclass]
        values = [key.split('"')[-2] for key in jsonlist for part in keys if part in key and 'DP' not in key]

        if binary and scorecheck[0] and len(values) == 4:
            if values[2] != binaryclass:
                values[2] = f'not_{binaryclass}'
        elif binary and not scorecheck[0] and len(values) == 3:
            if values[2] != binaryclass:
                values[2] = f'not_{binaryclass}'

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
        if minbinary and len(mostcommon) < 2:
            for amountrange in range(100, 0, -1):
                mostcommon = [f[0] for f in Counter(amountdict).most_common(mostfrequent) if f[1] > int(amountrange) / 100 * sum(amountdict.values())]
                if len(mostcommon) > 1:
                    break
            print(f'Amountcheck changed from {amountcheck[1]}% to {amountrange}% to get binary classes')
        minamount = min([f[1] for f in Counter(amountdict).most_common(mostfrequent) if f[0] == mostcommon[-1]])
    else:
        mostcommon = [f[0] for f in Counter(amountdict).most_common(mostfrequent)]
        minamount = min([f[1] for f in Counter(amountdict).most_common(mostfrequent) if f[0] == mostcommon[-1]])
    nclasses = len(mostcommon)

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

            if binary:
                if data[filterclass] == binaryclass:
                    lookup = binaryclass
                else:
                    lookup = f'not_{binaryclass}'
            else:
                lookup = data[filterclass]

            data[f'{filterclass}_class'] = str([index for index in mostcommon].index(lookup))
            outfile.write(json.dumps(data) + '\n')
            i += 1
    print(f'{i} lines written to filtered version \n{nclasses} classes: {[f for f in mostcommon]}')
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


def testfunc(path, imgpath):
    outfile = open(f'{path}subimage2.json', 'w')
    totallen = len(list(open(f'{path}subimage.json')))
    print(totallen)
    for i, f in enumerate(open(f'{path}subimage.json')):
        if i % 10000:
            print(i, end='\r')
        data = loads(f)

        image = data['image']

        with gzip.GzipFile(f'{imgpath}{image}', 'r') as fin:
            fullinfoimage = json.loads(fin.read().decode('utf-8'))

        ms2info = fullinfoimage['ms2']
        newms2info = [[mz, int] for mz in ms2info[0] for int in ms2info[1]]
        fullinfoimage['ms2'] = newms2info
        if data['ms2size'] == "[2, 0]":
            newms2size = [0]
        else:
            newms2size = str([f for f in np.array(newms2info).shape])

        data['ms2size'] = newms2size

        with gzip.GzipFile(f'{imgpath}{image}', 'w') as fout:
            fout.write(json.dumps(fullinfoimage).encode('utf-8'))

        outfile.write(json.dumps(data) + '\n')





if __name__ == '__main__':
    # Read datapath from config file
    with open('config.json') as json_file:
        data = json.load(json_file)

    path = f'{data["path"]}metadata/'
    imgpath = f'{data["path"]}images/'
    scorecheck = data['filterscore']
    scorecheck[0] = scorecheck[0] == 'True'
    amountcheck = data['filteramount']
    amountcheck[0] = amountcheck[0] == 'True'
    topxamount = data['topxamount']
    minbinary = data['minbinary'] == 'True'

    filterclass = sys.argv[1]
    if sys.argv[1].lower() == 'test':
        testfunc(path, imgpath)

    if sys.argv[1].lower() == 'combine':
        combine(path)

    elif sys.argv[1].lower() == 'accessions':
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
        boomfilter(path, scorecheck, amountcheck, size, scorepercentile, filterclass, topxamount, minbinary)
        stop = time.time()
        print(f'Creating filtered version complete - {round(stop-start,5)} seconds elapsed')

# python3 filehandler.py filter accessions
# python3 filehandler.py filter subimage PTM/Charge/Sequence/Length
