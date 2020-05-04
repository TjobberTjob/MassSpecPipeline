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


def get_size_and_score(path, add_score_filter):
    getsizes = [lines[re.search('\[', lines).span()[0]:re.search(']', lines).span()[1]] for lines in
                open(f'{path}subimage.json')]
    ms1size = max(set(getsizes), key=getsizes.count)

    if add_score_filter[0]:
        getscores = [float(line[9:-1]) for lines in open(f'{path}subimage.json') for line in lines.split(', "') if
                     'score' in line.lower() and 'dp' not in line.lower() and str(ms1size) in lines.lower()]
        getabovehere = np.percentile(getscores, add_score_filter[1])
    else:
        getabovehere = []
    return ms1size, getabovehere


def subimage_filter_r(path, add_score_filter, ms1size, getabovehere, filterclass):
    outfile = open(f'{path}subimage_filtered.json', 'w')
    i = 0
    for line in open(f'{path}subimage.json'):
        jsonlist = line.lower()[2:].split(', "')
        if add_score_filter[0]:
            lookupkeys = ['image', 'ms1size', filterclass, 'score']
        else:
            lookupkeys = ['image', 'ms1size', filterclass]
        values = [key.split('"')[-2] for key in jsonlist if key[0:re.search('"', key).span()[0]] in lookupkeys]

        if add_score_filter[0] and len(values) == 4 and values[1] == ms1size and float(values[3]) > getabovehere:
            data = loads(line)
            outfile.write(json.dumps(data) + '\n')
            i += 1
        elif not add_score_filter[0] and len(values) == 3 and values[1] == ms1size:
            data = loads(line)
            outfile.write(json.dumps(data) + '\n')
            i += 1

    print(f'{i} lines written to filtered version')


def subimage_filter_c(path, add_score_filter, ms1size, getabovehere, filterclass, min_amount_classes, max_amount_classes, min_amount_in_class):
    binary = False
    if len(sys.argv) > 2:
        binary = True
        binaryclass = str(sys.argv[2])

    outfile = open(f'{path}subimage_filtered.json', 'w')
    seen = defaultdict(list)
    for line in open(f'{path}subimage.json'):
        jsonlist = line[2:].split(', "')
        if add_score_filter[0]:
            lookupkeys = ['image', 'ms1size', filterclass, 'score']
        else:
            lookupkeys = ['image', 'ms1size', filterclass]
        values = [key.split('"')[-2] for key in jsonlist if key.lower()[0:re.search('"', key).span()[0]] in lookupkeys]

        if binary and add_score_filter[0] and len(values) == 4:
            if values[2] != binaryclass:
                values[2] = f'not_{binaryclass}'
        elif binary and not add_score_filter[0] and len(values) == 3:
            if values[2] != binaryclass:
                values[2] = f'not_{binaryclass}'

        if add_score_filter[0] and len(values) == 4 and values[1] == ms1size and float(values[3]) > getabovehere:
            seen[values[2]].append(values[0])
        elif not add_score_filter[0] and len(values) == 3 and values[1] == ms1size:
            seen[values[2]].append(values[0])

    amountdict = defaultdict()
    for classes in seen.keys():
        amountdict[classes] = len(seen[classes])

    if max_amount_classes == 'max':
        mostfrequent = len(seen.keys())
    else:
        mostfrequent = int(max_amount_classes)

    mostcommon = [f[0] for f in Counter(amountdict).most_common(mostfrequent) if f[1] > int(min_amount_in_class / 100 * sum(amountdict.values()))]
    if len(mostcommon) < min_amount_classes:
        for amountrange in range(min_amount_in_class, 0, -1):
            mostcommon = [f[0] for f in Counter(amountdict).most_common(mostfrequent) if f[1] > int(amountrange) / 100 * sum(amountdict.values())]
            if len(mostcommon) >= min_amount_classes:
                print(f'min_amount_in_class changed from {min_amount_in_class}% to {amountrange}% to get {len(mostcommon)} classes')
                break
        if len(mostcommon) < min_amount_classes:
            print(f'min_amount_classes isnt achieveable (Probably not enough classes in data). Only {len(mostcommon)} classes created')
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

            if binary:
                if data[filterclass] == binaryclass:
                    lookup = binaryclass
                else:
                    lookup = f'not_{binaryclass}'
            else:
                lookup = data[filterclass.capitalize()]

            data[f'{filterclass}_class'] = str([index for index in mostcommon].index(lookup))
            outfile.write(json.dumps(data) + '\n')
            i += 1
    print(f'{i} lines written to filtered version \n{len(mostcommon)} classes: {[f for f in mostcommon]}')
    outfile.close()


def accession_filter(path):
    outfile = open(f'{path}accessions_filtered.json', 'w')

    for line in open(f'{path}accessions.json', 'r'):
        data = json.loads(line)
        if 'allpeptides' in data and data['allpeptides'] and 'filetypes' in data and 'raw' in data['filetypes']:
            outfile.write(line)
            outfile.write(json.dumps(data) + '\n')
    outfile.close()


def combine_metadata(path):
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


def test_function(path, imgpath):
    for i in range(1000000):
        if i % 10000 == 0:
            print(i)
            a = 0

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
    add_score_filter = data['add_score_filter']
    add_score_filter[0] = add_score_filter[0] == 'True'
    min_amount_classes = data['min_amount_classes']
    max_amount_classes = data['max_amount_classes']
    min_amount_in_class = data['min_amount_in_class']

    filterclass = sys.argv[1].lower()
    if sys.argv[1].lower() == 'test':
        test_function(path, imgpath)

    if sys.argv[1].lower() == 'combine':
        combine_metadata(path)

    elif sys.argv[1].lower() == 'accessions':
        accession_filter(path)

    else:
        start = time.time()
        print('Getting sizes and scores', end='\r')
        output = get_size_and_score(path, add_score_filter)
        size = output[0]
        scorepercentile = output[1]
        stop = time.time()
        print(f'Getting sizes and scores - {round(stop - start, 5)} seconds elapsed')

        start = time.time()
        print('Creating filtered version', end='\r')
        if filterclass == 'm/z' or filterclass == 'Score':
            subimage_filter_r(path, add_score_filter, size, scorepercentile, filterclass)
        else:
            subimage_filter_c(path, add_score_filter, size, scorepercentile, filterclass, min_amount_classes, max_amount_classes, min_amount_in_class)
        stop = time.time()
        print(f'Creating filtered version complete - {round(stop - start, 5)} seconds elapsed')

# python3 filehandler.py filter accessions
# python3 filehandler.py filter subimage PTM/Charge/Sequence/Length
