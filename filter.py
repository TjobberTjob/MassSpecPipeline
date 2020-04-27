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


def getclass(word, string):
    if word == '"size"':
        output = string[re.search('\[', string).span()[0]:re.search(']', string).span()[1]]
    else:
        ab = [f for f in [m.start() for m in re.finditer('"', string)] if f > re.search(word, string).span()[1]]
        output = string[ab[0] + 1: ab[1]]
    return output


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


def clear(path):
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

    outfile = open(f'{path}subimage2.json', 'a')
    for line in open(f'{path}subimage.json', 'a'):
        data = json.loads(line)
        if data['size'] == ms1size:
            outfile.write(json.dumps(line) + '\n')
        else:
            os.remove(f'{datapath}/images/{data["image"]}')
    os.remove(f'{path}subimage.json')
    os.rename(f'{path}subimage2.json', f'{path}subimage.json')


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


def filtersequence(path, outfile, getabovehere, ms1size, scorecheck, amountcheck):
    mostfrequent = 10

    seen = defaultdict(list)
    for line in open(f'{path}subimage.json'):
        jsonlist = line.split(', "')
        if scorecheck[0]:
            keys = ['image', 'size', 'Sequence', 'Score']
        else:
            keys = ['image', 'size', 'Sequence']
        values = [key.split('"')[-2] for key in jsonlist for part in keys if part in key and 'DP' not in key]

        if scorecheck[0] and len(values) == 4 and values[1] == ms1size and float(values[3]) > getabovehere:
            seen[values[2]].append(values[0])
        elif not scorecheck[0] and len(values) == 3 and values[1] == ms1size:
            seen[values[2]].append(values[0])

    uniquesequences = [f for f in seen.keys()]
    amountdict = defaultdict()
    for f in uniquesequences:
        amountdict[f] = len(seen[f])
    if amountcheck[0]:
        mostcommon = [f[0] for f in Counter(amountdict).most_common(mostfrequent) if f[1] > amountcheck[1]/100 * sum(amountdict.values())]
        minamount = min([f[0] for f in Counter(amountdict).most_common(mostfrequent) if f[0] == mostcommon[-1]])
    else:
        mostcommon = [f[0] for f in Counter(amountdict).most_common(mostfrequent)]
        minamount = min([f[1] for f in Counter(amountdict).most_common(mostfrequent) if f[0] == mostcommon[-1]])

    Seen = defaultdict(list)
    for seq in seen:
        if seq in mostcommon:
            random.shuffle(seen[seq])
            Seen[f] = seen[seq][0:minamount]

    fullnamelist = list(chain.from_iterable([Seen[f] for f in Seen]))
    quickcheckdict = defaultdict(list)
    for names in fullnamelist:
        quickcheckdict[names.split('-')[-1][:-5]].append(f)

    i = 0
    for line in open(f'{path}subimage.json', 'r'):
        name = line.split(', "')[0][11:-1]

        if name in quickcheckdict[name.split('-')[-1][:-5]]:
            data = loads(line)
            data['Sequence_class'] = str([str(index) for index in Seen].index(data['Sequence']))
            outfile.write(json.dumps(data) + '\n')
            i += 1
    print(f'{i} lines written to filtered version')
    outfile.close()


def filterlength(path, outfile, getabovehere, ms1size):
    seen = [json.loads(line)['Length'] for line in open(f'{path}subimage.json') if
            'Length' in json.loads(line)]
    Seen = np.unique(seen)
    Seen = sorted(Seen)
    for line in open(f'{path}subimage.json', 'r'):
        data = json.loads(line)

        if 'size' in data and data['size'] == str(ms1size) and 'Length' in data and \
                data['Length'] in Seen and 'Score' in data and float(data['Score']) > getabovehere:
            data['Length_class'] = Seen.index(data['Length'])
            outfile.write(json.dumps(data) + '\n')
    outfile.close()


def filtercharge(path, outfile, getabovehere, ms1size, scorecheck, amountcheck):
    print('seperating data', end='\r')
    start = time.time()
    seen = defaultdict(list)
    for line in open(f'{path}subimage.json'):
        jsonlist = line.split(', "')
        if scorecheck[0]:
            keys = ['image', 'size', 'Charge', 'Score']
        else:
            keys = ['image', 'size', 'Charge']
        values = [key.split('"')[-2] for key in jsonlist for part in keys if part in key and 'DP' not in key]

        if scorecheck[0] and float(values[3]) > getabovehere and len(values) == 4 and values[1] == ms1size:
            seen[values[2]].append(values[0])
        elif not scorecheck[0] and len(values) == 3 and values[1] == ms1size:
            seen[values[2]].append(values[0])

    if amountcheck[0]:
        amounts = defaultdict(list)
        for f in seen:
            amounts[f] = len(seen[f])
        minamount = min(f for f in amounts.values() if f >= amountcheck[1]/100 * sum(amounts.values()))

        Seen = defaultdict(list)
        for f in seen:
            if len(seen[f]) >= minamount:
                random.shuffle(seen[f])
                Seen[f] = seen[f][0:minamount]
    else:
        Seen = seen

    fullnamelist = list(chain.from_iterable([Seen[f] for f in Seen]))
    quickcheckdict = defaultdict(list)
    for f in fullnamelist:
        quickcheckdict[f.split('-')[-1][:-5]].append(f)
    end = time.time()
    print(f'seperating data complete - {end - start} sec')

    print('writing to file', end='\r')
    start = time.time()
    i = 0
    for line in open(f'{path}subimage.json'):
        name = line.split(', "')[0][11:-1]
        if name in quickcheckdict[name.split('-')[-1][:-5]]:
            data = loads(line)
            data['Charge_class'] = str([str(index) for index in Seen].index(data['Charge']))
            outfile.write(json.dumps(data) + '\n')
            i += 1
    outfile.close()
    end = time.time()
    print(f'writing to file complete - {end - start} sec')
    print(f'Length of filtered file: {i}')


def filterptm(path, outfile, getabovehere, ms1size):
    seen = defaultdict(list)
    for line in open(f'{path}subimage.json'):
        if 'Modifications' in json.loads(line):
            if not json.loads(line)['Modifications'] == 'Unmodified':
                seen[0].append(json.loads(line)['image'])
            else:
                seen[1].append(json.loads(line)['image'])

    amounts = [len(seen[f]) for f in seen]
    Seen = defaultdict(list)
    for f in seen:
        random.shuffle(seen[f])
        Seen[f] = seen[f][0:min(amounts)]

    for line in open(f'{path}subimage.json', 'r'):
        data = json.loads(line)
        if 'size' in data and data['size'] == str(
                ms1size) and 'Modifications' in data:
            if data['image'] in Seen[0]:
                data['Modi_class'] = 0
                outfile.write(json.dumps(data) + '\n')
            elif data['image'] in Seen[1]:
                data['Modi_class'] = 1
                outfile.write(json.dumps(data) + '\n')
    outfile.close()


def filter(path, file):
    outfile = open(f'{path}{str(file)}_filtered.json', 'w')
    if file == 'subimage':
        if sys.argv[2] == 'dataframe':
            jsonlist = []
            for line in open(f'{path}subimage.json'):
                data = json.loads(line)
                jsonlist.append(data)
            df = pd.DataFrame(jsonlist)
            pd.DataFrame.to_csv(df, f'{path}subimage.csv')

        if sys.argv[2] == 'combine':
            combine(path)
            quit()

        elif sys.argv[2] == 'clear':
            clear(path)
            quit()

        with open('config.json') as json_file:
            data = json.load(json_file)

        scorecheck = data['filterscore']
        scorecheck[0] = scorecheck[0] == 'True'
        amountcheck = data['filteramount']
        amountcheck[0] = amountcheck[0] == 'True'

        print('Getting sizes and scores', end='\r')
        start = time.time()
        output = getsizeandscore(path, scorecheck)
        stop = time.time()
        print(f'Getting sizes and scores complete - {stop - start} seconds elapsed')
        size = output[0]
        scorepercentile = output[1]

        print('Creating filtered version', end='\r')
        start = time.time()
        if sys.argv[2] == 'Sequence':
            filtersequence(path, outfile, scorepercentile, size, scorecheck, amountcheck)

        elif sys.argv[2] == 'Length':
            filterlength(path, outfile, scorepercentile, size, scorecheck, amountcheck)

        elif sys.argv[2] == 'PTM':
            filterptm(path, outfile, scorepercentile, size)

        elif sys.argv[2] == 'Charge':
            filtercharge(path, outfile, scorepercentile, size, scorecheck, amountcheck)
        stop = time.time()
        print(f'Creating filtered version complete - {stop-start} seconds elapsed')


    elif file == 'accessions':
        if os.path.exists(f'{path}{str(file)}_filtered.json'):
            print('Removing old filtered version')
            os.remove(f'{path}{str(file)}_filtered.json')

        lines_seen = set()
        outfile = open(f'{path}accessions_filtered.json', 'w')

        for line in open(f'{path}{str(file)}.json', 'r'):
            data = json.loads(line)
            if 'allpeptides' in data and data['allpeptides'] and 'filetypes' in data and 'raw' in data[
                'filetypes'] and line not in lines_seen:
                outfile.write(line)
                lines_seen.add(line)

                outfile.write(json.dumps(data) + '\n')
        outfile.close()
    else:
        print('no bueno')
        quit()


if __name__ == '__main__':
    # Read datapath from config file
    with open('config.json') as json_file:
        data = json.load(json_file)

    datapath = f'{data["path"]}metadata/'

    filetofilter = sys.argv[1]
    begin = time.time()
    filter(datapath, filetofilter)
    finish = time.time()
    print(f'total time elapsed: {finish - begin}')
# python3 filehandler.py filter accessions
# python3 filehandler.py filter subimage PTM/Charge/Sequence/Length
