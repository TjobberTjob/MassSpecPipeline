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


def getclass(word, string):
    if word == '"size"':
        output = string[re.search('\[', string).span()[0]:re.search(']', string).span()[1]]
    else:
        ab = [f for f in [m.start() for m in re.finditer('"', string)] if f > re.search(word, string).span()[1]]
        output = string[ab[0]+1: ab[1]]
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
    getsizes = [lines[re.search('\[', lines).span()[0]:re.search(']', lines).span()[1]] for lines in open(f'{path}subimage.json')]

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
    print('Getting sizes and scores', end='\r')
    start = time.time()

    getsizes = [lines[re.search('\[', lines).span()[0]:re.search(']', lines).span()[1]] for lines in open(f'{path}subimage.json')]
    uniquesizes = np.unique(getsizes)

    sizedict = defaultdict()
    for sizes in uniquesizes:
        sizedict[str(sizes)] = getsizes.count(sizes)

    for sizes in Counter(sizedict).most_common():
        size = str(sizes[0]).replace(" ", "")[1:-1].split(',')
        if len(size) == 3:
            ms1size = sizes[0]
            break
    if scorecheck:
        getscores = [float(line[9:-1]) for lines in open(f'{path}subimage.json') for line in lines.split(', "') if 'score' in line.lower() and 'dp' not in line.lower() and str(ms1size) in lines.lower()]
        getabovehere = np.percentile(getscores, 1)
    else:
        getabovehere = 'Fuckthis'

    stop = time.time()
    print(f'Getting sizes and scores complete - {stop-start} sec')
    return ms1size, getabovehere


def filtersequence(path, outfile, getabovehere, ms1size):
    seen = [json.loads(line)['Sequence'] for line in open(f'{path}subimage.json') if 'Sequence' in json.loads(line)]
    Seen = np.unique(seen)
    a = {}
    for f in Seen:
        a[str(f)] = seen.count(f)
    seen = [f[0] for f in Counter(a).most_common(10)]

    amounts = [len(seen[f]) for f in Seen]
    Seen = defaultdict(list)
    for f in seen:
        random.shuffle(seen[f])
        Seen[f] = seen[f][0:min(amounts)]

    for line in open(f'{path}subimage.json', 'r'):
        data = json.loads(line)

        if 'size' in data and data['size'] == str(ms1size) and 'Sequence' in data and data[
            'Sequence'] in Seen:
            data['Seq_class'] = Seen.index(data['Sequence'])
            outfile.write(json.dumps(data) + '\n')
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


def filtercharge(path, outfile, getabovehere, ms1size, scorecheck):
    print('seperating data', end='\r')
    start = time.time()
    seen = defaultdict(list)
    for line in open(f'{path}subimage.json'):
        checklist = []
        if scorecheck:
            checklen = 4
        else:
            checklen = 3
        a = line.split(', "')
        for f in a:
            if len(checklist) == checklen:
                break
            elif 'image' in f:
                name = f[11:-1]
                checklist.append(True)
            elif 'size' in f:
                size = str(f[7:])
                checklist.append(True)
            elif 'Charge' in f:
                charge = int(f[-2:-1])
                checklist.append(True)
            elif 'Score' in f and 'DP' not in f and scorecheck:
                score = float(f[11:-1])
                checklist.append(True)

        if len(checklist) == checklen and size == ms1size:
            if scorecheck and score >= getabovehere:
                seen[charge].append(name)
            elif not scorecheck:
                seen[charge].append(name)


    amounts = defaultdict(list)
    for f in seen:
        amounts[f] = len(seen[f])
    minamount = min(f for f in amounts.values() if f >= 0.25 * sum(amounts.values()))

    Seen = defaultdict(list)
    for f in seen:
        if len(seen[f]) >= minamount:
            random.shuffle(seen[f])
            Seen[f] = seen[f][0:minamount]
    end = time.time()

    fullnamelist = list(chain.from_iterable([Seen[f] for f in Seen]))
    Seen2 = defaultdict(list)
    for f in fullnamelist:
        Seen2[f.split('-')[-1][:-5]].append(f)
    print(f'seperating data complete - {end - start} sec')


    print('writing to file', end='\r')
    start = time.time()
    i = 0
    for line in open(f'{path}subimage.json'):
        a = line.split(', "')
        name = a[0][11:-1]
        if name in Seen2[name.split('-')[-1][:-5]]:
            data = json.loads(line)
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

        scorecheck = False
        output = getsizeandscore(path, scorecheck)
        size = output[0]
        scorepercentile = output[1]

        print('Creating filtered version')
        if sys.argv[2] == 'Sequence':
            filtersequence(path, outfile, scorepercentile, size)
            quit()

        elif sys.argv[2] == 'Length':
            filterlength(path, outfile, scorepercentile, size)
            quit()

        elif sys.argv[2] == 'PTM':
            filterptm(path, outfile, scorepercentile, size)
            quit()

        elif sys.argv[2] == 'Charge':
            filtercharge(path, outfile, scorepercentile, size, scorecheck)
            quit()


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
