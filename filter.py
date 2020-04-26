import glob
import json
import os
import re
import time

import numpy as np
import random
import glob
import sys
from collections import defaultdict, Counter


def mostcommon(x, amount, path):
    seen = [json.loads(line)[x] for line in open(f'{path}subimage.json') if x in json.loads(line)]
    Seen = np.unique(seen)
    a = {}
    for f in Seen:
        a[str(f)] = seen.count(f)
    seen = [f[0] for f in Counter(a).most_common(amount)]

    amounts = [len(seen[f]) for f in Seen]
    Seen = defaultdict(list)
    for f in seen:
        random.shuffle(seen[f])
        Seen[f] = seen[f][0:min(amounts)]

    return Seen


def filter(path, file):
    if file == 'subimage':
        if sys.argv[2] == 'combine':

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
                    # os.remove(imagejson)
            quit()

        elif sys.argv[2] == 'clear':
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


        # GET MOST COMMON SIZES AND SCORE PERCENTILES #
        print('Getting sizes and scores')
        outfile = open(f'{path}{str(file)}_filtered.json', 'w')
        getsizes = [lines[re.search('\[', lines).span()[0]:re.search(']', lines).span()[1]] for lines in
                    open(f'{path}subimage.json')]
        # getsizes = [json.loads(lines)['size'] for lines in open(f'{path}subimage.json') if 'size' in json.loads(lines)]
        uniquesizes = np.unique(getsizes)

        sizedict = defaultdict()
        for sizes in uniquesizes:
            sizedict[str(sizes)] = getsizes.count(sizes)

        for sizes in Counter(sizedict).most_common():
            size = str(sizes[0]).replace(" ", "")[1:-1].split(',')
            if len(size) == 3:
                ms1size = sizes[0]
                break

        getscores = [float(line[9:-1]) for lines in open(f'{path}subimage.json') for line in lines.split(', "') if
                     'score' in line.lower() and 'dp' not in line.lower()]
        getabovehere = np.percentile(getscores, 60)
        ###############################################
        print('Creating filtered version')

        if sys.argv[2] == 'Sequence':
            Seen = mostcommon('Sequence', 10)

            for line in open(f'{path}{str(file)}.json', 'r'):
                data = json.loads(line)

                if 'size' in data and data['size'] == str(ms1size) and 'Sequence' in data and data[
                    'Sequence'] in Seen and line not in lines_seen:
                    data['Seq_class'] = Seen.index(data['Sequence'])
                    lines_seen.add(line)
                    outfile.write(json.dumps(data) + '\n')
            outfile.close()


        elif sys.argv[2] == 'Length':

            seen = [json.loads(line)['Length'] for line in open(f'{path}subimage.json') if
                    'Length' in json.loads(line)]
            Seen = np.unique(seen)
            Seen = sorted(Seen)
            for line in open(f'{path}{str(file)}.json', 'r'):
                data = json.loads(line)

                if 'size' in data and data['size'] == str(ms1size) and 'Length' in data and \
                        data['Length'] in Seen and 'Score' in data and float(data['Score']) > getabovehere:
                    data['Length_class'] = Seen.index(data['Length'])
                    lines_seen.add(line)
                    outfile.write(json.dumps(data) + '\n')
            outfile.close()


        elif sys.argv[2] == 'PTM':
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

            for line in open(f'{path}{str(file)}.json', 'r'):
                data = json.loads(line)
                if 'size' in data and data['size'] == str(
                        ms1size) and 'Modifications' in data and line not in lines_seen:
                    if data['image'] in Seen[0]:
                        data['Modi_class'] = 0
                        lines_seen.add(line)
                        outfile.write(json.dumps(data) + '\n')
                    elif data['image'] in Seen[1]:
                        data['Modi_class'] = 1
                        lines_seen.add(line)
                        outfile.write(json.dumps(data) + '\n')
            outfile.close()


        elif sys.argv[2] == 'Charge':

            seen = defaultdict(list)
            for line in open(f'{path}subimage.json'):
                a = line.split(', "')
                checklist = []
                for f in a:
                    if len(checklist) == 4:
                        break

                    if 'image' in f.lower():
                        name = f[11:-1]
                        checklist.append(True)
                    elif 'charge' in f.lower():
                        charge = int(f[10:-1])
                        checklist.append(True)
                    elif 'size' in f.lower():
                        size = f[7:]
                        checklist.append(True)
                    elif 'score' in f.lower() and 'dp' not in f.lower():
                        score = float(f[9:-1])
                        checklist.append(True)

                if score >= getabovehere and size == ms1size and len(checklist) == 4:
                    seen[charge].append(name)

            amounts = defaultdict(list)
            for f in seen:
                amounts[f] = len(seen[f])
            minamount = min(f for f in amounts.values() if f >= (0.1 * sum(amounts.values())))

            Seen = defaultdict(list)
            for f in seen:
                if len(seen[f]) >= minamount:
                    random.shuffle(seen[f])
                    Seen[f] = seen[f][0:minamount]

            start = time.time()
            i = 0
            # namesseen = []
            for line in open(f'{path}subimage.json'):
                a = line.split(', "')
                checklist = []
                for f in a:
                    if len(checklist) == 2:
                        break
                    if 'image' in f.lower():
                        name = f[11:-1]
                        checklist.append(True)
                    elif 'charge' in f.lower():
                        charge = int(f[10:-1])
                        checklist.append(True)

                if name in Seen[charge] and len(checklist) == 2:# and name not in namesseen:
                    outfile.write(line + '\n')
                    # namesseen.append(name)
                    i += 1
            outfile.close()
            print(f'Length of filtered file: {i}')
            end = time.time()
            print(end - start)

            start = time.time()
            i = 0
            # namesseen = []
            for line in open(f'{path}subimage.json'):
                data = json.loads(line)
                if 'Charge' in data and data['image'] in Seen[data['Charge']]
                    outfile.write(json.dumps(data) + '\n')
                    # namesseen.append(name)
                    i += 1
            outfile.close()
            print(f'Length of filtered file: {i}')
            end = time.time()
            print(end - start)


    elif file == 'accessions':
        if os.path.exists(f'{path}{str(file)}_filtered.json'):
            print('Removing old filtered version')
            os.remove(f'{path}{str(file)}_filtered.json')

        lines_seen = set()
        outfile = open(f'{path}{str(file)}_filtered.json', 'w')

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
    filter(datapath, filetofilter)

# python3 filehandler.py filter accessions
# python3 filehandler.py filter subimage PTM/Charge/Sequence/Length
