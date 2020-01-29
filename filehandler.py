import json
import os
import pickle
import sys
from collections import Counter, defaultdict
import random

import numpy as np


def debugger(path):
    f = open(f'{path}debugger.txt', 'r')
    debugger = f.readlines()
    seen = [line[14:-2] for line in debugger]
    Seen = np.unique(seen)

    for f in Seen:
        print(f'error: {f}', f'\noccourences: {len([g for g in debugger if g[14:-2] == f])}',
              f'\naccessions: {[g[2:11] for g in debugger if g[14:-2] == f]}\n')

    brokenAccessions = [f[2:11] for f in debugger if f[14:-2] == "KeyError('Sequence',)" or f[14:-2] == "ZeroDivisionError('integer division or modulo by zero',)" or f[14:-2] == "OSError(36, 'File name too long')" or f[14:-2] == "OSError(36, 'File name too long')" or f[14:-2] == "BadZipFile('File is not a zip file',)"]
    print(brokenAccessions)
    print(len(debugger))

    with open(f'{path}brokenlinks.txt', "wb") as pa:
        pickle.dump(brokenAccessions, pa)


def filter(path, file):
    if file == 'subimage':
        if os.path.exists(f'{path}subimage_filtered.json'):
            print('Removing old filtered version')
            os.remove(f'{path}{str(file)}_filtered.json')

        ####### Used to get only most abundant classes #######
        # seen = [json.loads(line)['Sequence'] for line in open(f'{path}subimage.json') if 'Sequence' in json.loads(line)]
        # Seen = np.unique(seen)
        # a = {}
        # for f in Seen:
        #     a[str(f)] = seen.count(f)
        # Seen = [f[0] for f in Counter(a).most_common(10)]

        seen = defaultdict(list)
        for line in open(f'{path}subimage.json'):
            if 'Modifications' in json.loads(line):
                if json.loads(line)['Modifications'] == 'Unmodified':
                    seen[0].append(json.loads(line)['image'])
                else:
                    seen[1].append(json.loads(line)['image'])
        amounts = [len(seen[f]) for f in seen]
        lowestamount = min(amounts)
        Seen = defaultdict(list)
        for f in seen:
            Seen[f].append(random.shuffle(seen[f])[0:lowestamount])
        ####### Used to get only most abundant classes #######

        lines_seen = set()
        outfile = open(f'{path}{str(file)}_filtered.json', 'w')

        for line in open(f'{path}{str(file)}.json', 'r'):
            data = json.loads(line)

            # if 'size' in data and data['size'] == [166, 66, 4] and 'Sequence' in data and data['Sequence'] in Seen:
            #     data['Seq_class'] = Seen.index(data['Sequence'])
            if 'size' in data and data['size'] == [166, 66, 4] and 'Modifications' in data and line not in lines_seen:
                if data['Modifications'] == 'Unmodified' and data['image'] in Seen[0]:
                    data['Modi_class'] = 0
                elif data['image'] in Seen[1]:
                    data['Modi_class'] = 1
                lines_seen.add(line)

                outfile.write(json.dumps(data) + '\n')
        outfile.close()

    elif file == 'accession':
        if os.path.exists(f'{path}accession_filtered.json'):
            print('Removing old filtered version')
            os.remove(f'{path}{str(file)}_filtered.json')

        lines_seen = set()
        outfile = open(f'{path}{str(file)}_filtered.json', 'w')

        for line in open(f'{path}{str(file)}.json', 'r'):
            data = json.loads(line)
            if 'allpeptides' in data and data['allpeptides'] and 'filetypes' in data and 'raw' in data['filetypes'] and line not in lines_seen:  ### FILTER HERE ###
                outfile.write(line)
                lines_seen.add(line)

                outfile.write(json.dumps(data) + '\n')
        outfile.close()
    else:
        print('nobueno')
        quit()


def moveserver(path, tarpath, ssh):
    os.chdir(path)
    os.system(f'tar -c . | ssh {ssh} "tar -xvf - -C  {tarpath}"')


if __name__ == '__main__':
    # Read datapath from config file
    with open('config.json') as json_file:
        data = json.load(json_file)

    datapath = f'{data["path"]}metadata/'

    if sys.argv[1] == 'debugger':
        debugger(datapath)

    if sys.argv[1] == 'filter':
        filetofilter = sys.argv[2]
        filter(datapath, filetofilter)

    elif sys.argv[1] == 'move':
        path = sys.argv[2]
        tarpath = sys.argv[3]
        ssh = sys.argv[4]
        moveserver(path, tarpath, ssh)

# python3 filehandler.py filter accessions
# python3 filehandler.py move /data/ProteomeToolsRaw/images/ /home/tochr15/images/ tochr15@yeast.imada.sdu.dk
# python3 filehandler.py move /data/ProteomeToolsRaw/metadata/ /home/tochr15/metadata/ tochr15@yeast.imada.sdu.dk
