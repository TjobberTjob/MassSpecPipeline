import json
import os
import sys
from collections import Counter
import numpy as np

def debugger(path):
    f = open(f'{path}debugger.txt', 'r')
    debugger = f.readlines()
    Seen = [line[14:-2] for line in debugger]
    Seen = np.unique(Seen)

    for f in Seen:
        print(f'error: {f}', f'\noccourences: {len([g[2:12] for g in debugger if g[14:-2] == f])}',
              f'\naccessions: {[g[2:12] for g in debugger if g[14:-2] == f]}\n')


def filter(path, file):
    if os.path.exists(f'{path}{str(file)}_filtered.json'):
        print('Removing old filtered version')
        os.remove(f'{path}{str(file)}_filtered.json')

    # Used to get only most abundant classes
    seen = [json.loads(line)['Sequence'] for line in open(f'{path}{str(file)}.json') if 'Sequence' in json.loads(line)]
    Seen = np.unique(seen)
    a = {}
    for f in Seen:
        a[str(f)] = seen.count(f)
    Seen = [f[0] for f in Counter(a).most_common(4)]

    # lines_seen = set()
    outfile = open(f'{path}{str(file)}_filtered.json', 'w')

    for line in open(f'{path}{str(file)}.json', 'r'):
        data = json.loads(line)

        # Filter for classification
        if 'size' in data and data['size'] == [166, 66, 4] and 'Sequence' in data and data['Sequence'] in Seen:
            data['Seq_class'] = Seen.index(data['Sequence'])
            outfile.write(json.dumps(data) + '\n')

        # # filter metadata for extractor
        # if 'allpeptides' in data and data['allpeptides'] and 'filetypes' in data and 'raw' in data['filetypes'] and line not in lines_seen:  ### FILTER HERE ###
        #     outfile.write(line)
        #     lines_seen.add(line)

    # outfile.close()


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
