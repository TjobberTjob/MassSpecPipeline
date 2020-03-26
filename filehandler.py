import json
import os
import pickle
import sys
from collections import Counter, defaultdict
import random
import numpy as np


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
                if json.loads(line)['Modifications'] == 'Oxidiation (M)':
                    seen[1].append(json.loads(line)['image'])
                else:
                    seen[0].append(json.loads(line)['image'])
        amounts = [len(seen[f]) for f in seen]
        Seen = defaultdict(list)
        for f in seen:
            random.shuffle(seen[f])
            Seen[f] = seen[f][0:min(amounts)]
        ####### Used to get only most abundant classes #######

        lines_seen = set()
        outfile = open(f'{path}{str(file)}_filtered.json', 'w')

        for line in open(f'{path}{str(file)}.json', 'r'):
            data = json.loads(line)

            if 'size' in data and data['size'] == [166, 66, 4] and 'Sequence' in data and data['Sequence'] in Seen and line not in lines_seen:
                # data['Seq_class'] = Seen.index(data['Sequence'])
                # lines_seen.add(line)
                # outfile.write(json.dumps(data) + '\n')
            if 'size' in data and data['size'] == [166, 66, 4] and 'Modifications' in data and line not in lines_seen:
                if data['image'] in Seen[0]:
                    data['Modi_class'] = 0
                    lines_seen.add(line)
                    outfile.write(json.dumps(data) + '\n')
                elif data['image'] in Seen[1]:
                    data['Modi_class'] = 1
                    lines_seen.add(line)
                    outfile.write(json.dumps(data) + '\n')
        outfile.close()

    elif file == 'accessions':
        if os.path.exists(f'{path}{str(file)}_filtered.json'):
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
        print('no bueno')
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
