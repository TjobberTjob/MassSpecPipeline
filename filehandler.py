import json
import os
import sys
from collections import Counter
import numpy as np


def filter(path, file):
    if os.path.exists(f'{path}{str(file)}_filtered.json'):
        print('Removing old filtered version')
        os.remove(f'{path}{str(file)}_filtered.json')

    # Used to get only most abundant classes
    seen = [json.loads(line)['Length'] for line in open(f'{path}{str(file)}.json') if 'Length' in json.loads(line)]
    Seen = np.unique(seen)
    a = {}
    for f in Seen:
        a[str(f)] = seen.count(f)
    Seen = [f[0] for f in Counter(a).most_common(4)]

    lines_seen = set()
    i = 0
    outfile = open(f'{path}{str(file)}_filtered.json', 'w')
    for line in open(f'{path}{str(file)}.json', 'r'):
        data = json.loads(line)

        if 'size' in data and data['size'] == [166, 66, 4] and 'Length' in data and data['Length'] in Seen:
            outfile.write(line)
            i += 1

        # if 'allpeptides' in data and data['allpeptides'] and 'filetypes' in data and 'raw' in data['filetypes'] and line not in lines_seen:  ### FILTER HERE ###
        #     outfile.write(line)
        #     lines_seen.add(line)
    outfile.close()
    print(i)


def moveserver(path, tarpath, ssh):
    os.chdir(path)
    os.system(f'tar -c . | ssh {ssh} "tar -xvf - -C  {tarpath}"')


if __name__ == '__main__':
    # Read datapath from config file
    with open('config.json') as json_file:
        data = json.load(json_file)

    datapath = f'{data["path"]}metadata/'

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
